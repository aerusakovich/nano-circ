#!/bin/bash
#SBATCH --job-name=nanocirc-linear
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=main
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

# Adds linear background reads to a dataset that already had circRNA reads
# simulated by run_read_simulation.sbatch, then reruns check-bsj.
#
# Required env vars:
#   DATASET_DIR       nanocirc-sim simulate output dir (has circRNAs.mature.fa, circRNA_metadata.tsv)
#   GENOME_FA         reference genome FASTA
#   TRANSCRIPTOME_FA  spliced transcriptome FASTA (scripts/build_transcriptome.py)
#   CONTROL_FASTQ     same real control FASTQ used for the circ run
#   N_LINEAR_READS    number of linear background reads to simulate (e.g. 1000000)
#   CIRC_OUTDIR       the circ run's OUTDIR (has control_* profile + pooled/combined_reads.fastq)
#   OUTDIR            where to write linear reads + merged pool + the rerun bsj_check
#
# Optional: PROJECT_ROOT, NANOSIM_SIF, NANOCIRC_SIM_SIF, NANOSIM_SEED (default 2)

set -euo pipefail

: "${DATASET_DIR:?must set DATASET_DIR}"
: "${GENOME_FA:?must set GENOME_FA}"
: "${TRANSCRIPTOME_FA:?must set TRANSCRIPTOME_FA}"
: "${CONTROL_FASTQ:?must set CONTROL_FASTQ}"
: "${N_LINEAR_READS:?must set N_LINEAR_READS}"
: "${CIRC_OUTDIR:?must set CIRC_OUTDIR}"
: "${OUTDIR:?must set OUTDIR}"
NANOSIM_SEED="${NANOSIM_SEED:-2}"

PROJECT_ROOT="${PROJECT_ROOT:-/nfs/sto-2/projects/cnt-000070/canine-genetics/nanocirc}"
NANOSIM_SIF="${NANOSIM_SIF:-$PROJECT_ROOT/simulator/new/containers/nanosim.sif}"
NANOCIRC_SIM_SIF="${NANOCIRC_SIM_SIF:-$PROJECT_ROOT/simulator/new/containers/nanocirc-sim.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUTDIR"
export APPTAINERENV_MPLCONFIGDIR=/tmp/mplcache

BIND="--bind /nfs/sto-2:/nfs/sto-2"
APPTAINER_NANOSIM="apptainer exec $BIND $NANOSIM_SIF"
APPTAINER_NANOCIRC="apptainer exec $BIND $NANOCIRC_SIM_SIF"

echo "=== Step 1: Quantify real linear transcript abundance from control ==="
$APPTAINER_NANOSIM read_analysis.py quantify \
    -e trans \
    -i "$CONTROL_FASTQ" \
    -rt "$TRANSCRIPTOME_FA" \
    -t "$THREADS" \
    -o "$OUTDIR/expression"

echo "=== Step 2: Simulate $N_LINEAR_READS linear background reads ==="
$APPTAINER_NANOSIM simulator.py transcriptome \
    -rg "$GENOME_FA" \
    -rt "$TRANSCRIPTOME_FA" \
    -c "$CIRC_OUTDIR/control" \
    -o "$OUTDIR/linear" \
    -n "$N_LINEAR_READS" \
    -b guppy \
    -e "$OUTDIR/expression_transcriptome_quantification.tsv" \
    --seed "$NANOSIM_SEED" \
    -t "$THREADS" \
    --fastq \
    --no_model_ir

echo "=== Step 3: Merge linear reads into the existing pooled FASTQ ==="
POOL_DIR="$CIRC_OUTDIR/pooled"
MERGED="$OUTDIR/combined_reads_with_linear.fastq"
cp "$POOL_DIR/combined_reads.fastq" "$MERGED"
for f in "$OUTDIR/linear_aligned_reads.fastq" "$OUTDIR/linear_unaligned_reads.fastq"; do
    if [ -f "$f" ]; then
        awk '{ if (NR % 4 == 1) { print substr($0,1,1) substr($0,2) "_linear" } else { print } }' "$f" >> "$MERGED"
    else
        echo "WARNING: expected NanoSim output missing: $f"
    fi
done
gzip -c "$MERGED" > "$MERGED.gz"
echo "Merged reads: $(( $(wc -l < "$MERGED") / 4 )) reads total"

echo "=== Step 4: Rerun mandatory BSJ check on full circ+linear mix ==="
$APPTAINER_NANOCIRC nanocirc-sim check-bsj \
    --fastq "$MERGED" \
    --mature-fasta "$DATASET_DIR/circRNAs.mature.fa" \
    --metadata "$DATASET_DIR/circRNA_metadata.tsv" \
    --outdir "$OUTDIR/bsj_check_full" \
    --threads "$THREADS"

echo "=== Done. Merged reads: $MERGED, BSJ check: $OUTDIR/bsj_check_full ==="
