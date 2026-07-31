#!/bin/bash
#SBATCH --job-name=nanocirc-linear-matched
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=main
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --exclusive

# Adds exactly N_LINEAR_READS linear background reads to a confirmed-only
# circular FASTQ, for a precise 1:1 circular:linear split. Normally
# self-chained by run_pool_select_and_filter.sbatch; can run standalone.
#
# Required env vars:
#   DATASET_DIR           nanocirc-sim simulate output dir
#   GENOME_FA
#   TRANSCRIPTOME_FA
#   CONTROL_FASTQ          same real control FASTQ used for the circ run
#   CIRC_OUTDIR            run_read_simulation.sbatch's OUTDIR (has control_* profile)
#   CONFIRMED_CIRC_FASTQ   confirmed-only circular FASTQ (run_topup_and_filter.sbatch's output)
#   N_LINEAR_READS         exact number of linear reads to simulate (== confirmed circular count)
#   OUTDIR                 where to write the final balanced dataset
#
# Optional: PROJECT_ROOT, NANOSIM_SIF, NANOCIRC_SIM_SIF, NANOSIM_SEED (default 3),
#   RUN_CIRCRNAS_BED / RUN_GROUND_TRUTH_TSV (copied alongside the output if set)

set -euo pipefail

: "${DATASET_DIR:?must set DATASET_DIR}"
: "${GENOME_FA:?must set GENOME_FA}"
: "${TRANSCRIPTOME_FA:?must set TRANSCRIPTOME_FA}"
: "${CONTROL_FASTQ:?must set CONTROL_FASTQ}"
: "${CIRC_OUTDIR:?must set CIRC_OUTDIR}"
: "${CONFIRMED_CIRC_FASTQ:?must set CONFIRMED_CIRC_FASTQ}"
: "${N_LINEAR_READS:?must set N_LINEAR_READS}"
: "${OUTDIR:?must set OUTDIR}"
NANOSIM_SEED="${NANOSIM_SEED:-3}"
RUN_CIRCRNAS_BED="${RUN_CIRCRNAS_BED:-}"
RUN_GROUND_TRUTH_TSV="${RUN_GROUND_TRUTH_TSV:-}"

PROJECT_ROOT="${PROJECT_ROOT:-/nfs/sto-2/projects/cnt-000070/canine-genetics/nanocirc}"
NANOSIM_SIF="${NANOSIM_SIF:-$PROJECT_ROOT/simulator/new/containers/nanosim.sif}"
NANOCIRC_SIM_SIF="${NANOCIRC_SIM_SIF:-$PROJECT_ROOT/simulator/new/containers/nanocirc-sim.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUTDIR"
export APPTAINERENV_MPLCONFIGDIR=/tmp/mplcache

BIND="--bind /nfs/sto-2:/nfs/sto-2"
APPTAINER_NANOSIM="apptainer exec $BIND $NANOSIM_SIF"
APPTAINER_NANOCIRC="apptainer exec $BIND $NANOCIRC_SIM_SIF"

n_confirmed=$(( $(wc -l < "$CONFIRMED_CIRC_FASTQ") / 4 ))
echo "=== Confirmed circular reads to match: $n_confirmed (requested N_LINEAR_READS=$N_LINEAR_READS) ==="

echo "=== Step 1: Quantify real linear transcript abundance from control ==="
$APPTAINER_NANOSIM read_analysis.py quantify \
    -e trans \
    -i "$CONTROL_FASTQ" \
    -rt "$TRANSCRIPTOME_FA" \
    -t "$THREADS" \
    -o "$OUTDIR/expression"

echo "=== Step 2: Simulate exactly $N_LINEAR_READS linear background reads ==="
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

echo "=== Step 3: Building the final 1:1 balanced dataset ==="
final_fastq="$OUTDIR/final_balanced_reads.fastq"
cp "$CONFIRMED_CIRC_FASTQ" "$final_fastq"
n_linear_added=0
for f in "$OUTDIR/linear_aligned_reads.fastq" "$OUTDIR/linear_unaligned_reads.fastq"; do
    if [ -f "$f" ]; then
        n=$(( $(wc -l < "$f") / 4 ))
        n_linear_added=$(( n_linear_added + n ))
        awk '{ if (NR % 4 == 1) { print substr($0,1,1) substr($0,2) "_linear" } else { print } }' "$f" >> "$final_fastq"
    fi
done
gzip -c "$final_fastq" > "$final_fastq.gz"

echo "=== Confirmed circular: $n_confirmed, linear added: $n_linear_added ==="
if [ "$n_linear_added" -ne "$n_confirmed" ]; then
    echo "WARNING: NanoSim produced $n_linear_added linear reads, not exactly the requested $N_LINEAR_READS " \
         "(it can occasionally simulate a handful fewer/more than requested). Final split is " \
         "$n_confirmed circular : $n_linear_added linear, not a perfect 1:1 match."
fi

echo "=== Step 4: Rechecking BSJ on the final balanced dataset (QC + final ground truth table) ==="
$APPTAINER_NANOCIRC nanocirc-sim check-bsj \
    --fastq "$final_fastq" \
    --mature-fasta "$DATASET_DIR/circRNAs.mature.fa" \
    --metadata "$DATASET_DIR/circRNA_metadata.tsv" \
    --outdir "$OUTDIR/bsj_check" \
    --threads "$THREADS"

if [ -n "$RUN_CIRCRNAS_BED" ] && [ -f "$RUN_CIRCRNAS_BED" ]; then
    cp "$RUN_CIRCRNAS_BED" "$OUTDIR/run_circrnas.bed"
fi
if [ -n "$RUN_GROUND_TRUTH_TSV" ] && [ -f "$RUN_GROUND_TRUTH_TSV" ]; then
    cp "$RUN_GROUND_TRUTH_TSV" "$OUTDIR/run_ground_truth.tsv"
fi

echo "=== Done. Final balanced dataset: $final_fastq ($((n_confirmed + n_linear_added)) reads total, " \
     "target 1:1 circular:linear) ==="
