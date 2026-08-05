#!/bin/bash
#SBATCH --job-name=nanocirc-readsim
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=main
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --exclusive

# --exclusive: concurrent NanoSim runs on the same node can die silently.

# Generates NanoSim long reads for a nanocirc-sim dataset, pools/labels
# them, then runs the mandatory BSJ check.
#
# Required env vars (pass via `sbatch --export=ALL,VAR=value,...`):
#   DATASET_DIR      output dir from `nanocirc-sim simulate` (has circRNAs.fa,
#                    circRNAs.mature.fa, abundances.tsv, circRNA_metadata.tsv)
#   GENOME_FA        reference genome FASTA
#   GTF              gene annotation GTF
#   TRANSCRIPTOME_FA spliced transcriptome FASTA (scripts/build_transcriptome.py)
#   CONTROL_FASTQ    real long-read FASTQ, for NanoSim's error-model characterization
#   N_READS          number of circRNA-supporting reads to simulate (e.g. 1000000)
#   OUTDIR           where to write reads + BSJ check results
#
# Optional:
#   PROJECT_ROOT      default: current directory, so run from the repo root
#                     or export PROJECT_ROOT before submitting
#   NANOSIM_SIF        default $PROJECT_ROOT/simulator/new/containers/nanosim.sif
#   NANOCIRC_SIM_SIF    default $PROJECT_ROOT/simulator/new/containers/nanocirc-sim.sif

set -euo pipefail

: "${DATASET_DIR:?must set DATASET_DIR}"
: "${GENOME_FA:?must set GENOME_FA}"
: "${GTF:?must set GTF}"
: "${TRANSCRIPTOME_FA:?must set TRANSCRIPTOME_FA}"
: "${CONTROL_FASTQ:?must set CONTROL_FASTQ}"
: "${N_READS:?must set N_READS}"
: "${OUTDIR:?must set OUTDIR}"
NANOSIM_SEED="${NANOSIM_SEED:-1}"

PROJECT_ROOT="$(cd "${PROJECT_ROOT:-.}" && pwd)"
NANOSIM_SIF="${NANOSIM_SIF:-$PROJECT_ROOT/simulator/new/containers/nanosim.sif}"
NANOCIRC_SIM_SIF="${NANOCIRC_SIM_SIF:-$PROJECT_ROOT/simulator/new/containers/nanocirc-sim.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUTDIR"
export MPLCONFIGDIR="$OUTDIR/.mplcache"
export APPTAINERENV_MPLCONFIGDIR=/tmp/mplcache
mkdir -p "$MPLCONFIGDIR"

BIND="--bind $PROJECT_ROOT:$PROJECT_ROOT"
APPTAINER_NANOSIM="apptainer exec $BIND $NANOSIM_SIF"
APPTAINER_NANOCIRC="apptainer exec $BIND $NANOCIRC_SIM_SIF"

echo "=== Step 1: NanoSim characterization (control=$CONTROL_FASTQ) ==="
$APPTAINER_NANOSIM read_analysis.py transcriptome \
    -i "$CONTROL_FASTQ" \
    -rg "$GENOME_FA" \
    -rt "$TRANSCRIPTOME_FA" \
    -a minimap2 \
    -o "$OUTDIR/control" \
    -t "$THREADS" \
    -annot "$GTF" \
    --no_intron_retention \
    --fastq

echo "=== Step 2: Simulate $N_READS circRNA-supporting reads ==="
$APPTAINER_NANOSIM simulator.py transcriptome \
    -rg "$GENOME_FA" \
    -rt "$DATASET_DIR/circRNAs.fa" \
    -c "$OUTDIR/control" \
    -o "$OUTDIR/circ" \
    -n "$N_READS" \
    -b guppy \
    -e "$DATASET_DIR/abundances.tsv" \
    --seed "$NANOSIM_SEED" \
    -t "$THREADS" \
    --fastq \
    --no_model_ir

echo "=== Step 3: Pool + label reads ==="
POOL_DIR="$OUTDIR/pooled"
mkdir -p "$POOL_DIR"
COMBINED="$POOL_DIR/combined_reads.fastq"
: > "$COMBINED"
for f in "$OUTDIR/circ_aligned_reads.fastq" "$OUTDIR/circ_unaligned_reads.fastq"; do
    if [ -f "$f" ]; then
        awk '{ if (NR % 4 == 1) { print substr($0,1,1) substr($0,2) "_circ" } else { print } }' "$f" >> "$COMBINED"
    else
        echo "WARNING: expected NanoSim output missing: $f"
    fi
done
gzip -c "$COMBINED" > "$COMBINED.gz"
echo "Pooled reads: $(wc -l < "$COMBINED") lines ($(( $(wc -l < "$COMBINED") / 4 )) reads)"

echo "=== Step 4: Mandatory BSJ check ==="
$APPTAINER_NANOCIRC nanocirc-sim check-bsj \
    --fastq "$COMBINED" \
    --mature-fasta "$DATASET_DIR/circRNAs.mature.fa" \
    --metadata "$DATASET_DIR/circRNA_metadata.tsv" \
    --outdir "$OUTDIR/bsj_check" \
    --threads "$THREADS"

echo "=== Done. Reads: $POOL_DIR, BSJ check: $OUTDIR/bsj_check ==="
