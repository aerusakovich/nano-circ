#!/bin/bash
#SBATCH --job-name=nanocirc-pool-select
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=main
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --exclusive

# Builds N_RUNS independent replicate datasets from a mega-pool: pool-select,
# iterative shortfall-weighted top-up for deficient circRNAs, then draw N_RUNS
# Poisson-sampled runs and self-chain each into run_linear_matched.sbatch. See README.
#
# Required env vars:
#   DATASET_DIR    nanocirc-sim simulate output dir
#   GENOME_FA
#   TRANSCRIPTOME_FA
#   CONTROL_FASTQ  same real control FASTQ used for the circ mega-pool run
#   CIRC_OUTDIR    run_read_simulation.sbatch's OUTDIR (has control_* profile,
#                  pooled/combined_reads.fastq, bsj_check/)
#   OUTDIR         where to write selection output + the per-run datasets
#
# Optional: PROJECT_ROOT (default: current directory, so run from the repo root
#   or export PROJECT_ROOT before submitting), NANOSIM_SIF, NANOCIRC_SIM_SIF,
#   FINAL_OUTDIR (default $OUTDIR/../final_balanced, each run written to
#   $FINAL_OUTDIR/run<N>/), SKIP_CHAIN=1 to not auto-submit the linear jobs,
#   TOPUP_YIELD_RATE (default 0.15; conservative confirmed-reads-per-simulated-
#   read estimate used to size each top-up round from the remaining shortfall;
#   lower = more margin, more reads requested per round),
#   MAX_TOPUP_ROUNDS (default 8; top-up repeats against the remaining
#   shortfall until it hits 0 or stops improving, capped here),
#   MIN_TOPUP_PROGRESS (default 0.02; stop early if a round improves the
#   shortfall by less than this fraction -- the rest is likely unreachable
#   for those circRNAs, e.g. too short/unalignable, not worth chasing),
#   N_RUNS (default 3), READS_PER_RUN (default 300000)

set -euo pipefail

: "${DATASET_DIR:?must set DATASET_DIR}"
: "${GENOME_FA:?must set GENOME_FA}"
: "${TRANSCRIPTOME_FA:?must set TRANSCRIPTOME_FA}"
: "${CONTROL_FASTQ:?must set CONTROL_FASTQ}"
: "${CIRC_OUTDIR:?must set CIRC_OUTDIR}"
: "${OUTDIR:?must set OUTDIR}"
FINAL_OUTDIR="${FINAL_OUTDIR:-$OUTDIR/../final_balanced}"
TOPUP_YIELD_RATE="${TOPUP_YIELD_RATE:-0.15}"
MAX_TOPUP_ROUNDS="${MAX_TOPUP_ROUNDS:-8}"
MIN_TOPUP_PROGRESS="${MIN_TOPUP_PROGRESS:-0.02}"
N_RUNS="${N_RUNS:-3}"
READS_PER_RUN="${READS_PER_RUN:-300000}"

PROJECT_ROOT="$(cd "${PROJECT_ROOT:-.}" && pwd)"
NANOSIM_SIF="${NANOSIM_SIF:-$PROJECT_ROOT/simulator/new/containers/nanosim.sif}"
NANOCIRC_SIM_SIF="${NANOCIRC_SIM_SIF:-$PROJECT_ROOT/simulator/new/containers/nanocirc-sim.sif}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
# NOT $(dirname "$0"): under sbatch, $0 is SLURM's transient spooled copy
# of *this* script (/var/spool/slurmd/job<id>/...), which does not contain
# run_linear_matched.sbatch alongside it. Self-chaining via that path
# fails with "Unable to open file". Use the known scripts/ location instead.
SCRIPTS_DIR="${SCRIPTS_DIR:-$PROJECT_ROOT/simulator/new/scripts}"

mkdir -p "$OUTDIR"
export APPTAINERENV_MPLCONFIGDIR=/tmp/mplcache

BIND="--bind $PROJECT_ROOT:$PROJECT_ROOT"
APPTAINER_NANOSIM="apptainer exec $BIND $NANOSIM_SIF"
APPTAINER_NANOCIRC="apptainer exec $BIND $NANOCIRC_SIM_SIF"

mega_pool="$CIRC_OUTDIR/pooled/combined_reads.fastq"
mega_bsj_reads="$CIRC_OUTDIR/bsj_check/bsj_check_reads.tsv"

echo "=== Step 1: Selecting up to n_reads_target confirmed reads per circRNA from the mega-pool ==="
current_bsj_reads="$mega_bsj_reads"
current_mega_fastq="$mega_pool"
current_dir="$OUTDIR/selection"
mkdir -p "$current_dir"
$APPTAINER_NANOCIRC nanocirc-sim pool-select \
    --bsj-reads "$current_bsj_reads" \
    --metadata "$DATASET_DIR/circRNA_metadata.tsv" \
    --mega-fastq "$current_mega_fastq" \
    --outdir "$current_dir" \
    --shortfall-profile-out "$current_dir/shortfall_profile.tsv"

total_shortfall=$(tail -n +2 "$current_dir/shortfall_profile.tsv" 2>/dev/null | awk -F'\t' '{s+=$3} END{print s+0}')
echo "=== Initial shortfall: $total_shortfall reads ==="

# shortfall_profile.tsv has TPM = each circRNA's own shortfall, so
# simulator.py's abundance-proportional sampling puts more reads where the
# gap is bigger instead of spreading the budget flat across every circRNA.
# Repeats against the *remaining* shortfall each round (not the original),
# since a flat TPM-proportional draw never lands exactly on target per locus.
round=0
while [ "$total_shortfall" -gt 0 ] && [ "$round" -lt "$MAX_TOPUP_ROUNDS" ]; do
    round=$(( round + 1 ))
    echo "=== Top-up round $round: total shortfall $total_shortfall reads, shortfall-weighted top-up pass ==="
    topup_dir="$OUTDIR/topup_round${round}"
    mkdir -p "$topup_dir"

    n_topup_reads=$(awk -v s="$total_shortfall" -v y="$TOPUP_YIELD_RATE" 'BEGIN{v=s/y; if (v<1) v=1; printf "%d", v}')
    echo "=== Requesting $n_topup_reads top-up reads (shortfall/${TOPUP_YIELD_RATE} yield estimate) ==="
    $APPTAINER_NANOSIM simulator.py transcriptome \
        -rg "$GENOME_FA" \
        -rt "$DATASET_DIR/circRNAs.fa" \
        -c "$CIRC_OUTDIR/control" \
        -o "$topup_dir/topup" \
        -n "$n_topup_reads" \
        -b guppy \
        -e "$current_dir/shortfall_profile.tsv" \
        --seed "$(( 200 + round ))" \
        -t "$THREADS" \
        --fastq \
        --no_model_ir

    topup_pool="$topup_dir/topup_pool.fastq"
    : > "$topup_pool"
    for f in "$topup_dir/topup_aligned_reads.fastq" "$topup_dir/topup_unaligned_reads.fastq"; do
        if [ -f "$f" ]; then
            awk '{ if (NR % 4 == 1) { print substr($0,1,1) substr($0,2) "_circ" } else { print } }' "$f" >> "$topup_pool"
        fi
    done

    $APPTAINER_NANOCIRC nanocirc-sim check-bsj \
        --fastq "$topup_pool" \
        --mature-fasta "$DATASET_DIR/circRNAs.mature.fa" \
        --metadata "$DATASET_DIR/circRNA_metadata.tsv" \
        --outdir "$topup_dir/bsj_check" \
        --threads "$THREADS"

    combined_mega="$topup_dir/mega_plus_topup.fastq"
    cat "$current_mega_fastq" "$topup_pool" > "$combined_mega"
    combined_bsj_reads="$topup_dir/combined_bsj_reads.tsv"
    head -1 "$current_bsj_reads" > "$combined_bsj_reads"
    tail -n +2 "$current_bsj_reads" >> "$combined_bsj_reads"
    tail -n +2 "$topup_dir/bsj_check/bsj_check_reads.tsv" >> "$combined_bsj_reads"

    # Both current_mega_fastq and current_bsj_reads (the previous round's files,
    # read above) live under the previous round's topup dir. Only safe to drop
    # it now that both have been fully read into this round's own files.
    if [ "$round" -gt 1 ]; then
        rm -rf "$OUTDIR/topup_round$((round - 1))"
    fi

    next_dir="$OUTDIR/selection_round${round}"
    mkdir -p "$next_dir"
    $APPTAINER_NANOCIRC nanocirc-sim pool-select \
        --bsj-reads "$combined_bsj_reads" \
        --metadata "$DATASET_DIR/circRNA_metadata.tsv" \
        --mega-fastq "$combined_mega" \
        --outdir "$next_dir" \
        --shortfall-profile-out "$next_dir/shortfall_profile.tsv"

    current_bsj_reads="$combined_bsj_reads"
    current_mega_fastq="$combined_mega"
    current_dir="$next_dir"
    prev_shortfall="$total_shortfall"
    total_shortfall=$(tail -n +2 "$current_dir/shortfall_profile.tsv" 2>/dev/null | awk -F'\t' '{s+=$3} END{print s+0}')
    echo "=== After round $round: remaining shortfall $total_shortfall reads (was $prev_shortfall) ==="

    if [ "$total_shortfall" -gt 0 ]; then
        stop_now=$(awk -v p="$prev_shortfall" -v c="$total_shortfall" -v m="$MIN_TOPUP_PROGRESS" \
            'BEGIN{ if (p<=0) { print 0 } else { print (((p-c)/p) < m) ? 1 : 0 } }')
        if [ "$stop_now" = "1" ]; then
            echo "=== Stopping top-up early: round $round barely improved the shortfall (below ${MIN_TOPUP_PROGRESS} fraction). Remaining $total_shortfall reads are likely unreachable for these circRNAs (e.g. too short/unalignable for this protocol), not worth chasing further. ==="
            break
        fi
    fi
done

if [ "$total_shortfall" -gt 0 ]; then
    echo "=== After $round top-up round(s): $total_shortfall reads still short across the catalog (reported, not hidden) ==="
else
    echo "=== Every circRNA reached its full n_reads_target after $round top-up round(s) ==="
fi

final_selected="$current_dir/selected_reads.fastq"
bsj_reads_for_subsample="$current_bsj_reads"
mega_fastq_for_subsample="$current_mega_fastq"

n_circ_full_pool=$(( $(wc -l < "$final_selected") / 4 ))
echo "=== Full proportional pool: $n_circ_full_pool confirmed circular reads across all circRNAs ==="
echo "$n_circ_full_pool" > "$OUTDIR/n_confirmed_circular_reads_full_pool.txt"
if [ "$n_circ_full_pool" -lt "$READS_PER_RUN" ]; then
    echo "WARNING: full pool ($n_circ_full_pool) is smaller than READS_PER_RUN ($READS_PER_RUN)." \
         "Every run below will fall short of $READS_PER_RUN circular reads (subsample-to caps at " \
         "whatever's available, it doesn't fail). Consider a larger mega-pool N_READS."
fi

echo "=== Step 4: Drawing $N_RUNS independent replicate runs of $READS_PER_RUN circular reads each ==="
for run in $(seq 1 "$N_RUNS"); do
    run_dir="$OUTDIR/run${run}"
    subsample_seed=$(( 1000 + run ))
    echo "--- Run $run: subsampling to $READS_PER_RUN reads (seed=$subsample_seed) ---"
    $APPTAINER_NANOCIRC nanocirc-sim pool-select \
        --bsj-reads "$bsj_reads_for_subsample" \
        --metadata "$DATASET_DIR/circRNA_metadata.tsv" \
        --mega-fastq "$mega_fastq_for_subsample" \
        --outdir "$run_dir" \
        --subsample-to "$READS_PER_RUN" \
        --seed "$subsample_seed" \
        --ground-truth "$DATASET_DIR/ground_truth.tsv"

    n_circ_run=$(( $(wc -l < "$run_dir/selected_reads.fastq") / 4 ))
    echo "--- Run $run: $n_circ_run circular reads selected (target $READS_PER_RUN) ---"

    if [ "${SKIP_CHAIN:-0}" = "1" ]; then
        echo "--- SKIP_CHAIN=1: not submitting run $run's linear-matching job ---"
    else
        linear_seed=$(( 300 + run ))
        echo "--- Submitting run_linear_matched.sbatch for run $run: exactly $n_circ_run linear reads (seed=$linear_seed) ---"
        sbatch --job-name="${SLURM_JOB_NAME}-run${run}-linear-matched" \
            --output="%x_%j.out" --error="%x_%j.err" \
            --export=ALL,\
DATASET_DIR="$DATASET_DIR",GENOME_FA="$GENOME_FA",TRANSCRIPTOME_FA="$TRANSCRIPTOME_FA",\
CONTROL_FASTQ="$CONTROL_FASTQ",CIRC_OUTDIR="$CIRC_OUTDIR",\
CONFIRMED_CIRC_FASTQ="$run_dir/selected_reads.fastq",N_LINEAR_READS="$n_circ_run",\
RUN_CIRCRNAS_BED="$run_dir/run_circrnas.bed",RUN_GROUND_TRUTH_TSV="$run_dir/run_ground_truth.tsv",\
NANOSIM_SEED="$linear_seed",OUTDIR="$FINAL_OUTDIR/run${run}" \
            "$SCRIPTS_DIR/run_linear_matched.sbatch"
    fi
done

echo "=== Done. $N_RUNS replicate run(s) submitted, $READS_PER_RUN circular + $READS_PER_RUN linear reads each, under $FINAL_OUTDIR/run<N>/ ==="
