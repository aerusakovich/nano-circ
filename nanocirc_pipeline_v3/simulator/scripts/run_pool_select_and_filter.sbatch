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
# optional top-up for deficient circRNAs, then draw N_RUNS Poisson-sampled
# runs and self-chain each into run_linear_matched.sbatch. See README.
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
# Optional: PROJECT_ROOT, NANOSIM_SIF, NANOCIRC_SIM_SIF,
#   FINAL_OUTDIR (default $OUTDIR/../final_balanced, each run written to
#   $FINAL_OUTDIR/run<N>/), SKIP_CHAIN=1 to not auto-submit the linear jobs,
#   TOPUP_READS_PER_DEFICIENT (default 500, single pass),
#   N_RUNS (default 3), READS_PER_RUN (default 300000)

set -euo pipefail

: "${DATASET_DIR:?must set DATASET_DIR}"
: "${GENOME_FA:?must set GENOME_FA}"
: "${TRANSCRIPTOME_FA:?must set TRANSCRIPTOME_FA}"
: "${CONTROL_FASTQ:?must set CONTROL_FASTQ}"
: "${CIRC_OUTDIR:?must set CIRC_OUTDIR}"
: "${OUTDIR:?must set OUTDIR}"
FINAL_OUTDIR="${FINAL_OUTDIR:-$OUTDIR/../final_balanced}"
TOPUP_READS_PER_DEFICIENT="${TOPUP_READS_PER_DEFICIENT:-500}"
N_RUNS="${N_RUNS:-3}"
READS_PER_RUN="${READS_PER_RUN:-300000}"

PROJECT_ROOT="${PROJECT_ROOT:-/nfs/sto-2/projects/cnt-000070/canine-genetics/nanocirc}"
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

BIND="--bind /nfs/sto-2:/nfs/sto-2"
APPTAINER_NANOSIM="apptainer exec $BIND $NANOSIM_SIF"
APPTAINER_NANOCIRC="apptainer exec $BIND $NANOCIRC_SIM_SIF"

mega_pool="$CIRC_OUTDIR/pooled/combined_reads.fastq"
mega_bsj_reads="$CIRC_OUTDIR/bsj_check/bsj_check_reads.tsv"

echo "=== Step 1: Selecting up to n_reads_target confirmed reads per circRNA from the mega-pool ==="
$APPTAINER_NANOCIRC nanocirc-sim pool-select \
    --bsj-reads "$mega_bsj_reads" \
    --metadata "$DATASET_DIR/circRNA_metadata.tsv" \
    --mega-fastq "$mega_pool" \
    --outdir "$OUTDIR/selection"

n_deficient=$(($(wc -l < "$OUTDIR/selection/deficient_circle_ids.txt" 2>/dev/null || echo 0)))
final_selected="$OUTDIR/selection/selected_reads.fastq"
bsj_reads_for_subsample="$mega_bsj_reads"
mega_fastq_for_subsample="$mega_pool"

if [ "$n_deficient" -gt 0 ]; then
    echo "=== Step 2: $n_deficient circRNA(s) had zero confirmed reads in the mega-pool: single top-up pass ==="
    topup_dir="$OUTDIR/topup"
    mkdir -p "$topup_dir"

    # A NanoSim expression profile listing only the deficient circRNAs, all
    # at equal TPM=1.0. NanoSim drops any row with TPM<=0, so this alone
    # is enough to focus the whole read budget on exactly this set.
    {
        echo -e "ID\tcount\tTPM"
        awk '{print $0"\t0\t1.0"}' "$OUTDIR/selection/deficient_circle_ids.txt"
    } > "$topup_dir/topup_profile.tsv"

    n_topup_reads=$(( n_deficient * TOPUP_READS_PER_DEFICIENT ))
    $APPTAINER_NANOSIM simulator.py transcriptome \
        -rg "$GENOME_FA" \
        -rt "$DATASET_DIR/circRNAs.fa" \
        -c "$CIRC_OUTDIR/control" \
        -o "$topup_dir/topup" \
        -n "$n_topup_reads" \
        -b guppy \
        -e "$topup_dir/topup_profile.tsv" \
        --seed 200 \
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

    echo "=== Step 3: Re-selecting, now including the top-up pool ==="
    combined_mega="$topup_dir/mega_plus_topup.fastq"
    cat "$mega_pool" "$topup_pool" > "$combined_mega"
    combined_bsj_reads="$topup_dir/combined_bsj_reads.tsv"
    head -1 "$mega_bsj_reads" > "$combined_bsj_reads"
    tail -n +2 "$mega_bsj_reads" >> "$combined_bsj_reads"
    tail -n +2 "$topup_dir/bsj_check/bsj_check_reads.tsv" >> "$combined_bsj_reads"

    $APPTAINER_NANOCIRC nanocirc-sim pool-select \
        --bsj-reads "$combined_bsj_reads" \
        --metadata "$DATASET_DIR/circRNA_metadata.tsv" \
        --mega-fastq "$combined_mega" \
        --outdir "$OUTDIR/selection_final"

    final_selected="$OUTDIR/selection_final/selected_reads.fastq"
    bsj_reads_for_subsample="$combined_bsj_reads"
    mega_fastq_for_subsample="$combined_mega"
    still_deficient=$(wc -l < "$OUTDIR/selection_final/deficient_circle_ids.txt" 2>/dev/null || echo 0)
    echo "=== After top-up: $still_deficient circRNA(s) still have zero confirmed reads (reported, not hidden) ==="
else
    echo "=== No deficient circRNAs: every circRNA got its full n_reads_target already ==="
fi

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
