#!/bin/bash
#SBATCH --job-name=circ_sim_10x
#SBATCH --output=circ_sim_10x_%j.log
#SBATCH --error=circ_sim_10x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=3T
#SBATCH -p bigmem


# Print job info
echo "Job started at $(date)"
echo "Running on node: $(hostname)"
echo "Allocated memory: 3TB"
echo "CPU cores available: $SLURM_CPUS_PER_TASK"

# Load necessary modules (adjust as needed for your environment)
source /local/env/envpython-3.7.1.sh
conda activate bed12

# Base paths
BASE_DIR="/scratch/aerusakovich/sim_ciri_long_jobim"

# Input files for simulation
GTF_FILE='/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/gencode.vM10.sorted.annotation.gtf'
GENOME_FILE='/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/GRCm38.p4.genome_corrected.fa'
CONTROL_FQ="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/rep1/CRR194180.fq"
TRANSCRIPTOME='/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/ref_transcriptome/gencode.vM10.transcripts.fa'

# Path to the circRNA simulator
CIRCRNA_SIMULATOR='/scratch/aerusakovich/~/new_git/simulation_circ_fa/circ_fa_generation_exoncount.py'

# Path to bedtools executable
BEDTOOLS_PATH='/home/genouest/cnrs_umr6290/aerusakovich/.conda/envs/bed12/bin/bedtools'
NANOSIM_VERSION='/cvmfs/singularity.galaxyproject.org/all/nanosim:3.1.0--hdfd78af_0'

# Create base output directory if it doesn't exist
mkdir -p $BASE_DIR

# Number of circRNAs to simulate in each run
NUM_CIRCRNAS=10000

# Number of parallel threads to use (adjust based on your resources)
# Using all available CPUs as specified in SLURM allocation
MAX_THREADS=$SLURM_CPUS_PER_TASK

# Function to run simulation
run_simulation() {
    local run_id=$1
    local output_dir="$BASE_DIR/run_${run_id}"
    
    echo "Starting simulation run ${run_id} at $(date)"
    
    # Create output directory
    mkdir -p $output_dir
    
    # Run circRNA simulation
    python $CIRCRNA_SIMULATOR \
        --gtf $GTF_FILE \
        --genome $GENOME_FILE \
        --output $output_dir \
        --count $NUM_CIRCRNAS \
        --bedtools $BEDTOOLS_PATH
    
    # Check if simulation was successful
    if [ $? -eq 0 ]; then
        echo "Simulation run ${run_id} completed successfully at $(date)"
    else
        echo "ERROR: Simulation run ${run_id} failed at $(date)"
    fi
}

# Run simulations in parallel using background jobs
echo "Starting $MAX_THREADS parallel simulation threads..."

# Track running jobs
declare -a pids

# Launch simulations in parallel
for i in {1..10}; do
    # If we're at max threads, wait for one to finish
    while [ ${#pids[@]} -ge $MAX_THREADS ]; do
        # Check each process to see if it's still running
        for j in "${!pids[@]}"; do
            if ! kill -0 ${pids[$j]} 2>/dev/null; then
                # Process has finished, remove it from the array
                unset pids[$j]
            fi
        done
        # Reindex array to remove gaps
        pids=("${pids[@]}")
        # Sleep briefly before checking again
        sleep 1
    done
    
    # Start a new simulation in the background
    run_simulation $i &
    
    # Store the process ID
    pids+=($!)
    
    echo "Launched simulation run $i (PID: ${pids[-1]})"
done

# Wait for all remaining background jobs to finish
echo "Waiting for all simulations to complete..."
wait

# Create a summary report of all runs
echo "Creating summary report..."
echo "Simulation Summary Report" > $BASE_DIR/summary_report.txt
echo "Generated on: $(date)" >> $BASE_DIR/summary_report.txt
echo "----------------------------------------" >> $BASE_DIR/summary_report.txt

for i in {1..10}; do
    output_dir="$BASE_DIR/run_${i}"
    meta_file="$output_dir/circRNA_metadata.tsv"
    
    if [ -f "$meta_file" ]; then
        num_circs=$(wc -l < "$meta_file")
        num_circs=$((num_circs - 1)) # Subtract header line
        echo "Run $i: Generated $num_circs circRNAs" >> $BASE_DIR/summary_report.txt
        
        # Count circRNA types
        echo " circRNA Types:" >> $BASE_DIR/summary_report.txt
        for type in eciRNA EIciRNA ciRNA intergenic; do
            type_count=$(grep -c "\t$type\t" "$meta_file")
            echo " - $type: $type_count" >> $BASE_DIR/summary_report.txt
        done
    else
        echo "Run $i: Metadata file not found" >> $BASE_DIR/summary_report.txt
    fi
    
    echo "----------------------------------------" >> $BASE_DIR/summary_report.txt
done

echo "All simulations completed at $(date)"
echo "See summary report at $BASE_DIR/summary_report.txt"

# Print resource usage statistics
echo "Resource usage statistics:"
echo "-------------------------"
sstat --format=MaxRSS,MaxVMSize,AveCPU,AveRSS,AveVMSize $SLURM_JOB_ID