#!/bin/bash
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr # Where to send mail
#SBATCH --job-name=circ-nick_bench # job name
##SBATCH --gres=gpu:1 -p gpu # if needed to specify a node (i.e GPU)
#SBATCH --output=circ-nick_bench_%j.out # output file name with job ID
#SBATCH --cpus-per-task=8 --mem=100G # nb cpu and mem

# Add resource tracking
#SBATCH --comment=resourcetracking

# Activate env
source /local/env/envpython-3.9.5.sh
conda activate long_read_circRNA

# Paths
ref_path="/projects/dog/anastasia/data/data_circNICK"
gene_anno_path="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/gencode.vM10.sorted.annotation.gtf.gz"
scripts_path="/projects/dog/anastasia/progs/circNICK/long_read_circRNA/scripts"
data_name="nanosim-ciri-long"
prog_name="circnick"
ref_name="GRCm38.p4"
annotation_name="union"

# Parameters
THREADS=8
species="mouse"

# Define the specific directory to process
target_dir="/scratch/aerusakovich/sim_ciri_long_jobim/combined"
run_folder=$(basename "$target_dir")
echo "Processing $run_folder"

# Create metrics directory
metrics_dir="${target_dir}/performance_metrics"
mkdir -p "$metrics_dir"

# Define input path
fastq_path="${target_dir}/pooled/combined_reads.fastq"

# Check if fastq file exists
if [ ! -f "$fastq_path" ]; then
    echo "Error: Input file not found: $fastq_path"
    exit 1
fi

# Create output directory
output_dir="${target_dir}/circnick_output"
mkdir -p "$output_dir"

# Create the compressed version if needed
compressed_path="${fastq_path%.fastq}.fq.gz"
echo "Checking for compressed file at $compressed_path"

# Check if the compressed file already exists
if [ -f "$compressed_path" ]; then
    echo "Compressed file already exists, skipping compression"
else
    echo "Compressing $fastq_path to $compressed_path"
    gzip -c "$fastq_path" > "$compressed_path"
    echo "Compression complete"
fi

# Record job start time and info
JOBID=$SLURM_JOB_ID
start_datetime=$(date "+%Y-%m-%d %H:%M:%S")
echo "Job started at: $start_datetime" > "${metrics_dir}/circnick_job_${JOBID}_metrics.txt"
echo "Tool: circnick" >> "${metrics_dir}/circnick_job_${JOBID}_metrics.txt"

# Run circNICK
echo "Running circ-nick on $data_name for $run_folder"
start_time=$(date +%s)

/projects/dog/anastasia/progs/circNICK/long_read_circRNA/long_read_circRNA run "$compressed_path" --species $species --reference-path "$ref_path" --script-path "$scripts_path" --output-path "$output_dir"

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "circNICK-lrs finished in $elapsed_time seconds for $run_folder"

# Record elapsed time in the metrics file
echo "Wall clock time: $elapsed_time seconds" >> "${metrics_dir}/circnick_job_${JOBID}_metrics.txt"

# After job completion, gather detailed SLURM metrics
end_datetime=$(date "+%Y-%m-%d %H:%M:%S")
echo "Job ended at: $end_datetime" >> "${metrics_dir}/circnick_job_${JOBID}_metrics.txt"

# Use sacct to get detailed job metrics and append to the metrics file
echo "Detailed SLURM metrics:" >> "${metrics_dir}/circnick_job_${JOBID}_metrics.txt"
sacct -j $JOBID --format=JobID,JobName,MaxRSS,MaxVMSize,NTasks,AllocCPUS,TotalCPU,CPUTime,Elapsed,Start,End -p >> "${metrics_dir}/circnick_job_${JOBID}_metrics.txt"

# Create a simple CSV file for easy plotting
echo "JobID,Tool,MaxRSS_KB,MaxCPU_Percent,ElapsedTime_Sec" > "${metrics_dir}/circnick_job_${JOBID}_summary.csv"
max_rss=$(sacct -j $JOBID --format=MaxRSS -p | grep -v "MaxRSS" | head -1 | tr -d "K|")
max_cpu_percent=$(echo "scale=2; $(sacct -j $JOBID --format=TotalCPU,Elapsed -p | grep -v "TotalCPU" | head -1 | awk -F"|" '{split($1,cpu,":");split($2,elapsed,":"); cpu_seconds=cpu[1]*3600+cpu[2]*60+cpu[3]; elapsed_seconds=elapsed[1]*3600+elapsed[2]*60+elapsed[3]; print cpu_seconds/elapsed_seconds*100}')" | bc)
elapsed_seconds=$(sacct -j $JOBID --format=Elapsed -p | grep -v "Elapsed" | head -1 | awk -F"|" '{split($1,t,":"); print t[1]*3600+t[2]*60+t[3]}')
echo "$JOBID,circnick,$max_rss,$max_cpu_percent,$elapsed_seconds" >> "${metrics_dir}/circnick_job_${JOBID}_summary.csv"

echo "Performance metrics saved to ${metrics_dir}/circnick_job_${JOBID}_metrics.txt"
echo "Summary data for plotting saved to ${metrics_dir}/circnick_job_${JOBID}_summary.csv"
echo "Processing completed"