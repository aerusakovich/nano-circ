#!/bin/bash
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr # Where to send mail
#SBATCH --job-name=isocirc_bench # job name
##SBATCH --gres=gpu:1 -p gpu # if needed to specify a node (i.e GPU)
#SBATCH --output=isocirc_bench_%j.out # output file name with job ID
#SBATCH --cpus-per-task=8 --mem=120G # nb cpu and mem

# Add resource tracking
#SBATCH --comment=resourcetracking

# Activate env
source /local/env/envpython-3.9.5.sh
conda activate isocirc

# Paths
data_name="nanosim-ciri-long"
prog_name="isocirc"
ref_name="GRCm38.p4"
annotation_name="union"
db_path="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/circ_anno/union_bed.bed"
ref_path="/projects/dog/anastasia/data/Requirements/GRCm38.p4.genome.fa"
gene_anno_path="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/gencode.vM10.annotation.gtf"

# Parameters
THREADS=8

# Define the specific directory to process
target_dir="/scratch/aerusakovich/sim_ciri_long_jobim/combined"
run_folder=$(basename "$target_dir")
echo "Processing $run_folder"

# Create metrics directory
metrics_dir="${target_dir}/performance_metrics"
mkdir -p "$metrics_dir"

# Define input and output paths
reads_path="${target_dir}/pooled/combined_reads.fastq"
output_dir="${target_dir}/isocirc_output"

# Check if input file exists
if [ ! -f "$reads_path" ]; then
    echo "Error: Input file not found: $reads_path"
    exit 1
fi

# Create output directory
mkdir -p "$output_dir"

# Record job start time and info
JOBID=$SLURM_JOB_ID
start_datetime=$(date "+%Y-%m-%d %H:%M:%S")
echo "Job started at: $start_datetime" > "${metrics_dir}/isocirc_job_${JOBID}_metrics.txt"
echo "Tool: isocirc" >> "${metrics_dir}/isocirc_job_${JOBID}_metrics.txt"

# Run isoCirc
start_time=$(date +%s)
echo "Running isocirc on $data_name for $run_folder"
cd /projects/dog/anastasia/progs/isoCirc/
python /projects/dog/anastasia/progs/isoCirc/isoCirc/isoCirc_pipeline/isocirc/isocirc "$reads_path" "$ref_path" "$gene_anno_path" "$db_path" "$output_dir" -t $THREADS

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Isocirc finished in $elapsed_time seconds for $run_folder"

# Record elapsed time in the metrics file
echo "Wall clock time: $elapsed_time seconds" >> "${metrics_dir}/isocirc_job_${JOBID}_metrics.txt"

# After job completion, gather detailed SLURM metrics
end_datetime=$(date "+%Y-%m-%d %H:%M:%S")
echo "Job ended at: $end_datetime" >> "${metrics_dir}/isocirc_job_${JOBID}_metrics.txt"

# Use sacct to get detailed job metrics and append to the metrics file
echo "Detailed SLURM metrics:" >> "${metrics_dir}/isocirc_job_${JOBID}_metrics.txt"
sacct -j $JOBID --format=JobID,JobName,MaxRSS,MaxVMSize,NTasks,AllocCPUS,TotalCPU,CPUTime,Elapsed,Start,End -p >> "${metrics_dir}/isocirc_job_${JOBID}_metrics.txt"

# Create a simple CSV file for easy plotting
echo "JobID,Tool,MaxRSS_KB,MaxCPU_Percent,ElapsedTime_Sec" > "${metrics_dir}/isocirc_job_${JOBID}_summary.csv"
max_rss=$(sacct -j $JOBID --format=MaxRSS -p | grep -v "MaxRSS" | head -1 | tr -d "K|")
max_cpu_percent=$(echo "scale=2; $(sacct -j $JOBID --format=TotalCPU,Elapsed -p | grep -v "TotalCPU" | head -1 | awk -F"|" '{split($1,cpu,":");split($2,elapsed,":"); cpu_seconds=cpu[1]*3600+cpu[2]*60+cpu[3]; elapsed_seconds=elapsed[1]*3600+elapsed[2]*60+elapsed[3]; print cpu_seconds/elapsed_seconds*100}')" | bc)
elapsed_seconds=$(sacct -j $JOBID --format=Elapsed -p | grep -v "Elapsed" | head -1 | awk -F"|" '{split($1,t,":"); print t[1]*3600+t[2]*60+t[3]}')
echo "$JOBID,isocirc,$max_rss,$max_cpu_percent,$elapsed_seconds" >> "${metrics_dir}/isocirc_job_${JOBID}_summary.csv"

echo "Performance metrics saved to ${metrics_dir}/isocirc_job_${JOBID}_metrics.txt"
echo "Summary data for plotting saved to ${metrics_dir}/isocirc_job_${JOBID}_summary.csv"
echo "Processing completed"