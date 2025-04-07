#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr
#SBATCH --job-name=ground_truth_map
#SBATCH --output=ground_truth_map_%j.out
#SBATCH --cpus-per-task=8 --mem=64G

# Add resource tracking
#SBATCH --comment=resourcetracking

# Reference genome
ref_path="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/GRCm38.p4.genome_corrected.fa"

# Target directory to process
target_dir="/scratch/aerusakovich/sim_ciri_long_jobim/combined"
run_folder=$(basename "$target_dir")

# Main output directory
MAIN_OUTPUT_DIR="${target_dir}/pooled/ground_truth"

# Create metrics directory
metrics_dir="${target_dir}/performance_metrics"
mkdir -p "$metrics_dir"

# Record job start time and info
JOBID=$SLURM_JOB_ID
start_datetime=$(date "+%Y-%m-%d %H:%M:%S")
echo "Job started at: $start_datetime" > "${metrics_dir}/ground_truth_job_${JOBID}_metrics.txt"
echo "Tool: ground_truth_mapper" >> "${metrics_dir}/ground_truth_job_${JOBID}_metrics.txt"

# Create main output directory
mkdir -p "$MAIN_OUTPUT_DIR"

# Create Python script as a temporary file
TEMP_PYTHON_SCRIPT=$(mktemp)

cat > "${TEMP_PYTHON_SCRIPT}" << 'PYTHON_SCRIPT'
#!/usr/bin/env python3

import os
import sys
import subprocess
import re
import argparse
from collections import defaultdict

def check_dependencies():
    """Check if required tools are installed."""
    required_tools = ['minimap2', 'samtools']
    missing_tools = []
    
    for tool in required_tools:
        try:
            subprocess.run([tool, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"✓ {tool} is installed")
        except FileNotFoundError:
            missing_tools.append(tool)
            print(f"✗ {tool} is NOT installed")
    
    if missing_tools:
        print("\nMissing dependencies. Please install them with:")
        print("module load minimap2 samtools")
        print("# OR")
        print("conda install -c bioconda minimap2 samtools")
        return False
    
    return True

def map_reads(fastq_file, ref_genome, output_prefix, threads=4):
    """
    Map reads to reference genome using minimap2.
    
    Args:
        fastq_file: Input FASTQ file
        ref_genome: Reference genome FASTA file
        output_prefix: Prefix for output files
        threads: Number of threads to use
        
    Returns:
        Path to the sorted BAM file
    """
    print(f"Mapping reads from {fastq_file} to {ref_genome}...")
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    
    # Output files
    sam_file = f"{output_prefix}.sam"
    bam_file = f"{output_prefix}.bam"
    sorted_bam = f"{output_prefix}.sorted.bam"
    
    # Map reads with minimap2
    # Use splice-aware mapping (-ax splice) for RNA-seq
    cmd_map = [
        'minimap2', '-ax', 'splice', 
        '-t', str(threads),
        '-p', '1.0',  # Only primary alignments (p=1.0 effectively means --secondary=no)
        '-N', '0',    # Don't output secondary alignments
        ref_genome, fastq_file
    ]
    
    with open(sam_file, 'w') as f:
        subprocess.run(cmd_map, stdout=f, check=True)
    
    # Convert SAM to BAM
    cmd_convert = ['samtools', 'view', '-bS', '-@', str(threads), '-o', bam_file, sam_file]
    subprocess.run(cmd_convert, check=True)
    
    # Sort BAM file
    cmd_sort = ['samtools', 'sort', '-@', str(threads), '-o', sorted_bam, bam_file]
    subprocess.run(cmd_sort, check=True)
    
    # Index BAM file
    cmd_index = ['samtools', 'index', sorted_bam]
    subprocess.run(cmd_index, check=True)
    
    # Clean up intermediate files
    os.remove(sam_file)
    os.remove(bam_file)
    
    print(f"Mapping complete. Sorted BAM file: {sorted_bam}")
    return sorted_bam

def parse_header(header):
    """
    Parse the FASTQ header to extract relevant information.
    Expected format includes _circ or _linear suffix after modification
    """
    # Remove the '@' symbol at the beginning
    header = header.strip()
    if header.startswith('@'):
        header = header[1:]
    
    # Check if header contains _circ or _linear
    if '_circ' in header:
        read_type = 'circ'
    elif '_linear' in header:
        read_type = 'linear'
    else:
        read_type = 'unknown'
    
    return {
        'type': read_type
    }

def identify_blocks_from_cigar(cigar_string, ref_start):
    """
    Parse CIGAR string to identify exon blocks.
    
    Args:
        cigar_string: CIGAR string from SAM/BAM file
        ref_start: Start position on reference (1-based)
        
    Returns:
        List of blocks (start, end) in 1-based coordinates
    """
    blocks = []
    current_pos = ref_start
    current_block_start = ref_start
    
    # Parse CIGAR string
    cigar_pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    for length, operation in cigar_pattern.findall(cigar_string):
        length = int(length)
        
        if operation == 'M' or operation == '=' or operation == 'X':  # Match or mismatch
            # These operations consume both query and reference
            current_pos += length
        elif operation == 'I' or operation == 'S':  # Insertion or soft clipping
            # These operations consume query but not reference
            pass
        elif operation == 'D':  # Deletion
            # Deletion consumes reference but not query
            current_pos += length
        elif operation == 'N':  # Intron (splice junction)
            # End current block
            blocks.append((current_block_start, current_pos - 1))
            
            # Skip the intron
            current_pos += length
            current_block_start = current_pos
    
    # Add the last block if there's one
    if current_block_start <= current_pos - 1:
        blocks.append((current_block_start, current_pos - 1))
    
    return blocks

def convert_blocks_to_bed12(blocks, chrom_start):
    """
    Convert blocks to BED12 block sizes and starts format.
    
    Args:
        blocks: List of (start, end) tuples in 1-based coordinates
        chrom_start: Start position on chromosome (0-based)
        
    Returns:
        block_count, block_sizes, block_starts
    """
    if not blocks:
        return 0, "", ""
    
    # Convert blocks to 0-based
    blocks_0based = [(start - 1, end - 1) for start, end in blocks]
    
    block_count = len(blocks_0based)
    
    block_sizes = []
    block_starts = []
    
    for start, end in blocks_0based:
        block_size = end - start + 1
        block_start = start - chrom_start
        
        block_sizes.append(str(block_size))
        block_starts.append(str(block_start))
    
    return block_count, ",".join(block_sizes), ",".join(block_starts)

def convert_to_bed12(qname, read_info, chrom, pos, cigar, mapq, flag):
    """
    Convert alignment to BED12 format.
    
    Args:
        qname: Read name
        read_info: Parsed read header info
        chrom: Chromosome name
        pos: 1-based position
        cigar: CIGAR string
        mapq: Mapping quality
        flag: SAM flag
        
    Returns:
        BED12 line or None if unmapped
    """
    # Skip unmapped reads
    if flag & 0x4:  # SAM flag for unmapped read
        return None
    
    # Identify blocks from CIGAR string
    blocks = identify_blocks_from_cigar(cigar, pos)
    
    # Skip if no valid blocks
    if not blocks:
        return None
    
    # Extract BED12 fields
    chrom_start = pos - 1  # Convert to 0-based
    chrom_end = blocks[-1][1] - 1  # End of last block (0-based)
    
    # Create name with read information
    name = f"{qname}"
    
    # Use mapping quality as score (max 1000)
    score = min(mapq * 25, 1000)  # Scale mapq (0-40) to 0-1000
    
    # Get strand from read info or determine from SAM flag
    strand = '-' if flag & 0x10 else '+'  # 0x10 flag means reverse strand
    
    # For BED12, thickStart and thickEnd define the coding region
    # For RNA-seq, we can use the same as chromStart and chromEnd
    thick_start = chrom_start
    thick_end = chrom_end
    
    # Color based on read type
    if read_info and read_info['type'] == 'circ':
        item_rgb = "255,0,0"  # Red for circular
    else:
        item_rgb = "0,0,255"  # Blue for linear
    
    # Convert blocks to BED12 format
    block_count, block_sizes, block_starts = convert_blocks_to_bed12(blocks, chrom_start)
    
    # Create BED12 line
    bed_line = f"{chrom}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\t"
    bed_line += f"{thick_start}\t{thick_end}\t{item_rgb}\t{block_count}\t{block_sizes}\t{block_starts}"
    
    return bed_line

def process_bam_to_bed12(bam_file, output_bed, read_type=None):
    """
    Process BAM file and convert alignments to BED12 format.
    
    Args:
        bam_file: Path to BAM file
        output_bed: Path to output BED12 file
        read_type: Optional filter for read type ('circ' or 'linear')
    """
    print(f"Converting {bam_file} to BED12 format...")
    
    # Statistics
    stats = {
        'total_reads': 0,
        'mapped_reads': 0,
        'unmapped_reads': 0,
        'circ_reads': 0,
        'linear_reads': 0,
        'multi_exon_reads': 0
    }
    
    # Run samtools view to get alignment information
    cmd = ['samtools', 'view', bam_file]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, check=True, text=True)
    
    with open(output_bed, 'w') as f_out:
        for line in result.stdout.splitlines():
            fields = line.split('\t')
            
            # Basic SAM fields
            qname = fields[0]       # Read name
            flag = int(fields[1])   # FLAG
            rname = fields[2]       # Reference name
            pos = int(fields[3])    # 1-based leftmost position
            mapq = int(fields[4])   # Mapping quality
            cigar = fields[5]       # CIGAR string
            
            stats['total_reads'] += 1
            
            # Parse read header to get additional information
            read_info = parse_header(qname)
            
            # Filter by read type if specified
            if read_type and read_info and read_info['type'] != read_type:
                continue
            
            # Update read type statistics
            if read_info:
                if read_info['type'] == 'circ':
                    stats['circ_reads'] += 1
                elif read_info['type'] == 'linear':
                    stats['linear_reads'] += 1
            
            # Skip unmapped reads
            if flag & 0x4:  # SAM flag for unmapped read
                stats['unmapped_reads'] += 1
                continue
            
            stats['mapped_reads'] += 1
            
            # Convert to BED12
            bed_line = convert_to_bed12(qname, read_info, rname, pos, cigar, mapq, flag)
            
            if bed_line:
                # Count multi-exon reads
                block_count = int(bed_line.split('\t')[9])
                if block_count > 1:
                    stats['multi_exon_reads'] += 1
                
                f_out.write(bed_line + '\n')
    
    # Print statistics
    print(f"Conversion complete for {bam_file}")
    print(f"Total reads: {stats['total_reads']}")
    percentage = 0
    if stats['total_reads'] > 0:
        percentage = stats['mapped_reads']/stats['total_reads']*100
    print(f"Mapped reads: {stats['mapped_reads']} ({percentage:.1f}%)")
    print(f"Unmapped reads: {stats['unmapped_reads']}")
    print(f"Circular reads: {stats['circ_reads']}")
    print(f"Linear reads: {stats['linear_reads']}")
    
    multi_exon_percentage = 0
    if stats['mapped_reads'] > 0:
        multi_exon_percentage = stats['multi_exon_reads']/stats['mapped_reads']*100
    print(f"Multi-exon reads: {stats['multi_exon_reads']} ({multi_exon_percentage:.1f}% of mapped)")
    
    return stats

def main():
    # Add argument parser
    parser = argparse.ArgumentParser(description='Map and convert combined RNA reads to BED12 format.')
    parser.add_argument('--combined_file', required=True, help='Path to combined RNA FASTQ file')
    parser.add_argument('--ref_path', required=True, help='Path to reference genome FASTA file')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads to use')
    parser.add_argument('--run_name', default='run', help='Run name to use in output files')
    args = parser.parse_args()
    
    # Check if required tools are installed
    if not check_dependencies():
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Output files
    mapping_prefix = os.path.join(args.output_dir, f"{args.run_name}_mapped")
    output_bed = os.path.join(args.output_dir, f"{args.run_name}_all.bed")
    circ_output_bed = os.path.join(args.output_dir, f"{args.run_name}_circular.bed")
    linear_output_bed = os.path.join(args.output_dir, f"{args.run_name}_linear.bed")
    stats_output = os.path.join(args.output_dir, f"{args.run_name}_mapping_stats.txt")
    
    # Process all reads
    print(f"\n=== Processing Combined RNA reads for {args.run_name} ===")
    bam_file = map_reads(args.combined_file, args.ref_path, mapping_prefix, threads=args.threads)
    
    # Process combined file for all reads
    all_stats = process_bam_to_bed12(bam_file, output_bed)
    
    # Process for circular reads only
    circ_stats = process_bam_to_bed12(bam_file, circ_output_bed, read_type="circ")
    
    # Process for linear reads only
    linear_stats = process_bam_to_bed12(bam_file, linear_output_bed, read_type="linear")
    
    # Write statistics to file
    with open(stats_output, 'w') as f:
        f.write(f"MAPPING STATISTICS FOR {args.run_name}\n")
        f.write("=================\n\n")
        
        f.write("All Reads Stats:\n")
        for key, value in all_stats.items():
            f.write(f"  {key}: {value}\n")
        
        f.write("\nCircular RNA Stats (filtered):\n")
        for key, value in circ_stats.items():
            f.write(f"  {key}: {value}\n")
        
        f.write("\nLinear RNA Stats (filtered):\n")
        for key, value in linear_stats.items():
            f.write(f"  {key}: {value}\n")
    
    print(f"\nStatistics written to: {stats_output}")
    print(f"Processing complete for {args.run_name}.")

if __name__ == "__main__":
    main()
PYTHON_SCRIPT

# Make the temporary Python script executable
chmod +x "${TEMP_PYTHON_SCRIPT}"

# Start timing
start_time=$(date +%s)

# Define the combined reads file
combined_file="${target_dir}/pooled/combined_reads.fastq"

# Check if the combined file exists
if [ ! -f "$combined_file" ]; then
    echo "Error: Combined reads file not found: $combined_file"
    exit 1
fi

# Run the Python script on the combined reads file
echo "Running ground truth mapping for ${run_folder}..."
python3 "${TEMP_PYTHON_SCRIPT}" \
    --combined_file "$combined_file" \
    --ref_path "$ref_path" \
    --output_dir "$MAIN_OUTPUT_DIR" \
    --threads 8 \
    --run_name "${run_folder}"

echo "Finished processing ${run_folder}"

# Clean up temporary Python script
rm "${TEMP_PYTHON_SCRIPT}"

# End timing
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Ground truth mapping finished in $elapsed_time seconds"

# Record elapsed time in the metrics file
echo "Wall clock time: $elapsed_time seconds" >> "${metrics_dir}/ground_truth_job_${JOBID}_metrics.txt"

# After job completion, gather detailed SLURM metrics
end_datetime=$(date "+%Y-%m-%d %H:%M:%S")
echo "Job ended at: $end_datetime" >> "${metrics_dir}/ground_truth_job_${JOBID}_metrics.txt"

# Use sacct to get detailed job metrics and append to the metrics file
echo "Detailed SLURM metrics:" >> "${metrics_dir}/ground_truth_job_${JOBID}_metrics.txt"
sacct -j $JOBID --format=JobID,JobName,MaxRSS,MaxVMSize,NTasks,AllocCPUS,TotalCPU,CPUTime,Elapsed,Start,End -p >> "${metrics_dir}/ground_truth_job_${JOBID}_metrics.txt"

# Create a simple CSV file for easy plotting
echo "JobID,Tool,MaxRSS_KB,MaxCPU_Percent,ElapsedTime_Sec" > "${metrics_dir}/ground_truth_job_${JOBID}_summary.csv"
max_rss=$(sacct -j $JOBID --format=MaxRSS -p | grep -v "MaxRSS" | head -1 | tr -d "K|")
max_cpu_percent=$(echo "scale=2; $(sacct -j $JOBID --format=TotalCPU,Elapsed -p | grep -v "TotalCPU" | head -1 | awk -F"|" '{split($1,cpu,":");split($2,elapsed,":"); cpu_seconds=cpu[1]*3600+cpu[2]*60+cpu[3]; elapsed_seconds=elapsed[1]*3600+elapsed[2]*60+elapsed[3]; print cpu_seconds/elapsed_seconds*100}')" | bc)
elapsed_seconds=$(sacct -j $JOBID --format=Elapsed -p | grep -v "Elapsed" | head -1 | awk -F"|" '{split($1,t,":"); print t[1]*3600+t[2]*60+t[3]}')
echo "$JOBID,ground_truth_mapping,$max_rss,$max_cpu_percent,$elapsed_seconds" >> "${metrics_dir}/ground_truth_job_${JOBID}_summary.csv"

echo "Performance metrics saved to ${metrics_dir}/ground_truth_job_${JOBID}_metrics.txt"
echo "Summary data for plotting saved to ${metrics_dir}/ground_truth_job_${JOBID}_summary.csv"
echo "Processing completed"