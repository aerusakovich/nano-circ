#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr
#SBATCH --job-name=pool_reads
#SBATCH --output=pool_reads_%j.out
#SBATCH --cpus-per-task=4 --mem=32G

# Define the specific directory to process
target_dir="/scratch/aerusakovich/sim_ciri_long_jobim/combined"
echo "Processing $(basename "$target_dir")"

# Create the pooled directory if it doesn't exist
pooled_dir="${target_dir}/pooled"
mkdir -p "$pooled_dir"

# Define the output file
output_file="${pooled_dir}/combined_reads.fastq"

# Skip if the output file already exists
if [ -f "$output_file" ]; then
  echo "Output file already exists, skipping"
  exit 0
fi

# Define the input files
circ_aligned="${target_dir}/circ_aligned_reads.fastq"
circ_unaligned="${target_dir}/circ_unaligned_reads.fastq"
linear_aligned="${target_dir}/linear_aligned_reads.fastq"
linear_unaligned="${target_dir}/linear_unaligned_reads.fastq"

# Check if all input files exist
missing_files=false
for file in "$circ_aligned" "$circ_unaligned" "$linear_aligned" "$linear_unaligned"; do
  if [ ! -f "$file" ]; then
    echo "Missing file: $file"
    missing_files=true
  fi
done

if [ "$missing_files" = true ]; then
  echo "Cannot proceed due to missing files"
  exit 1
fi

echo "Pooling files for $(basename "$target_dir")"

# Process each file and append to the output
# For circ_aligned reads
awk '{
  if (NR % 4 == 1) {
    # Modify header line by adding _circ suffix
    print substr($0, 1, 1) substr($0, 2) "_circ";
  } else {
    # Keep other lines unchanged
    print;
  }
}' "$circ_aligned" > "$output_file"

# For circ_unaligned reads
awk '{
  if (NR % 4 == 1) {
    # Modify header line by adding _circ suffix
    print substr($0, 1, 1) substr($0, 2) "_circ";
  } else {
    # Keep other lines unchanged
    print;
  }
}' "$circ_unaligned" >> "$output_file"

# For linear_aligned reads
awk '{
  if (NR % 4 == 1) {
    # Modify header line by adding _linear suffix
    print substr($0, 1, 1) substr($0, 2) "_linear";
  } else {
    # Keep other lines unchanged
    print;
  }
}' "$linear_aligned" >> "$output_file"

# For linear_unaligned reads
awk '{
  if (NR % 4 == 1) {
    # Modify header line by adding _linear suffix
    print substr($0, 1, 1) substr($0, 2) "_linear";
  } else {
    # Keep other lines unchanged
    print;
  }
}' "$linear_unaligned" >> "$output_file"

echo "Created pooled file: $output_file"

# Create a compressed version
echo "Compressing to ${output_file%.fastq}.fq.gz"
gzip -c "$output_file" > "${output_file%.fastq}.fq.gz"

echo "Completed processing $(basename "$target_dir")"