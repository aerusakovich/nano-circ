#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr
#SBATCH --job-name=cirilong_combined_bed12
#SBATCH --output=cirilong_combined_bed12_%j.out
#SBATCH --cpus-per-task=4 --mem=32G

# Base path for combined folder
base_path="/scratch/aerusakovich/sim_ciri_long_jobim/combined"

# Main output directory
MAIN_OUTPUT_DIR="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/ciri-long"

# Create the output directory if it doesn't exist
mkdir -p "$MAIN_OUTPUT_DIR"

# Create Python script as a temporary file
TEMP_PYTHON_SCRIPT=$(mktemp)

cat > "${TEMP_PYTHON_SCRIPT}" << 'PYTHON_SCRIPT'
#!/usr/bin/env python3

import os
import re
import argparse
import sys
try:
    from Bio import SeqIO
except ImportError:
    print("BioPython not found. Installing...")
    import subprocess
    subprocess.check_call(["pip", "install", "biopython", "--user"])
    from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Convert CIRI-long output to BED12 format')
    parser.add_argument('--fasta', required=True, help='Path to CIRI-long multifasta file')
    parser.add_argument('--info', required=True, help='Path to CIRI-long info file')
    parser.add_argument('--output', required=True, help='Output BED12 file path')
    return parser.parse_args()

def parse_ciri_info(info_file):
    """Parse CIRI-long info file to extract circRNA information."""
    circ_info = {}
    info_count = 0
    multi_exon_count = 0
    
    with open(info_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            info_count += 1
            chrom = fields[0]
            start = int(fields[3])  # 1-based in GFF
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            
            # Extract attributes
            attr_dict = {}
            attrs = fields[8].split(';')
            for attr in attrs:
                attr = attr.strip()
                if not attr:
                    continue
                if '"' in attr:
                    key, value = attr.split(' "', 1)
                    value = value.strip('"')
                    attr_dict[key.strip()] = value
                else:
                    try:
                        key, value = attr.split(' ', 1)
                        attr_dict[key.strip()] = value.strip()
                    except ValueError:
                        # Skip attributes that don't have a space separator
                        continue
            
            circ_id = attr_dict.get('circ_id', '').strip('"')
            if not circ_id:
                # Generate a circID if not found
                circ_id = f"{chrom}:{start}-{end}"
            
            # Parse isoform information
            isoform = attr_dict.get('isoform', '').strip('"')
            exons = []
            
            if isoform:
                exon_coords = isoform.split(',')
                for coord in exon_coords:
                    if '-' in coord:
                        try:
                            # Handle normal exon coordinates (start-end)
                            parts = coord.split('-')
                            if len(parts) == 2:
                                exon_start = int(parts[0])
                                exon_end = int(parts[1])
                                exons.append((exon_start, exon_end))
                        except ValueError as e:
                            print(f"Warning: Could not parse exon coordinate {coord}: {e}")
                            continue
            
            # If no exon information, treat as single exon
            if not exons:
                exons = [(start, end)]
            
            if len(exons) > 1:
                multi_exon_count += 1
                
            circ_info[circ_id] = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'score': score,
                'strand': strand,
                'exons': sorted(exons, key=lambda x: x[0]),
                'circ_id': circ_id,
                'gene_name': attr_dict.get('gene_name', '').strip('"'),
                'isoform': isoform
            }
    
    print(f"Parsed {info_count} entries from CIRI-long info file, extracted {len(circ_info)} unique circRNAs")
    print(f"Found {multi_exon_count} multi-exon circRNAs in info file")
    return circ_info

def load_fasta_ids(fasta_file):
    """Load all sequence IDs from the FASTA file for validation."""
    fasta_ids = set()
    circ_id_pattern = re.compile(r'chr\w+:\d+-\d+')
    
    for record in SeqIO.parse(fasta_file, 'fasta'):
        header = record.description
        # Try to find a circRNA ID pattern in the header
        matches = circ_id_pattern.findall(header)
        if matches:
            for match in matches:
                fasta_ids.add(match)
    
    print(f"Found {len(fasta_ids)} circRNA IDs in FASTA file")
    return fasta_ids

def create_bed12(circ_info, fasta_ids, output_file):
    """Create BED12 file from circRNA information."""
    found_in_fasta = 0
    missing_in_fasta = 0
    
    with open(output_file, 'w') as out:
        for circ_id, info in circ_info.items():
            # Check if this circRNA is in the FASTA file (for reporting only)
            if circ_id in fasta_ids:
                found_in_fasta += 1
            else:
                missing_in_fasta += 1
            
            chrom = info['chrom']
            start = info['start'] - 1  # Convert to 0-based for BED
            end = info['end']
            name = circ_id
            score = info['score']
            strand = info['strand']
            
            # BED12 requires these fields
            thick_start = start  # CDS start (or start if no CDS)
            thick_end = end  # CDS end (or end if no CDS)
            item_rgb = "0,0,0"  # Default color
            
            exons = info['exons']
            block_count = len(exons)
            
            # Calculate block sizes and starts (relative to start)
            block_sizes = []
            block_starts = []
            
            for exon_start, exon_end in exons:
                block_size = exon_end - exon_start + 1
                block_start = exon_start - 1 - start  # 0-based relative to start
                
                block_sizes.append(str(block_size))
                block_starts.append(str(block_start))
            
            # Join with commas
            block_sizes_str = ','.join(block_sizes)
            block_starts_str = ','.join(block_starts)
            
            # Write BED12 line
            bed_line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{item_rgb}\t{block_count}\t{block_sizes_str}\t{block_starts_str}\n"
            out.write(bed_line)
    
    print(f"CircRNAs found in both info and FASTA: {found_in_fasta}")
    print(f"CircRNAs in info but not in FASTA: {missing_in_fasta}")
    return len(circ_info)

def main():
    args = parse_args()
    
    # Parse CIRI-long info file (primary source)
    circ_info = parse_ciri_info(args.info)
    
    # Load FASTA IDs for validation/reporting
    fasta_ids = load_fasta_ids(args.fasta)
    
    # Create BED12 file
    count = create_bed12(circ_info, fasta_ids, args.output)
    
    print(f"Successfully created BED12 file with {count} entries: {args.output}")

if __name__ == "__main__":
    main()
PYTHON_SCRIPT

# Make the temporary Python script executable
chmod +x "${TEMP_PYTHON_SCRIPT}"

# Function to find CIRI-long output files
find_ciri_long_files() {
    local base_dir="$1"
    
    # Find all CIRI-long info files in the combined directory
    info_files=$(find "$base_dir" -name "CIRI-long.info")
    
    for info_file in $info_files; do
        # Get the directory containing the info file
        dir_path=$(dirname "$info_file")
        
        # Look for the matching FASTA file
        fasta_file="${dir_path}/CIRI-long.cand_circ.fa"
        
        if [ -f "$fasta_file" ]; then
            # Extract a suitable output folder name from the path
            rel_path=${dir_path#$base_path/}
            # Replace slashes with underscores for output folder
            output_name=$(echo "$rel_path" | tr '/' '_')
            
            echo "Processing CIRI-long output in: $dir_path"
            echo "Output will be saved as: $output_name"
            
            # Create output directory
            output_dir="${MAIN_OUTPUT_DIR}/${output_name}"
            mkdir -p "$output_dir"
            output_file="${output_dir}/CIRI-long.cand_circ.bed12"
            
            # Run the Python script
            echo "Running INFO-centric CIRI-long to BED12 conversion..."
            python3 "${TEMP_PYTHON_SCRIPT}" \
              --fasta "$fasta_file" \
              --info "$info_file" \
              --output "$output_file"
            
            # Check if the conversion was successful
            if [ $? -eq 0 ]; then
                echo "Conversion successful! BED12 file saved to: $output_file"
            
                # Count lines in output file
                LINE_COUNT=$(wc -l < "$output_file")
                echo "Created BED12 file with $LINE_COUNT entries"
                
                # Show a few examples from the output file
                echo "First few entries in the BED12 file:"
                head -n 3 "$output_file"
                
                # Show some multi-exon examples if they exist
                echo "Examples of multi-exon circRNAs:"
                grep -m 3 -P "\t[2-9]\t" "$output_file" || echo "No multi-exon circRNAs found"
                
                # Check if bedtools is available for validation
                if command -v bedtools &> /dev/null; then
                    echo "Validating BED12 format (first 10 lines) with bedtools sort..."
                    head -n 10 "$output_file" | bedtools sort -i - > /dev/null 2>&1
                    if [ $? -eq 0 ]; then
                        echo "BED12 format is valid."
                    else
                        echo "Warning: BED12 format validation failed."
                    fi
                fi
            else
                echo "Error occurred during conversion"
            fi
            
            echo "-----------------------------------------"
        else
            echo "Found info file but no matching FASTA file at: $info_file"
        fi
    done
}

# Process the combined directory
find_ciri_long_files "$base_path"

# Clean up temporary Python script
rm "${TEMP_PYTHON_SCRIPT}"

echo "All CIRI-long outputs have been processed"