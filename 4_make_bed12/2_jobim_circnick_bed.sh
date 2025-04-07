#!/bin/bash
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr # Where to send mail
#SBATCH --job-name=circnick_extract # job name
#SBATCH --output=circnick_extract_%j.out # output file name with job ID
#SBATCH --cpus-per-task=4 --mem=32G # nb cpu and mem

# Base path for input files (with trailing slash)
base_path="/scratch/aerusakovich/sim_ciri_long_jobim/combined/circnick_output/combined_reads/"

# Reference path
REF_PATH="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/GRCm38.p4.genome_corrected.fa"

# Main output directory
MAIN_OUTPUT_DIR="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/circnick"

conda activate bed12

# Create main output directory
mkdir -p "$MAIN_OUTPUT_DIR"

# Get list of candidate files
candidate_files=$(find "$base_path" -name "*circRNA_candidates.annotated.txt")

for CANDIDATES_FILE in $candidate_files; do
    # Extract run name from the path
    run_folder=$(echo "$CANDIDATES_FILE" | sed -E 's|.*/([^/]+)/circnick_output/.*|\1|')
    echo "Processing $run_folder"
    
    # Determine paths for other required files based on the candidate file
    base_dir=$(dirname "$CANDIDATES_FILE")
    ALL_EXONS_FILE=$(find "$base_dir" -name "*allExons.bed" | grep -v "notGencode")
    NOVEL_EXONS_FILE=$(find "$base_dir" -name "*allExons.notGencode.bed")
    
    # Check if required files exist
    if [ ! -f "$ALL_EXONS_FILE" ]; then
        echo "Skipping $run_folder: All exons file not found"
        continue
    fi
    
    # Create run-specific output directory
    OUTPUT_DIR="${MAIN_OUTPUT_DIR}/${run_folder}"
    mkdir -p "$OUTPUT_DIR"
    
    # Output files for this run
    DEBUG_LOG="$OUTPUT_DIR/debug.log"
    BED12_OUTPUT="$OUTPUT_DIR/comprehensive_circrna.bed12"
    FASTA_OUTPUT="$OUTPUT_DIR/comprehensive_circrna.fa"
    
    # Clear previous debug log
    > "$DEBUG_LOG"
    
    # Function to log debug information
    log_debug() {
        echo "$1" >> "$DEBUG_LOG"
        echo "$1"
    }
    
    # Combine all exon files for processing
    COMBINED_EXONS="$OUTPUT_DIR/all_combined_exons.bed"
    
    # Check if novel exons file exists, otherwise just use all exons file
    if [ -f "$NOVEL_EXONS_FILE" ]; then
        cat "$ALL_EXONS_FILE" "$NOVEL_EXONS_FILE" | sort -k1,1 -k2,2n > "$COMBINED_EXONS"
        log_debug "Combined exon files from both regular and novel sources"
    else
        cat "$ALL_EXONS_FILE" | sort -k1,1 -k2,2n > "$COMBINED_EXONS"
        log_debug "Using only regular exons (novel exons file not found)"
    fi
    
    # Python script for comprehensive extraction with improved exon handling
    log_debug "Running Python script to process circRNAs..."
    python3 << END > "$BED12_OUTPUT"
import sys

def parse_exons(exons_file):
    """Parse all exons into a chromosome-based dictionary"""
    exons_by_chrom = {}
    with open(exons_file, 'r') as f:
        for line in f:
            try:
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                strand = parts[5] if len(parts) > 5 else '+'
                
                if chrom not in exons_by_chrom:
                    exons_by_chrom[chrom] = []
                exons_by_chrom[chrom].append((start, end, strand))
            except (ValueError, IndexError):
                print(f"Skipping invalid exon: {line.strip()}", file=sys.stderr)
    
    # Sort exons for each chromosome
    for chrom in exons_by_chrom:
        exons_by_chrom[chrom].sort(key=lambda x: x[0])
    
    return exons_by_chrom

def find_matching_exons(chrom, start, end, strand, exons_by_chrom):
    """Find exons within the circRNA coordinates with matching strand"""
    matching_exons = []
    
    # If chromosome not in exons, return empty list
    if chrom not in exons_by_chrom:
        return matching_exons
    
    # Find exons within the circRNA coordinates with matching strand
    for ex_start, ex_end, ex_strand in exons_by_chrom[chrom]:
        # Strand must match for proper orientation
        if ex_strand != strand:
            continue
            
        # Exon fully within circRNA or overlapping
        if (ex_start >= start and ex_end <= end) or \
           (ex_start < start and ex_end > start) or \
           (ex_start < end and ex_end > end):
            # Handle partial overlaps by trimming to circRNA boundaries
            effective_start = max(start, ex_start)
            effective_end = min(end, ex_end)
            matching_exons.append((effective_start, effective_end))
    
    # If no exons found, create a single-block representation
    if not matching_exons:
        matching_exons.append((start, end))
    
    # Sort exons and merge overlapping ones
    matching_exons.sort(key=lambda x: x[0])
    merged_exons = []
    
    if matching_exons:
        current = matching_exons[0]
        for next_ex in matching_exons[1:]:
            if current[1] >= next_ex[0]:  # Overlapping
                current = (current[0], max(current[1], next_ex[1]))
            else:
                merged_exons.append(current)
                current = next_ex
        merged_exons.append(current)
    
    return merged_exons

def process_circrnas():
    # Parse all exons
    exons_by_chrom = parse_exons('$COMBINED_EXONS')
    
    # Read candidates and process
    with open('$CANDIDATES_FILE', 'r') as f:
        # Skip header
        next(f)
        
        line_num = 1
        for line in f:
            try:
                parts = line.strip().split('\t')
                if len(parts) < 7:  # Need at least 7 columns for required data
                    print(f"Line {line_num}: Not enough columns ({len(parts)}): {line.strip()}", file=sys.stderr)
                    line_num += 1
                    continue
                    
                circ_id = parts[0]
                chrom = parts[1]
                start = int(parts[2])
                end = int(parts[3])
                strand = parts[6]
                
                # Validate data
                if start >= end:
                    print(f"Line {line_num}: Invalid coordinates - start ({start}) >= end ({end})", file=sys.stderr)
                    line_num += 1
                    continue
                    
                if strand not in ['+', '-']:
                    print(f"Line {line_num}: Invalid strand: {strand}", file=sys.stderr)
                    line_num += 1
                    continue
                    
            except (ValueError, IndexError) as e:
                print(f"Line {line_num}: Error parsing: {e} - {line.strip()}", file=sys.stderr)
                line_num += 1
                continue
                
            line_num += 1
            
            # Find matching exons
            matching_exons = find_matching_exons(chrom, start, end, strand, exons_by_chrom)
            
            # Prepare exon sizes and starts
            exon_sizes = []
            exon_starts = []
            
            for ex_start, ex_end in matching_exons:
                # Calculate size and relative start
                exon_sizes.append(str(ex_end - ex_start))
                exon_starts.append(str(ex_start - start))
            
            # Validate exon data before creating BED12 entry
            if not matching_exons:
                print(f"Warning: No matching exons found for {circ_id}", file=sys.stderr)
                continue
                
            if not exon_sizes or not exon_starts:
                print(f"Warning: Empty exon sizes or starts for {circ_id}", file=sys.stderr)
                continue
                
            # Ensure all fields are properly formatted to avoid BED12 format issues
            # BED12 format: chrom, start, end, name, score, strand, thick_start, thick_end, rgb, block_count, block_sizes, block_starts
            bed_entry = [
                chrom, 
                str(start), 
                str(end), 
                circ_id, 
                '1000', 
                strand, 
                str(start), 
                str(end), 
                '0,0,0', 
                str(len(matching_exons)), 
                ','.join(exon_sizes) + ',', # Add trailing comma for BED12 format
                ','.join(exon_starts) + ',' # Add trailing comma for BED12 format
            ]
            
            # Final validation check
            if len(bed_entry) != 12:
                print(f"Error: Invalid BED12 entry (fields={len(bed_entry)}) for {circ_id}", file=sys.stderr)
                continue
                
            print('\t'.join(bed_entry))

# Catch and log any top-level exceptions
try:
    process_circrnas()
except Exception as e:
    print(f"Error processing circRNAs: {e}", file=sys.stderr)
END
    
    # Check BED12 file for format issues before running bedtools
    log_debug "Checking BED12 file for format issues..."
    line_number=1
    while IFS= read -r line; do
        field_count=$(echo "$line" | tr '\t' '\n' | wc -l)
        if [ "$field_count" -ne 12 ]; then
            log_debug "ERROR: Line $line_number has $field_count fields instead of 12:"
            log_debug "$line"
            # Output the problematic lines to a debug file
            echo "Line $line_number: $line" > "$OUTPUT_DIR/problematic_bed_entries.txt"
            head -n $((line_number + 5)) "$BED12_OUTPUT" | tail -n 10 >> "$OUTPUT_DIR/problematic_bed_entries.txt"
            break
        fi
        line_number=$((line_number + 1))
        # Print progress every 1000 lines
        if [ $((line_number % 1000)) -eq 0 ]; then
            log_debug "Checked $line_number lines..."
        fi
    done < "$BED12_OUTPUT"
    
    # Extract sequences using bedtools with the important -split option to respect exon structure
    log_debug "Extracting sequences with bedtools getfasta using -split option for proper exon handling..."
    if ! bedtools getfasta -split -name -s -fi "$REF_PATH" -bed "$BED12_OUTPUT" -fo "$FASTA_OUTPUT"; then
        log_debug "Error extracting FASTA sequences"
        log_debug "Inspect $OUTPUT_DIR/problematic_bed_entries.txt for details on problematic entries"
        # If error occurred, try to find problematic lines with awk
        awk '{if(NF != 12) print "Line " NR " has " NF " fields: " $0}' "$BED12_OUTPUT" > "$OUTPUT_DIR/all_problematic_entries.txt"
        log_debug "All problematic entries saved to $OUTPUT_DIR/all_problematic_entries.txt"
        continue
    fi
    
    # Log details
    log_debug "Total candidate circRNAs:"
    tail -n +2 "$CANDIDATES_FILE" | wc -l >> "$DEBUG_LOG"
    
    log_debug "Total circRNAs processed:"
    wc -l "$BED12_OUTPUT" >> "$DEBUG_LOG"
    
    log_debug "Unique chromosomes:"
    cut -f1 "$BED12_OUTPUT" | sort | uniq -c >> "$DEBUG_LOG"
    
    log_debug "FASTA file details:"
    wc -l "$FASTA_OUTPUT" >> "$DEBUG_LOG"
    
    # Verify sequences are properly extracted by checking for spliced sequences
    log_debug "Checking for proper spliced sequence extraction:"
    grep "^>" "$FASTA_OUTPUT" | head -n 5 >> "$DEBUG_LOG"
    
    # Additional FASTA verification
    log_debug "First few FASTA sequences:"
    head -n 10 "$FASTA_OUTPUT" >> "$DEBUG_LOG"
    
    # Display processing details
    echo "Processing complete for $run_folder:"
    echo "BED12 output: $BED12_OUTPUT"
    echo "FASTA output: $FASTA_OUTPUT"
    echo "Debug log: $DEBUG_LOG"
    
    # Show first few lines of debug log
    head -n 20 "$DEBUG_LOG"
done

echo "All runs have been processed"