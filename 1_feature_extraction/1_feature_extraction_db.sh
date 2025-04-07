#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr
#SBATCH --job-name=circRNA_pipeline
#SBATCH --output=circRNA_pipeline_%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G

# ===========================================================================
# Complete circRNA Analysis Pipeline Script
# This script performs the full circRNA analysis pipeline:
# 1. Processes databases and creates intersections
# 2. Analyzes circRNAs for splice sites and classification
# 3. Produces comprehensive CSV files with all metrics
# ===========================================================================

# --------------------------------
# Configuration
# --------------------------------

# Load required environment
source /local/env/envpython-3.7.1.sh
conda activate bed12

# Define paths
GENOME_REF="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/GRCm38.p4.genome_corrected.fa"
GTF_FILE="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/gencode.vM10.sorted.annotation.gtf"
BASE_DIR="/scratch/aerusakovich/database_study"
DATE=$(date +%Y%m%d)
RESULTS_DIR="$BASE_DIR/results/$DATE"
TEMP_DIR="$RESULTS_DIR/temp"
LOG_DIR="$BASE_DIR/logs"

# Create directory structure
mkdir -p "$RESULTS_DIR/circatlas" "$RESULTS_DIR/circbase" "$RESULTS_DIR/intersection" "$TEMP_DIR" "$LOG_DIR"

# Database paths
CIRCATLAS_BED="$BASE_DIR/database/$DATE/circatlas/circrna_mapped_bed12.bed"
CIRCBASE_BED="$BASE_DIR/database/$DATE/circbase/circrna_mapped_bed12.bed"
INTERSECTION_BED="$RESULTS_DIR/intersection/intersected_bed_files.bed"

# Command line arguments
PROCESS_ALL=true
PROCESS_CIRCATLAS=false
PROCESS_CIRCBASE=false
PROCESS_INTERSECTION=false
MAX_ENTRIES=0
TEST_MODE=false

# Parse command line options
while [[ $# -gt 0 ]]; do
    case $1 in
        --test)
            TEST_MODE=true
            MAX_ENTRIES=1000
            shift
            ;;
        --max)
            MAX_ENTRIES=$2
            shift 2
            ;;
        --circatlas-only)
            PROCESS_ALL=false
            PROCESS_CIRCATLAS=true
            shift
            ;;
        --circbase-only)
            PROCESS_ALL=false
            PROCESS_CIRCBASE=true
            shift
            ;;
        --intersection-only)
            PROCESS_ALL=false
            PROCESS_INTERSECTION=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--test] [--max NUMBER] [--circatlas-only] [--circbase-only] [--intersection-only]"
            exit 1
            ;;
    esac
done

if [ "$PROCESS_ALL" = true ]; then
    PROCESS_CIRCATLAS=true
    PROCESS_CIRCBASE=true
    PROCESS_INTERSECTION=true
fi

# --------------------------------
# Utility Functions
# --------------------------------

# Log function
log() {
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $1" >> "$LOG_DIR/pipeline_$DATE.log"
}

# Error handling function
error_exit() {
    log "ERROR: $1"
    exit 1
}

# Find bedtools executable
find_bedtools() {
    if command -v bedtools &> /dev/null; then
        echo "bedtools"
        return
    fi
    
    # Try to load module if needed
    if command -v module &> /dev/null; then
        module load bioinfo/bedtools-2.27.1 2>/dev/null || true
        if command -v bedtools &> /dev/null; then
            echo "bedtools"
            return
        fi
    fi
    
    # Try common locations
    for path in "/usr/bin/bedtools" "/usr/local/bin/bedtools" "/home/genouest/cnrs_umr6290/aerusakovich/.conda/envs/bed12/bin/bedtools"; do
        if [ -f "$path" ]; then
            echo "$path"
            return
        fi
    done
    
    error_exit "Could not find bedtools executable. Please install it or load the appropriate module."
}

# --------------------------------
# Python Scripts Generation
# --------------------------------

# Create the main analysis script
log "Creating analysis script..."
cat > "$TEMP_DIR/analyze_circrna.py" << 'ENDPYTHON'
#!/usr/bin/env python3
import argparse
import os
import logging
import pandas as pd
import numpy as np
import pysam
import tempfile
import shutil
import subprocess
from collections import defaultdict

def configure_logging(output_dir):
    """Set up logging configuration"""
    os.makedirs(output_dir, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(os.path.join(output_dir, 'circrna_analysis.log')),
            logging.StreamHandler()
        ]
    )

def parse_bed12_entries(bed_file, max_entries=0):
    """Parse BED12 entries efficiently"""
    entries = []
    with open(bed_file, 'r') as f:
        for i, line in enumerate(f):
            if max_entries > 0 and i >= max_entries:
                break
                
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
                
            try:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                score = fields[4]
                strand = fields[5]
                thick_start = int(fields[6])
                thick_end = int(fields[7])
                item_rgb = fields[8]
                block_count = int(fields[9])
                block_sizes = [int(x) for x in fields[10].strip(',').split(',') if x]
                block_starts = [int(x) for x in fields[11].strip(',').split(',') if x]
                
                entries.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                    'score': score,
                    'strand': strand,
                    'thick_start': thick_start,
                    'thick_end': thick_end,
                    'item_rgb': item_rgb,
                    'block_count': block_count,
                    'block_sizes': block_sizes,
                    'block_starts': block_starts,
                    'mature_length': sum(block_sizes),
                    'genomic_span': end - start
                })
            except Exception as e:
                logging.warning(f"Error parsing line {i+1}: {e}")
                
    return entries

def find_bedtools():
    """Find bedtools executable"""
    # Try to find bedtools in common locations
    bedtools_paths = [
        "bedtools",
        "/usr/bin/bedtools",
        "/usr/local/bin/bedtools",
        os.path.expanduser("~/.conda/envs/bed12/bin/bedtools")
    ]
    
    for path in bedtools_paths:
        try:
            subprocess.run([path, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return path
        except:
            continue
    
    # If we get here, we couldn't find bedtools
    logging.error("Could not find bedtools executable")
    return "bedtools"  # Default to system path even though it might fail

def run_bedtools_command(command):
    """Run a bedtools command and return the output"""
    try:
        result = subprocess.run(command, shell=True, check=True, 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Bedtools command failed: {e}")
        logging.error(f"Command: {command}")
        logging.error(f"Error: {e.stderr}")
        return ""

def classify_circrna(entry, bedtools_cmd, exons_gtf, genes_gtf, temp_dir):
    """Classify a single circRNA by genomic context"""
    try:
        chrom = entry['chrom']
        start = entry['start']
        end = entry['end']
        name = entry['name']
        strand = entry['strand']
        
        # Create temporary BED file for this entry
        temp_bed = os.path.join(temp_dir, f"{name.replace('|', '_')}.bed")
        with open(temp_bed, 'w') as f:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
        
        # Check if circRNA overlaps with any gene
        gene_cmd = f"{bedtools_cmd} intersect -a {temp_bed} -b {genes_gtf} -s -wa -wb"
        gene_result = run_bedtools_command(gene_cmd)
        
        # If no gene overlap, it's intergenic
        if not gene_result.strip():
            os.remove(temp_bed)
            return "intergenic"
        
        # Check overlap with exons using -split option
        exon_cmd = f"{bedtools_cmd} intersect -a {temp_bed} -b {exons_gtf} -s -wa -wb"
        exon_result = run_bedtools_command(exon_cmd)
        
        # Check if fully exonic
        full_exon_cmd = f"{bedtools_cmd} intersect -a {temp_bed} -b {exons_gtf} -s -f 1.0 -wa"
        full_exon_result = run_bedtools_command(full_exon_cmd)
        
        os.remove(temp_bed)
        
        # Determine type based on overlaps
        if full_exon_result.strip():
            return "eciRNA"  # Fully exonic
        elif exon_result.strip():
            return "EIciRNA"  # Mix of exonic and intronic
        else:
            return "ciRNA"  # Completely intronic within a gene
            
    except Exception as e:
        logging.warning(f"Error classifying circRNA {name}: {e}")
        return "unknown"

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                 'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                 'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def analyze_backsplice_junction(entry, genome_ref):
    """Analyze the backsplice junction of a circRNA"""
    try:
        chrom = entry['chrom']
        start = entry['start']
        end = entry['end']
        strand = entry['strand']
        block_starts = entry['block_starts']
        block_sizes = entry['block_sizes']
        block_count = entry['block_count']
        
        # Calculate junction positions
        if block_count == 1:
            # For single-exon circRNAs, backsplice junction connects the end and start of the same exon
            exon_start = start + block_starts[0]
            exon_end = exon_start + block_sizes[0]
            
            donor_region = (exon_end - 2, exon_end + 2)
            acceptor_region = (exon_start - 2, exon_start + 2)
        else:
            # For multi-exon circRNAs, backsplice junction connects the last and first exons
            first_exon_start = start + block_starts[0]
            last_exon_start = start + block_starts[-1]
            last_exon_end = last_exon_start + block_sizes[-1]
            
            donor_region = (last_exon_end - 2, last_exon_end + 2)
            acceptor_region = (first_exon_start - 2, first_exon_start + 2)
        
        # Fetch sequences
        donor_seq = genome_ref.fetch(chrom, donor_region[0], donor_region[1])
        acceptor_seq = genome_ref.fetch(chrom, acceptor_region[0], acceptor_region[1])
        
        # Adjust for strand
        if strand == '-':
            donor_seq = reverse_complement(donor_seq)
            acceptor_seq = reverse_complement(acceptor_seq)
            # Swap for negative strand
            donor_seq, acceptor_seq = acceptor_seq, donor_seq
        
        # Extract splice sites
        donor = donor_seq[2:4]
        acceptor = acceptor_seq[0:2]
        
        # Determine splice site type
        if (donor, acceptor) == ('GT', 'AG'):
            splice_type = "GT-AG"
            is_canonical = True
        elif (donor, acceptor) == ('GC', 'AG'):
            splice_type = "GC-AG"
            is_canonical = True
        elif (donor, acceptor) == ('AT', 'AC'):
            splice_type = "AT-AC"
            is_canonical = True
        else:
            splice_type = "Non-canonical"
            is_canonical = False
        
        return {
            'donor_site': donor,
            'acceptor_site': acceptor,
            'donor_context': donor_seq,
            'acceptor_context': acceptor_seq,
            'splice_site_type': splice_type,
            'is_canonical': is_canonical
        }
        
    except Exception as e:
        logging.warning(f"Error analyzing backsplice junction for {entry['name']}: {e}")
        return {
            'donor_site': '',
            'acceptor_site': '',
            'donor_context': '',
            'acceptor_context': '',
            'splice_site_type': 'Unknown',
            'is_canonical': False
        }

def extract_gtf_gene_info(entry, bedtools_cmd, gtf_file, temp_dir):
    """Extract gene information from GTF for a circRNA"""
    try:
        chrom = entry['chrom']
        start = entry['start']
        end = entry['end']
        name = entry['name']
        strand = entry['strand']
        
        # Create temporary BED file for this entry
        temp_bed = os.path.join(temp_dir, f"{name.replace('|', '_')}_gene.bed")
        with open(temp_bed, 'w') as f:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
        
        # Intersect with GTF
        intersect_cmd = f"{bedtools_cmd} intersect -a {temp_bed} -b {gtf_file} -s -wa -wb"
        intersect_result = run_bedtools_command(intersect_cmd)
        
        # Parse results
        gene_names = set()
        gene_types = set()
        transcript_types = set()
        
        for line in intersect_result.strip().split('\n'):
            if not line:
                continue
                
            fields = line.split('\t')
            if len(fields) < 15:  # BED(6) + GTF(9+)
                continue
                
            # Extract attributes from GTF
            attr_str = fields[14]
            attrs = {}
            
            for attr in attr_str.split(';'):
                if not attr.strip():
                    continue
                    
                try:
                    key, value = attr.strip().split(' ', 1)
                    attrs[key] = value.strip('"')
                except ValueError:
                    pass
            
            if 'gene_name' in attrs:
                gene_names.add(attrs['gene_name'])
            if 'gene_type' in attrs:
                gene_types.add(attrs['gene_type'])
            elif 'gene_biotype' in attrs:
                gene_types.add(attrs['gene_biotype'])
            if 'transcript_type' in attrs:
                transcript_types.add(attrs['transcript_type'])
            elif 'transcript_biotype' in attrs:
                transcript_types.add(attrs['transcript_biotype'])
        
        os.remove(temp_bed)
        
        # Default values if sets are empty
        if not gene_names:
            gene_names.add('Unknown')
        if not gene_types:
            gene_types.add('Unknown')
        if not transcript_types:
            transcript_types.add('Unknown')
            
        return {
            'gene_names': ','.join(gene_names),
            'gene_types': ','.join(gene_types),
            'transcript_types': ','.join(transcript_types)
        }
        
    except Exception as e:
        logging.warning(f"Error extracting gene info for {name}: {e}")
        return {
            'gene_names': 'Unknown',
            'gene_types': 'Unknown',
            'transcript_types': 'Unknown'
        }

def process_circrnas(bed_file, genome_ref, gtf_file, exons_gtf, genes_gtf, temp_dir, max_entries=0):
    """Process all circRNAs in a BED12 file"""
    # Find bedtools
    bedtools_cmd = find_bedtools()
    
    # Parse BED12 entries
    entries = parse_bed12_entries(bed_file, max_entries)
    logging.info(f"Parsed {len(entries)} entries from {bed_file}")
    
    # Process each entry
    results = []
    for i, entry in enumerate(entries):
        if i % 100 == 0:
            logging.info(f"Processing entry {i+1}/{len(entries)}")
        
        # Classify circRNA
        circrna_type = classify_circrna(entry, bedtools_cmd, exons_gtf, genes_gtf, temp_dir)
        
        # Analyze backsplice junction
        junction_info = analyze_backsplice_junction(entry, genome_ref)
        
        # Extract gene information
        gene_info = extract_gtf_gene_info(entry, bedtools_cmd, gtf_file, temp_dir)
        
        # Combine all information
        result = {
            'chrom': entry['chrom'],
            'start': entry['start'],
            'end': entry['end'],
            'name': entry['name'],
            'strand': entry['strand'],
            'block_count': entry['block_count'],
            'mature_length': entry['mature_length'],
            'genomic_span': entry['genomic_span'],
            'circrna_type': circrna_type,
            'gene_names': gene_info['gene_names'],
            'gene_types': gene_info['gene_types'],
            'transcript_types': gene_info['transcript_types'],
            'donor_site': junction_info['donor_site'],
            'acceptor_site': junction_info['acceptor_site'],
            'donor_context': junction_info['donor_context'],
            'acceptor_context': junction_info['acceptor_context'],
            'splice_site_type': junction_info['splice_site_type'],
            'is_canonical': junction_info['is_canonical']
        }
        
        results.append(result)
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Analyze circRNAs')
    parser.add_argument('--bed', required=True, help='Path to BED12 file with circRNAs')
    parser.add_argument('--genome', required=True, help='Path to reference genome FASTA')
    parser.add_argument('--gtf', required=True, help='Path to GTF annotation file')
    parser.add_argument('--exons', required=True, help='Path to exons GTF file')
    parser.add_argument('--genes', required=True, help='Path to genes GTF file')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--max', type=int, default=0, help='Maximum number of entries to process (0 = all)')
    args = parser.parse_args()
    
    # Configure logging
    configure_logging(args.output)
    
    # Create temporary directory
    temp_dir = os.path.join(args.output, 'temp')
    os.makedirs(temp_dir, exist_ok=True)
    
    try:
        # Load genome reference
        logging.info(f"Loading genome reference: {args.genome}")
        genome_ref = pysam.FastaFile(args.genome)
        
        # Process circRNAs
        logging.info(f"Processing circRNAs from: {args.bed}")
        results = process_circrnas(
            args.bed, genome_ref, args.gtf, args.exons, args.genes, temp_dir, args.max
        )
        
        # Convert results to DataFrame
        df = pd.DataFrame(results)
        
        # Calculate additional metrics
        if len(df) > 0:
            # Count by type
            type_counts = df['circrna_type'].value_counts()
            logging.info(f"circRNA type distribution: {type_counts.to_dict()}")
            
            # Count by splice site type
            splice_counts = df['splice_site_type'].value_counts()
            logging.info(f"Splice site distribution: {splice_counts.to_dict()}")
            
            # Generate summary statistics
            summary = {
                'total_circrnas': len(df),
                'mean_mature_length': df['mature_length'].mean(),
                'median_mature_length': df['mature_length'].median(),
                'mean_genomic_span': df['genomic_span'].mean(),
                'median_genomic_span': df['genomic_span'].median(),
                'mean_exon_count': df['block_count'].mean(),
                'canonical_count': df['is_canonical'].sum(),
                'multi_exon_count': (df['block_count'] > 1).sum(),
                'single_exon_count': (df['block_count'] == 1).sum(),
                'eciRNA_count': (df['circrna_type'] == 'eciRNA').sum(),
                'EIciRNA_count': (df['circrna_type'] == 'EIciRNA').sum(),
                'ciRNA_count': (df['circrna_type'] == 'ciRNA').sum(),
                'intergenic_count': (df['circrna_type'] == 'intergenic').sum()
            }
            
            # Save summary to file
            pd.DataFrame([summary]).to_csv(os.path.join(args.output, 'summary_stats.csv'), index=False)
        
        # Save full results to CSV
        csv_path = os.path.join(args.output, 'circRNA_results.csv')
        df.to_csv(csv_path, index=False)
        logging.info(f"Results saved to {csv_path}")
        
        # Clean up
        try:
            shutil.rmtree(temp_dir)
        except:
            logging.warning(f"Could not remove temporary directory: {temp_dir}")
        
    except Exception as e:
        logging.error(f"Error: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    main()
ENDPYTHON

# Make the script executable
chmod +x "$TEMP_DIR/analyze_circrna.py"

# --------------------------------
# Main Processing Pipeline
# --------------------------------

log "Starting circRNA analysis pipeline at $(date)"
log "Results will be saved to: $RESULTS_DIR"

# Find bedtools
BEDTOOLS=$(find_bedtools)
log "Using bedtools: $BEDTOOLS"

# Create filtered GTF files for faster processing
log "Extracting exons and genes from GTF..."
grep -P '\texon\t' "$GTF_FILE" > "$TEMP_DIR/exons.gtf"
grep -P '\tgene\t' "$GTF_FILE" > "$TEMP_DIR/genes.gtf"

# Show database sizes
log "Database sizes:"
log "CircAtlas: $(wc -l < "$CIRCATLAS_BED") entries"
log "CircBase: $(wc -l < "$CIRCBASE_BED") entries"

# Create test subsets if requested
if [[ "$TEST_MODE" == true ]]; then
    log "Creating test subsets with $MAX_ENTRIES entries..."
    head -n "$MAX_ENTRIES" "$CIRCATLAS_BED" > "$TEMP_DIR/circatlas_test.bed"
    head -n "$MAX_ENTRIES" "$CIRCBASE_BED" > "$TEMP_DIR/circbase_test.bed"
    CIRCATLAS_ANALYSIS_BED="$TEMP_DIR/circatlas_test.bed"
    CIRCBASE_ANALYSIS_BED="$TEMP_DIR/circbase_test.bed"
else
    CIRCATLAS_ANALYSIS_BED="$CIRCATLAS_BED"
    CIRCBASE_ANALYSIS_BED="$CIRCBASE_BED"
fi

# Process CircAtlas
if [[ "$PROCESS_CIRCATLAS" == true ]]; then
    log "Processing CircAtlas database..."
    
    # Run the analysis script
    python "$TEMP_DIR/analyze_circrna.py" \
        --bed "$CIRCATLAS_ANALYSIS_BED" \
        --genome "$GENOME_REF" \
        --gtf "$GTF_FILE" \
        --exons "$TEMP_DIR/exons.gtf" \
        --genes "$TEMP_DIR/genes.gtf" \
        --output "$RESULTS_DIR/circatlas" \
        --max "$MAX_ENTRIES"
    
    if [[ $? -ne 0 ]]; then
        log "Error processing CircAtlas"
    else
        log "CircAtlas processing complete"
    fi
fi

# Process CircBase
if [[ "$PROCESS_CIRCBASE" == true ]]; then
    log "Processing CircBase database..."
    
    # Run the analysis script
    python "$TEMP_DIR/analyze_circrna.py" \
        --bed "$CIRCBASE_ANALYSIS_BED" \
        --genome "$GENOME_REF" \
        --gtf "$GTF_FILE" \
        --exons "$TEMP_DIR/exons.gtf" \
        --genes "$TEMP_DIR/genes.gtf" \
        --output "$RESULTS_DIR/circbase" \
        --max "$MAX_ENTRIES"
    
    if [[ $? -ne 0 ]]; then
        log "Error processing CircBase"
    else
        log "CircBase processing complete"
    fi
fi

# Process Intersection
if [[ "$PROCESS_INTERSECTION" == true ]]; then
    log "Creating intersection file..."
    
    # Create sorted unique versions first to ensure accurate intersection
    sort -k1,1 -k2,2n -k3,3n "$CIRCATLAS_BED" | uniq > "$TEMP_DIR/circatlas_sorted.bed"
    sort -k1,1 -k2,2n -k3,3n "$CIRCBASE_BED" | uniq > "$TEMP_DIR/circbase_sorted.bed"
    
    # Create intersection with proper coordinate matching
    $BEDTOOLS intersect -a "$TEMP_DIR/circatlas_sorted.bed" \
                        -b "$TEMP_DIR/circbase_sorted.bed" \
                        -wa -f 1.0 -r > "$INTERSECTION_BED"
    
    log "Created intersection file with $(wc -l < "$INTERSECTION_BED") entries"
    
    # Create test subset if requested
    if [[ "$TEST_MODE" == true ]]; then
        head -n "$MAX_ENTRIES" "$INTERSECTION_BED" > "$TEMP_DIR/intersection_test.bed"
        INTERSECTION_ANALYSIS_BED="$TEMP_DIR/intersection_test.bed"
    else
        INTERSECTION_ANALYSIS_BED="$INTERSECTION_BED"
    fi
    
    # Process the intersection
    log "Processing intersection dataset..."
    
    python "$TEMP_DIR/analyze_circrna.py" \
        --bed "$INTERSECTION_ANALYSIS_BED" \
        --genome "$GENOME_REF" \
        --gtf "$GTF_FILE" \
        --exons "$TEMP_DIR/exons.gtf" \
        --genes "$TEMP_DIR/genes.gtf" \
        --output "$RESULTS_DIR/intersection" \
        --max "$MAX_ENTRIES"
    
    if [[ $? -ne 0 ]]; then
        log "Error processing intersection"
    else
        log "Intersection processing complete"
    fi
fi

# Create combined summary report
log "Creating combined summary report..."

if [ -f "$RESULTS_DIR/circatlas/summary_stats.csv" ] && \
   [ -f "$RESULTS_DIR/circbase/summary_stats.csv" ] && \
   [ -f "$RESULTS_DIR/intersection/summary_stats.csv" ]; then
    
    mkdir -p "$RESULTS_DIR/combined"
    
    # Create combined summary table
    echo "# Combined circRNA Database Summary" > "$RESULTS_DIR/combined/combined_summary.md"
    echo "## Database Statistics" >> "$RESULTS_DIR/combined/combined_summary.md"
    echo "| Metric | CircAtlas | CircBase | Intersection |" >> "$RESULTS_DIR/combined/combined_summary.md"
    echo "|--------|-----------|----------|--------------|" >> "$RESULTS_DIR/combined/combined_summary.md"
    
    # Extract metrics from summary files
    for metric in total_circrnas mean_mature_length median_mature_length \
                  mean_genomic_span median_genomic_span mean_exon_count \
                  canonical_count multi_exon_count single_exon_count \
                  eciRNA_count EIciRNA_count ciRNA_count intergenic_count; do
        
        atlas_val=$(grep -m1 "$metric" "$RESULTS_DIR/circatlas/summary_stats.csv" | cut -d, -f2)
        base_val=$(grep -m1 "$metric" "$RESULTS_DIR/circbase/summary_stats.csv" | cut -d, -f2)
        intersect_val=$(grep -m1 "$metric" "$RESULTS_DIR/intersection/summary_stats.csv" | cut -d, -f2)
        
        # Format the metric name for display
        display_metric=$(echo "$metric" | sed 's/_/ /g' | sed 's/\b\(.\)/\u\1/g')
        
        echo "| $display_metric | $atlas_val | $base_val | $intersect_val |" >> "$RESULTS_DIR/combined/combined_summary.md"
    done
    
    # Add dataset-specific metrics if needed
    
    log "Combined summary created at $RESULTS_DIR/combined/combined_summary.md"
else
    log "Cannot create combined summary: not all database results are available"
fi

# Clean up temporary files
if [[ "$TEST_MODE" != true ]]; then  # Keep temp files in test mode for debugging
    log "Cleaning up temporary files..."
    rm -rf "$TEMP_DIR"
fi

log "Analysis pipeline completed at $(date)"
log "Results are available in: $RESULTS_DIR"

# Display summary statistics
for db in circatlas circbase intersection; do
    if [ -f "$RESULTS_DIR/$db/summary_stats.csv" ]; then
        echo "================================================"
        echo "Summary Statistics for $db:"
        echo "================================================"
        cat "$RESULTS_DIR/$db/summary_stats.csv" | \
            sed 's/,/\t/g' | \
            column -t
        echo ""
    fi
done

exit 0