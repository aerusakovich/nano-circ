#!/usr/bin/env python3

import argparse
import os
import logging
import re
import subprocess
import numpy as np
import pandas as pd
import pybedtools
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

def configure_logging(output_dir):
    """Configure logging for the script."""
    os.makedirs(output_dir, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(os.path.join(output_dir, 'circRNA_analysis.log')),
            logging.StreamHandler()
        ]
    )

def rev_comp(seq):
    """Reverse complement of a DNA sequence."""
    return seq.translate(str.maketrans('ATCG', 'TAGC'))[::-1]

def check_canonical_splice_sites(seq, blocks, strand):
    """
    Check if circRNA has canonical splice site motifs.
    
    Args:
        seq (str): Full genomic sequence
        blocks (list): List of (start, size) tuples for exon blocks
        strand (str): Strand direction + or -
    
    Returns:
        tuple: (is_canonical, donor_site, acceptor_site)
    """
    if not seq or len(blocks) < 2:
        return False, "", ""
    
    try:
        if strand == '+':
            # For + strand:
            # Donor site: last 2nt of first exon + first 2nt of first intron
            first_exon_end = blocks[0][0] + blocks[0][1]
            donor = seq[first_exon_end-2:first_exon_end+2]
            
            # Acceptor site: last 2nt of last intron + first 2nt of last exon
            last_exon_start = blocks[-1][0]
            acceptor = seq[last_exon_start-2:last_exon_start+2]
        else:
            # For - strand: reverse complement and swap
            # Donor site: first 2nt of first intron + last 2nt of first exon (rev comp)
            last_exon_end = blocks[-1][0] + blocks[-1][1]
            donor = rev_comp(seq[last_exon_end-2:last_exon_end+2])
            
            # Acceptor site: first 2nt of last exon + last 2nt of last intron (rev comp)
            first_exon_start = blocks[0][0]
            acceptor = rev_comp(seq[first_exon_start-2:first_exon_start+2])
        
        # Check if splice sites are canonical
        is_canonical = False
        donor_site = donor[0:2] if len(donor) >= 2 else ""
        acceptor_site = acceptor[-2:] if len(acceptor) >= 2 else ""
        
        if (donor_site, acceptor_site) in [('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')]:
            is_canonical = True
            
        return is_canonical, donor_site, acceptor_site
    except Exception as e:
        logging.warning(f"Error checking splice sites: {e}")
        return False, "", ""

def create_visualizations(results_df, output_dir):
    """Generate visualizations for circRNA features."""
    # Create plots directory if not exists
    os.makedirs(os.path.join(output_dir, 'plots'), exist_ok=True)
    
    # Visualization helper function
    def plot_pie(data, title, filename):
        plt.figure(figsize=(10, 8))
        data_counts = data.value_counts()
        plt.pie(data_counts, labels=data_counts.index, autopct='%1.1f%%')
        plt.title(title)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'plots', filename))
        plt.close()
    
    # Visualizations based on available columns
    if 'circrna_type' in results_df.columns and not results_df['circrna_type'].isna().all():
        plot_pie(results_df['circrna_type'], 'CircRNA Types', 'circrna_types.png')
    
    if 'genomic_context' in results_df.columns and not results_df['genomic_context'].isna().all():
        plot_pie(results_df['genomic_context'], 'Genomic Context', 'genomic_context.png')
    
    if 'gene_biotype' in results_df.columns and not results_df['gene_biotype'].isna().all():
        plot_pie(results_df['gene_biotype'], 'Gene Biotypes', 'gene_biotypes.png')
    
    # Canonical splice site visualization
    if 'is_canonical' in results_df.columns:
        plot_pie(results_df['is_canonical'], 'Canonical Splice Sites', 'canonical_splice_sites.png')
    
    # Distributions by circRNA type
    if 'mature_length' in results_df.columns and 'circrna_type' in results_df.columns:
        # Length distribution for all circRNAs
        plt.figure(figsize=(10, 6))
        sns.histplot(results_df['mature_length'], kde=True, bins=50)
        plt.title('Mature Length Distribution')
        plt.xlabel('Length (bp)')
        plt.ylabel('Frequency')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'plots', 'mature_length_dist.png'))
        plt.close()
        
        # Length distribution by type
        valid_types = results_df[results_df['circrna_type'] != 'Unknown']
        if not valid_types.empty:
            plt.figure(figsize=(12, 8))
            sns.boxplot(x='circrna_type', y='mature_length', data=valid_types)
            plt.title('Mature Length Distribution by CircRNA Type')
            plt.xlabel('CircRNA Type')
            plt.ylabel('Length (bp)')
            plt.yscale('log')  # Log scale for better visualization
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, 'plots', 'length_by_type.png'))
            plt.close()
            
            # Summary statistics table for lengths by type
            length_stats = valid_types.groupby('circrna_type')['mature_length'].agg([
                'count', 'min', 'max', 'mean', 'median', 'std'
            ]).round(2)
            
            length_stats.to_csv(os.path.join(output_dir, 'length_stats_by_type.csv'))

def map_with_minimap2(fasta_file, reference_genome, output_dir, threads=4):
    """
    Map sequences to the reference genome using minimap2.
    
    Args:
        fasta_file: Path to FASTA file with sequences
        reference_genome: Path to reference genome FASTA
        output_dir: Directory to save outputs
        threads: Number of threads for minimap2
        
    Returns:
        str: Path to the SAM file with mapping results
    """
    # Check if minimap2 is installed
    try:
        subprocess.run(['minimap2', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except (subprocess.SubprocessError, FileNotFoundError):
        logging.error("minimap2 is not installed or not in PATH. Please install minimap2.")
        return None
    
    # Create output file path for SAM format (better for spliced alignments)
    sam_file = os.path.join(output_dir, 'minimap2_mapping.sam')
    
    # Set up minimap2 command for spliced alignment
    cmd = [
        'minimap2',
        '-a',              # Output SAM format
        '-x', 'splice',    # Splice-aware alignment mode
        '--cs',            # Output CS tag for CIGAR-like operations
        '-t', str(threads),
        '--secondary=no',  # Suppress secondary alignments
        reference_genome,
        fasta_file
    ]
    
    try:
        logging.info(f"Running minimap2: {' '.join(cmd)}")
        with open(sam_file, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        
        logging.info(f"Minimap2 mapping completed: {sam_file}")
        return sam_file
    
    except subprocess.SubprocessError as e:
        logging.error(f"minimap2 mapping failed: {e}")
        if e.stderr:
            logging.error(f"Error message: {e.stderr.decode('utf-8')}")
        return None

def convert_sam_to_bed12(sam_file, output_dir):
    """
    Convert SAM alignments to BED12 format for circRNA analysis.
    
    Args:
        sam_file: Path to SAM file from minimap2
        output_dir: Directory to save outputs
        
    Returns:
        str: Path to BED12 file
    """
    if not sam_file or not os.path.exists(sam_file):
        logging.error(f"SAM file does not exist: {sam_file}")
        return None
    
    bed12_file = os.path.join(output_dir, 'circrna_mapped_bed12.bed')
    
    try:
        # Read SAM file with pysam
        alignments = pysam.AlignmentFile(sam_file, "r")
        
        # Process each alignment to BED12 format
        with open(bed12_file, 'w') as out:
            for read in alignments:
                # Skip unmapped reads
                if read.is_unmapped:
                    continue
                
                # Extract information
                chrom = alignments.get_reference_name(read.reference_id)
                start = read.reference_start
                end = read.reference_end
                name = read.query_name
                score = read.mapping_quality
                strand = '-' if read.is_reverse else '+'
                
                # Extract blocks (exons) from CIGAR
                blocks = []
                position = start
                
                for op, length in read.cigartuples:
                    if op == 0:  # M - match or mismatch
                        blocks.append((position, length))
                        position += length
                    elif op == 1:  # I - insertion to reference
                        pass  # Insertions don't affect reference position
                    elif op == 2:  # D - deletion from reference
                        position += length
                    elif op in (3, 7, 8):  # N, = or X - Skipped region, match or mismatch
                        position += length
                    elif op == 4:  # S - soft clipping
                        pass  # Soft clips don't affect reference position
                
                # If no blocks, create a single block
                if not blocks:
                    blocks = [(start, end - start)]
                
                # Format BED12 entry
                block_count = len(blocks)
                block_sizes = ','.join(str(size) for pos, size in blocks) + ','
                block_starts = ','.join(str(pos - start) for pos, size in blocks) + ','
                
                # Write BED12 line
                bed12_line = [
                    chrom, str(start), str(end), name, str(score), strand,
                    str(start), str(end), "0", str(block_count), block_sizes, block_starts
                ]
                out.write('\t'.join(bed12_line) + '\n')
        
        logging.info(f"Created BED12 file: {bed12_file}")
        return bed12_file
    
    except Exception as e:
        logging.error(f"Failed to convert SAM to BED12: {e}")
        logging.exception("Stack trace:")
        return None

def classify_circrnas_with_bedtools(bed12_file, gtf_data, ref_genome, output_dir):
    """
    Classify circRNAs based on BedTools analysis using the specified logic.
    
    Args:
        bed12_file: Path to BED12 file with mappings
        gtf_data: GTF annotation as BedTool (filtered to exons)
        ref_genome: Reference genome as pysam.FastaFile
        output_dir: Output directory
    
    Returns:
        pd.DataFrame: Classification results
    """
    if not bed12_file or not os.path.exists(bed12_file):
        logging.error(f"BED12 file does not exist: {bed12_file}")
        return pd.DataFrame()
    
    try:
        circ_bed12 = pybedtools.BedTool(bed12_file)
        results = []

        for circ in circ_bed12:
            chrom = circ.chrom
            start = circ.start
            end = circ.end
            name = circ.name
            strand = circ.strand if circ.strand in ('+', '-') else '.'

            # Parse BED12 blocks
            if len(circ.fields) >= 11:
                block_sizes = [int(x) for x in circ.fields[10].rstrip(',').split(',') if x]
                block_starts = [int(x) for x in circ.fields[11].rstrip(',').split(',') if x]
                blocks = [(start + s, size) for s, size in zip(block_starts, block_sizes)]
                mature_length = sum(block_sizes) if block_sizes else end - start
            else:
                blocks = [(start, end - start)]
                mature_length = end - start

            # Splice site sequence check
            try:
                full_seq = ref_genome.fetch(chrom, start, end)
                is_canonical, donor, acceptor = check_canonical_splice_sites(full_seq, blocks, strand)
            except Exception as e:
                logging.warning(f"Error fetching sequence for {name}: {e}")
                full_seq, is_canonical, donor, acceptor = "", False, "", ""

            # Create region string for BedTool
            circ_region_str = f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand if strand in ('+', '-') else '+'}"
            circ_region = pybedtools.BedTool(circ_region_str, from_string=True)

            # Intersect logic with fallback if strand is unknown
            s_flag = strand in ('+', '-')
            genic_overlaps = gtf_data.intersect(circ_region, wa=True, s=s_flag)
            is_genic = len(genic_overlaps) > 0

            if not is_genic:
                row = {
                    'name': name,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'circrna_type': 'Unknown',
                    'genomic_context': 'intergenic',
                    'gene_biotype': 'NA',
                    'fraction_exonic': 0,
                    'is_canonical': is_canonical,
                    'donor_site': donor,
                    'acceptor_site': acceptor,
                    'mature_length': mature_length
                }
                results.append(row)
                continue

            exonic_overlaps = gtf_data.intersect(circ_region, wa=True, s=s_flag, split=True)
            has_exonic_overlap = len(exonic_overlaps) > 0

            exon_coverage = gtf_data.intersect(circ_region, wo=True, s=s_flag, split=True)
            total_overlap_bp = 0
            for overlap in exon_coverage:
                try:
                    total_overlap_bp += int(overlap[-1])
                except (IndexError, ValueError):
                    continue
            circrna_length = end - start
            fraction_exonic = total_overlap_bp / circrna_length if circrna_length > 0 else 0

            gene_biotypes = set()
            for feature in genic_overlaps:
                biotype = feature.attrs.get('gene_biotype', '') or feature.attrs.get('gene_type', 'unknown')
                if biotype:
                    gene_biotypes.add(biotype)
            gene_biotype = ','.join(sorted(gene_biotypes)) if gene_biotypes else 'unknown'

            if not has_exonic_overlap:
                circrna_type = 'ciRNA'
                genomic_context = 'intronic'
            elif fraction_exonic == 1.0:
                circrna_type = 'ecircRNA'
                genomic_context = 'exonic'
            else:
                circrna_type = 'EIcircRNA'
                genomic_context = 'exon-intron'

            row = {
                'name': name,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
                'circrna_type': circrna_type,
                'genomic_context': genomic_context,
                'gene_biotype': gene_biotype,
                'fraction_exonic': round(fraction_exonic, 3),
                'is_canonical': is_canonical,
                'donor_site': donor,
                'acceptor_site': acceptor,
                'mature_length': mature_length
            }
            results.append(row)

        return pd.DataFrame(results)

    except Exception as e:
        logging.error(f"Failed to classify circRNAs with BedTools: {e}")
        logging.exception("Stack trace:")
        return pd.DataFrame()

def process_circrna(fasta_record):
    """
    Process a single circRNA record for basic sequence properties.
    
    Args:
        fasta_record: SeqIO FASTA record
    
    Returns:
        dict: Basic circRNA features
    """
    try:
        record_id = fasta_record.id
        seq = str(fasta_record.seq)
        
        if not seq:
            return {
                'name': record_id,
                'sequence': seq,
                'mature_length': 0,
                'gc_content': 0
            }
        
        # Calculate GC content
        seq_lower = seq.lower()
        gc_content = (seq_lower.count('g') + seq_lower.count('c')) / len(seq) * 100
        
        # Create result with basic properties
        result = {
            'name': record_id,
            'sequence': seq,
            'mature_length': len(seq),
            'gc_content': gc_content
        }
        
        return result
        
    except Exception as e:
        logging.error(f"Error processing {fasta_record.id}: {e}")
        return None

def create_temp_fasta(sequences, output_dir):
    """Create a temporary FASTA file for the sequences."""
    temp_fasta = os.path.join(output_dir, 'temp_sequences.fasta')
    with open(temp_fasta, 'w') as f:
        for name, seq in sequences:
            if seq:  # Skip empty sequences
                f.write(f">{name}\n{seq}\n")
    return temp_fasta

def generate_summary(results_df, output_dir):
    """Generate summary report and type-specific files."""
    # Type-specific FASTA files
    if 'circrna_type' in results_df.columns and 'sequence' in results_df.columns:
        for circ_type in results_df['circrna_type'].unique():
            if pd.isna(circ_type) or circ_type == 'Unknown':
                continue
                
            type_df = results_df[results_df['circrna_type'] == circ_type]
            if len(type_df) > 0:
                type_fasta = os.path.join(output_dir, f'{circ_type}_sequences.fasta')
                with open(type_fasta, 'w') as out_file:
                    for _, row in type_df.iterrows():
                        if 'sequence' in row and pd.notna(row['sequence']):
                            out_file.write(f">{row['name']}\n{row['sequence']}\n")
                logging.info(f"Saved {len(type_df)} {circ_type} sequences to {type_fasta}")
    
    # Summary report
    summary_file = os.path.join(output_dir, 'analysis_summary.txt')
    with open(summary_file, 'w') as summary:
        total_circrnas = len(results_df)
        
        summary.write(f"CircRNA Analysis Summary\n")
        summary.write(f"======================\n\n")
        summary.write(f"Total circRNAs analyzed: {total_circrnas}\n\n")
        
        if 'circrna_type' in results_df.columns:
            type_counts = results_df['circrna_type'].value_counts()
            summary.write(f"CircRNA Types:\n")
            for type_name, count in type_counts.items():
                if total_circrnas > 0:
                    summary.write(f"  - {type_name}: {count} ({count/total_circrnas*100:.1f}%)\n")
                else:
                    summary.write(f"  - {type_name}: {count}\n")
        
        if 'genomic_context' in results_df.columns:
            context_counts = results_df['genomic_context'].value_counts()
            summary.write(f"\nGenomic Context:\n")
            for context, count in context_counts.items():
                if total_circrnas > 0:
                    summary.write(f"  - {context}: {count} ({count/total_circrnas*100:.1f}%)\n")
                else:
                    summary.write(f"  - {context}: {count}\n")
        
        if 'is_canonical' in results_df.columns:
            canonical_count = results_df['is_canonical'].sum()
            if total_circrnas > 0:
                canonical_pct = (canonical_count / total_circrnas) * 100
                summary.write(f"\nSplice Sites:\n")
                summary.write(f"  - Canonical: {canonical_count} ({canonical_pct:.1f}%)\n")
                summary.write(f"  - Non-canonical: {total_circrnas - canonical_count} ({100 - canonical_pct:.1f}%)\n")
        
        if 'mature_length' in results_df.columns:
            length_stats = results_df['mature_length'].describe()
            summary.write(f"\nLength Statistics (bp):\n")
            summary.write(f"  - Min: {length_stats['min']:.0f}\n")
            summary.write(f"  - Max: {length_stats['max']:.0f}\n")
            summary.write(f"  - Mean: {length_stats['mean']:.1f}\n")
            summary.write(f"  - Median: {length_stats['50%']:.0f}\n")
            summary.write(f"  - Std Dev: {length_stats['std']:.1f}\n")
    
    logging.info(f"Generated analysis summary: {summary_file}")

def main():
    """Main function to process circRNA sequences."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='CircRNA Sequence Analysis')
    parser.add_argument('--fasta', required=True, help='Input multi-FASTA file')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA')
    parser.add_argument('--gtf', required=True, help='Gene annotation GTF')
    parser.add_argument('--output', default='circRNA_output', help='Output directory')
    parser.add_argument('--workers', type=int, default=4, help='Number of threads for minimap2')
    args = parser.parse_args()
    
    # Configure logging
    configure_logging(args.output)
    
    try:
        # Load GTF data - filter for exon features only
        logging.info(f"Loading GTF annotation: {args.gtf}")
        gtf_data = pybedtools.BedTool(args.gtf).filter(lambda x: x[2] == 'exon')
        
        # Load reference genome
        logging.info(f"Loading reference genome: {args.genome}")
        ref_genome = pysam.FastaFile(args.genome)
        
        # Process FASTA records
        logging.info(f"Processing FASTA file: {args.fasta}")
        fasta_records = list(SeqIO.parse(args.fasta, "fasta"))
        logging.info(f"Found {len(fasta_records)} records in FASTA file")
        
        # Process each record for basic sequence analysis
        basic_results = []
        sequence_pairs = []  # For mapping
        
        for record in fasta_records:
            result = process_circrna(record)
            if result:
                basic_results.append(result)
                sequence_pairs.append((record.id, str(record.seq)))
        
        logging.info(f"Successfully processed {len(basic_results)} records for basic properties")
        
        # Convert basic sequence results to DataFrame
        basic_df = pd.DataFrame(basic_results)
        
        # Save basic results
        basic_csv = os.path.join(args.output, 'circRNA_basic_properties.csv')
        basic_df.to_csv(basic_csv, index=False)
        logging.info(f"Basic sequence properties saved to {basic_csv}")
        
        # Map sequences to genome with minimap2 (SAM output)
        logging.info("Mapping sequences to reference genome with minimap2")
        temp_fasta = create_temp_fasta(sequence_pairs, args.output)
        
        sam_file = map_with_minimap2(
            fasta_file=temp_fasta,
            reference_genome=args.genome,
            output_dir=args.output,
            threads=args.workers
        )
        
        if sam_file:
            # Convert SAM to BED12
            bed12_file = convert_sam_to_bed12(sam_file, args.output)
            
            if bed12_file:
                # Classify circRNAs using bedtools according to the specified logic
                classified_df = classify_circrnas_with_bedtools(
                    bed12_file, gtf_data, ref_genome, args.output
                )
                
                if not classified_df.empty:
                    # Merge with basic properties
                    final_df = pd.merge(basic_df, classified_df, on='name', how='left')
                    
                    # Save final results
                    final_csv = os.path.join(args.output, 'circRNA_complete_analysis.csv')
                    final_df.to_csv(final_csv, index=False)
                    logging.info(f"Complete analysis results saved to {final_csv}")
                    
                    # Create visualizations for final results
                    create_visualizations(final_df, args.output)
                    
                    # Generate summary
                    generate_summary(final_df, args.output)
                else:
                    logging.warning("No circRNAs were classified. Using basic results only.")
                    create_visualizations(basic_df, args.output)
            else:
                logging.warning("Failed to convert SAM to BED12. Using basic results only.")
                create_visualizations(basic_df, args.output)
        else:
            logging.warning("Mapping with minimap2 failed. Using basic results only.")
            create_visualizations(basic_df, args.output)
        
        logging.info(f"Analysis complete. Results saved to {args.output}")
        
    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        logging.exception("Stack trace:")
        raise
    finally:
        # Clean up temporary files
        if 'temp_fasta' in locals() and os.path.exists(temp_fasta):
            try:
                os.remove(temp_fasta)
            except:
                pass
        
        # Close reference genome
        if 'ref_genome' in locals():
            ref_genome.close()

if __name__ == "__main__":
    main()