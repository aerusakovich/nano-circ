#!/usr/bin/env python3
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
import argparse
import subprocess

# Tool colors
TOOL_COLORS = {
    'ciri_long': '#1f77b4',  # Blue
    'isocirc': '#2ca02c',    # Green
    'circnick': '#ff7f0e'    # Orange
}

def extract_circrna_types_using_bedtools(tool_file, circrna_db, output_dir, tool_name, overlap="0.95"):
    """
    Use bedtools to intersect tool predictions with the circrna database
    and extract circRNA types
    """
    if not os.path.exists(tool_file) or not os.path.exists(circrna_db):
        print(f"Error: Missing input file(s) for bedtools intersection")
        return {}
    
    try:
        # Create a temporary directory for bedtools output
        temp_dir = os.path.join(output_dir, "temp_bedtools")
        os.makedirs(temp_dir, exist_ok=True)
        
        # Output file for bedtools intersection
        intersect_file = os.path.join(temp_dir, f"{tool_name}_vs_circrna_db_{overlap}.bed")
        
        # Run bedtools intersect command with fraction 0.95
        bedtools_cmd = f"bedtools intersect -a {tool_file} -b {circrna_db} -f {overlap} -wb > {intersect_file}"
        print(f"Running: {bedtools_cmd}")
        ret_code = os.system(bedtools_cmd)
        
        if ret_code != 0:
            print(f"Error: bedtools command failed with return code {ret_code}")
            return {}
        
        # Parse the intersection file to extract circRNA types and IDs
        type_counts = defaultdict(int)
        circrna_types = set(['EIciRNA', 'eciRNA', 'ciRNA', 'intergenic'])
        
        # Store matched circRNA IDs with their TPM IDs
        matched_circrnas = {}
        
        with open(intersect_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # The type information is in the name field of the database (column 4+N)
                # where N is the number of columns in the tool file
                if len(parts) >= 16:  # Assuming at least BED6 for tool + BED6 for db
                    db_name_col = len(parts) // 2 + 3  # Estimate the position of the db name column
                    db_name = parts[db_name_col]
                    
                    # Get coordinates from the tool prediction
                    tool_chrom = parts[0]
                    tool_start = parts[1]
                    tool_end = parts[2]
                    tool_strand = parts[5]
                    tool_key = f"{tool_chrom}:{tool_start}-{tool_end}:{tool_strand}"
                    
                    # Extract type from name field
                    found_type = False
                    if '|' in db_name:
                        for part in db_name.split('|'):
                            if part in circrna_types:
                                type_counts[part] += 1
                                found_type = True
                                # Store the full ID for TPM matching
                                matched_circrnas[tool_key] = db_name
                                break
                    
                    if not found_type:
                        type_counts['Unknown'] += 1
        
        print(f"Found {sum(type_counts.values())} circRNAs with types using bedtools intersect -f {overlap}")
        for circ_type, count in sorted(type_counts.items(), key=lambda x: x[1], reverse=True):
            if circ_type != 'Unknown':
                print(f"  {circ_type}: {count}")
        
        return matched_circrnas, {k: v for k, v in type_counts.items() if k != 'Unknown'}
    
    except Exception as e:
        print(f"Error during bedtools intersection: {e}")
        import traceback
        traceback.print_exc()
        return {}, {}

def extract_read_names_from_fastq(fastq_file, output_dir):
    """
    Extract all read names from a FASTQ file using more efficient method
    """
    print(f"Extracting read names from {fastq_file}")
    
    if not os.path.exists(fastq_file):
        print(f"Error: FASTQ file {fastq_file} does not exist")
        return set()
    
    try:
        # Create a temporary file to store read names
        temp_file = os.path.join(output_dir, "temp_fastq_read_ids.txt")
        
        # Use awk for more efficient extraction (similar to bash script approach)
        awk_cmd = f"awk '/^@/ {{print substr($1, 2)}}' {fastq_file} > {temp_file}"
        print(f"Running command: {awk_cmd}")
        
        # Use subprocess to run the awk command
        process = subprocess.run(awk_cmd, shell=True, capture_output=True, text=True)
        
        if process.returncode != 0:
            print(f"Error running awk command: {process.stderr}")
            return set()
        
        # Read the names into a set
        fastq_read_names = set()
        with open(temp_file, 'r') as f:
            for line in f:
                fastq_read_names.add(line.strip())
        
        print(f"Extracted {len(fastq_read_names)} read names from FASTQ file")
        
        # Save a sample of read names for debugging
        debug_dir = os.path.join(output_dir, "debug")
        os.makedirs(debug_dir, exist_ok=True)
        
        with open(os.path.join(debug_dir, "fastq_read_ids_sample.txt"), 'w') as f:
            for name in list(fastq_read_names)[:10]:  # First 10 names
                f.write(f"{name}\n")
        
        return fastq_read_names
    
    except Exception as e:
        print(f"Error extracting read names: {e}")
        import traceback
        traceback.print_exc()
        return set()

def extract_read_info_from_bed(bed_file, output_dir):
    """
    Extract read information from a BED12 file
    """
    print(f"Extracting read info from {bed_file}")
    
    if not os.path.exists(bed_file):
        print(f"Error: BED file {bed_file} does not exist")
        return {}
    
    try:
        # Read the BED file
        reads_info = {}
        read_names = []
        with open(bed_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:  # At least BED6 format
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    read_name = parts[3]
                    strand = parts[5]
                    
                    # Store the coordinates and read name
                    key = f"{chrom}:{start}-{end}:{strand}"
                    reads_info[key] = read_name
                    read_names.append(read_name)
        
        # Save a sample of read names for debugging
        debug_dir = os.path.join(output_dir, "debug")
        os.makedirs(debug_dir, exist_ok=True)
        
        with open(os.path.join(debug_dir, "bed_read_names_sample.txt"), 'w') as f:
            for name in read_names[:10]:  # First 10 names
                f.write(f"{name}\n")
        
        print(f"Extracted information for {len(reads_info)} reads from BED file")
        return reads_info
    
    except Exception as e:
        print(f"Error extracting read info: {e}")
        import traceback
        traceback.print_exc()
        return {}

def create_expression_boxplot(abundances_file, circrna_db, tool_files, ground_truth_bed, fastq_file, output_dir):
    """
    Create a boxplot showing expression levels (TPM) of circRNAs found by each tool
    """
    print(f"Creating expression boxplot from {abundances_file}")
    
    if not os.path.exists(abundances_file):
        print(f"Error: Abundances file {abundances_file} does not exist")
        return False
    
    try:
        # Read the abundances file
        abundances_df = pd.read_csv(abundances_file, sep='\t')
        print(f"Read {len(abundances_df)} entries from abundances file")
        
        # Extract circRNA IDs and TPM values
        abundances_df['circRNA_id'] = abundances_df['circRNA_id'].astype(str)
        abundances_df['TPM'] = abundances_df['TPM'].astype(float)
        
        # Print first few rows to verify format
        print("Sample abundances entries:")
        print(abundances_df.head())
        
        # Extract read names from FASTQ file (if provided)
        fastq_read_names = set()
        if fastq_file and os.path.exists(fastq_file):
            fastq_read_names = extract_read_names_from_fastq(fastq_file, output_dir)
        
        # Extract read info from ground truth BED file (if provided)
        ground_truth_reads = {}
        if ground_truth_bed and os.path.exists(ground_truth_bed):
            ground_truth_reads = extract_read_info_from_bed(ground_truth_bed, output_dir)
        
        # Run bedtools intersection for each tool to get matched circRNAs
        tools_matched_circrnas = {}
        for tool_name, tool_file in tool_files.items():
            if os.path.exists(tool_file) and os.path.exists(circrna_db):
                print(f"Running bedtools intersection for {tool_name}...")
                matched_circrnas, type_counts = extract_circrna_types_using_bedtools(
                    tool_file, circrna_db, output_dir, tool_name, "0.95")
                
                if matched_circrnas:
                    tools_matched_circrnas[tool_name] = matched_circrnas
        
        # Create a dictionary to store TPM values by tool and tracking info for zero TPM circRNAs
        tpm_by_tool = {tool: [] for tool in tools_matched_circrnas.keys()}
        tpm_by_tool['ground_truth'] = abundances_df['TPM'].tolist()  # All ground truth TPMs
        
        # Store circRNAs with zero TPM
        zero_tpm_circrnas = {tool: [] for tool in tools_matched_circrnas.keys()}
        
        # For each tool, get the TPM values for matched circRNAs
        for tool, matched_circrnas in tools_matched_circrnas.items():
            for coord_key, circrna_id in matched_circrnas.items():
                # Find this circRNA in the abundances file
                matching_entries = abundances_df[abundances_df['circRNA_id'] == circrna_id]
                
                if not matching_entries.empty:
                    # Get TPM values
                    tpm_values = matching_entries['TPM'].tolist()
                    tpm_by_tool[tool].extend(tpm_values)
                    
                    # Check for zero TPM
                    for tpm in tpm_values:
                        if tpm == 0.0:
                            zero_tpm_circrnas[tool].append({
                                'coord_key': coord_key,
                                'circrna_id': circrna_id,
                                'tpm': 0.0
                            })
                else:
                    # Try partial matching if exact match fails
                    found = False
                    for idx, row in abundances_df.iterrows():
                        if circrna_id in row['circRNA_id']:
                            tpm_by_tool[tool].append(row['TPM'])
                            
                            # Check for zero TPM
                            if row['TPM'] == 0.0:
                                zero_tpm_circrnas[tool].append({
                                    'coord_key': coord_key,
                                    'circrna_id': circrna_id,
                                    'tpm': 0.0
                                })
                            
                            found = True
                            break
                    
                    # If still not found, consider it missing TPM data
                    if not found:
                        zero_tpm_circrnas[tool].append({
                            'coord_key': coord_key,
                            'circrna_id': circrna_id,
                            'tpm': None,
                            'reason': 'No TPM data found'
                        })
        
        # Create a DataFrame for plotting
        plot_data = []
        for tool, tpm_values in tpm_by_tool.items():
            for tpm in tpm_values:
                plot_data.append({'Tool': tool, 'TPM': tpm})
        
        plot_df = pd.DataFrame(plot_data)
        
        # Print some statistics
        print("TPM statistics by tool:")
        
        # Create summary dataframe to save to file
        summary_data = []
        
        for tool, tpm_values in tpm_by_tool.items():
            if tpm_values:
                min_tpm = min(tpm_values)
                max_tpm = max(tpm_values)
                median_tpm = np.median(tpm_values)
                mean_tpm = np.mean(tpm_values)
                count = len(tpm_values)
                
                summary_data.append({
                    'Tool': tool,
                    'Count': count,
                    'Min_TPM': min_tpm,
                    'Max_TPM': max_tpm,
                    'Median_TPM': median_tpm,
                    'Mean_TPM': mean_tpm
                })
                
                print(f"  {tool}: {count} circRNAs, " +
                      f"Min: {min_tpm:.2f}, " +
                      f"Max: {max_tpm:.2f}, " +
                      f"Median: {median_tpm:.2f}, " +
                      f"Mean: {mean_tpm:.2f}")
            else:
                summary_data.append({
                    'Tool': tool,
                    'Count': 0,
                    'Min_TPM': 0,
                    'Max_TPM': 0,
                    'Median_TPM': 0,
                    'Mean_TPM': 0
                })
                print(f"  {tool}: No data")
        
        # Save summary statistics to file
        summary_df = pd.DataFrame(summary_data)
        
        # Create a better formatted summary table
        tool_display_names = {
            "ciri_long": "CIRI-long",
            "isocirc": "IsoCirc",
            "circnick": "CircNick-lrs",
            "ground_truth": "Ground Truth"
        }
        
        # Copy summary_df but with better tool names
        formatted_summary_df = summary_df.copy()
        formatted_summary_df['Tool'] = formatted_summary_df['Tool'].map(
            lambda x: tool_display_names.get(x, x))
        
        # Save both versions (one with original names for programmatic use,
        # one with formatted names for human readability)
        summary_df.to_csv(f"{output_dir}/expression_summary_stats.tsv", sep='\t', index=False)
        formatted_summary_df.to_csv(f"{output_dir}/expression_summary_stats_formatted.tsv", 
                                   sep='\t', index=False)
        
        # Analyze zero TPM circRNAs
        with open(f"{output_dir}/zero_tpm_circrnas_analysis.txt", 'w') as f:
            f.write("Analysis of circRNAs with Zero TPM\n")
            f.write("=================================\n\n")
            
            for tool, zero_circs in zero_tpm_circrnas.items():
                if zero_circs:
                    display_name = tool_display_names.get(tool, tool)
                    f.write(f"{display_name}:\n")
                    f.write(f"  Found {len(zero_circs)} circRNAs with zero TPM\n\n")
                    
                    # Check each zero TPM circRNA against ground truth and FASTQ
                    for i, circ in enumerate(zero_circs, 1):
                        coord_key = circ['coord_key']
                        circrna_id = circ['circrna_id']
                        
                        f.write(f"  {i}. CircRNA: {coord_key}\n")
                        f.write(f"     ID: {circrna_id}\n")
                        
                        # Check if in ground truth BED
                        in_ground_truth = coord_key in ground_truth_reads
                        if in_ground_truth:
                            gt_read_name = ground_truth_reads[coord_key]
                            f.write(f"     Found in ground truth BED: Yes\n")
                            f.write(f"     Ground truth read name: {gt_read_name}\n")
                            
                            # Check if read name is in FASTQ - direct comparison
                            in_fastq = gt_read_name in fastq_read_names
                            f.write(f"     Read name found in FASTQ: {'Yes' if in_fastq else 'No'}\n")
                        else:
                            f.write(f"     Found in ground truth BED: No\n")
                            
                            # Try to find by similar coordinates
                            similar_coords = [k for k in ground_truth_reads.keys() 
                                            if coord_key.split(':')[0] in k]  # Same chromosome
                            
                            if similar_coords:
                                f.write(f"     Found {len(similar_coords)} entries on same chromosome in ground truth\n")
                                
                                # List first few similar coordinates
                                for j, similar_coord in enumerate(similar_coords[:3], 1):
                                    f.write(f"       {j}. Similar coordinate: {similar_coord}\n")
                                    f.write(f"          Read name: {ground_truth_reads[similar_coord]}\n")
                                
                                if len(similar_coords) > 3:
                                    f.write(f"       ... and {len(similar_coords) - 3} more\n")
                            else:
                                f.write(f"     No similar coordinates found in ground truth\n")
                        
                        f.write("\n")
                else:
                    display_name = tool_display_names.get(tool, tool)
                    f.write(f"{tool_display_names.get(tool, tool)}: No circRNAs with zero TPM\n\n")
        
        # Write minimum TPM values with full precision
        with open(f"{output_dir}/minimum_tpm_full_precision.txt", 'w') as f:
            f.write("Minimum TPM Values with Full Precision\n")
            f.write("====================================\n\n")
            
            for tool, tpm_values in tpm_by_tool.items():
                if tpm_values:
                    min_tpm = min(tpm_values)
                    display_name = tool_display_names.get(tool, tool)
                    
                    f.write(f"{display_name}:\n")
                    
                    if min_tpm == 0.0:
                        # Find the smallest non-zero value if any
                        non_zero_values = [v for v in tpm_values if v > 0]
                        if non_zero_values:
                            smallest_non_zero = min(non_zero_values)
                            f.write(f"  Minimum TPM: 0.0 (exactly zero)\n")
                            f.write(f"  Smallest non-zero TPM: {smallest_non_zero:.20f}\n")
                            # Count how many zeros
                            zero_count = sum(1 for v in tpm_values if v == 0.0)
                            f.write(f"  Number of circRNAs with exactly zero TPM: {zero_count}\n")
                        else:
                            f.write(f"  Minimum TPM: 0.0 (exactly zero for all values)\n")
                    else:
                        # For non-zero minimums, show full precision
                        f.write(f"  Minimum TPM: {min_tpm:.20f}\n")
                    
                    f.write("\n")
                else:
                    f.write(f"{tool_display_names.get(tool, tool)}: No data\n\n")
        
        
        with open(f"{output_dir}/expression_summary.txt", 'w') as f:
            f.write("CircRNA Expression Level Summary by Tool\n")
            f.write("=======================================\n\n")
            
            for idx, row in formatted_summary_df.iterrows():
                f.write(f"{row['Tool']} Statistics:\n")
                f.write(f"  Number of circRNAs detected: {row['Count']}\n")
                if row['Count'] > 0:
                    f.write(f"  Minimum expression level: {row['Min_TPM']:.4f} TPM\n")
                    f.write(f"  Maximum expression level: {row['Max_TPM']:.4f} TPM\n")
                    f.write(f"  Median expression level: {row['Median_TPM']:.4f} TPM\n")
                    f.write(f"  Mean expression level: {row['Mean_TPM']:.4f} TPM\n")
                else:
                    f.write("  No data available\n")
                f.write("\n")
        
        # Create boxplot for TPM values (log scale)
        plt.figure(figsize=(14, 10))
        
        # Use log scale for TPM values (add small value to avoid log(0))
        plot_df['TPM_log'] = np.log10(plot_df['TPM'] + 0.01)
        
        # Create the boxplot
        # Define tool order
        tool_order = ['ground_truth', 'ciri_long', 'isocirc', 'circnick']
        
        # Set the tool color palette
        tool_colors = {
            'ground_truth': '#808080',  # Gray
            'ciri_long': TOOL_COLORS['ciri_long'],
            'isocirc': TOOL_COLORS['isocirc'],
            'circnick': TOOL_COLORS['circnick']
        }
        
        # Create a color list in the right order for the plot
        color_list = [tool_colors[t] for t in tool_order if t in plot_df['Tool'].unique()]
        
        # Generate the boxplot
        ax = sns.boxplot(
            x='Tool', 
            y='TPM_log', 
            data=plot_df,
            order=[t for t in tool_order if t in plot_df['Tool'].unique()],
            palette=color_list,
            width=0.6,
            showmeans=True,
            meanprops={
                "marker":"D", 
                "markerfacecolor":"white", 
                "markeredgecolor":"black", 
                "markersize":"10"
            }
        )
        
        # Add text for median values
        for i, tool in enumerate([t for t in tool_order if t in plot_df['Tool'].unique()]):
            tool_data = plot_df[plot_df['Tool'] == tool]['TPM']
            if not tool_data.empty:
                median_tpm = tool_data.median()
                plt.text(i, plt.ylim()[1] * 0.9, f"Median: {median_tpm:.2f}", 
                        horizontalalignment='center', color='black', fontweight='bold')
        
        # Update x-tick labels
        plt.xticks(
            range(len([t for t in tool_order if t in plot_df['Tool'].unique()])),
            [tool_display_names.get(t, t) for t in tool_order if t in plot_df['Tool'].unique()],
            fontsize=14
        )
        
        plt.title("Distribution of CircRNA Expression Levels (TPM)", fontsize=18)
        plt.xlabel("Tool", fontsize=16)
        plt.ylabel("log10(TPM + 0.01)", fontsize=16)
        plt.yticks(fontsize=14)
        
        # Add grid lines for better readability
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/circrna_expression_by_tool.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Create a violin plot to better show the distribution
        plt.figure(figsize=(14, 10))
        
        ax = sns.violinplot(
            x='Tool', 
            y='TPM_log', 
            data=plot_df,
            order=[t for t in tool_order if t in plot_df['Tool'].unique()],
            palette=color_list,
            inner='box',
            scale='width'
        )
        
        # Add text for median values
        for i, tool in enumerate([t for t in tool_order if t in plot_df['Tool'].unique()]):
            tool_data = plot_df[plot_df['Tool'] == tool]['TPM']
            if not tool_data.empty:
                median_tpm = tool_data.median()
                plt.text(i, plt.ylim()[1] * 0.9, f"Median: {median_tpm:.2f}", 
                        horizontalalignment='center', color='black', fontweight='bold')
        
        # Update x-tick labels
        plt.xticks(
            range(len([t for t in tool_order if t in plot_df['Tool'].unique()])),
            [tool_display_names.get(t, t) for t in tool_order if t in plot_df['Tool'].unique()],
            fontsize=14
        )
        
        plt.title("Distribution of CircRNA Expression Levels (TPM)", fontsize=18)
        plt.xlabel("Tool", fontsize=16)
        plt.ylabel("log10(TPM + 0.01)", fontsize=16)
        plt.yticks(fontsize=14)
        
        # Add grid lines for better readability
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/circrna_expression_violin_by_tool.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        return True
    
    except Exception as e:
        print(f"Error creating expression boxplot: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    # Define paths to input files
    parser = argparse.ArgumentParser(description='Analyze circRNA expression data')
    parser.add_argument('--ground-truth', required=True, 
                       help='Path to ground truth BED file')
    parser.add_argument('--circrna-db', required=True, 
                       help='Path to circRNA database BED file')
    parser.add_argument('--abundances', required=True,
                        help='Path to the abundances.tsv file with TPM values')
    parser.add_argument('--isocirc', required=True, 
                        help='Path to IsoCirc output BED file')
    parser.add_argument('--ciri-long', required=True, 
                        help='Path to CIRI-long output BED12 file')
    parser.add_argument('--circnick', required=True, 
                        help='Path to CircNick output BED12 file')
    parser.add_argument('--output-dir', default='circrna_analysis_results', 
                        help='Output directory')
    parser.add_argument('--fastq', 
                        help='Path to the FASTQ file with reads', 
                        default='/scratch/aerusakovich/sim_ciri_long_jobim/combined/pooled/combined_reads.fastq')
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create debug directory
    debug_dir = os.path.join(args.output_dir, "debug")
    os.makedirs(debug_dir, exist_ok=True)
    
    # Define tool files
    tool_files = {
        'ciri_long': args.ciri_long,
        'isocirc': args.isocirc,
        'circnick': args.circnick
    }
    
    # Create expression plots and analyze zero TPM circRNAs
    create_expression_boxplot(
        args.abundances, 
        args.circrna_db, 
        tool_files, 
        args.ground_truth,  # Pass ground truth BED file
        args.fastq,         # Pass FASTQ file
        args.output_dir
    )
    
    print(f"Analysis complete. Results saved to {args.output_dir}/")

if __name__ == "__main__":
    main()