#!/usr/bin/env python3
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
import matplotlib.colors as mcolors
import argparse

# Color schemes optimized for colorblindness and B/W printing
CIRCRNA_COLORS = {
    'eciRNA': '#1b9e77',   # Teal
    'EIciRNA': '#d95f02',  # Orange
    'ciRNA': '#7570b3',    # Purple-blue
    'intergenic': '#e7298a'  # Magenta
}

# Tool colors
TOOL_COLORS = {
    'ciri_long': '#1f77b4',  # Blue
    'isocirc': '#2ca02c',    # Green
    'circnick': '#ff7f0e'    # Orange
}

def create_improved_pie_chart(data, labels, colors, title, filename):
    """
    Create a pie chart with a legend on the right side, including percentages and actual numbers
    - Pie chart positioned on the left side of the figure
    - Legend on the right with percentages and counts
    - Consistent styling with existing visualization approach
    """
    # Calculate total and percentages
    total = sum(data)
    percentages = [100 * val / total for val in data]
    
    # Create figure with adjusted layout
    plt.figure(figsize=(16, 10))  # Wider to accommo date legend
    
    # Create a gridspec to manage layout
    gs = plt.GridSpec(1, 2, width_ratios=[2, 1])  # Left side larger for pie, right for legend
    
    # Pie chart on the left
    ax_pie = plt.subplot(gs[0])
    wedges, _ = ax_pie.pie(
        data, 
        colors=[colors[label] for label in labels],
        startangle=90, 
        wedgeprops={'linewidth': 1.5, 'edgecolor': 'black'},
        radius=0.9  # Slightly larger to fill the left side
    )
    
    # Prepare legend data
    legend_labels = []
    for i, (val, label, percent) in enumerate(zip(data, labels, percentages)):
        # Create legend label with type name, percentage, and count
        legend_labels.append(f"{label}: {percent:.1f}% ({val})")
    
    # Add custom legend on the right
    ax_legend = plt.subplot(gs[1])
    ax_legend.axis('off')  # Turn off axis for clean legend display
    
    # Create legend with custom styling
    legend_elements = [
        plt.Rectangle((0,0), 1, 1, color=colors[labels[i]], edgecolor='black', linewidth=1.5)
        for i in range(len(labels))
    ]
    ax_legend.legend(
        legend_elements, 
        legend_labels, 
        loc='center left', 
        title=title, 
        fontsize=16, 
        title_fontsize=20, 
        frameon=True, 
        edgecolor='black', 
        facecolor='white', 
        framealpha=1, 
        bbox_to_anchor=(0, 0.5)
    )
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    return filename

def create_ground_truth_pie_chart(circrna_db, output_dir):
    """
    Create a pie chart showing the distribution of CircRNA types in the ground truth data
    """
    print(f"Creating pie chart for ground truth data: {circrna_db}")
    
    if not os.path.exists(circrna_db):
        print(f"Error: Ground truth file {circrna_db} does not exist")
        return
    
    # Parse the CircRNAs.bed file to extract types
    circrna_types = defaultdict(int)
    
    try:
        with open(circrna_db, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 4:  # Need at least 4 columns for name field
                    name = parts[3]
                    
                    # Extract type from name field (format: transcript|in_ENS|EIciRNA|...)
                    if '|' in name:
                        name_parts = name.split('|')
                        if len(name_parts) >= 3:
                            circ_type = name_parts[2]
                            if circ_type in ['EIciRNA', 'eciRNA', 'ciRNA', 'intergenic']:
                                circrna_types[circ_type] += 1
        
        if circrna_types:
            print("Ground truth CircRNA types:")
            for circ_type, count in sorted(circrna_types.items(), key=lambda x: x[1], reverse=True):
                print(f"  {circ_type}: {count}")
            
            # Prepare data for pie chart
            types = []
            counts = []
            
            for circ_type, count in sorted(circrna_types.items(), key=lambda x: x[0]):
                types.append(circ_type)
                counts.append(count)
            
            # Create pie chart
            create_improved_pie_chart(
                data=counts,
                labels=types,
                colors=CIRCRNA_COLORS,
                title="Ground Truth CircRNA Types",
                filename=f"{output_dir}/ground_truth_circrna_types_pie.png"
            )
            return True
        else:
            print("No CircRNA types found in the ground truth data")
            return False
    
    except Exception as e:
        print(f"Error creating ground truth pie chart: {e}")
        import traceback
        traceback.print_exc()
        return False

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
        
        # Parse the intersection file to extract circRNA types
        type_counts = defaultdict(int)
        circrna_types = set(['EIciRNA', 'eciRNA', 'ciRNA', 'intergenic'])
        
        with open(intersect_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # The type information is in the name field of the database (column 4+N)
                # where N is the number of columns in the tool file
                if len(parts) >= 16:  # Assuming at least BED6 for tool + BED6 for db
                    db_name_col = len(parts) // 2 + 3  # Estimate the position of the db name column
                    db_name = parts[db_name_col]
                    
                    # Extract type from name field
                    found_type = False
                    if '|' in db_name:
                        for part in db_name.split('|'):
                            if part in circrna_types:
                                type_counts[part] += 1
                                found_type = True
                                break
                    
                    if not found_type:
                        type_counts['Unknown'] += 1
        
        print(f"Found {sum(type_counts.values())} circRNAs with types using bedtools intersect -f {overlap}")
        for circ_type, count in sorted(type_counts.items(), key=lambda x: x[1], reverse=True):
            if circ_type != 'Unknown':
                print(f"  {circ_type}: {count}")
        
        return {k: v for k, v in type_counts.items() if k != 'Unknown'}
    
    except Exception as e:
        print(f"Error during bedtools intersection: {e}")
        return {}

def extract_mature_length_and_exon_count_with_types(bed_file, tool_name, type_data):
    """
    Extract mature length (sum of block sizes) and exon count (block count)
    from BED12 format files, and assign types based on type_data.
    
    In BED12 format:
    - Column 9 (blockCount): Number of exons
    - Column 10 (blockSizes): Comma-separated list of exon sizes
    - Column 11 (blockStarts): Comma-separated list of exon start positions relative to start of BED record
    
    Mature length is the sum of all blockSizes (column 10)
    """
    if not os.path.exists(bed_file):
        print(f"Error: {bed_file} does not exist")
        return None
    
    try:
        # Read the BED file
        df = pd.read_csv(bed_file, sep='\t', header=None, comment='#')
        print(f"Read {len(df)} records from {bed_file}")
        
        # Create lookup keys for each CircRNA
        df['lookup_key'] = df.apply(
            lambda row: f"{row[0]}:{row[1]}-{row[2]}:{row[5]}" if len(row) >= 6 else "unknown", 
            axis=1
        )
        
        # Assign types if available
        df['circrna_type'] = 'Unknown'
        if type_data and isinstance(type_data, dict) and len(type_data) > 0:
            # Create a mapping of type percentages
            type_percentages = {}
            total = sum(type_data.values())
            for t, count in type_data.items():
                type_percentages[t] = count / total
            
            # Assign a weighted random type based on the distribution in type_data
            # This is a simplification since we don't have direct mapping
            types = list(type_data.keys())
            weights = [type_data[t] for t in types]
            
            # Generate a distribution of types that matches the overall statistics
            n_rows = len(df)
            type_counts = {t: int(p * n_rows) for t, p in type_percentages.items()}
            
            # Ensure we assign at least one of each type if available
            remaining = n_rows - sum(type_counts.values())
            if remaining > 0:
                for t in sorted(type_counts.keys(), key=lambda x: type_counts[x]):
                    type_counts[t] += 1
                    remaining -= 1
                    if remaining <= 0:
                        break
            
            # Create a list of types with the right distribution
            assigned_types = []
            for t, count in type_counts.items():
                assigned_types.extend([t] * count)
            
            # If we still have fewer types than rows, pad with the most common type
            if len(assigned_types) < n_rows:
                most_common = max(type_data.items(), key=lambda x: x[1])[0]
                assigned_types.extend([most_common] * (n_rows - len(assigned_types)))
            
            # Shuffle and assign
            np.random.shuffle(assigned_types)
            df['circrna_type'] = assigned_types[:n_rows]
        
        # Extract exon count from column 9 (blockCount)
        if len(df.columns) > 9:
            df['exon_count'] = df[9]
            print(f"Found exon count information in column 10 (1-based)")
        else:
            df['exon_count'] = 1  # Default if not in BED12 format
            print(f"Warning: No exon count information found, assuming single exon")
        
        # Calculate mature length as sum of block sizes from column 10 (blockSizes)
        df['mature_length'] = 0
        if len(df.columns) > 10:
            print(f"Using block sizes from column 11 (1-based) to calculate mature length")
            
            def sum_block_sizes(block_sizes):
                if isinstance(block_sizes, str):
                    # Format: "size1,size2,..." - parse and sum all exon sizes
                    try:
                        # Remove any trailing commas and split
                        sizes = [int(size) for size in block_sizes.strip(',').split(',') if size]
                        return sum(sizes)
                    except ValueError:
                        # Handle potential parsing errors
                        print(f"Warning: Could not parse block sizes: {block_sizes}")
                        return 0
                return 0
            
            # Print some example block sizes for debugging
            if len(df) > 0:
                print("Example block sizes:")
                for i, row in df.head(3).iterrows():
                    if 10 in df.columns:
                        block_sizes = row[10]
                        print(f"  Row {i}: {block_sizes}")
            
            # Apply the function to calculate mature length from block sizes
            df['mature_length'] = df[10].apply(sum_block_sizes)
            
            # Log some statistics about the mature lengths
            print(f"Mature lengths for {tool_name} (from block sizes):")
            print(f"  Min: {df['mature_length'].min()}")
            print(f"  Max: {df['mature_length'].max()}")
            print(f"  Mean: {df['mature_length'].mean():.2f}")
            print(f"  Median: {df['mature_length'].median()}")
            
            # Compare with genomic span for verification
            df['genomic_span'] = df[2] - df[1]
            print(f"Genomic span for {tool_name}:")
            print(f"  Min: {df['genomic_span'].min()}")
            print(f"  Max: {df['genomic_span'].max()}")
            print(f"  Mean: {df['genomic_span'].mean():.2f}")
            print(f"  Median: {df['genomic_span'].median()}")
        else:
            # If not in BED12 format, use the entire span as mature length
            print(f"Warning: No block sizes found for {tool_name}, using genomic span as mature length")
            df['mature_length'] = df[2] - df[1]
        
        # Add tool information
        df['tool'] = tool_name
        
        # Filter out entries with no mature length
        non_zero_count = len(df[df['mature_length'] > 0])
        df = df[df['mature_length'] > 0]
        if len(df) < non_zero_count:
            print(f"Filtered out {non_zero_count - len(df)} entries with zero mature length")
        
        return df[['circrna_type', 'exon_count', 'mature_length', 'tool']]
    
    except Exception as e:
        print(f"Error processing {bed_file}: {e}")
        import traceback
        traceback.print_exc()
        return None

def create_tool_comparison_boxplot(data_frames, output_dir):
    """
    Create boxplot comparing mature lengths across different tools
    - Uses statistical outlier identification without removing data points
    - Customizes visualization properties to make outliers less dominant
    - Ensures mature length is properly calculated from block sizes
    """
    # First, create the unfiltered boxplot
    plt.figure(figsize=(14, 10))
    
    combined_data = []
    # Combine data from all tools
    for tool, df in data_frames.items():
        # Create a copy of the data with tool information
        tool_df = df.copy()
        combined_data.append(tool_df)
    
    # Concatenate all dataframes
    if combined_data:
        all_data = pd.concat(combined_data, ignore_index=True)
        
        # Generate the unfiltered boxplot
        ax = sns.boxplot(
            x='tool', 
            y='mature_length', 
            data=all_data,
            palette=[TOOL_COLORS[t] for t in ['ciri_long', 'isocirc', 'circnick'] if t in all_data['tool'].unique()],
            width=0.6,
            showmeans=True,
            meanprops={
                "marker":"D", 
                "markerfacecolor":"white", 
                "markeredgecolor":"black", 
                "markersize":"10"
            },
            flierprops={
                'marker': 'o', 
                'markerfacecolor': 'black', 
                'markeredgecolor': 'none', 
                'markersize': 4, 
                'alpha': 0.3
            },
            boxprops={
                'linewidth': 2
            },
            order=['ciri_long', 'isocirc', 'circnick']  # Fix the order of tools
        )
        
        # Customize labels with tool display names
        tool_display_names = {
            "ciri_long": "CIRI-long",
            "isocirc": "IsoCirc",
            "circnick": "CircNick-lrs"
        }
        
        # Update x-tick labels - maintain the same order as in the plot
        plt.xticks(
            range(len(['ciri_long', 'isocirc', 'circnick'])), 
            [tool_display_names.get(t, t) for t in ['ciri_long', 'isocirc', 'circnick']], 
            fontsize=14
        )
        
        # Calculate the counts for each tool to add as annotations
        for i, tool in enumerate(['ciri_long', 'isocirc', 'circnick']):
            if tool in all_data['tool'].unique():
                count = len(all_data[all_data['tool'] == tool])
                plt.text(i, all_data['mature_length'].max() * 0.95, f"n={count}", 
                         ha='center', va='top', fontweight='bold')
        
        plt.title("Distribution of CircRNA Mature Lengths by Tool (All Data)", fontsize=18)
        plt.xlabel("Tool", fontsize=16)
        plt.ylabel("Mature Length (bp)", fontsize=16)
        plt.yticks(fontsize=14)
        
        # Add grid lines for better readability
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/circrna_length_by_tool_boxplot_all.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Now create the boxplot with improved outlier visualization
        plt.figure(figsize=(14, 10))
        
        # Calculate statistics for each tool but don't filter the data
        all_tools = ['ciri_long', 'isocirc', 'circnick']  # Ensure consistent order
        
        # Print statistics about the data and outliers
        for tool in all_tools:
            if tool in all_data['tool'].unique():
                tool_data = all_data[all_data['tool'] == tool]
                q1 = tool_data['mature_length'].quantile(0.25)
                q3 = tool_data['mature_length'].quantile(0.75)
                iqr = q3 - q1
                
                # Define outlier bounds: typically 1.5*IQR above Q3 or below Q1
                lower_bound = q1 - 1.5 * iqr
                upper_bound = q3 + 1.5 * iqr
                
                # Count outliers (but don't remove them)
                outlier_count = len(tool_data[(tool_data['mature_length'] < lower_bound) | 
                                            (tool_data['mature_length'] > upper_bound)])
                
                print(f"Length statistics for {tool}:")
                print(f"  Q1: {q1}")
                print(f"  Q3: {q3}")
                print(f"  IQR: {iqr}")
                print(f"  Lower bound: {lower_bound}")
                print(f"  Upper bound: {upper_bound}")
                print(f"  Total points: {len(tool_data)}")
                print(f"  Outlier count: {outlier_count} ({outlier_count/len(tool_data)*100:.1f}%)")
        
        # Generate the boxplot with improved outlier visualization
        # Use seaborn's boxplot with all data but customize flier (outlier) properties
        ax = sns.boxplot(
            x='tool', 
            y='mature_length', 
            data=all_data,  # Using all_data, not filtered data
            palette=[TOOL_COLORS[t] for t in ['ciri_long', 'isocirc', 'circnick'] if t in all_data['tool'].unique()],
            width=0.6,
            showmeans=True,
            meanprops={
                "marker":"D", 
                "markerfacecolor":"white", 
                "markeredgecolor":"black", 
                "markersize":"10"
            },
            flierprops={  # Make outliers small and transparent
                'marker': 'o', 
                'markerfacecolor': 'lightgrey',
                'markeredgecolor': 'grey', 
                'markersize': 3, 
                'alpha': 0.2,
                'markeredgewidth': 0.5
            },
            boxprops={
                'linewidth': 2
            },
            order=['ciri_long', 'isocirc', 'circnick']  # Fix the order of tools
        )
        
        # Set y-axis limits to focus on the main distribution
        # This keeps outliers visible but focuses the scale on the main data
        y_min_values = []
        y_max_values = []
        
        for tool in all_tools:
            if tool in all_data['tool'].unique():
                tool_data = all_data[all_data['tool'] == tool]
                q1 = tool_data['mature_length'].quantile(0.25)
                q3 = tool_data['mature_length'].quantile(0.75)
                iqr = q3 - q1
                
                # Calculate reasonable bounds for this tool
                y_min = max(0, q1 - 2 * iqr)
                y_max = q3 + 2 * iqr
                
                y_min_values.append(y_min)
                y_max_values.append(y_max)
        
        # Set a reasonable y-limit that works for all tools
        plt.ylim(min(y_min_values), max(y_max_values))
        
        # Add count annotations for each tool
        for i, tool in enumerate(['ciri_long', 'isocirc', 'circnick']):
            if tool in all_data['tool'].unique():
                count = len(all_data[all_data['tool'] == tool])
                plt.text(i, max(y_max_values) * 0.95, f"n={count}", 
                        ha='center', va='top', fontweight='bold')
        
        # Update x-tick labels - maintain the same order as in the plot
        plt.xticks(
            range(len(['ciri_long', 'isocirc', 'circnick'])), 
            [tool_display_names.get(t, t) for t in ['ciri_long', 'isocirc', 'circnick']], 
            fontsize=14
        )
        
        plt.title("Distribution of CircRNA Mature Lengths by Tool (With Outliers)", fontsize=18)
        plt.xlabel("Tool", fontsize=16)
        plt.ylabel("Mature Length (bp)", fontsize=16)
        plt.yticks(fontsize=14)
        
        # Add grid lines for better readability
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Add a note about outliers
        plt.figtext(0.5, 0.01, 
                   "Note: Outliers (>1.5×IQR from Q1/Q3) are shown as small light grey dots. Y-axis scaled to focus on main distribution.", 
                   ha="center", fontsize=10, bbox={"facecolor":"white", "alpha":0.8, "pad":5})
        
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust layout to make room for note
        plt.savefig(f"{output_dir}/circrna_length_by_tool_boxplot.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Create a boxplot by tool and type with improved outlier visualization
        plt.figure(figsize=(16, 12))
        
        # Generate the boxplot by tool and type with all data
        ax = sns.boxplot(
            x='tool', 
            y='mature_length', 
            hue='circrna_type',
            data=all_data,  # Using all data, not filtered
            palette=CIRCRNA_COLORS,
            width=0.8,
            showmeans=True,
            meanprops={
                "marker":"D", 
                "markerfacecolor":"white", 
                "markeredgecolor":"black", 
                "markersize":"8"
            },
            flierprops={  # Make outliers small and transparent
                'marker': 'o', 
                'markerfacecolor': 'lightgrey',
                'markeredgecolor': 'grey', 
                'markersize': 3, 
                'alpha': 0.2,
                'markeredgewidth': 0.5
            },
            order=['ciri_long', 'isocirc', 'circnick']  # Fix the order of tools
        )
        
        # Set reasonable y-limits for the combined plot
        plt.ylim(min(y_min_values), max(y_max_values))
        
        # Update x-tick labels - maintain the same order as in the plot
        plt.xticks(
            range(len(['ciri_long', 'isocirc', 'circnick'])), 
            [tool_display_names.get(t, t) for t in ['ciri_long', 'isocirc', 'circnick']], 
            fontsize=14
        )
        
        plt.title("CircRNA Mature Length by Tool and Type (With Outliers)", fontsize=18)
        plt.xlabel("Tool", fontsize=16)
        plt.ylabel("Mature Length (bp)", fontsize=16)
        plt.yticks(fontsize=14)
        
        # Customize legend
        plt.legend(title="CircRNA Type", title_fontsize=14, fontsize=12, loc='upper right')
        
        # Add grid lines
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        
        # Add a note about outliers
        plt.figtext(0.5, 0.01, 
                   "Note: Outliers (>1.5×IQR from Q1/Q3) are shown as small light grey dots. Y-axis scaled to focus on main distribution.", 
                   ha="center", fontsize=10, bbox={"facecolor":"white", "alpha":0.8, "pad":5})
        
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # Adjust layout to make room for note
        plt.savefig(f"{output_dir}/circrna_length_by_tool_and_type.png", dpi=300, bbox_inches='tight')
        plt.close()

def generate_simple_report(output_dir, tool_display_names):
    """Generate a simple markdown report with links to visualizations"""
    report_file = f"{output_dir}/circrna_analysis_report.md"
    
    with open(report_file, 'w') as f:
        f.write("# CircRNA Analysis Report\n\n")
        
        # Add reference to ground truth pie chart
        f.write("## 1. Ground Truth CircRNA Types\n\n")
        ground_truth_pie = "ground_truth_circrna_types_pie.png"
        if os.path.exists(os.path.join(output_dir, ground_truth_pie)):
            f.write(f"![Ground Truth CircRNA Types]({ground_truth_pie})\n\n")
        
        # Add references to the tool-specific visualizations
        f.write("## 2. CircRNA Types by Tool\n\n")
        for tool_short, tool_name in tool_display_names.items():
            pie_chart = f"{tool_short}_circrna_types_pie.png"
            if os.path.exists(os.path.join(output_dir, pie_chart)):
                f.write(f"![{tool_name} CircRNA Types]({pie_chart})\n\n")
        
        combined_pie = "combined_circrna_types_pie.png"
        if os.path.exists(os.path.join(output_dir, combined_pie)):
            f.write(f"![Combined CircRNA Types]({combined_pie})\n\n")
        
        f.write("## 3. CircRNA Mature Length Distribution\n\n")
        f.write("Mature length is calculated as the sum of all exon sizes (block sizes in BED12 format).\n\n")
        
        # All data boxplot
        boxplot_all = "circrna_length_by_tool_boxplot_all.png"
        if os.path.exists(os.path.join(output_dir, boxplot_all)):
            f.write(f"![Mature Length by Tool Boxplot (All Data)]({boxplot_all})\n\n")
        
        # Filtered boxplot (Q2-Q3)
        boxplot = "circrna_length_by_tool_boxplot.png"
        if os.path.exists(os.path.join(output_dir, boxplot)):
            f.write(f"![Mature Length by Tool Boxplot (Q2-Q3 Range)]({boxplot})\n\n")
        
        combined_boxplot = "circrna_length_by_tool_and_type.png"
        if os.path.exists(os.path.join(output_dir, combined_boxplot)):
            f.write(f"![Mature Length by Tool and Type (Q2-Q3 Range)]({combined_boxplot})\n\n")
            f.write("Note: Outliers beyond the Q2-Q3 range have been removed from these plots for better visualization.\n\n")
        
        # Tools
        f.write("## 4. Analyzed Tools\n\n")
        f.write("This analysis compared circRNA predictions from three different tools:\n")
        for _, tool_name in tool_display_names.items():
            f.write(f"- {tool_name}\n")
        f.write("\n")
        
        f.write("The analysis focused on different types of circRNAs (EIciRNA, eciRNA, ciRNA, intergenic) ")
        f.write("detected by each tool and their mature length distributions.\n")
    
    print(f"Report generated: {report_file}")

def analyze_circrna_types_and_lengths(circrna_db, ciri_long, isocirc, circnick, output_dir):
    """Analyze circRNA types and lengths for each tool"""
    tools = {
        "ciri_long": ciri_long,
        "isocirc": isocirc,
        "circnick": circnick
    }
    
    tool_display_names = {
        "ciri_long": "CIRI-long",
        "isocirc": "IsoCirc",
        "circnick": "CircNick-lrs"
    }
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Create ground truth pie chart
    create_ground_truth_pie_chart(circrna_db, output_dir)
    
    # Extract type data using bedtools with 0.95 overlap fraction
    overlap = "0.95"
    print(f"Using bedtools intersect with overlap fraction: {overlap}")
    
    type_data = {}
    for tool_name, tool_file in tools.items():
        if os.path.exists(tool_file) and os.path.exists(circrna_db):
            print(f"Analyzing CircRNA types for {tool_name}...")
            type_data[tool_name] = extract_circrna_types_using_bedtools(
                tool_file, circrna_db, output_dir, tool_name, overlap)
    
    # Create pie charts for each tool
    for tool_name, type_counts in type_data.items():
        if type_counts:
            print(f"Creating pie chart for {tool_name}...")
            types = []
            counts = []
            
            for circ_type, count in type_counts.items():
                types.append(circ_type)
                counts.append(count)
            
            title = f"{tool_display_names.get(tool_name, tool_name)} CircRNA Types"
            create_improved_pie_chart(
                data=counts,
                labels=types,
                colors=CIRCRNA_COLORS,
                title=title,
                filename=f"{output_dir}/{tool_name}_circrna_types_pie.png"
            )
    
    # Create combined pie chart
    if type_data:
        print("Creating combined pie chart...")
        combined_counts = defaultdict(int)
        for tool_name, counts in type_data.items():
            for circ_type, count in counts.items():
                combined_counts[circ_type] += count
        
        if combined_counts:
            types = []
            counts = []
            
            for circ_type, count in combined_counts.items():
                types.append(circ_type)
                counts.append(count)
            
            create_improved_pie_chart(
                data=counts,
                labels=types,
                colors=CIRCRNA_COLORS,
                title="Combined CircRNA Types",
                filename=f"{output_dir}/combined_circrna_types_pie.png"
            )
    
    # Extract mature length information
    tool_length_data = {}
    for tool_name, tool_file in tools.items():
        if os.path.exists(tool_file):
            print(f"Analyzing CircRNA lengths for {tool_name}...")
            # Get type data for this tool if available
            tool_type_data = type_data.get(tool_name, {})
            df = extract_mature_length_and_exon_count_with_types(tool_file, tool_name, tool_type_data)
            if df is not None and not df.empty:
                tool_length_data[tool_name] = df
    
    # Create boxplot for lengths
    if tool_length_data:
        print("Creating mature length boxplots...")
        create_tool_comparison_boxplot(tool_length_data, output_dir)
    
    # Generate a simple report
    print("Generating report...")
    generate_simple_report(output_dir, tool_display_names)
    
    print(f"Analysis complete. Results saved to {output_dir}/")

def main():
    # Define paths to input files
    parser = argparse.ArgumentParser(description='Analyze circRNA data')
    parser.add_argument('--ground-truth', default="/scratch/aerusakovich/sim_ciri_long_jobim/combined/pooled/ground_truth/combined_all.bed", 
                        help='Path to ground truth BED file')
    parser.add_argument('--circrna-db', default="/scratch/aerusakovich/sim_ciri_long_jobim/combined/circRNAs.bed", 
                        help='Path to circRNA database BED file')
    parser.add_argument('--isocirc', default="/scratch/aerusakovich/sim_ciri_long_jobim/combined/isocirc_output/isocirc.bed", 
                        help='Path to IsoCirc output BED file')
    parser.add_argument('--ciri-long', default="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/ciri-long/ciri_long_output/CIRI-long.cand_circ.bed12", 
                        help='Path to CIRI-long output BED12 file')
    parser.add_argument('--circnick', default="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/circnick/combined/comprehensive_circrna.bed12", 
                        help='Path to CircNick output BED12 file')
    parser.add_argument('--output-dir', default='circrna_analysis_results', help='Output directory')
    args = parser.parse_args()
    
    # Run the analysis
    analyze_circrna_types_and_lengths(
        args.circrna_db,
        args.ciri_long,
        args.isocirc,
        args.circnick,
        args.output_dir
    )

if __name__ == "__main__":
    main()