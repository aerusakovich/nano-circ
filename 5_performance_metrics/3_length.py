#!/usr/bin/env python3
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse

# Tool colors with ground truth added
TOOL_COLORS = {
    'ciri_long': '#1f77b4',    # Blue
    'isocirc': '#2ca02c',      # Green
    'circnick': '#ff7f0e',     # Orange
    'ground_truth': '#808080'  # Gray
}

def extract_mature_length_from_bed(bed_file, tool_name):
    """
    Extract mature length from BED files.
    For BED12 format: sum of block sizes (column 10)
    For other formats: genomic span (end - start)
    """
    if not os.path.exists(bed_file):
        print(f"Error: {bed_file} does not exist")
        return None
    
    try:
        # Read the BED file
        df = pd.read_csv(bed_file, sep='\t', header=None, comment='#')
        print(f"Read {len(df)} records from {bed_file}")
        print(f"Number of columns: {len(df.columns)}")
        
        # Show first few rows for debugging
        if len(df) > 0:
            print("First few rows:")
            for i, row in df.head(2).iterrows():
                print(f"  Row {i}: {list(row)}")
        
        # Calculate mature length
        if len(df.columns) > 10:
            # BED12 format - sum of block sizes from column 10 (0-indexed: column 10)
            print(f"Using block sizes from column 11 (1-based) to calculate mature length for {tool_name}")
            
            def sum_block_sizes(block_sizes):
                if isinstance(block_sizes, str):
                    try:
                        # Handle both single values and comma-separated values
                        if ',' in block_sizes:
                            # Multi-exon: comma-separated values
                            sizes = [int(size) for size in block_sizes.strip(',').split(',') if size]
                            return sum(sizes)
                        else:
                            # Single exon: just one value
                            return int(block_sizes.strip())
                    except ValueError:
                        print(f"Warning: Could not parse block sizes: {block_sizes}")
                        return 0
                elif isinstance(block_sizes, (int, float)):
                    # Handle numeric values directly
                    return int(block_sizes)
                return 0
            
            # Calculate mature length from block sizes
            df['mature_length'] = df[10].apply(sum_block_sizes)
            
            # Print some statistics
            print(f"Mature lengths for {tool_name} (from block sizes):")
            print(f"  Min: {df['mature_length'].min()}")
            print(f"  Max: {df['mature_length'].max()}")
            print(f"  Mean: {df['mature_length'].mean():.2f}")
            print(f"  Median: {df['mature_length'].median()}")
            
        else:
            # Standard BED format - use genomic span
            print(f"Using genomic span (end - start) as mature length for {tool_name}")
            df['mature_length'] = df[2] - df[1]
        
        # Add tool information
        df['tool'] = tool_name
        
        # Filter out entries with no mature length
        df = df[df['mature_length'] > 0]
        print(f"Retained {len(df)} entries with positive mature length")
        
        return df[['mature_length', 'tool']]
    
    except Exception as e:
        print(f"Error processing {bed_file}: {e}")
        import traceback
        traceback.print_exc()
        return None

def create_length_boxplot(data_frames, output_dir):
    """
    Create boxplot comparing mature lengths across different tools including ground truth
    Applied aesthetic specifications:
    - Title: fontsize=20
    - X-axis label: fontsize=16  
    - Y-axis label: fontsize=16
    - Tick labels (X and Y): fontsize=14
    """
    # Combine data from all tools
    combined_data = []
    for tool, df in data_frames.items():
        combined_data.append(df.copy())
    
    if not combined_data:
        print("No data available for plotting")
        return
    
    # Concatenate all dataframes
    all_data = pd.concat(combined_data, ignore_index=True)
    
    # Define tool order including ground truth
    tool_order = ['ground_truth', 'ciri_long', 'isocirc', 'circnick']
    available_tools = [t for t in tool_order if t in all_data['tool'].unique()]
    
    # Tool display names
    tool_display_names = {
        "ground_truth": "Ground Truth",
        "ciri_long": "CIRI-long",
        "isocirc": "IsoCirc", 
        "circnick": "CircNick-LRS"
    }
    
    # Create the main boxplot
    plt.figure(figsize=(14, 10))
    
    ax = sns.boxplot(
        x='tool', 
        y='mature_length', 
        data=all_data,
        palette=[TOOL_COLORS[t] for t in available_tools],
        width=0.6,
        showmeans=True,
        meanprops={
            "marker": "D", 
            "markerfacecolor": "white", 
            "markeredgecolor": "black", 
            "markersize": 10
        },
        flierprops={
            'marker': 'o', 
            'markerfacecolor': 'lightgrey',
            'markeredgecolor': 'grey', 
            'markersize': 3, 
            'alpha': 0.2,
            'markeredgewidth': 0.5
        },
        boxprops={'linewidth': 2},
        order=available_tools,
        hue='tool',
        legend=False
    )
    
    # Calculate statistics and set reasonable y-limits
    y_min_values = []
    y_max_values = []
    
    for tool in available_tools:
        tool_data = all_data[all_data['tool'] == tool]
        q1 = tool_data['mature_length'].quantile(0.25)
        q3 = tool_data['mature_length'].quantile(0.75)
        iqr = q3 - q1
        
        # Calculate reasonable bounds
        y_min = max(0, q1 - 2 * iqr)
        y_max = q3 + 2 * iqr
        
        y_min_values.append(y_min)
        y_max_values.append(y_max)
        
        # Print statistics
        print(f"Length statistics for {tool}:")
        print(f"  Count: {len(tool_data)}")
        print(f"  Mean: {tool_data['mature_length'].mean():.2f}")
        print(f"  Median: {tool_data['mature_length'].median():.2f}")
        print(f"  Q1: {q1:.2f}")
        print(f"  Q3: {q3:.2f}")
    
    # Set y-axis limits to focus on main distribution
    plt.ylim(min(y_min_values), max(y_max_values))
    
    # Add count annotations for each tool
    for i, tool in enumerate(available_tools):
        count = len(all_data[all_data['tool'] == tool])
        plt.text(i, max(y_max_values) * 0.95, f"n={count}", 
                ha='center', va='top', fontweight='bold', fontsize=16)
    
    # Apply aesthetic specifications
    plt.title("CircRNA Mature Length Distribution by Tool", fontsize=24, fontweight='bold')
    plt.xlabel("Tool", fontsize=20)
    plt.ylabel("Mature Length (bp)", fontsize=20)
    
    # Update x-tick labels with display names
    plt.xticks(
        range(len(available_tools)), 
        [tool_display_names.get(t, t) for t in available_tools], 
        fontsize=18
    )
    plt.yticks(fontsize=18)
    
    # Add grid lines for better readability
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add subtitle with reads count information (italic style)
    plt.figtext(0.5, 0.92, 
               "Analysis based on simulated reads (mature length = sum of exon sizes)", 
               ha="center", fontsize=16, style='italic')
    
    # Add note about outliers
    plt.figtext(0.5, 0.02, 
               "Note: Outliers (>1.5×IQR from Q1/Q3) shown as small light grey dots. Y-axis scaled to focus on main distribution.", 
               ha="center", fontsize=14, 
               bbox={"facecolor": "lightyellow", "alpha": 0.8, "pad": 5})
    
    plt.tight_layout(rect=[0, 0.06, 1, 0.90])  # Adjust layout for subtitle and note
    plt.savefig(f"{output_dir}/circrna_length_by_tool_with_ground_truth.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a version without y-axis limits for full data view
    plt.figure(figsize=(14, 10))
    
    ax = sns.boxplot(
        x='tool', 
        y='mature_length', 
        data=all_data,
        palette=[TOOL_COLORS[t] for t in available_tools],
        width=0.6,
        showmeans=True,
        meanprops={
            "marker": "D", 
            "markerfacecolor": "white", 
            "markeredgecolor": "black", 
            "markersize": 10
        },
        flierprops={
            'marker': 'o', 
            'markerfacecolor': 'black', 
            'markeredgecolor': 'none', 
            'markersize': 4, 
            'alpha': 0.3
        },
        boxprops={'linewidth': 2},
        order=available_tools,
        hue='tool',
        legend=False
    )
    
    # Add count annotations
    for i, tool in enumerate(available_tools):
        count = len(all_data[all_data['tool'] == tool])
        plt.text(i, all_data['mature_length'].max() * 0.95, f"n={count}", 
                ha='center', va='top', fontweight='bold', fontsize=16)
    
    # Apply aesthetic specifications
    plt.title("CircRNA Mature Length Distribution by Tool (Full Range)", fontsize=24, fontweight='bold')
    plt.xlabel("Tool", fontsize=20)
    plt.ylabel("Mature Length (bp)", fontsize=20)
    
    # Update x-tick labels with display names
    plt.xticks(
        range(len(available_tools)), 
        [tool_display_names.get(t, t) for t in available_tools], 
        fontsize=18
    )
    plt.yticks(fontsize=18)
    
    # Add grid lines
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add subtitle
    plt.figtext(0.5, 0.92, 
               "Analysis based on simulated reads (mature length = sum of exon sizes)", 
               ha="center", fontsize=16, style='italic')
    
    plt.tight_layout(rect=[0, 0.02, 1, 0.90])
    plt.savefig(f"{output_dir}/circrna_length_by_tool_with_ground_truth_full.png", dpi=300, bbox_inches='tight')
    plt.close()

def generate_length_report(output_dir, tool_counts):
    """Generate a simple markdown report for length analysis"""
    report_file = f"{output_dir}/circrna_length_analysis_report.md"
    
    with open(report_file, 'w') as f:
        f.write("# CircRNA Mature Length Analysis Report\n\n")
        
        f.write("## Analysis Overview\n\n")
        f.write("This report compares the mature length distributions of circRNAs detected by different tools.\n")
        f.write("Mature length is calculated as:\n")
        f.write("- **BED12 format**: Sum of all exon sizes (block sizes)\n")
        f.write("- **Standard BED format**: Genomic span (end - start)\n\n")
        
        f.write("## Tool Comparison\n\n")
        f.write("| Tool | Count | Color |\n")
        f.write("|------|-------|-------|\n")
        
        tool_display_names = {
            "ground_truth": "Ground Truth",
            "ciri_long": "CIRI-long", 
            "isocirc": "IsoCirc",
            "circnick": "CircNick-LRS"
        }
        
        for tool, count in tool_counts.items():
            display_name = tool_display_names.get(tool, tool)
            color = TOOL_COLORS.get(tool, '#000000')
            f.write(f"| {display_name} | {count} | {color} |\n")
        
        f.write("\n## Visualizations\n\n")
        
        # Main plot
        main_plot = "circrna_length_by_tool_with_ground_truth.png"
        if os.path.exists(os.path.join(output_dir, main_plot)):
            f.write(f"![CircRNA Length Distribution (Focused)]({main_plot})\n\n")
            f.write("*Y-axis scaled to focus on main distribution (outliers visible but compressed)*\n\n")
        
        # Full range plot  
        full_plot = "circrna_length_by_tool_with_ground_truth_full.png"
        if os.path.exists(os.path.join(output_dir, full_plot)):
            f.write(f"![CircRNA Length Distribution (Full Range)]({full_plot})\n\n")
            f.write("*Y-axis shows full data range including all outliers*\n\n")
        
        f.write("## Notes\n\n")
        f.write("- **Diamond markers**: Mean values\n")
        f.write("- **Box plots**: Show median, Q1, Q3, and whiskers (1.5×IQR)\n") 
        f.write("- **Outliers**: Data points beyond 1.5×IQR from Q1/Q3\n")
        f.write("- **Grid lines**: Horizontal guidelines for easier reading\n")
    
    print(f"Length analysis report generated: {report_file}")

def analyze_circrna_lengths(ground_truth, ciri_long, isocirc, circnick, output_dir):
    """Main function to analyze circRNA lengths across all tools including ground truth"""
    
    tools = {
        "ground_truth": ground_truth,
        "ciri_long": ciri_long,
        "isocirc": isocirc,
        "circnick": circnick
    }
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Extract length data for each tool
    tool_length_data = {}
    tool_counts = {}
    
    for tool_name, tool_file in tools.items():
        if os.path.exists(tool_file):
            print(f"\nAnalyzing CircRNA lengths for {tool_name}...")
            df = extract_mature_length_from_bed(tool_file, tool_name)
            if df is not None and not df.empty:
                tool_length_data[tool_name] = df
                tool_counts[tool_name] = len(df)
                print(f"Successfully processed {len(df)} circRNAs for {tool_name}")
            else:
                print(f"No valid data found for {tool_name}")
        else:
            print(f"Warning: File not found for {tool_name}: {tool_file}")
    
    # Create boxplots if we have data
    if tool_length_data:
        print(f"\nCreating length comparison plots with {len(tool_length_data)} tools...")
        create_length_boxplot(tool_length_data, output_dir)
        
        # Generate report
        generate_length_report(output_dir, tool_counts)
        
        print(f"\nLength analysis complete. Results saved to {output_dir}/")
        print("\nGenerated files:")
        print(f"- circrna_length_by_tool_with_ground_truth.png")
        print(f"- circrna_length_by_tool_with_ground_truth_full.png") 
        print(f"- circrna_length_analysis_report.md")
    else:
        print("Error: No valid data found for any tool")

def main():
    parser = argparse.ArgumentParser(description='Analyze circRNA mature lengths across tools')
    parser.add_argument('--ground-truth', 
                        default="/scratch/aerusakovich/sim_ciri_long_jobim/combined/circRNAs.bed", 
                        help='Path to ground truth BED file')
    parser.add_argument('--ciri-long', 
                        default="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/ciri-long/ciri_long_output/CIRI-long.cand_circ.bed12", 
                        help='Path to CIRI-long output BED12 file')
    parser.add_argument('--isocirc', 
                        default="/scratch/aerusakovich/sim_ciri_long_jobim/combined/isocirc_output/isocirc.bed", 
                        help='Path to IsoCirc output BED file')
    parser.add_argument('--circnick', 
                        default="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/circnick/combined/comprehensive_circrna.bed12", 
                        help='Path to CircNick output BED12 file')
    parser.add_argument('--output-dir', 
                        default='circrna_length_analysis', 
                        help='Output directory')
    
    args = parser.parse_args()
    
    # Run the length analysis
    analyze_circrna_lengths(
        args.ground_truth,
        args.ciri_long, 
        args.isocirc,
        args.circnick,
        args.output_dir
    )

if __name__ == "__main__":
    main()