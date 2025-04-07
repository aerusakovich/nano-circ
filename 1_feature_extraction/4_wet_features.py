#!/usr/bin/env python3
"""
Feature Extraction Script for Nanopore Sequencing Data
Extracts key parameters from real sequencing data for use in simulation:
- Read length distributions
- Rolling circle features (period sizes, copy numbers)
- Tandem repeat characteristics

This script focuses specifically on extracting features from real data
that can be used to parameterize in silico read simulators.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from collections import defaultdict
import math


def set_style():
    """Set consistent styling for all visualizations"""
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams['legend.fontsize'] = 14
    
    # Better contrast and readability
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['grid.linewidth'] = 0.8
    plt.rcParams['lines.linewidth'] = 2.0
    
    # Use a slightly off-white background
    sns.set_style("whitegrid", {'axes.facecolor': '#F9F9F9'})


def create_histogram(data, bins=None, binwidth=None, title=None, xlabel=None, 
                     xlim=None, ylim=None, color='#1f78b4', filename=None,
                     log_scale=False):
    """Create a histogram with consistent styling"""
    plt.figure(figsize=(10, 6))
    
    # Create bins if binwidth specified
    if binwidth is not None:
        if xlim is None:
            min_val = min(data)
            max_val = max(data)
        else:
            min_val, max_val = xlim
        bins = np.arange(min_val, max_val + binwidth, binwidth)
    
    # Create the histogram
    plt.hist(data, bins=bins, color=color, edgecolor='black', linewidth=1, alpha=0.8)
    
    # Add title and labels
    if title:
        plt.title(title, fontsize=18, pad=20, fontweight='bold')
    if xlabel:
        plt.xlabel(xlabel, fontsize=16, labelpad=15)
    plt.ylabel('Count', fontsize=16, labelpad=15)
    
    # Set axis limits if provided
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    
    # Set log scale if requested (useful for nanopore read lengths)
    if log_scale:
        plt.xscale('log')
    
    # Add grid for better readability
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Save the figure
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Created histogram: {filename}")
    
    plt.close()
    return filename


def create_boxplot(data, categories, title=None, xlabel=None, ylabel=None, 
                   color='#1f78b4', filename=None):
    """Create a boxplot with consistent styling"""
    plt.figure(figsize=(12, 8))
    
    # Create boxplot
    boxprops = dict(linewidth=2)
    whiskerprops = dict(linewidth=2)
    capprops = dict(linewidth=2)
    medianprops = dict(linewidth=2, color='black')
    meanprops = dict(marker='D', markerfacecolor='white', markeredgecolor='black', markersize=10)
    
    # Create a list of data for each category
    boxplot_data = [data[data['category'] == cat]['value'].values for cat in categories]
    
    # Create boxplot
    bp = plt.boxplot(boxplot_data, patch_artist=True, notch=True, showfliers=True, 
                     meanprops=meanprops, showmeans=True, boxprops=boxprops, 
                     whiskerprops=whiskerprops, capprops=capprops, medianprops=medianprops)
    
    # Color the boxes
    for box in bp['boxes']:
        box.set(facecolor=color, alpha=0.7)
    
    # Set x-axis labels
    plt.xticks(range(1, len(categories) + 1), categories, rotation=45, ha='right')
    
    # Add title and labels
    if title:
        plt.title(title, fontsize=18, pad=20, fontweight='bold')
    if xlabel:
        plt.xlabel(xlabel, fontsize=16, labelpad=15)
    if ylabel:
        plt.ylabel(ylabel, fontsize=16, labelpad=15)
    
    # Add grid for better readability
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Save the figure
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Created boxplot: {filename}")
    
    plt.close()
    return filename


def create_stacked_barplot(df, x_col, y_col, stack_col, title=None, xlabel=None, ylabel=None, 
                          filename=None, color_palette='viridis'):
    """Create a stacked barplot for repeat length vs copy number distribution"""
    plt.figure(figsize=(14, 8))
    
    # Create pivot table for stacked barplot
    pivot_df = df.pivot_table(values=y_col, index=x_col, columns=stack_col, aggfunc='count', fill_value=0)
    
    # Get color palette
    colors = plt.cm.get_cmap(color_palette, len(pivot_df.columns))
    
    # Create stacked barplot
    pivot_df.plot(kind='bar', stacked=True, ax=plt.gca(), colormap=color_palette, 
                  edgecolor='black', linewidth=0.8)
    
    # Add title and labels
    if title:
        plt.title(title, fontsize=18, pad=20, fontweight='bold')
    if xlabel:
        plt.xlabel(xlabel, fontsize=16, labelpad=15)
    if ylabel:
        plt.ylabel(ylabel, fontsize=16, labelpad=15)
    
    # Improve legend
    plt.legend(title='Copy Number', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
    
    # Add grid for better readability
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Save the figure
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Created stacked barplot: {filename}")
    
    plt.close()
    return filename


def calculate_read_statistics(bed_file):
    """Calculate statistics from a BED file"""
    try:
        df = pd.read_csv(bed_file, sep='\t', header=None)
        df.columns = ['chr', 'start', 'end'] + [f'col{i+4}' for i in range(len(df.columns)-3)]
        
        # Calculate lengths
        df['length'] = df['end'] - df['start']
        
        stats = {
            'count': len(df),
            'mean_length': df['length'].mean(),
            'median_length': df['length'].median(),
            'min_length': df['length'].min(),
            'max_length': df['length'].max(),
            'lengths': df['length'].values
        }
        
        return stats
    except Exception as e:
        print(f"Error processing BED file {bed_file}: {e}")
        return None


def parse_trf_output(trf_file):
    """Parse Tandem Repeats Finder output file
    
    Extracts information about tandem repeats, including:
    - Period sizes (k-mer lengths)
    - Number of copies (possible rolling circles)
    - Repeat lengths
    - Match percentages
    """
    try:
        data = []
        current_read = None
        
        with open(trf_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('@'):
                    # This is a read header
                    read_info = line.split()
                    current_read = read_info[0][1:]  # Remove @ from read ID
                elif current_read and not line.startswith('Sequence:'):
                    # This is a repeat data line
                    try:
                        fields = line.split()
                        if len(fields) >= 15:  # Make sure we have the minimum required fields
                            repeat_data = {
                                'ReadID': current_read,
                                'Start': int(fields[0]),
                                'End': int(fields[1]),
                                'PeriodSize': int(fields[2]),
                                'NumCopies': float(fields[3]),
                                'ConsensusSize': int(fields[4]),
                                'MatchPercent': float(fields[5]),
                                'IndelPercent': float(fields[6]),
                                'Score': float(fields[7]),
                                'PercA': float(fields[8]),
                                'PercC': float(fields[9]),
                                'PercG': float(fields[10]),
                                'PercT': float(fields[11]),
                                'Entropy': float(fields[12]),
                                'Consensus': fields[13],
                                'Repeat': fields[14] if len(fields) > 14 else ''
                            }
                            
                            # Calculate repeat length
                            repeat_data['RepeatLength'] = repeat_data['End'] - repeat_data['Start'] + 1
                            
                            data.append(repeat_data)
                    except (ValueError, IndexError) as e:
                        # Skip malformed lines
                        continue
        
        # Convert to DataFrame
        df = pd.DataFrame(data)
        
        # Calculate summary statistics per read
        if len(df) > 0:
            read_summary = df.groupby('ReadID').agg({
                'RepeatLength': ['mean', 'sum', 'count'],
                'PeriodSize': ['mean', 'median', 'max'],
                'NumCopies': ['mean', 'max']
            })
            
            # Flatten multi-index columns
            read_summary.columns = ['_'.join(col).strip() for col in read_summary.columns.values]
            
            # Rename for clarity
            read_summary = read_summary.rename(columns={
                'RepeatLength_count': 'TotalRepeats',
                'RepeatLength_mean': 'AvgRepeatLength',
                'RepeatLength_sum': 'TotalRepeatLength',
                'PeriodSize_mean': 'AvgPeriodSize',
                'PeriodSize_median': 'MedianPeriodSize',
                'PeriodSize_max': 'MaxPeriodSize',
                'NumCopies_mean': 'AvgNumCopies',
                'NumCopies_max': 'MaxNumCopies'
            })
            
            return {
                'repeat_data': df,
                'read_summary': read_summary
            }
        else:
            print(f"Warning: No valid repeat data found in {trf_file}")
            return None
            
    except Exception as e:
        print(f"Error processing TRF file {trf_file}: {e}")
        return None


def analyze_read_lengths(real_bed_file, output_dir):
    """Analyze read length distributions from real data"""
    print("\nAnalyzing Read Length Distributions...")
    
    # Calculate statistics
    real_stats = calculate_read_statistics(real_bed_file)
    
    if not real_stats:
        print("Error: Could not process real data file")
        return
    
    # Print summary statistics
    print(f"Real Dataset: {real_stats['count']} reads")
    print(f"  Mean Length: {real_stats['mean_length']:.1f} bp")
    print(f"  Median Length: {real_stats['median_length']:.1f} bp")
    print(f"  Min Length: {real_stats['min_length']} bp")
    print(f"  Max Length: {real_stats['max_length']} bp")
    
    # Create histograms - both linear and log scale
    create_histogram(
        data=real_stats['lengths'],
        binwidth=100,
        title='Read Length Distribution',
        xlabel='Read Length (bp)',
        xlim=(0, min(50000, real_stats['max_length'])),
        color='navy',
        filename=os.path.join(output_dir, 'read_length_distribution.png')
    )
    
    create_histogram(
        data=real_stats['lengths'],
        bins=np.logspace(np.log10(max(1, real_stats['min_length'])), 
                         np.log10(real_stats['max_length']), 50),
        title='Read Length Distribution (Log Scale)',
        xlabel='Read Length (bp)',
        color='navy',
        log_scale=True,
        filename=os.path.join(output_dir, 'read_length_distribution_log.png')
    )
    
    return real_stats


def analyze_rolling_circles(trf_file, output_dir):
    """Analyze rolling circle features from TRF output"""
    print("\nAnalyzing Rolling Circle Features...")
    
    # Parse TRF file
    trf_data = parse_trf_output(trf_file)
    
    if not trf_data:
        print("Error: Could not process TRF file")
        return
    
    # Print summary statistics
    print(f"Dataset: {len(trf_data['repeat_data'])} repeats in {len(trf_data['read_summary'])} reads")
    print(f"  Mean Repeats per Read: {trf_data['read_summary']['TotalRepeats'].mean():.1f}")
    print(f"  Mean Repeat Length: {trf_data['repeat_data']['RepeatLength'].mean():.1f} bp")
    print(f"  Mean Period Size: {trf_data['repeat_data']['PeriodSize'].mean():.1f} bp")
    print(f"  Mean Copy Number: {trf_data['repeat_data']['NumCopies'].mean():.1f}")
    
    # Create histograms for key metrics
    
    # Number of copies per repeat (rolling circle indicator)
    create_histogram(
        data=trf_data['repeat_data']['NumCopies'],
        bins=np.arange(1, min(20, trf_data['repeat_data']['NumCopies'].max()) + 0.5, 0.5),
        title='Number of Copies per Repeat',
        xlabel='Number of Copies',
        color='forestgreen',
        filename=os.path.join(output_dir, 'num_copies_distribution.png')
    )
    
    # Period size distribution (k-mer sizes)
    create_histogram(
        data=trf_data['repeat_data']['PeriodSize'],
        bins=np.arange(1, min(100, trf_data['repeat_data']['PeriodSize'].max()) + 5, 5),
        title='Period Size Distribution',
        xlabel='Period Size (bp)',
        color='forestgreen',
        filename=os.path.join(output_dir, 'period_size_distribution.png')
    )
    
    # Total repeats per read
    create_histogram(
        data=trf_data['read_summary']['TotalRepeats'],
        bins=np.arange(1, min(50, trf_data['read_summary']['TotalRepeats'].max()) + 1, 1),
        title='Repeats per Read Distribution',
        xlabel='Number of Repeats',
        color='forestgreen',
        filename=os.path.join(output_dir, 'repeats_per_read_distribution.png')
    )
    
    # Create length vs copy number analysis
    analyze_length_vs_copies(trf_data['repeat_data'], output_dir)
    
    return trf_data


def analyze_length_vs_copies(repeat_data, output_dir):
    """Analyze the relationship between repeat length and copy number (rolling circles)
    with fixed binning to ensure 1-5 repeats are properly displayed"""
    print("\nAnalyzing Repeat Length vs Copy Number Relationship...")
    
    # Create length categories
    length_bins = [0, 50, 100, 300, 600, 1000, 1500, 4000, float('inf')]
    length_labels = ['1-50', '51-100', '101-300', '301-600', '601-1000', '1001-1500', '1501-4000', '4001+']
    
    # Create copy number categories - FIXED to properly include 1-5 range
    # Starting at 1 since copy numbers less than 1 are biologically impossible
    copy_bins = [1, 6, 12, 17, 21, 26, 31, float('inf')]
    copy_labels = ['1-5', '6-11', '12-16', '17-20', '21-25', '26-30', '31+']
    
    # Create a copy of the data for analysis
    df = repeat_data.copy()
    
    # Add length category - using right=False to ensure proper binning
    df['LengthCategory'] = pd.cut(df['RepeatLength'], bins=length_bins, labels=length_labels, right=False)
    
    # Add copy number category - using right=False to ensure proper binning
    df['CopyCategory'] = pd.cut(df['NumCopies'], bins=copy_bins, labels=copy_labels, right=False)
    
    # Create stacked barplot of length vs copy number
    create_stacked_barplot(
        df=df,
        x_col='LengthCategory',
        y_col='RepeatLength',
        stack_col='CopyCategory',
        title='Repeat Length vs Copy Number Distribution',
        xlabel='Repeat Length (bp)',
        ylabel='Count',
        filename=os.path.join(output_dir, 'length_vs_copies_stacked_fixed.png'),
        color_palette='viridis'
    )
    
    # Create boxplot data
    boxplot_data = []
    for length_cat in length_labels:
        for copy_cat in copy_labels:
            subset = df[(df['LengthCategory'] == length_cat) & (df['CopyCategory'] == copy_cat)]
            if len(subset) > 0:
                for _, row in subset.iterrows():
                    boxplot_data.append({
                        'length_category': length_cat,
                        'copy_category': copy_cat,
                        'value': row['NumCopies']
                    })
    
    boxplot_df = pd.DataFrame(boxplot_data)
    
    if len(boxplot_df) > 0:
        # Create boxplot of copy numbers by length category
        create_boxplot(
            data=pd.DataFrame({
                'category': boxplot_df['length_category'],
                'value': boxplot_df['value']
            }),
            categories=length_labels,
            title='Copy Number Distribution by Repeat Length',
            xlabel='Repeat Length (bp)',
            ylabel='Number of Copies',
            color='forestgreen',
            filename=os.path.join(output_dir, 'length_vs_copies_boxplot_fixed.png')
        )
    
    # Create heatmap with fixed categories
    # First ensure we have all combinations present in the pivot table
    all_combinations = [(length, copy) for length in length_labels for copy in copy_labels]
    empty_df = pd.DataFrame(all_combinations, columns=['LengthCategory', 'CopyCategory'])
    
    # Merge with actual data counts
    pivot_table = pd.crosstab(df['LengthCategory'], df['CopyCategory'])
    
    # Make sure all categories are represented (even if zero)
    for length in length_labels:
        if length not in pivot_table.index:
            pivot_table.loc[length] = 0
    
    for copy in copy_labels:
        if copy not in pivot_table.columns:
            pivot_table[copy] = 0
    
    # Sort indices to match original order
    pivot_table = pivot_table.reindex(index=length_labels, columns=copy_labels, fill_value=0)
    
    # Create the heatmap
    plt.figure(figsize=(14, 10))
    sns.heatmap(pivot_table, annot=True, cmap='YlGnBu', fmt='d', linewidths=.5)
    plt.title('Repeat Length vs Copy Number Heatmap', fontsize=18, pad=20, fontweight='bold')
    plt.xlabel('Number of Copies', fontsize=16, labelpad=15)
    plt.ylabel('Repeat Length (bp)', fontsize=16, labelpad=15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'length_vs_copies_heatmap_fixed.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Created fixed heatmap: {os.path.join(output_dir, 'length_vs_copies_heatmap_fixed.png')}")

    # Create circRNA-specific analysis for the 100-4000 bp range (typical circRNA size)
    circrna_df = df[(df['RepeatLength'] >= 100) & (df['RepeatLength'] <= 4000)]

    if len(circrna_df) > 0:
        print(f"\nPotential circRNA Analysis (100-4000 bp):")
        print(f"  Total potential circRNAs: {len(circrna_df)}")
        print(f"  Mean copy number: {circrna_df['NumCopies'].mean():.2f}")
        print(f"  Max copy number: {circrna_df['NumCopies'].max():.2f}")
        
        # Create histogram for circRNA copy numbers
        create_histogram(
            data=circrna_df['NumCopies'],
            bins=np.arange(1, min(30, circrna_df['NumCopies'].max()) + 1, 1),
            title='Copy Number Distribution for Potential circRNAs (100-4000 bp)',
            xlabel='Number of Copies',
            color='darkred',
            filename=os.path.join(output_dir, 'circrna_copies_distribution_fixed.png')
        )
        
        # Create length distribution for potential circRNAs
        create_histogram(
            data=circrna_df['RepeatLength'],
            binwidth=100,
            title='Length Distribution for Potential circRNAs (100-4000 bp)',
            xlabel='Repeat Length (bp)',
            color='darkred',
            filename=os.path.join(output_dir, 'circrna_length_distribution_fixed.png')
        )
        
        # Create scatter plot of length vs copies for circRNAs
        plt.figure(figsize=(12, 8))
        plt.scatter(circrna_df['RepeatLength'], circrna_df['NumCopies'], 
                alpha=0.7, c='darkred', edgecolor='black', s=50)
        plt.title('Repeat Length vs Copy Number for Potential circRNAs', fontsize=18, pad=20, fontweight='bold')
        plt.xlabel('Repeat Length (bp)', fontsize=16, labelpad=15)
        plt.ylabel('Number of Copies', fontsize=16, labelpad=15)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'circrna_length_vs_copies_fixed.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Created scatter plot: {os.path.join(output_dir, 'circrna_length_vs_copies_fixed.png')}")
        
    return pivot_table

def export_parameters(analysis_results, output_dir):
    print("\nExporting parameters for simulation...")
    
    params_file = os.path.join(output_dir, 'simulation_parameters.txt')
    
    with open(params_file, 'w') as f:
        f.write("# Nanopore Sequencing Simulation Parameters\n")
        f.write(f"# Generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Write read length parameters
        if 'read_lengths' in analysis_results:
            read_lengths = analysis_results['read_lengths']
            f.write("## Read Length Parameters\n")
            f.write(f"mean_read_length = {read_lengths['mean_length']:.2f}\n")
            f.write(f"median_read_length = {read_lengths['median_length']:.2f}\n")
            f.write(f"min_read_length = {read_lengths['min_length']}\n")
            f.write(f"max_read_length = {read_lengths['max_length']}\n\n")
        
        # Write rolling circle parameters
        if 'rolling_circles' in analysis_results:
            rc = analysis_results['rolling_circles']
            f.write("## Rolling Circle Parameters\n")
            f.write(f"mean_repeats_per_read = {rc['read_summary']['TotalRepeats'].mean():.2f}\n")
            f.write(f"mean_repeat_length = {rc['repeat_data']['RepeatLength'].mean():.2f}\n")
            f.write(f"mean_period_size = {rc['repeat_data']['PeriodSize'].mean():.2f}\n")
            f.write(f"mean_copy_number = {rc['repeat_data']['NumCopies'].mean():.2f}\n")
            f.write(f"max_copy_number = {rc['repeat_data']['NumCopies'].max():.2f}\n\n")
            
            # Add circRNA-specific parameters
            circrna_df = rc['repeat_data'][(rc['repeat_data']['RepeatLength'] >= 100) & 
                                          (rc['repeat_data']['RepeatLength'] <= 4000)]
            if len(circrna_df) > 0:
                f.write("## circRNA-specific Parameters (100-4000 bp)\n")
                f.write(f"circrna_count = {len(circrna_df)}\n")
                f.write(f"mean_circrna_length = {circrna_df['RepeatLength'].mean():.2f}\n")
                f.write(f"mean_circrna_copies = {circrna_df['NumCopies'].mean():.2f}\n")
                f.write(f"max_circrna_copies = {circrna_df['NumCopies'].max():.2f}\n")
                
                # Add distribution of copy numbers for simulation
                copy_counts = circrna_df['NumCopies'].value_counts().sort_index()
                f.write("\n# Copy number distribution for circRNAs (for simulation)\n")
                for copies, count in copy_counts.items():
                    f.write(f"circrna_copies_{int(copies) if copies.is_integer() else copies} = {count}\n")
    
    print(f"Parameters exported to {params_file}")
    return params_file
def main(args):
    """Main function to run the feature extraction pipeline"""
    print("=" * 80)
    print("Nanopore Sequencing Feature Extraction")
    print("=" * 80)

    # Set up output directory
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Set consistent styling for all plots
    set_style()

    # Store all analysis results
    analysis_results = {}

    # Analyze read length distributions
    if args.bed_file:
        read_length_results = analyze_read_lengths(
            args.bed_file,
            output_dir
        )
        
        if read_length_results:
            analysis_results['read_lengths'] = read_length_results

    # Analyze rolling circle features
    if args.trf_file:
        rolling_circle_results = analyze_rolling_circles(
            args.trf_file,
            output_dir
        )
        
        if rolling_circle_results:
            analysis_results['rolling_circles'] = rolling_circle_results

    # Export parameters for simulation
    if analysis_results:
        export_parameters(analysis_results, output_dir)

    print("\nCompleted feature extraction from real data.")
    print(f"All results saved to {output_dir}")
    print("=" * 80)
if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Extract features from real nanopore data for simulation')
    parser.add_argument('--bed_file', help='BED file with mapped reads for length analysis')
    parser.add_argument('--trf_file', help='Tandem Repeats Finder output file for rolling circle analysis')
    parser.add_argument('--output_dir', default='./feature_extraction_results',
                        help='Directory to save results (default: ./feature_extraction_results)')
    args = parser.parse_args()
    
    if not (args.bed_file or args.trf_file):
        parser.error("At least one of --bed_file or --trf_file must be provided")

    main(args)