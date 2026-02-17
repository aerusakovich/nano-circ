#!/usr/bin/env python3
"""
Complete visualization script for circRNA data with all improvements:
1. Properly positioned labels on pie charts with slight adjustment for ciRNA
2. Adjusted Y-axis on boxplot to match data range
3. Increased tick marks on Y-axis (numeric format)
4. Consistent styling across all visualizations
"""
import argparse
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter

# Set consistent styling for all plots
def set_style():
    # Set up larger fonts for all plots (increased sizes)
    mpl.rcParams['font.size'] = 18         # Increased from 16
    mpl.rcParams['axes.titlesize'] = 24    # Increased from 20
    mpl.rcParams['axes.labelsize'] = 22    # Increased from 18
    mpl.rcParams['xtick.labelsize'] = 18   # Increased from 16
    mpl.rcParams['ytick.labelsize'] = 18   # Increased from 16
    mpl.rcParams['legend.fontsize'] = 18   # Increased from 16
    mpl.rcParams['figure.titlesize'] = 26  # Increased from 22
    
    # Use a slightly off-white background that's easier on the eyes
    sns.set_style("whitegrid", {'axes.facecolor': '#F9F9F9'})
    
    # Increase line widths for better visibility
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['grid.linewidth'] = 0.8
    mpl.rcParams['lines.linewidth'] = 2.0
    mpl.rcParams['patch.linewidth'] = 1.5
    
    # Increase marker sizes
    mpl.rcParams['lines.markersize'] = 10
    
    # Ensure high-contrast text
    mpl.rcParams['text.color'] = 'black'
    mpl.rcParams['axes.labelcolor'] = 'black'
    mpl.rcParams['xtick.color'] = 'black'
    mpl.rcParams['ytick.color'] = 'black'
    
    # Use a slightly off-white background that's easier on the eyes
    sns.set_style("whitegrid", {'axes.facecolor': '#F9F9F9'})
    
    # Increase line widths for better visibility
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['grid.linewidth'] = 0.8
    mpl.rcParams['lines.linewidth'] = 2.0
    mpl.rcParams['patch.linewidth'] = 1.5
    
    # Increase marker sizes
    mpl.rcParams['lines.markersize'] = 10
    
    # Ensure high-contrast text
    mpl.rcParams['text.color'] = 'black'
    mpl.rcParams['axes.labelcolor'] = 'black'
    mpl.rcParams['xtick.color'] = 'black'
    mpl.rcParams['ytick.color'] = 'black'

# Color schemes optimized for colorblindness and B/W printing
CIRCRNA_COLORS = {
    'eciRNA': '#1b9e77',    # Teal
    'EIciRNA': '#d95f02',   # Orange
    'ciRNA': '#7570b3',     # Purple-blue
    'intergenic': '#e7298a' # Magenta
}

SPLICE_COLORS = {
    'GT-AG': '#66c2a5',      # Teal
    'GC-AG': '#fc8d62',      # Orange
    'AT-AC': '#8da0cb',      # Blue
    'Non-canonical': '#e78ac3', # Pink
    'Unknown': '#a6d854'      # Green
}

EXON_COUNT_COLORS = {
    1: '#8dd3c7',    # Mint
    2: '#ffffb3',    # Light yellow
    3: '#bebada',    # Lavender
    4: '#fb8072',    # Salmon
    '5+': '#80b1d3'  # Light blue
}

# Labels for plot legends
TYPE_LABELS = {
    'eciRNA': 'eciRNA (Exonic)',
    'EIciRNA': 'EIciRNA (Exon-Intron)',
    'ciRNA': 'ciRNA (Intronic)',
    'intergenic': 'Intergenic'
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
    plt.figure(figsize=(16, 10))  # Wider to accommodate legend
    
    # Create a gridspec to manage layout
    gs = plt.GridSpec(1, 2, width_ratios=[2, 1])  # Left side larger for pie, right for legend
    
    # Pie chart on the left
    ax_pie = plt.subplot(gs[0])
    wedges, _ = ax_pie.pie(
        data,
        colors=colors,
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
        plt.Rectangle((0,0), 1, 1, color=colors[i], edgecolor='black', linewidth=1.5) 
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

def create_improved_boxplot(df, output_dir):
    """Create boxplot with Y-axis adjusted to match actual data range"""
    plt.figure(figsize=(14, 10))
    
    # Define order of categories for consistency
    category_order = ['eciRNA', 'EIciRNA', 'ciRNA', 'intergenic']
    
    # Generate boxplot with custom colors
    ax = sns.boxplot(
        x='circrna_type', 
        y='mature_length', 
        data=df,
        order=category_order,
        palette=[CIRCRNA_COLORS[t] for t in category_order if t in CIRCRNA_COLORS],
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
        }
    )
    
    # Replace x-axis labels with more descriptive ones
    x_labels = [TYPE_LABELS.get(cat, cat) for cat in category_order if cat in df['circrna_type'].unique()]
    ax.set_xticklabels(x_labels, fontweight='bold')
    
    # Set labels and title
    plt.xlabel('circRNA Type', fontsize=18, labelpad=15, fontweight='bold')
    plt.ylabel('Mature Length (nt)', fontsize=18, labelpad=15, fontweight='bold')
    plt.title('Mature Length Distribution by circRNA Type', fontsize=22, pad=20, fontweight='bold')
    
    # Calculate the max length with some padding for the plot range
    max_length = df['mature_length'].max()
    max_y = max(5000, max_length * 1.1)  # Don't go below 5000, but don't go too much higher than max
    
    # Use log scale with appropriate limits
    plt.yscale('log')
    plt.ylim(50, max_y)  # Set minimum to 50 to avoid too much empty space at bottom
    
    # Determine appropriate tick locations based on data range
    if max_y <= 2000:
        yticks = [100, 200, 500, 1000, 2000]
    elif max_y <= 5000:
        yticks = [100, 200, 500, 1000, 2000, 5000]
    else:
        yticks = [100, 200, 500, 1000, 2000, 5000, 10000]
    
    # Keep only ticks that are within the plot range
    yticks = [y for y in yticks if y <= max_y]
    
    plt.yticks(yticks)
    
    # Format Y-axis with regular numbers, not scientific notation
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.get_major_formatter().set_scientific(False)
    
    # Add grid lines for easier reading of values
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Adjust layout and save
    plt.tight_layout()
    filename = os.path.join(output_dir, 'mature_length_by_type_boxplot.png')
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    return filename

def create_stacked_bar(data, x_categories, stack_categories, colors, title, xlabel, ylabel, filename, 
                      add_percentages=True, legend_title=None, count_note=None):
    """Generic function to create a stacked bar chart with consistent styling"""
    plt.figure(figsize=(14, 10))
    
    # Plot stacked bar chart
    ax = data.plot(
        kind='bar', 
        stacked=True, 
        figsize=(14, 10),
        color=colors,
        width=0.7,
        edgecolor='black',
        linewidth=1
    )
    
    # Add percentage labels on bars if requested
    if add_percentages:
        for i, (idx, row) in enumerate(data.iterrows()):
            cumulative = 0
            for j, val in enumerate(row):
                if val > 5:  # Only add text if percentage is large enough
                    ax.text(
                        i, 
                        cumulative + val/2, 
                        f"{val:.1f}%", 
                        ha='center', 
                        va='center',
                        fontsize=14,
                        fontweight='bold',
                        color='black'
                    )
                cumulative += val
    
    # Set x-tick labels
    if x_categories:
        plt.xticks(
            range(len(data.index)),
            [x_categories.get(t, t) for t in data.index],
            rotation=45, 
            ha='right', 
            fontsize=16,
            fontweight='bold'
        )
    
    # Set labels and title
    plt.xlabel(xlabel, fontsize=18, labelpad=15, fontweight='bold')
    plt.ylabel(ylabel, fontsize=18, labelpad=15, fontweight='bold')
    plt.title(title, fontsize=22, pad=20, fontweight='bold')
    
    # Create legend
    legend_kwargs = {
        'fontsize': 14,
        'bbox_to_anchor': (1.02, 1),
        'loc': 'upper left',
        'frameon': True,
        'framealpha': 1,
        'edgecolor': 'black'
    }
    
    if legend_title:
        legend_kwargs['title'] = legend_title
        legend_kwargs['title_fontsize'] = 16
    
    plt.legend(**legend_kwargs)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    return filename

def create_stacked_bar_exon_count(df, output_dir):
    """Create stacked bar chart showing exon count distribution by circRNA type"""
    # Group exon counts 5 and above
    df['exon_group'] = df['block_count'].apply(lambda x: x if x < 5 else 5)
    df['exon_group'] = df['exon_group'].map({1: 1, 2: 2, 3: 3, 4: 4, 5: '5+'})
    
    # Calculate percentages of exon counts within each circRNA type
    cross_tab = pd.crosstab(
        df['circrna_type'], 
        df['exon_group'], 
        normalize='index'
    ) * 100
    
    # Ensure consistent order
    circ_type_order = ['eciRNA', 'EIciRNA', 'ciRNA', 'intergenic']
    exon_order = [1, 2, 3, 4, '5+']
    
    # Reindex to ensure all categories are present (fill missing with 0)
    cross_tab = cross_tab.reindex(columns=exon_order, fill_value=0)
    cross_tab = cross_tab.reindex(circ_type_order, fill_value=0)
    
    # Create custom exon count color map
    exon_colors = [EXON_COUNT_COLORS.get(e, '#a6761d') for e in exon_order]
    
    # Create the stacked bar chart
    filename = os.path.join(output_dir, 'exon_count_by_type_stacked.png')
    return create_stacked_bar(
        data=cross_tab,
        x_categories=TYPE_LABELS,
        stack_categories=None,
        colors=exon_colors,
        title='Exon Count Distribution by circRNA Type',
        xlabel='circRNA Type',
        ylabel='Percentage',
        filename=filename,
        add_percentages=True,
        legend_title='Number of Exons'
    )

def create_splice_site_by_type_stacked_bar(df, output_dir):
    """Create stacked bar chart showing splice site distribution by circRNA type"""
    if 'splice_site_type' not in df.columns:
        print("Warning: 'splice_site_type' column not found, skipping splice site stacked bar chart")
        return None
    
    # Calculate percentages of splice site types within each circRNA type
    cross_tab = pd.crosstab(
        df['circrna_type'], 
        df['splice_site_type'], 
        normalize='index'
    ) * 100
    
    # Ensure consistent order
    circ_type_order = ['eciRNA', 'EIciRNA', 'ciRNA', 'intergenic']
    splice_order = ['GT-AG', 'GC-AG', 'AT-AC', 'Non-canonical', 'Unknown']
    
    # Keep only columns that exist in the data
    splice_order = [s for s in splice_order if s in cross_tab.columns]
    
    # Reindex to ensure consistent ordering
    cross_tab = cross_tab.reindex(columns=splice_order, fill_value=0)
    cross_tab = cross_tab.reindex(circ_type_order, fill_value=0)
    
    # Create custom color map for splice sites
    splice_colors = [SPLICE_COLORS.get(s, '#a6761d') for s in splice_order]
    
    # Get type counts for the note at bottom
    type_counts = df['circrna_type'].value_counts()
    count_note = "Total counts: " + ", ".join([
        f"{TYPE_LABELS[t].split(' ')[0]}: {type_counts.get(t, 0)}" 
        for t in circ_type_order if t in df['circrna_type'].unique()
    ])
    
    # Create the stacked bar chart
    filename = os.path.join(output_dir, 'splice_site_by_type_stacked.png')
    return create_stacked_bar(
        data=cross_tab,
        x_categories=TYPE_LABELS,
        stack_categories=None,
        colors=splice_colors,
        title='Splice Site Distribution by circRNA Type',
        xlabel='circRNA Type',
        ylabel='Percentage',
        filename=filename,
        add_percentages=True,
        legend_title='Splice Site Type',
        count_note=count_note
    )

def generate_visualizations(csv_file, output_dir):
    """Generate all visualizations from the results CSV file"""
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Set consistent style
    set_style()
    
    # Load data
    print(f"Loading data from {csv_file}")
    df = pd.read_csv(csv_file)
    
    # Track all generated files
    generated_files = []
    
    try:
        # 1. circRNA Type Distribution - Improved pie chart with optimally positioned labels
        print("Creating circRNA type pie chart with optimized labels")
        type_counts = df['circrna_type'].value_counts()
        labels = [TYPE_LABELS.get(idx, idx) for idx in type_counts.index]
        colors = [CIRCRNA_COLORS.get(idx, '#a6761d') for idx in type_counts.index]
        
        filename = os.path.join(output_dir, 'circrna_type_pie.png')
        create_improved_pie_chart(
            type_counts.values, 
            labels, 
            colors, 
            'circRNA Type Distribution', 
            filename
        )
        generated_files.append(filename)
        
        # 2. Mature Length by circRNA Type - Improved Box Plot with better Y-axis
        print("Creating improved mature length boxplot")
        filename = create_improved_boxplot(df, output_dir)
        generated_files.append(filename)
        
        # 3. Exon Count by circRNA Type - Stacked Bar Plot
        print("Creating exon count stacked bar chart")
        filename = create_stacked_bar_exon_count(df, output_dir)
        generated_files.append(filename)
        
        # 4. Splice Site Type Distribution - Improved pie chart
        if 'splice_site_type' in df.columns:
            print("Creating splice site type pie chart with optimized labels")
            splice_counts = df['splice_site_type'].value_counts()
            labels = splice_counts.index
            colors = [SPLICE_COLORS.get(idx, '#a6761d') for idx in splice_counts.index]
            
            filename = os.path.join(output_dir, 'splice_site_pie.png')
            create_improved_pie_chart(
                splice_counts.values, 
                labels, 
                colors, 
                'Splice Site Type Distribution', 
                filename
            )
            generated_files.append(filename)
            
            # 5. Splice Site by circRNA Type - Stacked Bar Plot
            print("Creating splice site by circRNA type stacked bar chart")
            filename = create_splice_site_by_type_stacked_bar(df, output_dir)
            if filename:
                generated_files.append(filename)
        
        print(f"Successfully created {len(generated_files)} visualizations")
        
        # Generate HTML report for easy viewing
        html_report = create_html_report(generated_files, df, output_dir)
        print(f"Generated HTML report: {html_report}")
        
        return generated_files
        
    except Exception as e:
        print(f"Error generating visualizations: {e}")
        import traceback
        traceback.print_exc()
        return []

def create_html_report(image_files, df, output_dir):
    """Create an HTML report with all visualizations"""
    # Calculate summary statistics
    summary = {
        'total_circrnas': len(df),
        'circrna_types': df['circrna_type'].value_counts().to_dict(),
        'canonical_count': df['is_canonical'].sum() if 'is_canonical' in df.columns else 'N/A',
        'multi_exon_count': (df['block_count'] > 1).sum(),
        'single_exon_count': (df['block_count'] == 1).sum(),
    }
    
    # Calculate median lengths by type
    for circ_type in df['circrna_type'].unique():
        summary[f'{circ_type}_median_length'] = df[df['circrna_type'] == circ_type]['mature_length'].median()
    
    # Create HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>circRNA Analysis Visualization Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                line-height: 1.6;
                color: #333;
                max-width: 1200px;
                margin: 0 auto;
                padding: 20px;
            }}
            h1, h2, h3 {{
                color: #2c3e50;
            }}
            .summary-table {{
                border-collapse: collapse;
                width: 100%;
                margin: 20px 0;
            }}
            .summary-table th, .summary-table td {{
                border: 1px solid #ddd;
                padding: 12px;
                text-align: left;
            }}
            .summary-table th {{
                background-color: #f2f2f2;
            }}
            .visualization {{
                margin: 30px 0;
                text-align: center;
            }}
            .visualization img {{
                max-width: 100%;
                height: auto;
                border: 1px solid #ddd;
                box-shadow: 0 4px 8px rgba(0,0,0,0.1);
            }}
            .caption {{
                margin-top: 10px;
                font-style: italic;
                color: #7f8c8d;
            }}
        </style>
    </head>
    <body>
        <h1>circRNA Analysis Visualization Report</h1>
        
        <h2>Summary Statistics</h2>
        <table class="summary-table">
            <tr>
                <th>Metric</th>
                <th>Value</th>
            </tr>
            <tr>
                <td>Total circRNAs analyzed</td>
                <td>{summary['total_circrnas']}</td>
            </tr>
    """
    
    # Add circRNA type counts
    for circ_type, count in summary['circrna_types'].items():
        type_label = TYPE_LABELS.get(circ_type, circ_type)
        percentage = (count / summary['total_circrnas']) * 100
        html_content += f"""
            <tr>
                <td>{type_label}</td>
                <td>{count} ({percentage:.1f}%)</td>
            </tr>
        """
    
    # Add exon counts
    html_content += f"""
            <tr>
                <td>Single-exon circRNAs</td>
                <td>{summary['single_exon_count']} ({summary['single_exon_count']/summary['total_circrnas']*100:.1f}%)</td>
            </tr>
            <tr>
                <td>Multi-exon circRNAs</td>
                <td>{summary['multi_exon_count']} ({summary['multi_exon_count']/summary['total_circrnas']*100:.1f}%)</td>
            </tr>
    """
    
    # Add canonical splice site info if available
    if summary['canonical_count'] != 'N/A':
        canonical_percent = (summary['canonical_count'] / summary['total_circrnas']) * 100
        html_content += f"""
            <tr>
                <td>Canonical splice sites</td>
                <td>{summary['canonical_count']} ({canonical_percent:.1f}%)</td>
            </tr>
        """
    
    # Add median lengths by type
    html_content += """
            <tr>
                <th colspan="2">Median Mature Lengths</th>
            </tr>
    """
    
    for circ_type in df['circrna_type'].unique():
        type_label = TYPE_LABELS.get(circ_type, circ_type)
        median_length = summary.get(f'{circ_type}_median_length', 'N/A')
        html_content += f"""
            <tr>
                <td>{type_label}</td>
                <td>{median_length:.1f} nt</td>
            </tr>
        """
    
    html_content += """
        </table>
        
        <h2>Visualizations</h2>
    """
    
    # Add all visualization images
    for img_file in image_files:
        img_name = os.path.basename(img_file)
        caption = img_name.replace('_', ' ').replace('.png', '')
        caption = ' '.join(word.capitalize() for word in caption.split())
        
        html_content += f"""
        <div class="visualization">
            <img src="{img_name}" alt="{caption}">
            <p class="caption">{caption}</p>
        </div>
        """
    
    html_content += """
    </body>
    </html>
    """
    
    # Write HTML to file
    html_path = os.path.join(output_dir, 'visualization_report.html')
    with open(html_path, 'w') as f:
        f.write(html_content)
    
    return html_path

def main():
    parser = argparse.ArgumentParser(description='Create complete set of circRNA visualizations')
    parser.add_argument('--csv', required=True, help='Path to circRNA results CSV file')
    parser.add_argument('--output', required=True, help='Output directory for visualizations')
    args = parser.parse_args()
    
    if not os.path.exists(args.csv):
        print(f"Error: CSV file not found at {args.csv}")
        return 1
    
    # Generate visualizations
    generated_files = generate_visualizations(args.csv, args.output)
    
    if generated_files:
        print("Visualization generation completed successfully.")
        return 0
    else:
        print("Visualization generation failed.")
        return 1

if __name__ == "__main__":
    exit(main())