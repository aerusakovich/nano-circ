#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import numpy as np
import re

# Tool colors 
TOOL_COLORS = {
    'ciri-long': '#0072B2',  # Blue
    'isocirc': '#009E73',    # Green
    'circnick': '#D55E00',   # Orange-red
    'circnick-lrs': '#D55E00'  # Orange-red (alternative naming)
}

# Tool display names
TOOL_DISPLAY_NAMES = {
    'ciri-long': 'CIRI-long',
    'isocirc': 'IsoCirc',
    'circnick': 'CircNick-lrs',
    'circnick-lrs': 'CircNick-lrs'
}

# Tool markers for the scatter plot
TOOL_MARKERS = {
    'ciri-long': 'o',      # Circle
    'isocirc': 's',        # Square
    'circnick': '^',       # Triangle
    'circnick-lrs': '^'    # Triangle (alternative naming)
}

def parse_memory_file(memory_file):
    """
    Parse the memory file which has a non-standard format
    """
    memory_rows = []
    
    try:
        with open(memory_file, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            # Use regex to extract tool name and memory value
            match = re.match(r'([\w-]+)[,\s]+([\d\.]+K)', line)
            if match:
                tool, memory = match.groups()
                memory_rows.append({'Tool': tool.strip(), 'MaxCPU_K': memory.strip()})
        
        # Create DataFrame
        memory_df = pd.DataFrame(memory_rows)
        return memory_df
        
    except Exception as e:
        print(f"Error parsing memory file: {e}")
        return None

def create_combined_performance_plot(time_csv, memory_csv, output_dir, total_reads=400000, use_log_scale=True):
    """
    Create a dot plot showing the relationship between running time and RAM usage
    for each cirRNA analysis tool.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # Read the time data
        print(f"Reading performance data from: {time_csv}")
        time_df = pd.read_csv(time_csv)
        
        # Process time data
        time_df['Tool'] = time_df['Tool'].str.strip()
        time_df['Tool_Lower'] = time_df['Tool'].str.lower()
        time_df['DisplayName'] = time_df['Tool_Lower'].map(lambda x: TOOL_DISPLAY_NAMES.get(x, x.capitalize()))
        time_df['ElapsedTime_Min'] = time_df['ElapsedTime_Sec'] / 60
        
        # Read memory data using custom parser
        print(f"Reading memory data from: {memory_csv}")
        memory_df = parse_memory_file(memory_csv)
        
        if memory_df is None or memory_df.empty:
            print("ERROR: Failed to parse memory data")
            return False
            
        print("Memory data:")
        print(memory_df)
        
        # Process memory data
        memory_df['Tool_Lower'] = memory_df['Tool'].str.lower()
        
        # Convert memory KB to GB
        memory_df['MaxCPU_GB'] = memory_df['MaxCPU_K'].apply(
            lambda x: float(str(x).replace('K', '')) / 1024 / 1024)
        
        # Create a final dataframe with all tools and their metrics
        result_data = []
        
        # Process each tool in time data
        for _, time_row in time_df.iterrows():
            tool = time_row['Tool_Lower']
            display_name = time_row['DisplayName']
            elapsed_time_min = time_row['ElapsedTime_Min']
            
            # Find matching memory data
            memory_row = None
            
            # Try exact match
            memory_matches = memory_df[memory_df['Tool_Lower'] == tool]
            if not memory_matches.empty:
                memory_row = memory_matches.iloc[0]
            else:
                # Try alternative names for CircNick
                if tool == 'circnick':
                    alt_matches = memory_df[memory_df['Tool_Lower'] == 'circnick-lrs']
                    if not alt_matches.empty:
                        memory_row = alt_matches.iloc[0]
                elif tool == 'circnick-lrs':
                    alt_matches = memory_df[memory_df['Tool_Lower'] == 'circnick']
                    if not alt_matches.empty:
                        memory_row = alt_matches.iloc[0]
            
            if memory_row is not None:
                max_cpu_gb = memory_row['MaxCPU_GB']
                
                # Add to results
                result_data.append({
                    'Tool': tool,
                    'DisplayName': display_name,
                    'ElapsedTime_Min': elapsed_time_min,
                    'MaxCPU_GB': max_cpu_gb
                })
            else:
                print(f"WARNING: No memory data found for tool: {tool}")
        
        # Process tools in memory data that might not be in time data
        for _, memory_row in memory_df.iterrows():
            tool = memory_row['Tool_Lower']
            
            # Check if we already added this tool
            if any(row['Tool'] == tool for row in result_data):
                continue
                
            # Try to find matching time data
            time_row = None
            
            # Try exact match
            time_matches = time_df[time_df['Tool_Lower'] == tool]
            if not time_matches.empty:
                time_row = time_matches.iloc[0]
            else:
                # Try alternative names for CircNick
                if tool == 'circnick-lrs':
                    alt_matches = time_df[time_df['Tool_Lower'] == 'circnick']
                    if not alt_matches.empty:
                        time_row = alt_matches.iloc[0]
                elif tool == 'circnick':
                    alt_matches = time_df[time_df['Tool_Lower'] == 'circnick-lrs']
                    if not alt_matches.empty:
                        time_row = alt_matches.iloc[0]
            
            if time_row is not None:
                display_name = time_row['DisplayName']
                elapsed_time_min = time_row['ElapsedTime_Min']
                max_cpu_gb = memory_row['MaxCPU_GB']
                
                # Add to results
                result_data.append({
                    'Tool': tool,
                    'DisplayName': display_name,
                    'ElapsedTime_Min': elapsed_time_min,
                    'MaxCPU_GB': max_cpu_gb
                })
                
        # Convert to DataFrame
        result_df = pd.DataFrame(result_data)
        
        # Verify we have data
        print("Final performance data:")
        print(result_df[['DisplayName', 'ElapsedTime_Min', 'MaxCPU_GB']])
        
        if result_df.empty:
            print("ERROR: No matching data found")
            return False
        
        # Create a more square plot with similar style to benchmark plots
        # New dimensions: 10x7 (width x height) - slightly taller for better proportions
        plt.figure(figsize=(10, 7))
        
        # Create scatter plot with custom markers
        for _, row in result_df.iterrows():
            tool = row['Tool']
            marker = TOOL_MARKERS.get(tool, 'o')
            color = TOOL_COLORS.get(tool, '#333333')
            
            plt.scatter(
                row['ElapsedTime_Min'], 
                row['MaxCPU_GB'], 
                s=300,  # Marker size
                color=color, 
                marker=marker, 
                label=row['DisplayName'],
                edgecolors='black',
                linewidths=1.5,
                alpha=0.8,
                zorder=3
            )
            
            # Add tool name as text label with benchmark-like styling
            # For horizontal layout, adjust label positions to avoid overlap
            if 'ciri' in tool:
                x_offset = 0
                y_offset = -15  # Place below point
                ha = 'center'
                va = 'top'
            elif 'circ' in tool:
                x_offset = 0
                y_offset = 15  # Place above point
                ha = 'center'
                va = 'bottom'
            else:  # isocirc
                x_offset = 15  # Place to the right
                y_offset = 0
                ha = 'left'
                va = 'center'
                
            plt.annotate(
                row['DisplayName'],
                (row['ElapsedTime_Min'], row['MaxCPU_GB']),
                xytext=(x_offset, y_offset),
                textcoords='offset points',
                fontsize=12,
                fontweight='bold',
                va=va,
                ha=ha,
                color=color
            )
        
        # Add styling and labels
        plt.title('CircRNA Analysis Tools: Time vs Memory Performance', fontsize=20)
        plt.xlabel('Elapsed Time (minutes)', fontsize=16)
        plt.ylabel('Memory Usage (GB)', fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7, zorder=0)
        
        # Add subtitle with number of reads
        plt.figtext(0.5, 0.01, f'Based on analysis of {total_reads:,} reads', 
                   ha='center', fontsize=12, style='italic')
        
        # Calculate 10 minutes position as a fraction of the x-axis
        x_max_val = result_df['ElapsedTime_Min'].max() * 1.1
        ten_min_fraction = min(10 / x_max_val, 0.3)  # Adjusted for horizontal layout
        
        # Create the ideal corner box - positioned around 10 minutes mark
        # Adjusted for horizontal layout
        plt.annotate(
            'Ideal\n(faster & lower memory)',
            xy=(ten_min_fraction, 0.15),  # Adjusted y position for horizontal layout
            xycoords='axes fraction',
            bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", ec="gold", alpha=0.8),
            fontsize=14,  # Adjusted for horizontal layout
            ha='center',
            va='bottom'
        )
        
        # Keep the original arrow position but adjust for horizontal layout
        plt.arrow(0.05, 0.1, -0.03, -0.05, head_width=0.01, head_length=0.01, 
                 fc='black', ec='black', transform=plt.gca().transAxes, zorder=5)
        
        # Set axis limits
        x_min = 0
        x_max = result_df['ElapsedTime_Min'].max() * 1.1
        plt.xlim(x_min, x_max)
        
        # Use log scale if needed
        memory_ratio = result_df['MaxCPU_GB'].max() / result_df['MaxCPU_GB'].min()
        if use_log_scale and memory_ratio > 10:
            plt.yscale('log')
            print(f"Using logarithmic scale for memory axis (ratio: {memory_ratio:.1f}x)")
            plt.grid(True, which="both", linestyle='--', alpha=0.6, zorder=0)
            plt.ylabel('Memory Usage (GB, log scale)', fontsize=16)
        else:
            y_min = 0
            y_max = result_df['MaxCPU_GB'].max() * 1.1
            plt.ylim(y_min, y_max)
        
        # Save the plot
        plt.tight_layout(rect=[0, 0.03, 1, 1])
        output_file = os.path.join(output_dir, 'time_vs_memory_performance.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Combined performance plot saved to: {output_file}")
        return True
    
    except Exception as e:
        print(f"Error creating combined performance plot: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Create combined performance plot')
    parser.add_argument('--time-csv', 
                        default=None,
                        help='Path to CSV file with time performance metrics')
    parser.add_argument('--memory-csv', 
                        default=None,
                        help='Path to CSV file with memory metrics')
    parser.add_argument('--output-dir', 
                        default='circrna_analysis_results',
                        help='Output directory for plots')
    parser.add_argument('--total-reads', 
                        type=int, 
                        default=400000,
                        help='Total number of reads processed (default: 400,000)')
    parser.add_argument('--use-log-scale',
                        action='store_true',
                        default=True,
                        help='Use logarithmic scale for memory axis')
    args = parser.parse_args()
    
    # Create the combined plot
    create_combined_performance_plot(args.time_csv, args.memory_csv, args.output_dir, args.total_reads, args.use_log_scale)

if __name__ == "__main__":
    main()