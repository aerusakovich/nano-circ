process COLLECT_PERFORMANCE {
    tag "performance_summary"
    
    publishDir "${params.outdir}/performance", mode: 'copy'
    
    cpus = 1
    memory = '4 GB'
    time = '30m'
    
    input:
    path trace_file
    
    output:
    path "performance_time_individual.csv", emit: time_csv
    path "performance_memory_individual.csv", emit: memory_csv
    
    script:
    """
    #!/bin/bash
    
    set -e
    
    python3 << 'EOFPYTHON'
import pandas as pd
import re
from pathlib import Path

print("=" * 60)
print("COLLECTING PERFORMANCE DATA")
print("=" * 60)

# Read trace file
trace_df = pd.read_csv('${trace_file}', sep='\\t')

# Filter for tool processes
tool_processes = trace_df[trace_df['name'].str.contains('RUN_CIRILONG|RUN_ISOCIRC|RUN_CIRCNICK', case=False, na=False)]

print(f"\\nFound {len(tool_processes)} tool processes in trace file")

# Read timing data from performance_metrics.txt files
print("\\n" + "=" * 60)
print("READING TIMING DATA")
print("=" * 60)

tool_outputs_dir = Path("${params.outdir}/tool_outputs")

time_data = []
tools = ['cirilong', 'isocirc', 'circnick']
tool_name_map = {
    'cirilong': 'ciri-long',
    'isocirc': 'isocirc',
    'circnick': 'circnick'
}

# Build list of all runs we expect
expected_runs = []

for tool_dir in tools:
    tool_path = tool_outputs_dir / tool_dir
    if not tool_path.exists():
        continue
    
    for run_dir in tool_path.glob('run_*'):
        run_name = run_dir.name
        perf_file = run_dir / 'performance_metrics.txt'
        
        if not perf_file.exists():
            continue
        
        with open(perf_file, 'r') as f:
            lines = f.readlines()
        
        metrics = {}
        for line in lines:
            if ':' in line:
                key, value = line.strip().split(':', 1)
                metrics[key.strip()] = value.strip()
        
        tool_name = tool_name_map[tool_dir]
        elapsed_sec = int(metrics.get('Elapsed_Seconds', 0))
        threads = int(metrics.get('Threads', ${params.threads}))
        
        time_data.append({
            'Tool': tool_name,
            'Run': run_name,
            'ElapsedTime_Sec': elapsed_sec,
            'CPUs': threads,
            'Status': 'COMPLETED'
        })
        
        expected_runs.append((tool_name, run_name, elapsed_sec))
        print(f"  {tool_name} {run_name}: {elapsed_sec} seconds")

time_df = pd.DataFrame(time_data)
print(f"\\nCollected timing data for {len(time_df)} runs")

# Read memory data from trace file
print("\\n" + "=" * 60)
print("READING MEMORY DATA FROM TRACE FILE")
print("=" * 60)

memory_data = []

# For each expected run, find matching trace entry
for exp_tool, exp_run, exp_time in expected_runs:
    # Find matching trace entry by tool and approximate timing
    matched = False
    
    for _, row in tool_processes.iterrows():
        name = row['name']
        
        # Check if tool matches
        trace_tool = None
        if 'CIRILONG' in name.upper():
            trace_tool = 'ciri-long'
        elif 'ISOCIRC' in name.upper():
            trace_tool = 'isocirc'
        elif 'CIRCNICK' in name.upper():
            trace_tool = 'circnick'
        
        if trace_tool != exp_tool:
            continue
        
        # Check if run name matches (it's in parentheses in the name)
        if f'({exp_run})' in name:
            # Found matching trace entry - extract memory
            peak_rss = row.get('peak_rss', '')
            
            if pd.notna(peak_rss) and peak_rss and peak_rss != '-':
                # Parse memory value (format: "283.8 GB" or similar)
                peak_rss_str = str(peak_rss).strip()
                
                # Extract number and unit
                match = re.match(r'([0-9.]+)\\s*([KMGT]?B?)', peak_rss_str)
                if match:
                    value = float(match.group(1))
                    unit = match.group(2).replace('B', '')  # Remove 'B' if present
                    
                    # Convert to KB
                    if unit == 'M':
                        value_kb = value * 1024
                    elif unit == 'G':
                        value_kb = value * 1024 * 1024
                    elif unit == 'T':
                        value_kb = value * 1024 * 1024 * 1024
                    elif unit == 'K' or unit == '':
                        value_kb = value
                    else:
                        value_kb = value
                    
                    memory_data.append({
                        'Tool': exp_tool,
                        'Run': exp_run,
                        'MaxCPU_K': f"{int(value_kb)}K"
                    })
                    print(f"  {exp_tool} {exp_run}: {int(value_kb):,} KB ({value_kb/1024/1024:.1f} GB)")
                    matched = True
                    break
    
    if not matched:
        print(f"  WARNING: No memory data found for {exp_tool} {exp_run}")

memory_df = pd.DataFrame(memory_data)
print(f"\\nCollected memory data for {len(memory_df)} runs")

# Save data
time_df.to_csv('performance_time_individual.csv', index=False)

with open('performance_memory_individual.csv', 'w') as f:
    for _, row in memory_df.iterrows():
        f.write(f"{row['Tool']}, {row['Run']}, {row['MaxCPU_K']}\\n")

print("\\n" + "=" * 60)
print("PERFORMANCE DATA COLLECTION COMPLETE")
print("=" * 60)
print(f"Time data: {len(time_df)} runs")
print(f"Memory data: {len(memory_df)} runs")

EOFPYTHON
    """
}

process VISUALIZE_PERFORMANCE {
    tag "performance_plots"
    
    publishDir "${params.outdir}/performance", mode: 'copy'
    
    cpus = 1
    memory = '4 GB'
    time = '30m'
    
    input:
    path time_csv
    path memory_csv
    
    output:
    path "time_vs_memory_performance.png"
    path "performance_summary.txt"
    
    script:
    def total_reads = params.total_reads
    def circrna_reads = params.circrna_reads
    def linear_reads = params.linear_reads
    def threads = params.threads
    """
    #!/bin/bash
    
    # Activate conda environment with matplotlib
    source /local/env/envpython-3.9.5.sh
    conda activate bed12
    
    # Run Python script
    python3 << 'EOF'
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import re

# Tool colors 
TOOL_COLORS = {
    'ciri-long': '#0072B2',
    'isocirc': '#009E73',
    'circnick': '#D55E00',
    'circnick-lrs': '#D55E00'
}

# Tool display names
TOOL_DISPLAY_NAMES = {
    'ciri-long': 'CIRI-long',
    'isocirc': 'IsoCirc',
    'circnick': 'CircNick-lrs',
    'circnick-lrs': 'CircNick-lrs'
}

# Tool markers
TOOL_MARKERS = {
    'ciri-long': 'o',
    'isocirc': 's',
    'circnick': '^',
    'circnick-lrs': '^'
}

def parse_memory_file(memory_file):
    memory_rows = []
    with open(memory_file, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        # Parse format: Tool, Run, MaxCPU_K
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 3:
            tool = parts[0]
            run = parts[1]
            memory = parts[2]
            memory_rows.append({'Tool': tool, 'Run': run, 'MaxCPU_K': memory})
    
    return pd.DataFrame(memory_rows)

# Read INDIVIDUAL run data
print("Reading time data from: ${time_csv}")
time_df = pd.read_csv('${time_csv}')

print("Reading memory data from: ${memory_csv}")
memory_df = parse_memory_file('${memory_csv}')

print("\\nTime data (individual runs):")
print(time_df)

print("\\nMemory data (individual runs):")
print(memory_df)

if memory_df.empty:
    print("ERROR: No memory data available!")
    open('time_vs_memory_performance.png', 'w').close()
    with open('performance_summary.txt', 'w') as f:
        f.write("No memory data available\\n")
    exit(1)

# Process time data
time_df['Tool'] = time_df['Tool'].str.strip()
time_df['Tool_Lower'] = time_df['Tool'].str.lower()
time_df['DisplayName'] = time_df['Tool_Lower'].map(lambda x: TOOL_DISPLAY_NAMES.get(x, x.capitalize()))
time_df['ElapsedTime_Min'] = time_df['ElapsedTime_Sec'] / 60

# Calculate mean and std for each tool (across ALL runs)
time_stats = time_df.groupby('Tool_Lower').agg({
    'ElapsedTime_Min': ['mean', 'std', 'count'],
    'DisplayName': 'first'
}).reset_index()

time_stats.columns = ['Tool_Lower', 'ElapsedTime_Min_Mean', 'ElapsedTime_Min_Std', 'Run_Count', 'DisplayName']
time_stats['ElapsedTime_Min_Std'] = time_stats['ElapsedTime_Min_Std'].fillna(0)

print("\\nTime statistics (mean ± std):")
for _, row in time_stats.iterrows():
    print("  {}: {:.3f} ± {:.3f} min (n={})".format(
        row['DisplayName'], row['ElapsedTime_Min_Mean'], 
        row['ElapsedTime_Min_Std'], int(row['Run_Count'])))

# Process memory data
memory_df['Tool'] = memory_df['Tool'].str.strip()
memory_df['Tool_Lower'] = memory_df['Tool'].str.lower()
memory_df['MaxCPU_GB'] = memory_df['MaxCPU_K'].apply(
    lambda x: float(str(x).replace('K', '')) / 1024 / 1024)

# Calculate mean and std for memory across runs
memory_stats = memory_df.groupby('Tool_Lower').agg({
    'MaxCPU_GB': ['mean', 'std', 'count']
}).reset_index()

memory_stats.columns = ['Tool_Lower', 'MaxCPU_GB_Mean', 'MaxCPU_GB_Std', 'Memory_Run_Count']
memory_stats['MaxCPU_GB_Std'] = memory_stats['MaxCPU_GB_Std'].fillna(0)

print("\\nMemory statistics (mean ± std):")
for _, row in memory_stats.iterrows():
    print("  {}: {:.2f} ± {:.2f} GB (n={})".format(
        row['Tool_Lower'], row['MaxCPU_GB_Mean'], 
        row['MaxCPU_GB_Std'], int(row['Memory_Run_Count'])))

# Merge time and memory statistics
result_df = time_stats.merge(memory_stats, on='Tool_Lower', how='inner')

print("\\nFinal merged data:")
print(result_df[['DisplayName', 'Run_Count', 'ElapsedTime_Min_Mean', 'ElapsedTime_Min_Std', 
                 'MaxCPU_GB_Mean', 'MaxCPU_GB_Std']])

if result_df.empty:
    print("ERROR: No data to plot")
    open('time_vs_memory_performance.png', 'w').close()
    with open('performance_summary.txt', 'w') as f:
        f.write("No performance data available\\n")
    exit(0)

# Create plot
fig, ax = plt.subplots(figsize=(10, 7))

for _, row in result_df.iterrows():
    tool = row['Tool_Lower']
    marker = TOOL_MARKERS.get(tool, 'o')
    color = TOOL_COLORS.get(tool, '#333333')
    
    # Plot main point with error bars on BOTH axes
    ax.errorbar(
        row['ElapsedTime_Min_Mean'], 
        row['MaxCPU_GB_Mean'],
        xerr=row['ElapsedTime_Min_Std'] if row['ElapsedTime_Min_Std'] > 0 else None,
        yerr=row['MaxCPU_GB_Std'] if row['MaxCPU_GB_Std'] > 0 else None,
        fmt=marker,
        markersize=15,
        color=color,
        label=row['DisplayName'],
        markeredgecolor='black',
        markeredgewidth=1.5,
        capsize=5,
        capthick=2,
        elinewidth=2,
        alpha=0.8,
        zorder=3
    )
    
    # Position labels
    if 'ciri' in tool:
        x_offset, y_offset = 0, -25
        ha, va = 'center', 'top'
    elif 'circ' in tool:
        x_offset, y_offset = 0, 25
        ha, va = 'center', 'bottom'
    else:
        x_offset, y_offset = 25, 0
        ha, va = 'left', 'center'
        
    label_text = "{} (n={})".format(row['DisplayName'], int(row['Run_Count']))
    ax.annotate(
        label_text,
        (row['ElapsedTime_Min_Mean'], row['MaxCPU_GB_Mean']),
        xytext=(x_offset, y_offset),
        textcoords='offset points',
        fontsize=11,
        fontweight='bold',
        va=va,
        ha=ha,
        color=color
    )

ax.set_title('CircRNA Analysis Tools: Time vs Memory Performance', fontsize=20)
ax.set_xlabel('Elapsed Time (minutes)', fontsize=16)
ax.set_ylabel('Memory Usage (GB)', fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.grid(True, linestyle='--', alpha=0.7, zorder=0)

fig.text(0.5, 0.02, 
         'Total reads: {:,} ({:,} circRNAs + {:,} linear) | Mean ± SD across runs'.format(
             ${total_reads}, ${circrna_reads}, ${linear_reads}), 
         ha='center', fontsize=11, style='italic')

# Ideal corner annotation
x_max_val = (result_df['ElapsedTime_Min_Mean'] + result_df['ElapsedTime_Min_Std']).max() * 1.1
if x_max_val > 0:
    ten_min_fraction = min(10 / x_max_val, 0.3)
else:
    ten_min_fraction = 0.3

ax.annotate(
    'Ideal\\n(faster & lower memory)',
    xy=(ten_min_fraction, 0.15),
    xycoords='axes fraction',
    bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", ec="gold", alpha=0.8),
    fontsize=14,
    ha='center',
    va='bottom'
)

ax.arrow(0.05, 0.1, -0.03, -0.05, head_width=0.01, head_length=0.01, 
         fc='black', ec='black', transform=ax.transAxes, zorder=5)

# Set limits
x_min = 0
x_max_with_err = (result_df['ElapsedTime_Min_Mean'] + result_df['ElapsedTime_Min_Std']).max()
x_max = max(x_max_with_err * 1.15, result_df['ElapsedTime_Min_Mean'].max() * 1.2)
ax.set_xlim(x_min, x_max)

# Check if we should use log scale for memory
memory_ratio = result_df['MaxCPU_GB_Mean'].max() / result_df['MaxCPU_GB_Mean'].min()
if memory_ratio > 10:
    ax.set_yscale('log')
    ax.set_ylabel('Memory Usage (GB, log scale)', fontsize=16)
    ax.grid(True, which="both", linestyle='--', alpha=0.6, zorder=0)
    y_min = result_df['MaxCPU_GB_Mean'].min() * 0.5
    y_max = (result_df['MaxCPU_GB_Mean'] + result_df['MaxCPU_GB_Std']).max() * 1.5
    ax.set_ylim(y_min, y_max)
else:
    y_min_with_err = (result_df['MaxCPU_GB_Mean'] - result_df['MaxCPU_GB_Std']).min()
    y_max_with_err = (result_df['MaxCPU_GB_Mean'] + result_df['MaxCPU_GB_Std']).max()
    y_min = max(0, y_min_with_err * 0.9)
    y_max = y_max_with_err * 1.15
    ax.set_ylim(y_min, y_max)

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig('time_vs_memory_performance.png', dpi=300, bbox_inches='tight')
plt.close()

# Create summary text
with open('performance_summary.txt', 'w') as f:
    f.write("CircRNA Analysis Tool Performance Summary\\n")
    f.write("=" * 60 + "\\n\\n")
    f.write("Dataset Information:\\n")
    f.write("  Total reads: {:,}\\n".format(${total_reads}))
    f.write("  - circRNAs: {:,}\\n".format(${circrna_reads}))
    f.write("  - Linear RNAs: {:,}\\n".format(${linear_reads}))
    f.write("  CPU cores per tool: {}\\n\\n".format(${threads}))
    f.write("Performance Metrics (mean ± standard deviation):\\n")
    f.write("-" * 60 + "\\n\\n")
    
    for _, row in result_df.iterrows():
        f.write("{}:\\n".format(row['DisplayName']))
        f.write("  Number of runs: {}\\n".format(int(row['Run_Count'])))
        f.write("  Elapsed Time: {:.2f} ± {:.2f} minutes\\n".format(
            row['ElapsedTime_Min_Mean'], row['ElapsedTime_Min_Std']))
        f.write("  Memory Usage: {:.2f} ± {:.2f} GB\\n".format(
            row['MaxCPU_GB_Mean'], row['MaxCPU_GB_Std']))
        
        time_cv = (row['ElapsedTime_Min_Std'] / row['ElapsedTime_Min_Mean'] * 100) if row['ElapsedTime_Min_Mean'] > 0 else 0
        memory_cv = (row['MaxCPU_GB_Std'] / row['MaxCPU_GB_Mean'] * 100) if row['MaxCPU_GB_Mean'] > 0 else 0
        
        f.write("  Time variability (CV): {:.1f}%\\n".format(time_cv))
        f.write("  Memory variability (CV): {:.1f}%\\n\\n".format(memory_cv))

print("\\nPerformance visualization complete!")
EOF
    """
}

process VISUALIZE_THROUGHPUT {
    tag "throughput_plots"
    
    publishDir "${params.outdir}/performance", mode: 'copy'
    
    cpus = 1
    memory = '4 GB'
    time = '30m'
    
    input:
    path time_csv
    path memory_csv
    
    output:
    path "throughput_vs_memory_performance.png"
    path "throughput_summary.txt"
    
    script:
    def total_reads = params.total_reads
    def circrna_reads = params.circrna_reads
    def linear_reads = params.linear_reads
    def threads = params.threads
    """
    #!/bin/bash
    
    # Activate conda environment with matplotlib
    source /local/env/envpython-3.9.5.sh
    conda activate bed12
    
    # Run Python script
    python3 << 'EOF'
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import re

# Tool colors 
TOOL_COLORS = {
    'ciri-long': '#0072B2',
    'isocirc': '#009E73',
    'circnick': '#D55E00',
    'circnick-lrs': '#D55E00'
}

# Tool display names
TOOL_DISPLAY_NAMES = {
    'ciri-long': 'CIRI-long',
    'isocirc': 'IsoCirc',
    'circnick': 'CircNick-lrs',
    'circnick-lrs': 'CircNick-lrs'
}

# Tool markers
TOOL_MARKERS = {
    'ciri-long': 'o',
    'isocirc': 's',
    'circnick': '^',
    'circnick-lrs': '^'
}

def parse_memory_file(memory_file):
    memory_rows = []
    with open(memory_file, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 3:
            tool = parts[0]
            run = parts[1]
            memory = parts[2]
            memory_rows.append({'Tool': tool, 'Run': run, 'MaxCPU_K': memory})
    
    return pd.DataFrame(memory_rows)

# Read data
print("Reading time data from: ${time_csv}")
time_df = pd.read_csv('${time_csv}')

print("Reading memory data from: ${memory_csv}")
memory_df = parse_memory_file('${memory_csv}')

if memory_df.empty:
    print("ERROR: No memory data available!")
    open('throughput_vs_memory_performance.png', 'w').close()
    with open('throughput_summary.txt', 'w') as f:
        f.write("No memory data available\\n")
    exit(1)

# Process time data and calculate throughput (reads per hour)
time_df['Tool'] = time_df['Tool'].str.strip()
time_df['Tool_Lower'] = time_df['Tool'].str.lower()
time_df['DisplayName'] = time_df['Tool_Lower'].map(lambda x: TOOL_DISPLAY_NAMES.get(x, x.capitalize()))
time_df['ElapsedTime_Hours'] = time_df['ElapsedTime_Sec'] / 3600.0
time_df['Throughput_ReadsPerHour'] = ${total_reads} / time_df['ElapsedTime_Hours']

# Calculate mean and std for throughput
time_stats = time_df.groupby('Tool_Lower').agg({
    'Throughput_ReadsPerHour': ['mean', 'std', 'count'],
    'DisplayName': 'first'
}).reset_index()

time_stats.columns = ['Tool_Lower', 'Throughput_Mean', 'Throughput_Std', 'Run_Count', 'DisplayName']
time_stats['Throughput_Std'] = time_stats['Throughput_Std'].fillna(0)

print("\\nThroughput statistics (mean ± std):")
for _, row in time_stats.iterrows():
    print("  {}: {:.0f} ± {:.0f} reads/hour (n={})".format(
        row['DisplayName'], row['Throughput_Mean'], 
        row['Throughput_Std'], int(row['Run_Count'])))

# Process memory data
memory_df['Tool'] = memory_df['Tool'].str.strip()
memory_df['Tool_Lower'] = memory_df['Tool'].str.lower()
memory_df['MaxCPU_GB'] = memory_df['MaxCPU_K'].apply(
    lambda x: float(str(x).replace('K', '')) / 1024 / 1024)

# Calculate mean and std for memory
memory_stats = memory_df.groupby('Tool_Lower').agg({
    'MaxCPU_GB': ['mean', 'std', 'count']
}).reset_index()

memory_stats.columns = ['Tool_Lower', 'MaxCPU_GB_Mean', 'MaxCPU_GB_Std', 'Memory_Run_Count']
memory_stats['MaxCPU_GB_Std'] = memory_stats['MaxCPU_GB_Std'].fillna(0)

# Merge statistics
result_df = time_stats.merge(memory_stats, on='Tool_Lower', how='inner')

print("\\nFinal merged data:")
print(result_df[['DisplayName', 'Run_Count', 'Throughput_Mean', 'Throughput_Std', 
                 'MaxCPU_GB_Mean', 'MaxCPU_GB_Std']])

if result_df.empty:
    print("ERROR: No data to plot")
    open('throughput_vs_memory_performance.png', 'w').close()
    with open('throughput_summary.txt', 'w') as f:
        f.write("No performance data available\\n")
    exit(0)

# Create plot
fig, ax = plt.subplots(figsize=(10, 7))

for _, row in result_df.iterrows():
    tool = row['Tool_Lower']
    marker = TOOL_MARKERS.get(tool, 'o')
    color = TOOL_COLORS.get(tool, '#333333')
    
    # Plot with error bars on both axes
    ax.errorbar(
        row['Throughput_Mean'], 
        row['MaxCPU_GB_Mean'],
        xerr=row['Throughput_Std'] if row['Throughput_Std'] > 0 else None,
        yerr=row['MaxCPU_GB_Std'] if row['MaxCPU_GB_Std'] > 0 else None,
        fmt=marker,
        markersize=15,
        color=color,
        label=row['DisplayName'],
        markeredgecolor='black',
        markeredgewidth=1.5,
        capsize=5,
        capthick=2,
        elinewidth=2,
        alpha=0.8,
        zorder=3
    )
    
    # Position labels
    if 'ciri' in tool:
        x_offset, y_offset = 0, -25
        ha, va = 'center', 'top'
    elif 'circ' in tool:
        x_offset, y_offset = 0, 25
        ha, va = 'center', 'bottom'
    else:
        x_offset, y_offset = 25, 0
        ha, va = 'left', 'center'
        
    label_text = "{} (n={})".format(row['DisplayName'], int(row['Run_Count']))
    ax.annotate(
        label_text,
        (row['Throughput_Mean'], row['MaxCPU_GB_Mean']),
        xytext=(x_offset, y_offset),
        textcoords='offset points',
        fontsize=11,
        fontweight='bold',
        va=va,
        ha=ha,
        color=color
    )

ax.set_title('CircRNA Analysis Tools: Throughput vs Memory Performance', fontsize=20)
ax.set_xlabel('Throughput (reads/hour)', fontsize=16)
ax.set_ylabel('Memory Usage (GB)', fontsize=16)

# Format x-axis to show values in millions with 'M' suffix
def millions_formatter(x, pos):
    return '{:.1f}M'.format(x / 1e6)

ax.xaxis.set_major_formatter(ticker.FuncFormatter(millions_formatter))

ax.tick_params(axis='both', which='major', labelsize=14)
ax.grid(True, linestyle='--', alpha=0.7, zorder=0)

fig.text(0.5, 0.02, 
         'Total reads: {:,} ({:,} circRNAs + {:,} linear) | Mean ± SD across runs'.format(
             ${total_reads}, ${circrna_reads}, ${linear_reads}), 
         ha='center', fontsize=11, style='italic')

# Ideal corner annotation
ax.annotate(
    'Ideal\\n(high throughput\\nlow memory)',
    xy=(0.92, 0.12),
    xycoords='axes fraction',
    bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", ec="gold", alpha=0.8),
    fontsize=12,
    ha='center',
    va='bottom'
)

# Arrow pointing to the corner
ax.annotate(
    '',
    xy=(0.97, 0.05),
    xytext=(0.75, 0.25),
    xycoords='axes fraction',
    arrowprops=dict(
        arrowstyle='->',
        lw=2.5,
        color='black'
    ),
    zorder=5
)

# Set limits with some padding
x_min = 0
x_max_with_err = (result_df['Throughput_Mean'] + result_df['Throughput_Std']).max()
x_max = x_max_with_err * 1.15
ax.set_xlim(x_min, x_max)

# Memory axis
memory_ratio = result_df['MaxCPU_GB_Mean'].max() / result_df['MaxCPU_GB_Mean'].min()
if memory_ratio > 10:
    ax.set_yscale('log')
    ax.set_ylabel('Memory Usage (GB, log scale)', fontsize=16)
    ax.grid(True, which="both", linestyle='--', alpha=0.6, zorder=0)
    y_min = result_df['MaxCPU_GB_Mean'].min() * 0.5
    y_max = (result_df['MaxCPU_GB_Mean'] + result_df['MaxCPU_GB_Std']).max() * 1.5
    ax.set_ylim(y_min, y_max)
else:
    y_min_with_err = (result_df['MaxCPU_GB_Mean'] - result_df['MaxCPU_GB_Std']).min()
    y_max_with_err = (result_df['MaxCPU_GB_Mean'] + result_df['MaxCPU_GB_Std']).max()
    y_min = max(0, y_min_with_err * 0.9)
    y_max = y_max_with_err * 1.15
    ax.set_ylim(y_min, y_max)

plt.tight_layout(rect=[0, 0.04, 1, 1])
plt.savefig('throughput_vs_memory_performance.png', dpi=300, bbox_inches='tight')
plt.close()

# Create summary text
with open('throughput_summary.txt', 'w') as f:
    f.write("CircRNA Analysis Tool Throughput Summary\\n")
    f.write("=" * 60 + "\\n\\n")
    f.write("Dataset Information:\\n")
    f.write("  Total reads: {:,}\\n".format(${total_reads}))
    f.write("  - circRNAs: {:,}\\n".format(${circrna_reads}))
    f.write("  - Linear RNAs: {:,}\\n".format(${linear_reads}))
    f.write("  CPU cores per tool: {}\\n\\n".format(${threads}))
    f.write("Throughput Metrics (mean ± standard deviation):\\n")
    f.write("-" * 60 + "\\n\\n")
    
    for _, row in result_df.iterrows():
        f.write("{}:\\n".format(row['DisplayName']))
        f.write("  Number of runs: {}\\n".format(int(row['Run_Count'])))
        f.write("  Throughput: {:.0f} ± {:.0f} reads/hour ({:.2f}M ± {:.2f}M)\\n".format(
            row['Throughput_Mean'], row['Throughput_Std'],
            row['Throughput_Mean']/1e6, row['Throughput_Std']/1e6))
        f.write("  Memory Usage: {:.2f} ± {:.2f} GB\\n".format(
            row['MaxCPU_GB_Mean'], row['MaxCPU_GB_Std']))
        
        throughput_cv = (row['Throughput_Std'] / row['Throughput_Mean'] * 100) if row['Throughput_Mean'] > 0 else 0
        memory_cv = (row['MaxCPU_GB_Std'] / row['MaxCPU_GB_Mean'] * 100) if row['MaxCPU_GB_Mean'] > 0 else 0
        
        f.write("  Throughput variability (CV): {:.1f}%\\n".format(throughput_cv))
        f.write("  Memory variability (CV): {:.1f}%\\n\\n".format(memory_cv))

print("\\nThroughput visualization complete!")
print("Plot saved: throughput_vs_memory_performance.png")
print("Summary saved: throughput_summary.txt")
EOF
    """
}