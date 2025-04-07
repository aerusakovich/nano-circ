#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr
#SBATCH --job-name=circRNA_benchmark
#SBATCH --output=circRNA_benchmark_%j.out
#SBATCH --cpus-per-task=4 --mem=16G

# Complete circRNA Tools Benchmark Script
# This script evaluates the performance of CIRI-long, isocirc, and circnick
# against ground truth data using reciprocal overlap with and without -split option.

set -e

# =====================================================
# CONFIGURATION
# =====================================================

# Tool output paths
ISOCIRC="/scratch/aerusakovich/sim_ciri_long_jobim/combined/isocirc_output/isocirc.bed"
CIRI_LONG="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/ciri-long/ciri_long_output/CIRI-long.cand_circ.bed12"
CIRCNICK="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/circnick/combined/comprehensive_circrna.bed12"

# Ground truth path - directly using the original circRNAs file
BENCHMARK_DIR="/scratch/aerusakovich/benchmark"
GROUND_TRUTH="/scratch/aerusakovich/sim_ciri_long_jobim/combined/circRNAs.bed"

# Output directory
RESULTS_DIR="${BENCHMARK_DIR}/results_benchmark_split_comp_$(date +%Y%m%d)"

# Define overlap fractions to test
FRACTIONS=("0.25" "0.5" "0.75")

# =====================================================
# SETUP
# =====================================================

# Create directory structure
mkdir -p ${RESULTS_DIR}/{intersections,metrics,plots,ground_truth,debug}/{split,nosplit}

# Function to print section header
print_section() {
    echo
    echo "============================================================"
    echo "  $1"
    echo "============================================================"
    echo
}

# Check for required tools
check_dependencies() {
    print_section "CHECKING DEPENDENCIES"
    
    if ! command -v bedtools &> /dev/null; then
        echo "Error: bedtools is not installed or not in PATH."
        echo "Please install bedtools using: mamba install -c bioconda bedtools"
        exit 1
    fi
    
    if ! command -v R &> /dev/null; then
        echo "Error: R is not installed or not in PATH."
        echo "Please install R using: mamba install -c conda-forge r-base r-essentials"
        exit 1
    fi
    
    echo "All dependencies are installed."
}

# Create R package check script
create_r_check_script() {
    cat > "${RESULTS_DIR}/check_packages.R" << 'EOFR'
#!/usr/bin/env Rscript

required_packages <- c("ggplot2", "dplyr", "tidyr", "gridExtra", "RColorBrewer", "patchwork", "scales")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing R packages:", paste(missing_packages, collapse=", "), "\n")
  install.packages(missing_packages, repos="https://cloud.r-project.org")
} else {
  cat("All required R packages are installed.\n")
}

# Make sure they're all loaded properly
for (pkg in required_packages) {
  library(pkg, character.only = TRUE)
}
EOFR

    chmod +x "${RESULTS_DIR}/check_packages.R"
}

# =====================================================
# GROUND TRUTH PROCESSING
# =====================================================

# Verify ground truth file
verify_ground_truth() {
    print_section "VERIFYING GROUND TRUTH DATA"
    
    echo "Ground truth circRNAs file: ${GROUND_TRUTH}"
    
    if [ ! -f "${GROUND_TRUTH}" ]; then
        echo "ERROR: Ground truth file not found: ${GROUND_TRUTH}"
        exit 1
    fi
    
    # Count original circRNAs
    GROUND_TRUTH_COUNT=$(wc -l < "${GROUND_TRUTH}")
    echo "Total ground truth circRNAs: ${GROUND_TRUTH_COUNT}"
    
    # Copy to results directory for reference
    cp "${GROUND_TRUTH}" "${RESULTS_DIR}/ground_truth/ground_truth.bed"
    
    # Check tool files existence
    check_tool_files
}

# Check if tool files exist
check_tool_files() {
    print_section "CHECKING TOOL FILES"
    
    echo "CIRI-long file: ${CIRI_LONG}"
    if [ ! -f "${CIRI_LONG}" ]; then
        echo "WARNING: CIRI-long file not found. This tool will be skipped."
    fi
    
    echo "isocirc file: ${ISOCIRC}"
    if [ ! -f "${ISOCIRC}" ]; then
        echo "WARNING: isocirc file not found. This tool will be skipped."
    fi
    
    echo "circnick file: ${CIRCNICK}"
    if [ ! -f "${CIRCNICK}" ]; then
        echo "WARNING: circnick file not found. This tool will be skipped."
    fi
}

# =====================================================
# BENCHMARKING FUNCTIONS
# =====================================================

# Run benchmark for a specific tool, fraction and split option
run_benchmark() {
    local tool_name=$1
    local tool_bed=$2
    local fraction=$3
    local use_split=$4  # "true" or "false"
    
    # Set split option suffix for directory structure
    local split_dir="split"
    local split_suffix="with_split"
    local split_option="-split"
    
    if [ "${use_split}" = "false" ]; then
        split_dir="nosplit"
        split_suffix="no_split"
        split_option=""
    fi
    
    print_section "RUNNING BENCHMARK FOR ${tool_name} (f=${fraction}, ${split_suffix})"
    
    # Prepare output directories
    local out_dir="${RESULTS_DIR}/intersections/${split_dir}/${tool_name}_f${fraction}"
    mkdir -p "${out_dir}"
    
    # Extract unique circRNA coordinates from tool predictions
    echo "Extracting unique circRNA coordinates from ${tool_name} predictions..."
    local unique_tool_circrnas="${out_dir}/unique_circrnas.bed"
    cut -f 1-3 "${tool_bed}" | sort -k1,1 -k2,2n -k3,3n | uniq > "${unique_tool_circrnas}"
    local total_predictions=$(wc -l < "${unique_tool_circrnas}")
    echo "Total unique circRNA predictions by ${tool_name}: ${total_predictions}"
    
    # Find true positives - circRNAs correctly identified by the tool
    echo "Finding true positives..."
    local tp_file="${out_dir}/TP.bed"
    if [ "${use_split}" = "true" ]; then
        bedtools intersect -a "${unique_tool_circrnas}" -b "${GROUND_TRUTH}" -f ${fraction} -r -wa -u -split > "${tp_file}"
    else
        bedtools intersect -a "${unique_tool_circrnas}" -b "${GROUND_TRUTH}" -f ${fraction} -r -wa -u > "${tp_file}"
    fi
    local TP=$(wc -l < "${tp_file}")
    echo "True positives (TP): ${TP}"
    
    # Find false positives - predicted circRNAs that don't match any ground truth
    echo "Finding false positives..."
    local fp_file="${out_dir}/FP.bed"
    if [ "${use_split}" = "true" ]; then
        bedtools intersect -a "${unique_tool_circrnas}" -b "${GROUND_TRUTH}" -f ${fraction} -r -v -split > "${fp_file}"
    else
        bedtools intersect -a "${unique_tool_circrnas}" -b "${GROUND_TRUTH}" -f ${fraction} -r -v > "${fp_file}"
    fi
    local FP=$(wc -l < "${fp_file}")
    echo "False positives (FP): ${FP}"
    
    # Find false negatives - ground truth circRNAs missed by the tool
    echo "Finding false negatives..."
    local fn_file="${out_dir}/FN.bed"
    if [ "${use_split}" = "true" ]; then
        bedtools intersect -a "${GROUND_TRUTH}" -b "${unique_tool_circrnas}" -f ${fraction} -r -v -split > "${fn_file}"
    else
        bedtools intersect -a "${GROUND_TRUTH}" -b "${unique_tool_circrnas}" -f ${fraction} -r -v > "${fn_file}"
    fi
    local FN=$(wc -l < "${fn_file}")
    echo "False negatives (FN): ${FN}"
    
    # Verification of counts
    echo "Verifying counts:"
    echo "TP(${TP}) + FP(${FP}) = $((TP + FP)), which should equal total predictions (${total_predictions})"
    if [ $((TP + FP)) -ne ${total_predictions} ]; then
        echo "WARNING: TP + FP doesn't match total predictions. This might indicate an issue with the analysis."
    fi
    
    local ground_truth_count=$(wc -l < "${GROUND_TRUTH}")
    echo "TP(${TP}) + FN(${FN}) = $((TP + FN)), which should equal ground truth count (${ground_truth_count})"
    if [ $((TP + FN)) -ne ${ground_truth_count} ]; then
        echo "WARNING: TP + FN doesn't match ground truth count. This might indicate an issue with the analysis."
    fi
    
    # Calculate metrics
    calculate_metrics "${tool_name}" "${fraction}" "${use_split}" ${TP} ${FP} ${FN}
}

# Calculate performance metrics
calculate_metrics() {
    local tool_name=$1
    local fraction=$2
    local use_split=$3
    local TP=$4
    local FP=$5
    local FN=$6
    
    # Set split option suffix for directory structure
    local split_dir="split"
    
    if [ "${use_split}" = "false" ]; then
        split_dir="nosplit"
    fi
    
    # Calculate precision
    local precision=0
    if [ $((TP + FP)) -ne 0 ]; then
        precision=$(echo "scale=4; ${TP} / (${TP} + ${FP})" | bc)
    fi
    
    # Calculate recall (sensitivity)
    local recall=0
    if [ $((TP + FN)) -ne 0 ]; then
        recall=$(echo "scale=4; ${TP} / (${TP} + ${FN})" | bc)
    fi
    
    # Calculate F1 score
    local f1_score=0
    if [ $(echo "${precision} + ${recall} > 0" | bc) -eq 1 ]; then
        f1_score=$(echo "scale=4; 2 * (${precision} * ${recall}) / (${precision} + ${recall})" | bc)
    fi
    
    # Print metrics
    echo "Precision: ${precision}"
    echo "Recall: ${recall}"
    echo "F1 Score: ${f1_score}"
    
    # Save to tool-specific metrics file
    mkdir -p "${RESULTS_DIR}/metrics/${split_dir}"
    echo -e "${tool_name}\t${fraction}\t${TP}\t${FP}\t${FN}\t${precision}\t${recall}\t${f1_score}" >> "${RESULTS_DIR}/metrics/${split_dir}/${tool_name}.tsv"
    
    # Append to combined metrics and contingency files
    echo -e "${tool_name}\t${fraction}\t${use_split}\t${TP}\t${FP}\t${FN}\t${precision}\t${recall}\t${f1_score}" >> "${RESULTS_DIR}/all_metrics.tsv"
    echo -e "${tool_name}\t${fraction}\t${use_split}\t${TP}\t${FP}\t${FN}" >> "${RESULTS_DIR}/contingency_table.tsv"
}

# =====================================================
# VISUALIZATION SCRIPTS
# =====================================================

# Create visualization scripts
create_visualization_scripts() {
    print_section "CREATING VISUALIZATION SCRIPTS"
    
    # Create metrics comparison dot plot script
    cat > "${RESULTS_DIR}/plot_metrics_comparison.R" << 'EOF_METRICS'
#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(scales)
})

# Define output directory
output_dir <- commandArgs(trailingOnly = TRUE)[1]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat(paste0("Output will be saved to: ", output_dir, "\n"))

# Read metrics data
metrics_data <- read.delim("all_metrics.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Convert columns to proper types
metrics_data$Tool <- factor(metrics_data$Tool, levels=c("cirilong", "isocirc", "circnick"))
metrics_data$Fraction <- as.factor(metrics_data$Fraction)
metrics_data$use_split <- as.logical(metrics_data$use_split)

# Create nice tool names
metrics_data <- metrics_data %>%
  mutate(ToolName = case_when(
    Tool == "cirilong" ~ "CIRI-long",
    Tool == "isocirc" ~ "IsoCirc",
    Tool == "circnick" ~ "CircNick-LRS",
    TRUE ~ Tool
  ),
  SplitOption = ifelse(use_split, "Exon-aware", "Not exon-aware"))

# Create tool colors
tool_colors <- c(
  "CIRI-long" = "#0072B2",   # Blue
  "IsoCirc" = "#009E73",     # Green
  "CircNick-LRS" = "#D55E00" # Orange-red
)

# Create summary table
metrics_summary <- metrics_data %>%
  group_by(Tool, ToolName, Fraction, use_split, SplitOption) %>%
  summarize(
    Precision = mean(Precision, na.rm = TRUE),
    Recall = mean(Recall, na.rm = TRUE),
    F1_Score = mean(F1_Score, na.rm = TRUE),
    .groups = "drop"
  )

# Create dot plot data
dot_plot_data <- metrics_summary %>%
  select(ToolName, Fraction, use_split, SplitOption, Precision, Recall, F1_Score) %>%
  pivot_longer(cols = c(Precision, Recall, F1_Score),
               names_to = "Metric", values_to = "Value")

# Make sure Metric is in the correct order
dot_plot_data$Metric <- factor(dot_plot_data$Metric, 
                              levels = c("Precision", "Recall", "F1_Score"))

# Create dot plot with split/no-split comparison using alpha
split_comparison_plot <- ggplot(dot_plot_data, 
                  aes(x = Metric, y = Value, color = ToolName, shape = Fraction, alpha = SplitOption)) +
  geom_point(size = 4, stroke = 1.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = tool_colors) +
  scale_shape_manual(values = c("0.25" = 16, "0.5" = 17, "0.75" = 15)) +
  scale_alpha_manual(values = c("Exon-aware" = 1.0, "Not exon-aware" = 0.3)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Performance Metrics Comparison: Exon-aware vs. Not exon-aware",
    subtitle = "Points with lower opacity = Not exon-aware",
    x = "Metric",
    y = "Value",
    color = "Tool",
    shape = "Overlap Fraction",
    alpha = "Exon Awareness"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Save the split comparison dot plot
ggsave(file.path(output_dir, "split_nosplit_comparison.png"), split_comparison_plot, width = 12, height = 8, dpi = 300)

# Create faceted version to see split/no-split side by side
facet_plot <- ggplot(dot_plot_data, 
                   aes(x = Metric, y = Value, color = ToolName, shape = Fraction)) +
  geom_point(size = 4, stroke = 1.5, position = position_dodge(width = 0.6)) +
  facet_wrap(~ SplitOption) +
  scale_color_manual(values = tool_colors) +
  scale_shape_manual(values = c("0.25" = 16, "0.5" = 17, "0.75" = 15)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Performance Metrics: Exon-aware vs. Not exon-aware",
    x = "Metric",
    y = "Value",
    color = "Tool",
    shape = "Overlap Fraction"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold")
  )

# Save the faceted split/no-split comparison
ggsave(file.path(output_dir, "split_nosplit_facet.png"), facet_plot, width = 12, height = 8, dpi = 300)

# Create a delta plot (showing the difference between split and no-split)
delta_data <- metrics_summary %>%
  select(Tool, ToolName, Fraction, use_split, Precision, Recall, F1_Score) %>%
  pivot_wider(names_from = use_split, 
              values_from = c(Precision, Recall, F1_Score),
              names_glue = "{.value}_{use_split}") %>%
  mutate(
    Precision_diff = Precision_TRUE - Precision_FALSE,
    Recall_diff = Recall_TRUE - Recall_FALSE,
    F1_Score_diff = F1_Score_TRUE - F1_Score_FALSE
  ) %>%
  select(Tool, ToolName, Fraction, Precision_diff, Recall_diff, F1_Score_diff) %>%
  pivot_longer(cols = c(Precision_diff, Recall_diff, F1_Score_diff),
               names_to = "Metric", values_to = "Difference") %>%
  mutate(Metric = gsub("_diff", "", Metric))

# Create delta plot
delta_plot <- ggplot(delta_data, 
                    aes(x = Metric, y = Difference, fill = ToolName)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  facet_wrap(~ Fraction, labeller = labeller(Fraction = function(x) paste0("Overlap Fraction = ", x))) +
  scale_fill_manual(values = tool_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Impact of Exon Awareness on Performance Metrics",
    subtitle = "Difference between exon-aware and non-exon-aware (positive = exon-aware is better)",
    x = "Metric",
    y = "Difference (exon-aware - not exon-aware)",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12, face = "bold")
  )

# Save the delta plot
ggsave(file.path(output_dir, "split_impact_delta.png"), delta_plot, width = 12, height = 8, dpi = 300)

# Write summary tables to CSV
write.csv(metrics_summary, file.path(output_dir, "split_comparison_metrics.csv"), row.names = FALSE)
write.csv(delta_data, file.path(output_dir, "split_impact_metrics.csv"), row.names = FALSE)

cat("Split/No-split comparison visualization completed.\n")
EOF_METRICS

    chmod +x "${RESULTS_DIR}/plot_metrics_comparison.R"

    # Create contingency table comparison script
    cat > "${RESULTS_DIR}/plot_contingency_comparison.R" << 'EOF_CONTINGENCY'
#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(scales)
})

# Define output directory
output_dir <- commandArgs(trailingOnly = TRUE)[1]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat(paste0("Output will be saved to: ", output_dir, "\n"))

# Read contingency table data
contingency_data <- read.delim("contingency_table.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Convert columns to proper types
contingency_data$Tool <- factor(contingency_data$Tool, levels=c("cirilong", "isocirc", "circnick"))
contingency_data$Fraction <- as.factor(contingency_data$Fraction)
contingency_data$use_split <- as.logical(contingency_data$use_split)

# Create nice tool names
contingency_data <- contingency_data %>%
  mutate(ToolName = case_when(
    Tool == "cirilong" ~ "CIRI-long",
    Tool == "isocirc" ~ "IsoCirc",
    Tool == "circnick" ~ "CircNick-LRS",
    TRUE ~ Tool
  ),
  SplitOption = ifelse(use_split, "With -split", "Without -split"))

# Create summary table
summary_table <- contingency_data %>%
  group_by(Tool, ToolName, Fraction, use_split, SplitOption) %>%
  summarize(
    TP = sum(TP),
    FP = sum(FP),
    FN = sum(FN),
    .groups = "drop"
  )

# Create heatmap data
heatmap_data <- summary_table %>%
  select(ToolName, Fraction, SplitOption, TP, FP, FN) %>%
  pivot_longer(cols = c(TP, FP, FN), names_to = "Metric", values_to = "Count") %>%
  arrange(ToolName, Fraction, SplitOption, Metric)

# Make sure Metric is in the correct order
heatmap_data$Metric <- factor(heatmap_data$Metric, levels = c("TP", "FP", "FN"))

# Create a contingency heatmap faceted by split option
contingency_heatmap <- ggplot(heatmap_data, 
                             aes(x = Fraction, y = ToolName, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_grid(SplitOption ~ Metric) +
  scale_fill_gradient(low = "#FFF5F0", high = "#67000D",
                    trans = "log10",
                    breaks = c(10, 100, 1000, 10000),
                    labels = comma) +
  geom_text(aes(label = format(Count, big.mark = ",")), color = "black", size = 3) +
  labs(
    title = "circRNA Detection Performance: Exon-aware vs. Not exon-aware",
    subtitle = "TP = True Positives, FP = False Positives, FN = False Negatives",
    x = "Reciprocal Overlap Fraction",
    y = "Tool",
    fill = "Count (log scale)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold")
  )

# Save the contingency heatmap
ggsave(file.path(output_dir, "contingency_heatmap_comparison.png"), contingency_heatmap, width = 14, height = 10, dpi = 300)

# Calculate delta contingency (differences between split and no-split)
delta_contingency <- summary_table %>%
  select(Tool, ToolName, Fraction, use_split, TP, FP, FN) %>%
  pivot_wider(names_from = use_split, 
              values_from = c(TP, FP, FN),
              names_glue = "{.value}_{use_split}") %>%
  mutate(
    TP_diff = TP_TRUE - TP_FALSE,
    FP_diff = FP_TRUE - FP_FALSE,
    FN_diff = FN_TRUE - FN_FALSE
  ) %>%
  select(Tool, ToolName, Fraction, TP_diff, FP_diff, FN_diff) %>%
  pivot_longer(cols = c(TP_diff, FP_diff, FN_diff),
               names_to = "Metric", values_to = "Difference") %>%
  mutate(Metric = gsub("_diff", "", Metric))

# Create delta contingency bar chart
delta_contingency_chart <- ggplot(delta_contingency, 
                                aes(x = Metric, y = Difference, fill = ToolName)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  facet_wrap(~ Fraction, scales = "free_y") +
  scale_fill_manual(values = c(
    "CIRI-long" = "#0072B2", 
    "IsoCirc" = "#009E73", 
    "CircNick-LRS" = "#D55E00"
  )) +
  labs(
    title = "Impact of Exon Awareness on Contingency Counts",
    subtitle = "Difference between exon-aware and non-exon-aware counts",
    x = "Metric",
    y = "Difference (exon-aware - not exon-aware)",
    fill = "Tool"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 12, face = "bold")
  )

# Save the delta contingency chart
ggsave(file.path(output_dir, "contingency_delta_chart.png"), delta_contingency_chart, width = 12, height = 8, dpi = 300)

# Write summary tables to CSV
write.csv(summary_table, file.path(output_dir, "contingency_split_comparison.csv"), row.names = FALSE)
write.csv(delta_contingency, file.path(output_dir, "contingency_delta.csv"), row.names = FALSE)

cat("Contingency comparison visualization completed.\n")
EOF_CONTINGENCY

    chmod +x "${RESULTS_DIR}/plot_contingency_comparison.R"
}

# Run visualization
run_visualization() {
    print_section "RUNNING VISUALIZATION"
    
    cd "${RESULTS_DIR}"
    
    # Check R packages
    echo "Checking required R packages..."
    Rscript "${RESULTS_DIR}/check_packages.R"
    
    # Run visualization scripts
    echo "Running metrics comparison visualization..."
    Rscript "${RESULTS_DIR}/plot_metrics_comparison.R" "${RESULTS_DIR}/plots"
    
    echo "Running contingency comparison visualization..."
    Rscript "${RESULTS_DIR}/plot_contingency_comparison.R" "${RESULTS_DIR}/plots"
    
    echo "Visualizations completed."
}

# =====================================================
# MAIN EXECUTION
# =====================================================

# Create metrics headers
create_metrics_headers() {
    echo -e "Tool\tFraction\tuse_split\tTP\tFP\tFN\tPrecision\tRecall\tF1_Score" > "${RESULTS_DIR}/all_metrics.tsv"
    echo -e "Tool\tFraction\tuse_split\tTP\tFP\tFN" > "${RESULTS_DIR}/contingency_table.tsv"
}

# Main function
main() {
    # Setup
    check_dependencies
    create_r_check_script
    create_metrics_headers
    
    # Verify ground truth
    verify_ground_truth
    
    # Create visualization scripts
    create_visualization_scripts
    
    # Run benchmarks for each tool, fraction, and split option
    for FRAC in "${FRACTIONS[@]}"; do
        echo "Using overlap fraction: ${FRAC}"
        
        if [ -f "${CIRI_LONG}" ]; then
            # Run with -split option
            run_benchmark "cirilong" "${CIRI_LONG}" "${FRAC}" "true"
            # Run without -split option
            run_benchmark "cirilong" "${CIRI_LONG}" "${FRAC}" "false"
        fi
        
        if [ -f "${ISOCIRC}" ]; then
            # Run with -split option
            run_benchmark "isocirc" "${ISOCIRC}" "${FRAC}" "true"
            # Run without -split option
            run_benchmark "isocirc" "${ISOCIRC}" "${FRAC}" "false"
        fi
        
        if [ -f "${CIRCNICK}" ]; then
            # Run with -split option
            run_benchmark "circnick" "${CIRCNICK}" "${FRAC}" "true"
            # Run without -split option
            run_benchmark "circnick" "${CIRCNICK}" "${FRAC}" "false"
        fi
    done
    
    # Generate visualizations
    run_visualization
    
    # Print summary
    print_section "BENCHMARK SUMMARY"
    echo "Results are available in: ${RESULTS_DIR}"
    echo
    echo "Ground truth information:"
    echo "- Original circRNAs (ground truth): ${GROUND_TRUTH} (${GROUND_TRUTH_COUNT} entries)"
    echo
    echo "Key visualizations:"
    echo "- Exon-aware vs Not exon-aware comparison (alpha): ${RESULTS_DIR}/plots/split_nosplit_comparison.png"
    echo "- Exon-aware vs Not exon-aware comparison (faceted): ${RESULTS_DIR}/plots/split_nosplit_facet.png"
    echo "- Exon awareness impact delta: ${RESULTS_DIR}/plots/split_impact_delta.png"
    echo "- Contingency heatmap comparison: ${RESULTS_DIR}/plots/contingency_heatmap_comparison.png"
    echo "- Contingency delta chart: ${RESULTS_DIR}/plots/contingency_delta_chart.png"
    echo
    echo "Benchmark completed successfully at $(date)"
}

# Execute main function
main