#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr
#SBATCH --job-name=combinatorial_circrna_benchmark
#SBATCH --output=combinatorial_circrna_benchmark_%j.out
#SBATCH --cpus-per-task=4 --mem=16G

# Combinatorial circRNA Tools Benchmark Script
# This script evaluates the performance of tool combinations:
# CIRI-long+IsoCirc, CIRI-long+CircNick, IsoCirc+CircNick, and union of all tools
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
RESULTS_DIR="${BENCHMARK_DIR}/results_combinatorial_benchmark_$(date +%Y%m%d)"

# Define overlap fractions to test
FRACTIONS=("0.25" "0.5" "0.75")

# =====================================================
# SETUP
# =====================================================

# Create directory structure
mkdir -p ${RESULTS_DIR}/{intersections,metrics,plots,ground_truth,debug,combined_tools}/{split,nosplit}

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
# GROUND TRUTH AND TOOL PROCESSING
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
        echo "ERROR: CIRI-long file not found. Required for combinations."
        exit 1
    fi
    
    echo "isocirc file: ${ISOCIRC}"
    if [ ! -f "${ISOCIRC}" ]; then
        echo "ERROR: isocirc file not found. Required for combinations."
        exit 1
    fi
    
    echo "circnick file: ${CIRCNICK}"
    if [ ! -f "${CIRCNICK}" ]; then
        echo "ERROR: circnick file not found. Required for combinations."
        exit 1
    fi
    
    # All files must exist for combinations
    echo "All tool files found. Proceeding with combinations."
}

# Extract unique circRNA coordinates from all tool predictions
prepare_individual_tools() {
    print_section "PREPARING INDIVIDUAL TOOL PREDICTIONS"
    
    # Create directory for individual tool BED files
    mkdir -p "${RESULTS_DIR}/tools_individual"
    
    echo "Extracting unique circRNA coordinates from CIRI-long predictions..."
    cut -f 1-3 "${CIRI_LONG}" | sort -k1,1 -k2,2n -k3,3n | uniq > "${RESULTS_DIR}/tools_individual/cirilong.bed"
    ciri_count=$(wc -l < "${RESULTS_DIR}/tools_individual/cirilong.bed")
    echo "CIRI-long unique circRNA predictions: ${ciri_count}"
    
    echo "Extracting unique circRNA coordinates from IsoCirc predictions..."
    cut -f 1-3 "${ISOCIRC}" | sort -k1,1 -k2,2n -k3,3n | uniq > "${RESULTS_DIR}/tools_individual/isocirc.bed"
    iso_count=$(wc -l < "${RESULTS_DIR}/tools_individual/isocirc.bed")
    echo "IsoCirc unique circRNA predictions: ${iso_count}"
    
    echo "Extracting unique circRNA coordinates from CircNick predictions..."
    cut -f 1-3 "${CIRCNICK}" | sort -k1,1 -k2,2n -k3,3n | uniq > "${RESULTS_DIR}/tools_individual/circnick.bed"
    nick_count=$(wc -l < "${RESULTS_DIR}/tools_individual/circnick.bed")
    echo "CircNick unique circRNA predictions: ${nick_count}"
}

# Generate tool combinations
generate_tool_combinations() {
    print_section "GENERATING TOOL COMBINATIONS"
    
    ciri="${RESULTS_DIR}/tools_individual/cirilong.bed"
    iso="${RESULTS_DIR}/tools_individual/isocirc.bed"
    nick="${RESULTS_DIR}/tools_individual/circnick.bed"
    
    mkdir -p "${RESULTS_DIR}/tools_combined"
    
    # CIRI-long + IsoCirc
    echo "Creating CIRI-long + IsoCirc combination..."
    bedtools intersect -a "${ciri}" -b "${iso}" -wa | sort -k1,1 -k2,2n -k3,3n | uniq > "${RESULTS_DIR}/tools_combined/ciri_iso.bed"
    ciri_iso_count=$(wc -l < "${RESULTS_DIR}/tools_combined/ciri_iso.bed")
    echo "CIRI-long + IsoCirc common circRNAs: ${ciri_iso_count}"
    
    # CIRI-long + CircNick
    echo "Creating CIRI-long + CircNick combination..."
    bedtools intersect -a "${ciri}" -b "${nick}" -wa | sort -k1,1 -k2,2n -k3,3n | uniq > "${RESULTS_DIR}/tools_combined/ciri_nick.bed"
    ciri_nick_count=$(wc -l < "${RESULTS_DIR}/tools_combined/ciri_nick.bed")
    echo "CIRI-long + CircNick common circRNAs: ${ciri_nick_count}"
    
    # IsoCirc + CircNick
    echo "Creating IsoCirc + CircNick combination..."
    bedtools intersect -a "${iso}" -b "${nick}" -wa | sort -k1,1 -k2,2n -k3,3n | uniq > "${RESULTS_DIR}/tools_combined/iso_nick.bed"
    iso_nick_count=$(wc -l < "${RESULTS_DIR}/tools_combined/iso_nick.bed")
    echo "IsoCirc + CircNick common circRNAs: ${iso_nick_count}"
    
    # Create the union of all tools
    echo "Creating union of all tool predictions..."
    cat "${ciri}" "${iso}" "${nick}" | sort -k1,1 -k2,2n -k3,3n | uniq > "${RESULTS_DIR}/tools_combined/all_union.bed"
    union_count=$(wc -l < "${RESULTS_DIR}/tools_combined/all_union.bed")
    echo "Union of all tool predictions: ${union_count}"
}

# =====================================================
# BENCHMARKING FUNCTIONS
# =====================================================

# Run benchmark for a specific tool combination, fraction and split option
run_benchmark() {
    local comb_name=$1
    local comb_bed=$2
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
    
    print_section "RUNNING BENCHMARK FOR ${comb_name} (f=${fraction}, ${split_suffix})"
    
    # Prepare output directories
    local out_dir="${RESULTS_DIR}/intersections/${split_dir}/${comb_name}_f${fraction}"
    mkdir -p "${out_dir}"
    
    echo "Using combined predictions from: ${comb_bed}"
    local total_predictions=$(wc -l < "${comb_bed}")
    echo "Total circRNA predictions in combination: ${total_predictions}"
    
    # Find true positives - circRNAs correctly identified by the combination
    echo "Finding true positives..."
    local tp_file="${out_dir}/TP.bed"
    if [ "${use_split}" = "true" ]; then
        bedtools intersect -a "${comb_bed}" -b "${GROUND_TRUTH}" -f ${fraction} -r -wa -u -split > "${tp_file}"
    else
        bedtools intersect -a "${comb_bed}" -b "${GROUND_TRUTH}" -f ${fraction} -r -wa -u > "${tp_file}"
    fi
    local TP=$(wc -l < "${tp_file}")
    echo "True positives (TP): ${TP}"
    
    # Find false positives - predicted circRNAs that don't match any ground truth
    echo "Finding false positives..."
    local fp_file="${out_dir}/FP.bed"
    if [ "${use_split}" = "true" ]; then
        bedtools intersect -a "${comb_bed}" -b "${GROUND_TRUTH}" -f ${fraction} -r -v -split > "${fp_file}"
    else
        bedtools intersect -a "${comb_bed}" -b "${GROUND_TRUTH}" -f ${fraction} -r -v > "${fp_file}"
    fi
    local FP=$(wc -l < "${fp_file}")
    echo "False positives (FP): ${FP}"
    
    # Find false negatives - ground truth circRNAs missed by the combination
    echo "Finding false negatives..."
    local fn_file="${out_dir}/FN.bed"
    if [ "${use_split}" = "true" ]; then
        bedtools intersect -a "${GROUND_TRUTH}" -b "${comb_bed}" -f ${fraction} -r -v -split > "${fn_file}"
    else
        bedtools intersect -a "${GROUND_TRUTH}" -b "${comb_bed}" -f ${fraction} -r -v > "${fn_file}"
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
    calculate_metrics "${comb_name}" "${fraction}" "${use_split}" ${TP} ${FP} ${FN}
}

# Calculate performance metrics
calculate_metrics() {
    local comb_name=$1
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
    
    # Save to combination-specific metrics file
    mkdir -p "${RESULTS_DIR}/metrics/${split_dir}"
    echo -e "${comb_name}\t${fraction}\t${TP}\t${FP}\t${FN}\t${precision}\t${recall}\t${f1_score}" >> "${RESULTS_DIR}/metrics/${split_dir}/${comb_name}.tsv"
    
    # Append to combined metrics file
    echo -e "${comb_name}\t${fraction}\t${use_split}\t${TP}\t${FP}\t${FN}\t${precision}\t${recall}\t${f1_score}" >> "${RESULTS_DIR}/all_metrics.tsv"
}

# =====================================================
# VISUALIZATION SCRIPTS
# =====================================================

# Create visualization scripts
create_visualization_scripts() {
    print_section "CREATING VISUALIZATION SCRIPTS"
    
    # Create metrics comparison dot plot script
    cat > "${RESULTS_DIR}/plot_combination_metrics.R" << 'EOF_METRICS'
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
colnames(metrics_data) <- c("Combination", "Fraction", "use_split", "TP", "FP", "FN", "Precision", "Recall", "F1_Score")

# Convert columns to proper types
metrics_data$Combination <- factor(metrics_data$Combination)
metrics_data$Fraction <- as.factor(metrics_data$Fraction)
metrics_data$use_split <- as.logical(metrics_data$use_split)

# Create nice combination names
metrics_data <- metrics_data %>%
  mutate(CombName = case_when(
    Combination == "ciri_iso" ~ "CIRI-long + IsoCirc",
    Combination == "ciri_nick" ~ "CIRI-long + CircNick",
    Combination == "iso_nick" ~ "IsoCirc + CircNick",
    Combination == "all_union" ~ "Union of All Tools",
    TRUE ~ Combination
  ),
  SplitOption = ifelse(use_split, "Exon-aware", "Not exon-aware"))

# Create color palette for combinations
comb_colors <- c(
  "CIRI-long + IsoCirc" = "#4575B4",      # Blue
  "CIRI-long + CircNick" = "#91BFDB",     # Light Blue
  "IsoCirc + CircNick" = "#FC8D59",       # Orange
  "Union of All Tools" = "#D73027"        # Red
)

# Create summary table
metrics_summary <- metrics_data %>%
  group_by(Combination, CombName, Fraction, use_split, SplitOption) %>%
  summarize(
    Precision = mean(Precision, na.rm = TRUE),
    Recall = mean(Recall, na.rm = TRUE),
    F1_Score = mean(F1_Score, na.rm = TRUE),
    .groups = "drop"
  )

# Create dot plot data
dot_plot_data <- metrics_summary %>%
  select(CombName, Fraction, use_split, SplitOption, Precision, Recall, F1_Score) %>%
  pivot_longer(cols = c(Precision, Recall, F1_Score),
               names_to = "Metric", values_to = "Value")

# Make sure Metric is in the correct order
dot_plot_data$Metric <- factor(dot_plot_data$Metric, 
                              levels = c("Precision", "Recall", "F1_Score"))

# Create dot plot with split/no-split comparison using alpha
split_comparison_plot <- ggplot(dot_plot_data, 
                  aes(x = Metric, y = Value, color = CombName, shape = Fraction, alpha = SplitOption)) +
  geom_point(size = 4, stroke = 1.5, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = comb_colors) +
  scale_shape_manual(values = c("0.25" = 16, "0.5" = 17, "0.75" = 15)) +
  scale_alpha_manual(values = c("Exon-aware" = 1.0, "Not exon-aware" = 0.3)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Performance Metrics: Tool Combinations",
    subtitle = "Points with lower opacity = Not exon-aware",
    x = "Metric",
    y = "Value",
    color = "Tool Combination",
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
ggsave(file.path(output_dir, "combinations_performance.png"), split_comparison_plot, width = 14, height = 10, dpi = 300)

# Create faceted version to see split/no-split side by side
facet_plot <- ggplot(dot_plot_data, 
                   aes(x = Metric, y = Value, color = CombName, shape = Fraction)) +
  geom_point(size = 4, stroke = 1.5, position = position_dodge(width = 0.6)) +
  facet_wrap(~ SplitOption) +
  scale_color_manual(values = comb_colors) +
  scale_shape_manual(values = c("0.25" = 16, "0.5" = 17, "0.75" = 15)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Performance Metrics: Tool Combinations",
    x = "Metric",
    y = "Value",
    color = "Tool Combination",
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
ggsave(file.path(output_dir, "combinations_split_facet.png"), facet_plot, width = 14, height = 10, dpi = 300)

# Create a plot faceted by overlap fraction
fraction_facet_plot <- ggplot(dot_plot_data, 
                            aes(x = Metric, y = Value, color = CombName, shape = SplitOption)) +
  geom_point(size = 4, stroke = 1.5, position = position_dodge(width = 0.6)) +
  facet_wrap(~ Fraction, labeller = labeller(Fraction = function(x) paste0("Overlap Fraction: ", x))) +
  scale_color_manual(values = comb_colors) +
  scale_shape_manual(values = c("Exon-aware" = 16, "Not exon-aware" = 4)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Tool Combinations Performance by Overlap Fraction",
    x = "Metric",
    y = "Value",
    color = "Tool Combination",
    shape = "Exon Awareness"
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

# Save the fraction facet plot
ggsave(file.path(output_dir, "combinations_fraction_facet.png"), fraction_facet_plot, width = 14, height = 10, dpi = 300)

# Create a plot faceted by metric
metric_facet_plot <- ggplot(dot_plot_data, 
                          aes(x = Fraction, y = Value, color = CombName, shape = SplitOption, group = interaction(CombName, SplitOption))) +
  geom_point(size = 3.5) +
  geom_line(position = position_dodge(width = 0.3)) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_color_manual(values = comb_colors) +
  scale_shape_manual(values = c("Exon-aware" = 16, "Not exon-aware" = 4)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Tool Combinations Performance by Metric",
    x = "Overlap Fraction",
    y = "Value",
    color = "Tool Combination",
    shape = "Exon Awareness"
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

# Save the metric facet plot
ggsave(file.path(output_dir, "combinations_metric_facet.png"), metric_facet_plot, width = 14, height = 10, dpi = 300)

# Write summary table to CSV
write.csv(metrics_summary, file.path(output_dir, "combinations_metrics_summary.csv"), row.names = FALSE)

cat("Tool combinations visualization completed.\n")
EOF_METRICS

    chmod +x "${RESULTS_DIR}/plot_combination_metrics.R"
}

# Run visualization
run_visualization() {
    print_section "RUNNING VISUALIZATION"
    
    cd "${RESULTS_DIR}"
    
    # Check R packages
    echo "Checking required R packages..."
    Rscript "${RESULTS_DIR}/check_packages.R"
    
    # Run visualization scripts
    echo "Running tool combinations visualization..."
    Rscript "${RESULTS_DIR}/plot_combination_metrics.R" "${RESULTS_DIR}/plots"
    
    echo "Visualizations completed."
}

# =====================================================
# MAIN EXECUTION
# =====================================================

# Create metrics headers
create_metrics_headers() {
    echo -e "Combination\tFraction\tuse_split\tTP\tFP\tFN\tPrecision\tRecall\tF1_Score" > "${RESULTS_DIR}/all_metrics.tsv"
}

# Main function
main() {
    # Setup
    check_dependencies
    create_r_check_script
    create_metrics_headers
    
    # Verify ground truth
    verify_ground_truth
    
    # Prepare individual tool predictions
    prepare_individual_tools
    
    # Generate tool combinations
    generate_tool_combinations
    
    # Create visualization scripts
    create_visualization_scripts
    
    # Run benchmarks for each combination, fraction, and split option
    for FRAC in "${FRACTIONS[@]}"; do
        echo "Using overlap fraction: ${FRAC}"
        
        # CIRI-long + IsoCirc
        run_benchmark "ciri_iso" "${RESULTS_DIR}/tools_combined/ciri_iso.bed" "${FRAC}" "true"
        run_benchmark "ciri_iso" "${RESULTS_DIR}/tools_combined/ciri_iso.bed" "${FRAC}" "false"
        
        # CIRI-long + CircNick
        run_benchmark "ciri_nick" "${RESULTS_DIR}/tools_combined/ciri_nick.bed" "${FRAC}" "true"
        run_benchmark "ciri_nick" "${RESULTS_DIR}/tools_combined/ciri_nick.bed" "${FRAC}" "false"
        
        # IsoCirc + CircNick
        run_benchmark "iso_nick" "${RESULTS_DIR}/tools_combined/iso_nick.bed" "${FRAC}" "true"
        run_benchmark "iso_nick" "${RESULTS_DIR}/tools_combined/iso_nick.bed" "${FRAC}" "false"
        
        # Union of all tools
        run_benchmark "all_union" "${RESULTS_DIR}/tools_combined/all_union.bed" "${FRAC}" "true"
        run_benchmark "all_union" "${RESULTS_DIR}/tools_combined/all_union.bed" "${FRAC}" "false"
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
    echo "Tool combination information:"
    echo "- CIRI-long + IsoCirc: $(wc -l < ${RESULTS_DIR}/tools_combined/ciri_iso.bed) circRNAs"
    echo "- CIRI-long + CircNick: $(wc -l < ${RESULTS_DIR}/tools_combined/ciri_nick.bed) circRNAs"
    echo "- IsoCirc + CircNick: $(wc -l < ${RESULTS_DIR}/tools_combined/iso_nick.bed) circRNAs"
    echo "- Union of all tools: $(wc -l < ${RESULTS_DIR}/tools_combined/all_union.bed) circRNAs (combined predictions)"
    echo
    echo "Key visualizations:"
    echo "- Tool combinations performance: ${RESULTS_DIR}/plots/combinations_performance.png"
    echo "- Split/no-split faceted comparison: ${RESULTS_DIR}/plots/combinations_split_facet.png"
    echo "- Performance by fraction: ${RESULTS_DIR}/plots/combinations_fraction_facet.png"
    echo "- Performance by metric: ${RESULTS_DIR}/plots/combinations_metric_facet.png"
    echo
    echo "Benchmark completed successfully at $(date)"
}

# Execute main function
main