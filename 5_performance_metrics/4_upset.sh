#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr
#SBATCH --job-name=circrna_overlap_split
#SBATCH --output=circrna_overlap_split_%j.out
#SBATCH --cpus-per-task=4 --mem=16G

# CircRNA Isoform Comparison Script (Combined version with -split option)
# This script compares circRNA isoform predictions from three tools:
# 1. Ciri-long
# 2. isocirc
# 3. circnick
# Uses bedtools intersect with -split option and different overlap fractions and generates upset plots

# Set paths to input files
ISOCIRC="/scratch/aerusakovich/sim_ciri_long_jobim/combined/isocirc_output/isocirc.bed"
CIRI_LONG="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/ciri-long/ciri_long_output/CIRI-long.cand_circ.bed12"
CIRCNICK="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/circnick/combined/comprehensive_circrna.bed12"

# Create output directory
OUTPUT_DIR="circRNA_comparison_results_split"
mkdir -p $OUTPUT_DIR

# Create directory for bedtools results
BEDTOOLS_DIR="${OUTPUT_DIR}/bedtools_results"
mkdir -p $BEDTOOLS_DIR

# Array of overlap fractions to test
FRACTIONS=("0.95")

echo "Starting comparison of circRNA isoform predictions using -split option..."

# Print the paths for debugging
echo "Checking files:"
echo "CIRI-long: $CIRI_LONG"
echo "IsoCirc: $ISOCIRC"
echo "CircNick: $CIRCNICK"

# Skip if any of the required files don't exist
if [ ! -f "$CIRI_LONG" ]; then
    echo "ERROR: Missing CIRI-long file at $CIRI_LONG"
    exit 1
fi

if [ ! -f "$ISOCIRC" ]; then
    echo "ERROR: Missing IsoCirc file at $ISOCIRC"
    exit 1
fi

if [ ! -f "$CIRCNICK" ]; then
    echo "ERROR: Missing CircNick file at $CIRCNICK"
    exit 1
fi

# Create run-specific directories
RUN_BEDTOOLS_DIR="${BEDTOOLS_DIR}/combined"
mkdir -p $RUN_BEDTOOLS_DIR

# Step 1: Use bedtools to compare isoforms with different overlap fractions
for FRACTION in "${FRACTIONS[@]}"; do
    echo "Processing with overlap fraction: $FRACTION"
    
    # Create directory for this fraction
    FRAC_DIR="${RUN_BEDTOOLS_DIR}/f${FRACTION}"
    mkdir -p $FRAC_DIR
    
    # For BED12 format, we'll use the full files without simplification
    # This is important for the -split option to work correctly
    echo "Using complete BED files for comparison with -split option..."
    
    # Create unique identifier BED files but retain full BED12 structure for intersect
    # Add a unique ID column based on coordinates
    awk '{print $0"\t"$1"_"$2"_"$3}' $CIRI_LONG > "${FRAC_DIR}/ciri_full.bed"
    awk '{print $0"\t"$1"_"$2"_"$3}' $ISOCIRC > "${FRAC_DIR}/isocirc_full.bed"
    awk '{print $0"\t"$1"_"$2"_"$3}' $CIRCNICK > "${FRAC_DIR}/circnick_full.bed"
    
    # Extract only coordinates for ID purposes
    awk '{print $1"_"$2"_"$3}' $CIRI_LONG > "${FRAC_DIR}/ciri_unique_ids.txt"
    awk '{print $1"_"$2"_"$3}' $ISOCIRC > "${FRAC_DIR}/isocirc_unique_ids.txt"
    awk '{print $1"_"$2"_"$3}' $CIRCNICK > "${FRAC_DIR}/circnick_unique_ids.txt"
    
    # Compare tool pairs using -split option to handle exon structure
    echo "Comparing Ciri-long and isocirc with -split option..."
    bedtools intersect -a "${FRAC_DIR}/ciri_full.bed" -b "${FRAC_DIR}/isocirc_full.bed" -f $FRACTION -wa -wb -split > "${FRAC_DIR}/ciri_isocirc_intersect.bed"
    
    echo "Comparing Ciri-long and circnick with -split option..."
    bedtools intersect -a "${FRAC_DIR}/ciri_full.bed" -b "${FRAC_DIR}/circnick_full.bed" -f $FRACTION -wa -wb -split > "${FRAC_DIR}/ciri_circnick_intersect.bed"
    
    echo "Comparing isocirc and circnick with -split option..."
    bedtools intersect -a "${FRAC_DIR}/isocirc_full.bed" -b "${FRAC_DIR}/circnick_full.bed" -f $FRACTION -wa -wb -split > "${FRAC_DIR}/isocirc_circnick_intersect.bed"
    
    # Find common isoforms (triple intersection)
    echo "Finding common isoforms between all three tools with -split option..."
    
    # Get IDs from ciri that intersect with isocirc
    awk '{print $(NF-3)"_"$(NF-2)"_"$(NF-1)}' "${FRAC_DIR}/ciri_isocirc_intersect.bed" > "${FRAC_DIR}/ciri_isocirc_ids.txt"
    
    # Find which of these IDs also intersect with circnick
    while read -r id; do
        grep "$id" "${FRAC_DIR}/ciri_circnick_intersect.bed" >> "${FRAC_DIR}/common_isoforms.bed"
    done < "${FRAC_DIR}/ciri_isocirc_ids.txt"
    
    # Extract IDs for common isoforms
    if [ -f "${FRAC_DIR}/common_isoforms.bed" ]; then
        awk '{print $(NF-3)"_"$(NF-2)"_"$(NF-1)}' "${FRAC_DIR}/common_isoforms.bed" | sort -u > "${FRAC_DIR}/common_ids.txt"
    else
        touch "${FRAC_DIR}/common_ids.txt" # Create empty file if no common isoforms
    fi
    
    # Get unique ids from intersection files
    awk '{print $(NF-3)"_"$(NF-2)"_"$(NF-1)}' "${FRAC_DIR}/ciri_isocirc_intersect.bed" | sort -u > "${FRAC_DIR}/ciri_isocirc_ids.txt"
    awk '{print $(NF-3)"_"$(NF-2)"_"$(NF-1)}' "${FRAC_DIR}/ciri_circnick_intersect.bed" | sort -u > "${FRAC_DIR}/ciri_circnick_ids.txt"
    awk '{print $(NF-3)"_"$(NF-2)"_"$(NF-1)}' "${FRAC_DIR}/isocirc_circnick_intersect.bed" | sort -u > "${FRAC_DIR}/isocirc_circnick_ids.txt"
    
    # Get all unique coordinates across all tools
    cat "${FRAC_DIR}/ciri_unique_ids.txt" "${FRAC_DIR}/isocirc_unique_ids.txt" "${FRAC_DIR}/circnick_unique_ids.txt" | sort -u > "${FRAC_DIR}/all_coords_ids.txt"
    
    # Get counts
    CIRI_COUNT=$(wc -l < "${FRAC_DIR}/ciri_unique_ids.txt")
    ISOCIRC_COUNT=$(wc -l < "${FRAC_DIR}/isocirc_unique_ids.txt")
    CIRCNICK_COUNT=$(wc -l < "${FRAC_DIR}/circnick_unique_ids.txt")
    CIRI_ISOCIRC_COUNT=$(wc -l < "${FRAC_DIR}/ciri_isocirc_ids.txt")
    CIRI_CIRCNICK_COUNT=$(wc -l < "${FRAC_DIR}/ciri_circnick_ids.txt")
    ISOCIRC_CIRCNICK_COUNT=$(wc -l < "${FRAC_DIR}/isocirc_circnick_ids.txt")
    COMMON_COUNT=$(wc -l < "${FRAC_DIR}/common_ids.txt")
    
    # Write counts to file for R processing
    echo "Tool,Count" > "${FRAC_DIR}/isoform_counts.csv"
    echo "Ciri-long,$CIRI_COUNT" >> "${FRAC_DIR}/isoform_counts.csv"
    echo "isocirc,$ISOCIRC_COUNT" >> "${FRAC_DIR}/isoform_counts.csv"
    echo "circnick,$CIRCNICK_COUNT" >> "${FRAC_DIR}/isoform_counts.csv"
    
    echo "Intersection,Count" > "${FRAC_DIR}/intersection_counts.csv"
    echo "Ciri-long_isocirc,$CIRI_ISOCIRC_COUNT" >> "${FRAC_DIR}/intersection_counts.csv"
    echo "Ciri-long_circnick,$CIRI_CIRCNICK_COUNT" >> "${FRAC_DIR}/intersection_counts.csv"
    echo "isocirc_circnick,$ISOCIRC_CIRCNICK_COUNT" >> "${FRAC_DIR}/intersection_counts.csv"
    echo "Common,$COMMON_COUNT" >> "${FRAC_DIR}/intersection_counts.csv"
    
    # Create presence/absence matrix for UpSetR using efficient AWK approach
    echo "Creating presence/absence matrix for UpSetR using AWK..."
    
    awk 'BEGIN {print "coords,CIRI-long,isoCIRC,circNICK-lrs"} 
    FILENAME == ARGV[1] {a[$1]=1; ids[$1]=1}
    FILENAME == ARGV[2] {b[$1]=1; ids[$1]=1}
    FILENAME == ARGV[3] {c[$1]=1; ids[$1]=1}
    END {
        n = asorti(ids, sorted);
        for (i=1; i<=n; i++) {
            id = sorted[i];
            printf "%s,%d,%d,%d\n", id, (id in a), (id in b), (id in c)
        }
    }' "${FRAC_DIR}/ciri_unique_ids.txt" "${FRAC_DIR}/isocirc_unique_ids.txt" "${FRAC_DIR}/circnick_unique_ids.txt" > "${FRAC_DIR}/presence_absence.csv"
    
    # Clean up temporary files to save space
    echo "Cleaning up temporary files..."
    rm "${FRAC_DIR}/ciri_full.bed" "${FRAC_DIR}/isocirc_full.bed" "${FRAC_DIR}/circnick_full.bed"
    
    echo "Completed processing with overlap fraction: $FRACTION"
    
    # Print summary information to console
    echo "Summary for fraction $FRACTION with -split option:"
    echo "------------------------------"
    echo "Ciri-long: $CIRI_COUNT isoforms"
    echo "isocirc: $ISOCIRC_COUNT isoforms"
    echo "circnick: $CIRCNICK_COUNT isoforms"
    echo "Ciri-long ∩ isocirc: $CIRI_ISOCIRC_COUNT isoforms"
    echo "Ciri-long ∩ circnick: $CIRI_CIRCNICK_COUNT isoforms"
    echo "isocirc ∩ circnick: $ISOCIRC_CIRCNICK_COUNT isoforms"
    echo "Common (all three tools): $COMMON_COUNT isoforms"
    
    # Calculate percentages of overlap
    TOTAL_UNIQUE=$(($(wc -l < "${FRAC_DIR}/all_coords_ids.txt")))
    COMMON_PERCENT=$(echo "scale=2; ($COMMON_COUNT / $TOTAL_UNIQUE) * 100" | bc)
    echo "Percentage shared by all tools: $COMMON_PERCENT%"
    echo "------------------------------"
done

# Step 2: Create and execute R script for analysis and visualization
echo "Creating R script for analysis..."
cat > "${OUTPUT_DIR}/analyze_results.R" << 'EOL'
#!/usr/bin/env Rscript

# R script to analyze circRNA isoform prediction results and create upset plots
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
  if (!require("UpSetR")) {
    install.packages("UpSetR", repos="https://cloud.r-project.org")
    library(UpSetR)
  }
  if (!require("ggplot2")) {
    install.packages("ggplot2", repos="https://cloud.r-project.org")
    library(ggplot2)
  }
  if (!require("dplyr")) {
    install.packages("dplyr", repos="https://cloud.r-project.org")
    library(dplyr)
  }
  if (!require("tidyr")) {
    install.packages("tidyr", repos="https://cloud.r-project.org")
    library(tidyr)
  }
  if (!require("readr")) {
    install.packages("readr", repos="https://cloud.r-project.org")
    library(readr)
  }
})

# Set output directory for plots
plot_dir <- "/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/results/upset_split"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
cat(paste0("Output will be saved to: ", plot_dir, "\n"))

# Process results for each overlap fraction
fractions <- c("0.95")

# Process data for the combined folder
base_dir <- "circRNA_comparison_results_split/bedtools_results"

# Create upset plots and summary statistics for each fraction
for (fraction in fractions) {
  # Read presence/absence data
  data_path <- paste0(base_dir, "/combined/f", fraction, "/presence_absence.csv")
  
  cat(paste0("Processing overlap fraction: ", fraction, "\n"))
  
  if (file.exists(data_path)) {
    presence_data <- read_csv(data_path, col_types = cols())
    
    # Create upset plot
    cat("  Creating UpSet plot...\n")
    png_path <- paste0(plot_dir, "/upset_plot_f", fraction, "_split.png")
    png(file = png_path, width = 1600, height = 1200, res = 150)
    
    # Ensure data is properly formatted for UpSetR
    binary_data <- as.data.frame(presence_data)
    rownames(binary_data) <- binary_data$coords
    binary_data$coords <- NULL
    
    # Print column names for debugging
    cat("Column names in the data: ", paste(colnames(binary_data), collapse=", "), "\n")
    
    # Create a custom color vector matching what's shown in the example
    # isoCIRC: green (#00A86B)
    # CIRI-long: blue (#4682B4)
    # circNICK-lrs: orange (#FF8C00)
    
    print(upset(
      binary_data,
      nsets = 3,
      order.by = "freq",
      sets = c("CIRI-long", "isoCIRC", "circNICK-lrs"),
      sets.bar.color = c("#4682B4", "#00A86B", "#FF8C00"),  # Blue, Green, Orange
      main.bar.color = "black",
      matrix.color = "black",
      point.size = 4,
      line.size = 1.5,
      text.scale = c(2, 2, 1.5, 1.5, 2, 1.5),  # Increase all text elements
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size",
      shade.color = "black"  # Add shading with black outline
    ))
    dev.off()
    cat(paste0("  UpSet plot saved to: ", png_path, "\n"))
    
    # Create a Venn diagram like visualization
    # (We'll use the upset plot data to create a proportional representation)
    cat("  Creating intersection barplot...\n")
    ciri_only <- sum(binary_data$`CIRI-long` == 1 & binary_data$isoCIRC == 0 & binary_data$`circNICK-lrs` == 0)
    iso_only <- sum(binary_data$`CIRI-long` == 0 & binary_data$isoCIRC == 1 & binary_data$`circNICK-lrs` == 0)
    circ_only <- sum(binary_data$`CIRI-long` == 0 & binary_data$isoCIRC == 0 & binary_data$`circNICK-lrs` == 1)
    ciri_iso <- sum(binary_data$`CIRI-long` == 1 & binary_data$isoCIRC == 1 & binary_data$`circNICK-lrs` == 0)
    ciri_circ <- sum(binary_data$`CIRI-long` == 1 & binary_data$isoCIRC == 0 & binary_data$`circNICK-lrs` == 1)
    iso_circ <- sum(binary_data$`CIRI-long` == 0 & binary_data$isoCIRC == 1 & binary_data$`circNICK-lrs` == 1)
    all_three <- sum(binary_data$`CIRI-long` == 1 & binary_data$isoCIRC == 1 & binary_data$`circNICK-lrs` == 1)
    
    # Create a data frame for the bar plot
    bar_data <- data.frame(
      Category = c("CIRI-long only", "IsoCirc only", "CircNick only", 
                  "CIRI-long & IsoCirc", "CIRI-long & CircNick", "IsoCirc & CircNick", 
                  "All three tools"),
      Count = c(ciri_only, iso_only, circ_only, ciri_iso, ciri_circ, iso_circ, all_three),
      stringsAsFactors = FALSE
    )
    
    # Arrange in descending order
    bar_data <- bar_data[order(-bar_data$Count),]
    
    # Create a bar plot with black outlines
    png_path <- paste0(plot_dir, "/intersection_barplot_f", fraction, "_split.png")
    png(file = png_path, width = 1800, height = 1200, res = 150)
    p <- ggplot(bar_data, aes(x = reorder(Category, -Count), y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue", color = "black", size = 0.5) +  # Added black outline
      theme_minimal() +
      labs(
        title = paste0("CircRNA Isoform Intersections (overlap fraction = ", fraction, ", with -split)"),
        x = "",
        y = "Number of isoforms"
      ) +
      theme(
        text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        title = element_text(size = 20)
      )
    print(p)
    dev.off()
    cat(paste0("  Intersection barplot saved to: ", png_path, "\n"))
    
    # Calculate total and overlap percentages
    total_predictions <- nrow(presence_data)
    ciri_total <- sum(binary_data$`CIRI-long`)
    iso_total <- sum(binary_data$isoCIRC)
    circ_total <- sum(binary_data$`circNICK-lrs`)
    
    # Write summary to file
    cat("  Creating summary report...\n")
    summary_file <- paste0("/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/results/upset_split/summary_f", fraction, "_split.txt")
    writeLines(
      c(
        paste0("Summary of results with overlap fraction = ", fraction, " and -split option"),
        "----------------------------------------------------",
        paste0("Total unique circRNA isoforms across all tools: ", total_predictions),
        paste0("Isoforms shared by all three tools: ", all_three, " (", round((all_three / total_predictions) * 100, 2), "%)"),
        "\nIndividual tool statistics:",
        paste0("Ciri-long: ", ciri_total, " isoforms"),
        paste0("isocirc: ", iso_total, " isoforms"),
        paste0("circnick: ", circ_total, " isoforms"),
        "\nPairwise intersections:",
        paste0("Ciri-long ∩ isocirc: ", ciri_iso + all_three, " isoforms"),
        paste0("Ciri-long ∩ circnick: ", ciri_circ + all_three, " isoforms"),
        paste0("isocirc ∩ circnick: ", iso_circ + all_three, " isoforms"),
        "\nExclusive counts:",
        paste0("Unique to Ciri-long: ", ciri_only, " isoforms"),
        paste0("Unique to isocirc: ", iso_only, " isoforms"),
        paste0("Unique to circnick: ", circ_only, " isoforms")
      ),
      summary_file
    )
    cat(paste0("  Summary report saved to: ", summary_file, "\n"))
    
  } else {
    cat(paste0("Data file not found for overlap fraction: ", fraction, "\n"))
  }
}

# Create tool comparison plots
for (fraction in fractions) {
  # Read the counts for a pie chart
  data_path <- paste0(base_dir, "/combined/f", fraction, "/isoform_counts.csv")
  
  if (file.exists(data_path)) {
    tool_counts <- read_csv(data_path, col_types = cols())
    
    cat(paste0("Creating tool comparison plots for fraction ", fraction, "...\n"))
    
    # Create a pie chart of tool contributions
    png_path <- paste0(plot_dir, "/tool_comparison_pie_f", fraction, "_split.png")
    png(file = png_path, width = 1400, height = 1400, res = 150)
    
    # Rename the tools for consistency
    tool_counts$Tool <- gsub("Ciri-long", "CIRI-long", tool_counts$Tool)
    tool_counts$Tool <- gsub("circnick", "circNICK-lrs", tool_counts$Tool)
    tool_counts$Tool <- gsub("isocirc", "isoCIRC", tool_counts$Tool)
    
    pie_data <- tool_counts
    pie_data$percentage <- pie_data$Count / sum(pie_data$Count) * 100
    pie_data$label <- paste0(pie_data$Tool, "\n", pie_data$Count, " (", round(pie_data$percentage, 1), "%)")
    
    p <- ggplot(pie_data, aes(x = "", y = Count, fill = Tool)) +
      geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +  # Added black outline
      coord_polar("y", start = 0) +
      scale_fill_manual(values = c("CIRI-long" = "#4682B4", "isoCIRC" = "#00A86B", "circNICK-lrs" = "#FF8C00")) +
      theme_void() +
      labs(title = paste0("Tool Contributions (overlap fraction = ", fraction, ", with -split)")) +
      geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 6) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 22, hjust = 0.5)
      )
    
    print(p)
    dev.off()
    cat(paste0("  Pie chart saved to: ", png_path, "\n"))
    
    # Create a simple bar chart
    png_path <- paste0(plot_dir, "/tool_comparison_bar_f", fraction, "_split.png")
    png(file = png_path, width = 1400, height = 1000, res = 150)
    
    # Tool names should already be updated from pie chart code
    
    p <- ggplot(tool_counts, aes(x = reorder(Tool, -Count), y = Count, fill = Tool)) +
      geom_bar(stat = "identity", color = "black", size = 0.5) +  # Added black outline
      scale_fill_manual(values = c("CIRI-long" = "#4682B4", "isoCIRC" = "#00A86B", "circNICK-lrs" = "#FF8C00")) +
      theme_minimal() +
      labs(
        title = paste0("Number of CircRNA Isoforms per Tool (overlap fraction = ", fraction, ", with -split)"),
        x = "",
        y = "Number of isoforms"
      ) +
      theme(
        legend.position = "none", 
        text = element_text(size = 18),
        axis.text = element_text(size = 16),
        title = element_text(size = 20)
      ) +
      geom_text(aes(label = Count), vjust = -0.5, size = 6)
    
    print(p)
    dev.off()
    cat(paste0("  Bar chart saved to: ", png_path, "\n"))
  }
}

# Print completion message
cat("Analysis complete! Results are saved in /scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/results/upset_split/\n")
EOL

# Make the R script executable
chmod +x "${OUTPUT_DIR}/analyze_results.R"

# Execute the R script
echo "Running R analysis script..."
Rscript "${OUTPUT_DIR}/analyze_results.R"

# Generate summary table using presence_absence.csv for consistency with UpSet plot
FRACTION_DIR="${RUN_BEDTOOLS_DIR}/f0.95"
PRESENCE_FILE="${FRACTION_DIR}/presence_absence.csv"
SUMMARY_OUTPUT="${FRACTION_DIR}/final_summary_counts.tsv"

if [ -f "$PRESENCE_FILE" ]; then
    echo -e "Combination\tCount" > "$SUMMARY_OUTPUT"
    awk -F',' '
    NR > 1 {
        combo = ($2==1 ? "CIRI-long" : "") \
              (($2==1 && $3==1) ? "_" : "") ($3==1 ? "isoCIRC" : "") \
              ((($2==1 || $3==1) && $4==1) ? "_" : "") ($4==1 ? "circNICK-lrs" : "")
        if (combo == "") combo = "None"
        count[combo]++
    }
    END {
        for (k in count) print k "\t" count[k]
    }' "$PRESENCE_FILE" >> "$SUMMARY_OUTPUT"

    echo "Saved presence-based summary to: $SUMMARY_OUTPUT"
else
    echo "WARNING: presence_absence.csv not found. Skipping summary table generation."
fi


echo "CircRNA isoform comparison with -split option complete!"
echo "Results have been saved to: /scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/results/upset_split/"