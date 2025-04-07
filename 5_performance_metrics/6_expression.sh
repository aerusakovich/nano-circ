#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr
#SBATCH --job-name=circrna_analysis
#SBATCH --output=circrna_analysis_%j.out
#SBATCH --cpus-per-task=4 
#SBATCH --mem=16G

# Unified CircRNA Analysis Script
# This script combines bedtools analysis with Python visualization
# to analyze circRNA predictions from different tools and generate comprehensive reports

# Set paths to input files
GROUND_TRUTH="/scratch/aerusakovich/sim_ciri_long_jobim/combined/pooled/ground_truth/combined_all.bed"
CIRCRNA_DB="/scratch/aerusakovich/sim_ciri_long_jobim/combined/circRNAs.bed"
ISOCIRC="/scratch/aerusakovich/sim_ciri_long_jobim/combined/isocirc_output/isocirc.bed"
CIRI_LONG="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/ciri-long/ciri_long_output/CIRI-long.cand_circ.bed12"
CIRCNICK="/scratch/aerusakovich/sim_ciri_long_jobim/all/bed12/circnick/combined/comprehensive_circrna.bed12"

# Create output directory structure
OUTPUT_DIR="circRNA_analysis_results"
BEDTOOLS_DIR="${OUTPUT_DIR}/bedtools_results"
STATS_DIR="${OUTPUT_DIR}/stats"
VIS_DIR="${OUTPUT_DIR}/visualizations"
mkdir -p $BEDTOOLS_DIR $STATS_DIR $VIS_DIR

echo "Starting comprehensive circRNA analysis..."

# Check if all required files exist
for FILE in "$GROUND_TRUTH" "$CIRCRNA_DB" "$ISOCIRC" "$CIRI_LONG" "$CIRCNICK"; do
    if [ ! -f "$FILE" ]; then
        echo "ERROR: Missing file at $FILE"
        exit 1
    fi
done

# Check required tools and packages
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools is required but not installed. Aborting."; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "Error: Python 3 is required but not installed. Aborting."; exit 1; }

# Try to check for required Python packages
python3 -c "import pandas, matplotlib, seaborn, numpy" >/dev/null 2>&1 || { 
    echo "Warning: One or more required Python packages may be missing.";
    echo "Required packages: pandas, matplotlib, seaborn, numpy, matplotlib_venn, upsetplot";
    echo "Continuing, but script may fail later...";
}

# Extract circRNA types from the circRNAs.bed file
echo "Extracting circRNA types from database..."
cut -f4 "$CIRCRNA_DB" | grep -o 'EIciRNA\|eciRNA\|ciRNA\|intergenic' | sort | uniq -c > "${STATS_DIR}/circrna_types_in_db.txt"

# Function to perform bedtools intersect and count circRNA types
analyze_tool_predictions() {
    local TOOL_NAME=$1
    local TOOL_FILE=$2
    local OVERLAP=$3
    
    echo "Analyzing ${TOOL_NAME} predictions with overlap fraction ${OVERLAP}..."
    
    # Intersect tool predictions with ground truth
    bedtools intersect -a "$TOOL_FILE" -b "$GROUND_TRUTH" -f "$OVERLAP" -wa > "${BEDTOOLS_DIR}/${TOOL_NAME}_vs_ground_truth_${OVERLAP}.bed"
    
    # Intersect tool predictions with circRNA database to identify types
    bedtools intersect -a "$TOOL_FILE" -b "$CIRCRNA_DB" -f "$OVERLAP" -wb > "${BEDTOOLS_DIR}/${TOOL_NAME}_vs_circrna_db_${OVERLAP}.bed"
    
    # Count circRNA types detected by the tool using awk
    awk -F'|' '{
        for(i=1;i<=NF;i++) {
            if($i=="EIciRNA" || $i=="eciRNA" || $i=="ciRNA" || $i=="intergenic")
                count[$i]++
        }
    } END {
        for(type in count)
            print count[type], type
    }' "${BEDTOOLS_DIR}/${TOOL_NAME}_vs_circrna_db_${OVERLAP}.bed" | sort -nr > "${STATS_DIR}/${TOOL_NAME}_circrna_types_${OVERLAP}.txt"
    
    # Count total predictions and matches
    local TOTAL_PREDICTIONS=$(wc -l < "$TOOL_FILE")
    local MATCHED_TO_GROUND=$(wc -l < "${BEDTOOLS_DIR}/${TOOL_NAME}_vs_ground_truth_${OVERLAP}.bed")
    
    # Use awk to count non-empty lines
    local MATCHED_TO_DB=$(awk 'END{print NR}' "${BEDTOOLS_DIR}/${TOOL_NAME}_vs_circrna_db_${OVERLAP}.bed")
    
    echo "${TOOL_NAME},${TOTAL_PREDICTIONS},${MATCHED_TO_GROUND},${MATCHED_TO_DB}" >> "${STATS_DIR}/summary_stats_${OVERLAP}.csv"
}

# Array of overlap fractions to test
FRACTIONS=("0.95")

# Initialize summary files with headers
for OVERLAP in "${FRACTIONS[@]}"; do
    echo "Tool,Total_Predictions,Matched_to_Ground_Truth,Matched_to_DB" > "${STATS_DIR}/summary_stats_${OVERLAP}.csv"
done

# Process each tool with different overlap fractions
for OVERLAP in "${FRACTIONS[@]}"; do
    analyze_tool_predictions "ciri_long" "$CIRI_LONG" "$OVERLAP"
    analyze_tool_predictions "isocirc" "$ISOCIRC" "$OVERLAP"
    analyze_tool_predictions "circnick" "$CIRCNICK" "$OVERLAP"
done

