#!/bin/bash
#SBATCH --job-name=circRNA_seq_analysis # Job name
#SBATCH --output=/home/genouest/cnrs_umr6290/aerusakovich/logs/circRNA_seq_analysis_%j.out # Output log file
#SBATCH --error=/home/genouest/cnrs_umr6290/aerusakovich/logs/circRNA_seq_analysis_%j.err # Error log file
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr # Where to send mail
#SBATCH --ntasks=1 # Run on a single CPU
#SBATCH --cpus-per-task=8 # Number of CPU cores per task
#SBATCH --mem=60G # Memory per node
#SBATCH --time=24:00:00 # Time limit hrs:min:sec

# Logging function
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}

# Load required environment
log "Loading environment"
source /local/env/envpython-3.7.1.sh
conda activate nanopore_env

# Define base directories and files
BASE_OUTPUT_DIR="/scratch/aerusakovich/circRNA_analysis/database/$(date +%Y%m%d)"
SCRIPT_PATH="/scratch/aerusakovich/~/new_git/database_study/database_study.py"
LIFTOVER_SCRIPT="/scratch/aerusakovich/~/new_git/database_study/liftover_script.py"  # Path to the liftOver script

# Reference genome and annotation files
GENOME_REF="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/GRCm38.p4.genome_corrected.fa"
GTF_FILE="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/gencode.vM10.sorted.annotation.gtf"

# Define database sources
declare -A DATABASES=(
    ["circatlas"]="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/db/circatlas/mouse_sequence_v3.fasta"
    ["circbase"]="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/db/circbase/mouse_mm9_circRNAs_putative_spliced_sequence.fa"
)

# Define which databases need conversion from mm9 to mm10
declare -A NEEDS_LIFTOVER=(
    ["circatlas"]=0  # 0 = false
    ["circbase"]=1   # 1 = true, this is mm9 and needs conversion
)

# Create base output directory
mkdir -p "$BASE_OUTPUT_DIR"
log "Created base output directory: $BASE_OUTPUT_DIR"

# Process each database
for DB_NAME in "${!DATABASES[@]}"; do
    # Create database-specific output directory
    DB_OUTPUT_DIR="${BASE_OUTPUT_DIR}/${DB_NAME}"
    mkdir -p "$DB_OUTPUT_DIR"
    log "Starting analysis for ${DB_NAME}"
    
    FASTA_TO_USE="${DATABASES[$DB_NAME]}"
    
    # Check if this database needs liftOver from mm9 to mm10
    if [ "${NEEDS_LIFTOVER[$DB_NAME]}" -eq 1 ]; then
        log "This database (${DB_NAME}) uses mm9 coordinates. Converting to mm10 first..."
        
        # Create liftOver directory
        LIFTOVER_DIR="${DB_OUTPUT_DIR}/liftover"
        mkdir -p "$LIFTOVER_DIR"
        
        # Run the liftOver script
        python "$LIFTOVER_SCRIPT" --fasta "$FASTA_TO_USE" --output "$LIFTOVER_DIR"
        
        # Check if conversion was successful
        if [ $? -eq 0 ] && [ -f "${LIFTOVER_DIR}/mm10_circrnas.fasta" ]; then
            log "Successfully converted to mm10 coordinates"
            FASTA_TO_USE="${LIFTOVER_DIR}/mm10_circrnas.fasta"
        else
            log "ERROR: Failed to convert coordinates. Using original file."
        fi
    fi
    
    # Run the analysis script with the appropriate FASTA file
    log "Running analysis with FASTA: $FASTA_TO_USE"
    python "$SCRIPT_PATH" \
        --fasta "$FASTA_TO_USE" \
        --genome "$GENOME_REF" \
        --gtf "$GTF_FILE" \
        --output "$DB_OUTPUT_DIR" \
        --workers 8
    
    # Check script exit status
    if [ $? -eq 0 ]; then
        log "Completed analysis for ${DB_NAME} successfully"
    else
        log "ERROR: Analysis failed for ${DB_NAME}"
    fi
done

# Final log
log "CircRNA sequence analysis complete"