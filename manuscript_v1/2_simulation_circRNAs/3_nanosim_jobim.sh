#!/bin/bash
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=anastasia.rusakovich@univ-rennes1.fr # Where to send mail
#SBATCH --job-name=Nanosim # Job name
##SBATCH --gres=gpu:1 -p gpu # If needed to specify a node (i.e GPU)
#SBATCH --output=/scratch/aerusakovich/Nanosim_isocirc.out # Log file
#SBATCH --cpus-per-task=16 --mem=10G # Number of CPUs and memory

# Define file paths
GTF_FILE='/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/gencode.vM10.sorted.annotation.gtf'
GENOME_FILE='/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/GRCm38.p4.genome_corrected.fa'
CONTROL_FQ="/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/rep1/CRR194180.fq"
TRANSCRIPTOME='/home/genouest/cnrs_umr6290/aerusakovich/benchmark/data/ciri-long/files/m38/ref_transcriptome/gencode.vM10.transcripts.fa'
NANOSIM_VERSION='/cvmfs/singularity.galaxyproject.org/all/nanosim:3.1.0--hdfd78af_0'
THREADS=16

# Initialize Python environment
echo "Initializing Python environment"
source /local/env/envpython-3.7.1.sh
conda activate nanosim

# Get all run folders
BASE_DIR="/scratch/aerusakovich/sim_ciri_long_jobim"
RUN_DIRS=$(find ${BASE_DIR} -maxdepth 1 -type d -name "run_*")

# Process each Run folder
for OUTPUT in ${RUN_DIRS}; do
    RUN_NAME=$(basename ${OUTPUT})
    echo "Processing ${RUN_NAME}"
    
    echo "Step1. Characterization stage for ${RUN_NAME}"
    singularity exec -B /home/genouest/cnrs_umr6290/aerusakovich/ -B /scratch $NANOSIM_VERSION read_analysis.py transcriptome \
    -i ${CONTROL_FQ} \
    -rg ${GENOME_FILE} \
    -rt ${TRANSCRIPTOME} \
    -a minimap2 \
    -o ${OUTPUT}/control \
    -t ${THREADS} \
    -annot ${GTF_FILE} \
    --no_intron_retention
    
    echo "Step2. Simulation circRNA reads for ${RUN_NAME}"
    singularity exec -B /home/genouest/cnrs_umr6290/aerusakovich/ -B /scratch $NANOSIM_VERSION simulator.py transcriptome \
    -rg ${GENOME_FILE} \
    -rt ${OUTPUT}/circRNAs.fa \
    -c ${OUTPUT}/control \
    -o ${OUTPUT}/circ \
    -n 200000 \
    -b guppy \
    -r cDNA_1D \
    -e ${OUTPUT}/abundances.tsv \
    --fastq \
    --no_model_ir
    
    echo "Step3. Quantification of linear genes for ${RUN_NAME}"
    singularity exec -B /home/genouest/cnrs_umr6290/aerusakovich/ -B /scratch $NANOSIM_VERSION read_analysis.py quantify \
    -e trans \
    -i ${CONTROL_FQ} \
    -rt ${TRANSCRIPTOME} \
    -t ${THREADS} \
    -o ${OUTPUT}/expression
    
    echo "Step4. Simulate linear mRNA reads for ${RUN_NAME}"
    singularity exec -B /home/genouest/cnrs_umr6290/aerusakovich/ -B /scratch $NANOSIM_VERSION simulator.py transcriptome \
    -rg ${GENOME_FILE} \
    -rt ${TRANSCRIPTOME} \
    -c ${OUTPUT}/control \
    -o ${OUTPUT}/linear \
    -n 200000 \
    -b guppy \
    -r cDNA_1D \
    -e ${OUTPUT}/expression_transcriptome_quantification.tsv \
    --fastq \
    --no_model_ir
    
    echo "Completed processing ${RUN_NAME}"
done

echo "All runs completed"