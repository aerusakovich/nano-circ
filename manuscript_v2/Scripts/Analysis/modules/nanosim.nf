nextflow.enable.dsl=2

/*
 * NanoSim module (with pooling)
 *
 * Input:
 *   sim_runs       : (run_name, run_dir) where run_dir = params.simulated_data_dir/run_*
 *   genome_fasta   : reference genome FASTA
 *   gtf_file       : gene annotation GTF
 *   transcriptome  : transcriptome FASTA
 *   control_fastq  : real control FASTQ used for NanoSim characterization
 *
 * Params (from params.yaml):
 *   simulated_data_dir
 *   base_dir
 *   genome_fasta
 *   gtf_file
 *   transcriptome
 *   control_fastq
 *   nanosim_container
 *   circrna_reads
 *   linear_reads
 *
 * Output:
 *   - run_* directories under params.base_dir, each containing:
 *       circ/   linear/   control/   expression/   ...
 *       pooled/combined_reads.fastq
 *       pooled/combined_reads.fq.gz
 */

process NANOSIM_SIMULATE {

    tag "${run_name}"
    cpus 16
    memory '150 GB'

    // Publish each run_* dir (with NanoSim + pooled reads) into params.base_dir
    publishDir "${params.base_dir}", mode: 'copy', pattern: "run_*"

    input:
    // (run_name, run_dir_in_simulated_data_dir)
    tuple val(run_name), path(run_dir)

    // reference files
    path genome_fasta
    path gtf_file
    path transcriptome
    path control_fastq

    output:
    // Each run_* directory produced by this process
    path "run_*", emit: runs

    script:
    """
    #!/bin/bash
    source /local/env/envpython-3.7.1.sh || true
    conda activate nanosim || true
    set -euo pipefail

    echo "==================================="
    echo "NanoSim simulation for ${run_name}"
    echo "==================================="

    # Input directory with circRNAs.fa and abundances.tsv
    INPUT_DIR="${run_dir}"

    # Output directory name for this run (will be published as run_*)
    # We'll name it run_<run_name> to align with the rest of the pipeline.
    OUTPUT_DIR="${run_name}"
    mkdir -p "\${OUTPUT_DIR}"

    # Reference files
    GTF_FILE="${gtf_file}"
    GENOME_FILE="${genome_fasta}"
    TRANSCRIPTOME_FILE="${transcriptome}"
    CONTROL_FQ="${control_fastq}"

    # NanoSim container (from params)
    NANOSIM_CONTAINER="${params.nanosim_container}"

    echo "Using:"
    echo "  INPUT_DIR         = \${INPUT_DIR}"
    echo "  OUTPUT_DIR        = \${OUTPUT_DIR}"
    echo "  GENOME_FILE       = \${GENOME_FILE}"
    echo "  GTF_FILE          = \${GTF_FILE}"
    echo "  TRANSCRIPTOME     = \${TRANSCRIPTOME_FILE}"
    echo "  CONTROL_FQ        = \${CONTROL_FQ}"
    echo "  circrna_reads     = ${params.circrna_reads}"
    echo "  linear_reads      = ${params.linear_reads}"
    echo "  NanoSim container = \${NANOSIM_CONTAINER}"
    echo "  Threads (task.cpus) = ${task.cpus}"

    # Sanity checks for input
    if [ ! -f "\${INPUT_DIR}/circRNAs.fa" ]; then
        echo "ERROR: circRNAs.fa not found in \${INPUT_DIR}"
        exit 1
    fi

    if [ ! -f "\${INPUT_DIR}/abundances.tsv" ]; then
        echo "ERROR: abundances.tsv not found in \${INPUT_DIR}"
        exit 1
    fi



    # ----------------------------------------------------
    # STEP 1: NanoSim characterization (transcriptome)
    # ----------------------------------------------------
    echo "Step 1. NanoSim characterization for ${run_name}"
    singularity exec -B /home/genouest/cnrs_umr6290/aerusakovich/ -B /scratch -B /projects "\${NANOSIM_CONTAINER}" read_analysis.py transcriptome \\
        -i "\${CONTROL_FQ}" \\
        -rg "\${GENOME_FILE}" \\
        -rt "\${TRANSCRIPTOME_FILE}" \\
        -a minimap2 \\
        -o "\${OUTPUT_DIR}/control" \\
        -t ${task.cpus} \\
        -annot "\${GTF_FILE}" \\
        --no_intron_retention \\
        --fastq

    # ----------------------------------------------------
    # STEP 2: Simulate circRNA reads
    # ----------------------------------------------------
    echo "Step 2. Simulate circRNA reads for ${run_name}"
    singularity exec -B /home/genouest/cnrs_umr6290/aerusakovich/ -B /scratch -B /projects "\${NANOSIM_CONTAINER}" simulator.py transcriptome \\
        -rg "\${GENOME_FILE}" \\
        -rt "\${INPUT_DIR}/circRNAs.fa" \\
        -c "\${OUTPUT_DIR}/control" \\
        -o "\${OUTPUT_DIR}/circ" \\
        -n ${params.circrna_reads} \\
        -b guppy \\
        -e "\${INPUT_DIR}/abundances.tsv" \\
        --fastq \\
        --no_model_ir

    # ----------------------------------------------------
    # STEP 3: Quantify linear transcripts
    # ----------------------------------------------------
    echo "Step 3. NanoSim quantification of linear transcripts for ${run_name}"
    singularity exec -B /home/genouest/cnrs_umr6290/aerusakovich/ -B /scratch -B /projects "\${NANOSIM_CONTAINER}" read_analysis.py quantify \\
        -e trans \\
        -i "\${CONTROL_FQ}" \\
        -rt "\${TRANSCRIPTOME_FILE}" \\
        -t ${task.cpus} \\
        -o "\${OUTPUT_DIR}/expression"

    # ----------------------------------------------------
    # STEP 4: Simulate linear mRNA reads
    # ----------------------------------------------------
    echo "Step 4. Simulate linear mRNA reads for ${run_name}"
    singularity exec -B /home/genouest/cnrs_umr6290/aerusakovich/ -B /scratch -B /projects "\${NANOSIM_CONTAINER}" simulator.py transcriptome \\
        -rg "\${GENOME_FILE}" \\
        -rt "\${TRANSCRIPTOME_FILE}" \\
        -c "\${OUTPUT_DIR}/control" \\
        -o "\${OUTPUT_DIR}/linear" \\
        -n ${params.linear_reads} \\
        -b guppy \\
        -e "\${OUTPUT_DIR}/expression_transcriptome_quantification.tsv" \\
        --fastq \\
        --no_model_ir

    echo "NanoSim simulation finished for ${run_name}"
    echo "Now pooling circ/linear NanoSim outputs..."

    # ----------------------------------------------------
    # STEP 5: Pool NanoSim-generated reads
    # ----------------------------------------------------

    RUN_DIR="\${OUTPUT_DIR}"

    POOL_DIR="\${RUN_DIR}/pooled"
    mkdir -p "\${POOL_DIR}"

    OUTPUT_FILE="\${POOL_DIR}/combined_reads.fastq"

    circ_aligned="\${RUN_DIR}/circ_aligned_reads.fastq"
    circ_unaligned="\${RUN_DIR}/circ_unaligned_reads.fastq"
    linear_aligned="\${RUN_DIR}/linear_aligned_reads.fastq"
    linear_unaligned="\${RUN_DIR}/linear_unaligned_reads.fastq"

    missing=false
    for f in "\${circ_aligned}" "\${circ_unaligned}" "\${linear_aligned}" "\${linear_unaligned}"; do
        if [ ! -f "\${f}" ]; then
            echo "ERROR: Expected NanoSim output file not found for ${run_name}: \${f}"
            missing=true
        fi
    done

    if [ "\${missing}" = true ]; then
        echo "Cannot pool reads for ${run_name} because some NanoSim FASTQ files are missing."
        exit 1
    fi

    echo "Pooling NanoSim files for ${run_name} into \${OUTPUT_FILE}"

    # circ_aligned: add _circ suffix to header
    awk '{
        if (NR % 4 == 1) {
            print substr(\$0, 1, 1) substr(\$0, 2) "_circ";
        } else {
            print;
        }
    }' "\${circ_aligned}" > "\${OUTPUT_FILE}"

    # circ_unaligned: add _circ suffix to header
    awk '{
        if (NR % 4 == 1) {
            print substr(\$0, 1, 1) substr(\$0, 2) "_circ";
        } else {
            print;
        }
    }' "\${circ_unaligned}" >> "\${OUTPUT_FILE}"

    # linear_aligned: add _linear suffix to header
    awk '{
        if (NR % 4 == 1) {
            print substr(\$0, 1, 1) substr(\$0, 2) "_linear";
        } else {
            print;
        }
    }' "\${linear_aligned}" >> "\${OUTPUT_FILE}"

    # linear_unaligned: add _linear suffix to header
    awk '{
        if (NR % 4 == 1) {
            print substr(\$0, 1, 1) substr(\$0, 2) "_linear";
        } else {
            print;
        }
    }' "\${linear_unaligned}" >> "\${OUTPUT_FILE}"

    echo "Created pooled FASTQ: \${OUTPUT_FILE}"

    echo "Compressing pooled FASTQ to \${OUTPUT_FILE%.fastq}.fq.gz"
    gzip -c "\${OUTPUT_FILE}" > "\${OUTPUT_FILE%.fastq}.fq.gz"

    echo "Completed NanoSim + pooling for ${run_name}"
    echo "Final run directory: \${OUTPUT_DIR}"
    """
}


workflow RUN_NANOSIM {

    take:
    sim_runs
    genome
    gtf
    transcriptome
    control_fastq

    main:
    NANOSIM_SIMULATE(
        sim_runs,
        genome,
        gtf,
        transcriptome,
        control_fastq
    )

    emit:
    runs = NANOSIM_SIMULATE.out.runs
}
