process RUN_ISOCIRC {
    tag "${run_name}"
    
    publishDir "${params.outdir}/tool_outputs/isocirc/${run_name}", mode: 'copy'

    input:
    tuple val(run_name), path(fastq), path(ground_truth), path(abundances), path(run_dir)
    path genome
    path gtf
    path circrna_db
    path container
    
    output:
    tuple val(run_name), path("isocirc_output"), path("performance_metrics.txt"), emit: results
    
    script:
    """
    #!/bin/bash
    
    set -e
    
    START_TIME=\$(date +%s)
    START_DATETIME=\$(date "+%Y-%m-%d %H:%M:%S")
    
    # Get absolute paths
    ABS_FASTQ=\$(realpath ${fastq})
    ABS_GENOME=\$(realpath ${genome})
    ABS_GTF=\$(realpath ${gtf})
    ABS_DB=\$(realpath ${circrna_db})
    ABS_OUTPUT=\$(realpath .)/isocirc_output
    
    mkdir -p isocirc_output
    
    # Run IsoCirc using Singularity container with bind mounts and proper environment
    singularity exec \\
        --bind \$(pwd):\$(pwd) \\
        --bind /scratch:/scratch \\
        --bind /projects:/projects \\
        --env LC_ALL=C \\
        --env LANG=C \\
        --cleanenv \\
        ${container} isocirc \\
        "\${ABS_FASTQ}" \\
        "\${ABS_GENOME}" \\
        "\${ABS_GTF}" \\
        "\${ABS_DB}" \\
        "\${ABS_OUTPUT}" \\
        -t ${params.threads}
    
    END_TIME=\$(date +%s)
    END_DATETIME=\$(date "+%Y-%m-%d %H:%M:%S")
    ELAPSED=\$((END_TIME - START_TIME))
    
    cat > performance_metrics.txt <<PERF
Run: ${run_name}
Tool: isocirc
Start: \${START_DATETIME}
End: \${END_DATETIME}
Elapsed_Seconds: \${ELAPSED}
Threads: ${params.threads}
PERF
    
    echo "IsoCirc completed for ${run_name} in \${ELAPSED} seconds"
    """
}