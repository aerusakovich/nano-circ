process RUN_CIRILONG {
    tag "${run_name}"
    
    publishDir "${params.outdir}/tool_outputs/cirilong/${run_name}", mode: 'copy'
    
    input:
    tuple val(run_name), path(fastq), path(ground_truth), path(abundances), path(run_dir)
    path genome
    path gtf
    path circrna_db
    path container
    
    output:
    tuple val(run_name), path("cirilong_output"), path("performance_metrics.txt"), emit: results
    
    script:
    """
    #!/bin/bash
    set -e
    
    START_TIME=\$(date +%s)
    START_DATETIME=\$(date "+%Y-%m-%d %H:%M:%S")
    
    # Get absolute paths for all inputs
    ABS_FASTQ=\$(realpath ${fastq})
    ABS_GENOME=\$(realpath ${genome})
    ABS_GTF=\$(realpath ${gtf})
    ABS_DB=\$(realpath ${circrna_db})
    ABS_WORKDIR=\$(pwd)
    
    mkdir -p cirilong_output
    
    singularity exec --writable-tmpfs \\
        -B \${ABS_WORKDIR}  \\
        -B /scratch \\
        -B /home \\
        -B /projects \\
        ${container} \\
        CIRI-long call \\
            -i "\${ABS_FASTQ}" \\
            -r "\${ABS_GENOME}" \\
            -o "\${ABS_WORKDIR}/cirilong_output" \\
            -t ${params.threads} \\
            -a "\${ABS_GTF}" \\
            -c "\${ABS_DB}"
    
    sample_name=\$(basename ${fastq} .fastq)
    echo "\${sample_name} \${ABS_WORKDIR}/cirilong_output/CIRI-long.cand_circ.fa" > cirilong_output/list.txt
    
    singularity exec --writable-tmpfs \\
        -B \${ABS_WORKDIR} \\
        -B /scratch \\
        -B /projects \\
        -B /home \\
        ${container} \\
        CIRI-long collapse \\
            -i "\${ABS_WORKDIR}/cirilong_output/list.txt" \\
            -o "\${ABS_WORKDIR}/cirilong_output" \\
            -r "\${ABS_GENOME}" \\
            -a "\${ABS_GTF}"
    
    END_TIME=\$(date +%s)
    END_DATETIME=\$(date "+%Y-%m-%d %H:%M:%S")
    ELAPSED=\$((END_TIME - START_TIME))
    
    cat > performance_metrics.txt <<PERF
Run: ${run_name}
Tool: cirilong
Start: \${START_DATETIME}
End: \${END_DATETIME}
Elapsed_Seconds: \${ELAPSED}
Threads: ${params.threads}
PERF
    
    echo "CIRI-long completed for ${run_name} in \${ELAPSED} seconds"
    """
}