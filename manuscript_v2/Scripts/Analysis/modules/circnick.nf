process RUN_CIRCNICK {
    tag "${run_name}"
    
    publishDir "${params.outdir}/tool_outputs/circnick/${run_name}", mode: 'copy'
    
    input:
    tuple val(run_name), path(fastq), path(ground_truth), path(abundances), path(run_dir)
    path genome
    path gtf
    
    output:
    tuple val(run_name), path("circnick_output"), path("performance_metrics.txt"), emit: results
    
    script:
    """
    #!/bin/bash
    
    # Activate conda environment
    source /local/env/envpython-3.9.5.sh || true
    conda activate long_read_circRNA || true
    
    set -e

    START_TIME=\$(date +%s)
    START_DATETIME=\$(date "+%Y-%m-%d %H:%M:%S")
    
    mkdir -p circnick_output
    
    # Compress FASTQ if needed
    compressed_fastq="${fastq.baseName}.fq.gz"
    if [ ! -f "\${compressed_fastq}" ]; then
        gzip -c ${fastq} > \${compressed_fastq}
    fi
    
    # Run CircNick
    ${params.circnick_path}/long_read_circRNA run \\
        \${compressed_fastq} \\
        --species ${params.species} \\
        --reference-path ${params.circnick_ref_path} \\
        --script-path ${params.circnick_scripts} \\
        --output-path circnick_output
    
    END_TIME=\$(date +%s)
    END_DATETIME=\$(date "+%Y-%m-%d %H:%M:%S")
    ELAPSED=\$((END_TIME - START_TIME))
    
    cat > performance_metrics.txt <<PERF
Run: ${run_name}
Tool: circnick
Start: \${START_DATETIME}
End: \${END_DATETIME}
Elapsed_Seconds: \${ELAPSED}
Threads: ${params.threads}
PERF
    
    echo "CircNick completed for ${run_name} in \${ELAPSED} seconds"
    """
}