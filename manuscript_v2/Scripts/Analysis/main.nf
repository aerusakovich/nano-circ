#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Help message
def helpMessage() {
    log.info"""
    ====================================
    circRNA Analysis Pipeline
    ====================================
    
    Usage:
    nextflow run main.nf -params-file params.yml
    
    Required arguments:
      --base_dir          Base directory containing run_* folders
      --genome_fasta      Reference genome FASTA file
      --gtf_file          Gene annotation GTF file
      --circrna_db        CircRNA database BED file
      --transcriptome     Transcriptome FASTA file
      --simulated_data_dir Directory with original circRNA definitions
    
    Optional arguments:
      --outdir            Output directory (default: results)
      --threads           Number of threads (default: 8)
      --run_cirilong      Run CIRI-long (default: true)
      --run_isocirc       Run IsoCirc (default: true)
      --run_circnick      Run CircNick (default: true)
      --track_performance Track performance (default: true)
      --cirilong_container Path to CIRI-long singularity container
      --isocirc_container Path to IsoCirc singularity container
      --overlap_fractions Overlap fractions for UpSet analysis (default: 0.95)
      --min_reads_per_circ Min reads to keep circRNA in ground truth (default: 5)
      --mapq_threshold    MAPQ threshold for mapping (default: 0)
      --initial_tile      Initial tiling factor (default: 1)
      --rescue_tile       Rescue tiling factor (default: 5)
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.base_dir) {
    error "Missing required parameter: --base_dir"
}
if (!params.genome_fasta) {
    error "Missing required parameter: --genome_fasta"
}
if (!params.gtf_file) {
    error "Missing required parameter: --gtf_file"
}
if (!params.circrna_db) {
    error "Missing required parameter: --circrna_db"
}
if (!params.simulated_data_dir) {
    error "Missing required parameter: --simulated_data_dir"
}

// Import modules
include { RUN_CIRILONG } from './modules/cirilong.nf'
include { RUN_ISOCIRC } from './modules/isocirc.nf'
include { RUN_CIRCNICK } from './modules/circnick.nf'
include { CONVERT_TO_BED12_CIRILONG } from './modules/bed12_conversion.nf'
include { CONVERT_TO_BED12_CIRCNICK } from './modules/bed12_conversion.nf'
include { CONVERT_TO_BED12_GROUNDTRUTH } from './modules/bed12_conversion.nf'
include { COLLECT_PERFORMANCE } from './modules/performance.nf'
include { VISUALIZE_PERFORMANCE } from './modules/performance.nf'
include { VISUALIZE_THROUGHPUT } from './modules/performance.nf'


workflow {
    
    // Reference files
    genome = file(params.genome_fasta)
    gtf = file(params.gtf_file)
    circrna_db = file(params.circrna_db)
    
    // Create ground truth BED12 files FIRST (filtered and validated)
    Channel
        .fromPath("${params.base_dir}/run_*", type: 'dir')
        .map { dir -> tuple(dir.name, dir) }
        .set { runs_for_ground_truth }

    ground_truth_bed12 = CONVERT_TO_BED12_GROUNDTRUTH(runs_for_ground_truth, genome)
    
    // Create input channel from run directories for tool execution
    run_dirs = Channel
        .fromPath("${params.base_dir}/run_*", type: 'dir')
        .map { dir ->
            def run_name = dir.name
            def fastq = file("${dir}/pooled/combined_reads.fastq")
            def ground_truth = file("${dir}/circRNAs.bed")
            def abundances = file("${dir}/abundances.tsv")
            
            if (!fastq.exists()) {
                log.warn "Missing FASTQ file for ${run_name}: ${fastq}"
                return null
            }
            
            tuple(run_name, fastq, ground_truth, abundances, dir)
        }
        .filter { it != null }
    
    // Run CIRI-long
    if (params.run_cirilong && params.cirilong_container) {
        container = file(params.cirilong_container)
        cirilong_results = RUN_CIRILONG(
            run_dirs,
            genome,
            gtf,
            circrna_db,
            container
        )
        
        cirilong_bed12 = CONVERT_TO_BED12_CIRILONG(cirilong_results)
    } else {
        cirilong_results = Channel.empty()
        cirilong_bed12 = Channel.empty()
    }
    
    // Run IsoCirc
    if (params.run_isocirc && params.isocirc_container) {
        isocirc_container = file(params.isocirc_container)
        isocirc_results = RUN_ISOCIRC(
            run_dirs,
            genome,
            gtf,
            circrna_db,
            isocirc_container
        )
    } else {
        isocirc_results = Channel.empty()
    }
    
    // Run CircNick
    if (params.run_circnick) {
        circnick_results = RUN_CIRCNICK(
            run_dirs,
            genome,
            gtf
        )
        
        circnick_bed12 = CONVERT_TO_BED12_CIRCNICK(circnick_results)
    } else {
        circnick_results = Channel.empty()
        circnick_bed12 = Channel.empty()
    }
    
    // Collect performance metrics
    if (params.track_performance) {
        // Wait for all tools to complete
        all_results = Channel.empty()
            .mix(
                cirilong_results.ifEmpty(Channel.empty()),
                isocirc_results.ifEmpty(Channel.empty()),
                circnick_results.ifEmpty(Channel.empty())
            )
            .collect()
        
        // Process trace file after all tools complete
        all_results
            .map { file("${params.outdir}/pipeline_info/trace.txt") }
            .set { trace_input }
        
        perf_data = COLLECT_PERFORMANCE(trace_input)
        
        // Create both visualizations
        VISUALIZE_PERFORMANCE(perf_data.time_csv, perf_data.memory_csv)
        VISUALIZE_THROUGHPUT(perf_data.time_csv, perf_data.memory_csv)
    }
}

workflow.onComplete {
    log.info """
    ====================================
    Pipeline completed!
    ====================================
    Status:    ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration:  ${workflow.duration}
    Output:    ${params.outdir}
    
    Ground truth filtering applied:
      - Min reads per circRNA: ${params.min_reads_per_circ ?: 5}
      - MAPQ threshold: ${params.mapq_threshold ?: 0}
      - Initial tiling: ${params.initial_tile ?: 1}x
      - Rescue tiling: ${params.rescue_tile ?: 5}x
    """.stripIndent()
}