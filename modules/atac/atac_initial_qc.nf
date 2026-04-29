process ATAC_INITIAL_QC {
    label 'process_high'
    publishDir "${params.outdir}/atac/initial_qc", mode: 'copy'

    input:
    path fragment_files
    path metadata

    output:
    path "atac_initial_qc.h5ads", emit: anndataset
    path "sample_statistics.csv", emit: sample_stats
    path "qc_plots/*", emit: qc_plots, optional: true
    // JSON summary actually written by atac_initial_qc.py
    path "atac_initial_qc_summary.json", emit: summary, optional: true
    // Per-sample h5ad files (written in working dir, e.g. L1_Donor_0_july.h5ad etc.)
    path "*.h5ad", emit: individual_samples

    script:
    def frag_list = fragment_files.collect { "'${it}'" }.join(' ')
    
    def min_counts_arg  = params.atac.min_counts  ? "--min_counts ${params.atac.min_counts}"     : ""
    def min_tsse_arg    = params.atac.min_tsse    ? "--min_tsse ${params.atac.min_tsse}"         : (params.atac.initial_min_tsse ? "--min_tsse ${params.atac.initial_min_tsse}" : "")
    def max_counts_arg  = params.atac.max_counts  ? "--max_counts ${params.atac.max_counts}"     : ""
    def batch_key_arg   = params.atac.batch_key   ? "--batch_key ${params.atac.batch_key}"       : ""
    def tempdir_arg     = params.atac.tempdir     ? "--tempdir ${params.atac.tempdir}"           : ""
    
    """
    mkdir -p qc_plots

    bash ${projectDir}/bin/atac_pipeline_wrapper.sh \
      --fragment_files ${frag_list} \
      --metadata ${metadata} \
      --species ${params.species} \
      --genome_build ${params.scprinter.genome} \
      --output_dir . \
      --min_fragments ${params.atac.min_fragments} \
      ${min_counts_arg} \
      ${min_tsse_arg} \
      ${max_counts_arg} \
      --n_features ${params.atac.n_features} \
      --batch_correction ${params.atac.batch_correction} \
      ${batch_key_arg} \
      --clustering_resolutions ${params.atac.clustering_resolutions.join(' ')} \
      --peak_fdr ${params.atac.peak_fdr} \
      ${tempdir_arg} \
      --n_jobs ${task.cpus}
    """
}

process ATAC_MAKE_THRESHOLDS {
    tag "ATAC_Thresholds"
    label 'process_low'
    publishDir "${params.outdir}/atac/initial_qc", mode: 'copy'

    input:
    path stats_file  // Expects sample_statistics.csv from ATAC_INITIAL_QC

    output:
    path "sample_thresholds.json", emit: thresholds_file

    script:
    """
    python ${projectDir}/bin/make_atac_thresholds_from_csv.py \\
        --stats ${stats_file} \\
        --min_cells 100 \\
        --output sample_thresholds.json
    """
}
