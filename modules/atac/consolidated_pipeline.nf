process ATAC_FINAL_PIPELINE {
    tag "ATAC_Final"
    label 'process_high'
    publishDir "${params.outdir}/atac/final", mode: 'copy'

    input:
    path fragment_files
    path metadata
    path thresholds_file

    output:
    path "atac_complete.h5ads", emit: anndataset
    path "peak_matrix.h5ad", emit: peak_matrix
    path "qc_plots/*", emit: qc_plots
    path "atac_pipeline_summary.json", emit: summary
    path "*.h5ad", emit: individual_samples
    path "cluster_avg_gene_scores.csv", emit: cluster_scores

    script:
    // FIX-29b: Pass local GTF annotation to avoid internet download on HPC compute nodes.
    // snap.genome.hg38 lazily downloads gencode GFF3 from ftp.ebi.ac.uk which fails
    // on offline compute nodes. The --gtf flag uses the local gencode GTF instead.
    def gtf_path = params.species == 'human' ?
        (params.gtf_human_full ?: '') :
        (params.gtf_mouse_full ?: '')
    def gtf_arg = gtf_path ? "--gtf ${gtf_path}" : ""

    """
    python \${PROJECTDIR:-${projectDir}}/bin/atac_consolidated_pipeline.py \\
        --fragment_files ${fragment_files} \\
        --metadata ${metadata} \\
        --species ${params.species} \\
        --output_dir . \\
        --thresholds_file ${thresholds_file} \\
        --min_fragments ${params.atac.min_fragments} \\
        --n_features ${params.atac.n_features} \\
        --batch_correction ${params.atac.batch_correction} \\
        --batch_key sample \\
        --clustering_resolutions ${params.atac.clustering_resolutions.join(' ')} \\
        --peak_fdr ${params.atac.peak_fdr} \\
        --tempdir ${params.tempdir} \\
        --n_jobs ${task.cpus} \\
        ${gtf_arg}
    """
}
