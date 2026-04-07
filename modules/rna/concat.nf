process CONCAT_BATCHES {
    tag "concatenating_all"
    label 'process_high'
    publishDir "${params.outdir}/rna/concatenated", mode: 'copy'
    
    input:
    path h5ad_files
    
    output:
    path "concatenated.h5ad", emit: concatenated
    path "concat_metrics.csv", emit: metrics
    path "concat_plots/*", emit: plots, optional: true
    
    script:
    """
    # v2: fixed suffix stripping regex and removed legacy batch inference
    python ${baseDir}/bin/concat_batches.py \
        --inputs ${h5ad_files} \
        --output concatenated.h5ad \
        --metrics concat_metrics.csv \
        --plots_dir concat_plots
    """
}
