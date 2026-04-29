process ATAC_MOTIF_ANALYSIS {
    label 'process_high'
    publishDir "${params.outdir}/atac/motifs", mode: 'copy'

    input:
    path peak_matrix
    path annotated_data

    output:
    path "motif_enrichment_results.csv", emit: motifs

    script:
    """
    python3 ${baseDir}/bin/atac_motif_analysis.py \\
        --peak_matrix ${peak_matrix} \\
        --species ${params.species} \\
        --groupby 'cell_type' \\
        --max_fdr 0.0001 \\
        --output_dir .
    """
}
