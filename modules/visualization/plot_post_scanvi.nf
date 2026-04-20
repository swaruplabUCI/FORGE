process PLOT_POST_SCANVI {
    label 'process_medium'
    publishDir "${params.outdir}/rna/post_integration_plots", mode: 'copy'

    input:
    path(h5ad)
    val(cell_type_key)

    output:
    path("*.png"), emit: plots
    path("annotated_with_celltype.h5ad"), emit: annotated_updated
    path("hdwgcna_cell_types.txt"), emit: cell_types_file
    path("hdwgcna_cell_type_counts.csv"), emit: cell_type_counts

    script:
    """
    plot_postscanvi.py \\
        --input ${h5ad} \\
        --output_dir . \\
        --cell_type_key ${cell_type_key} \\
        --hdwgcna_min_cells ${params.hdwgcna.min_cells} \\
        --tissue_type ${params.celltypist.tissue_type ?: 'pbmc'}
    """
}
