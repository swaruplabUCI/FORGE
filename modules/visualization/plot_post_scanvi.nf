process PLOT_POST_SCANVI {
    publishDir "${params.outdir}/rna/post_integration_plots", mode: 'copy'
    
    input:
    path(h5ad)
    
    output:
    path("*.png"), emit: plots
    path("annotated_with_scanvi_clustering.h5ad"), emit: annotated_updated  // ← Explicit filename
    path("hdwgcna_cell_types.txt"), emit: cell_types_file
    path("hdwgcna_cell_type_counts.csv"), emit: cell_type_counts
    
    script:
    """
    plot_postscanvi.py \\
        --input ${h5ad} \\
        --output_dir . \\
        --hdwgcna_min_cells ${params.hdwgcna.min_cells}
    """
}
