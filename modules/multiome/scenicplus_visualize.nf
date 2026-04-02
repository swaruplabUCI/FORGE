// modules/multiome/scenicplus_visualize.nf
//
// SCENIC+ downstream visualization: RSS heatmaps, AUCell UMAPs,
// eRegulon correlation, top regulons per cell type, TF target summary.

nextflow.enable.dsl=2

process SCENICPLUS_VISUALIZE {
    label 'process_high'

    tag "scenicplus_viz"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/scenicplus/plots", mode: 'copy', overwrite: true

    input:
    path scplus_mudata       // scplusmdata.h5mu from SCENICPLUS_RUN
    path aucell_direct       // AUCell_direct.h5mu
    path eregulon_tsv        // eRegulon_direct.tsv
    val  cell_type_key       // Cell type column name

    output:
    path "eregulon_rss_heatmap.png",         emit: rss_heatmap
    path "eregulon_rss.csv",                 emit: rss_csv
    path "eregulon_aucell_umap.png",         emit: aucell_umap,        optional: true
    path "eregulon_aucell_umap_celltypes.png", emit: celltype_umap,    optional: true
    path "eregulon_correlation.png",         emit: correlation,        optional: true
    path "eregulon_top_per_celltype.png",    emit: top_per_celltype,   optional: true
    path "tf_target_summary.png",            emit: tf_summary,         optional: true
    path "*.pdf",                            emit: pdfs,               optional: true

    script:
    """
    export HDF5_USE_FILE_LOCKING=FALSE

    python ${projectDir}/bin/plot_scenicplus_results.py \\
        --scplus-mudata '${scplus_mudata}' \\
        --aucell-direct '${aucell_direct}' \\
        --eregulon-tsv '${eregulon_tsv}' \\
        --cell-type-key '${cell_type_key}' \\
        --top-n 10 \\
        --dpi 200 \\
        --outdir .
    """
}
