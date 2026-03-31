// modules/multiome/scenicplus_grn_viz.nf
//
// SCENIC+ GRN network visualization with graph-tool.
// Generates static PDF + PNG graphs: holistic GRN, per-cell-type GRNs,
// and TF co-regulation network.

nextflow.enable.dsl=2

process SCENICPLUS_GRN_VIZ {

    tag "scenicplus_grn"
    publishDir "${params.outdir}/scenicplus/grn_networks", mode: 'copy', overwrite: true

    input:
    path eregulon_tsv        // eRegulon_direct.tsv
    path aucell_direct       // AUCell_direct.h5mu
    path rss_csv             // eregulon_rss.csv from SCENICPLUS_VISUALIZE
    path tf_to_gene_adj      // tf_to_gene_adj.tsv from SCENICPLUS_RUN
    val  cell_type_key       // Cell type column name

    output:
    path "*.pdf",                               emit: pdfs,               optional: true
    path "*.png",                               emit: pngs,               optional: true
    path "per_celltype/*.pdf",                  emit: per_ct_pdfs,        optional: true
    path "per_celltype/*.png",                  emit: per_ct_pngs,        optional: true

    script:
    """
    export HDF5_USE_FILE_LOCKING=FALSE

    # graph-tool requires the gt conda env; fall back to system python if missing
    GT_PYTHON="/opt/miniforge3/envs/gt/bin/python"
    if [ ! -x "\${GT_PYTHON}" ]; then
        echo "WARNING: graph-tool conda env not found at \${GT_PYTHON}, falling back to system python"
        GT_PYTHON="python"
    fi

    \${GT_PYTHON} ${projectDir}/bin/plot_eregulon_networks.py \\
        --eregulon-tsv '${eregulon_tsv}' \\
        --aucell-direct '${aucell_direct}' \\
        --rss-csv '${rss_csv}' \\
        --tf-to-gene-adj '${tf_to_gene_adj}' \\
        --cell-type-key '${cell_type_key}' \\
        --top-n 5 \\
        --outdir .
    """
}
