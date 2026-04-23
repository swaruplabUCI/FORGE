// modules/scprint/plot_tf_gene_network.nf
//
// PLOT_TF_GENE_NETWORK — graph-tool visualization of the SDas TF-gene
// regulatory network (downstream of BUILD_TF_GENE_NETWORK).
//
// Renders three views (PDF + PNG): per-cell-type GRNs, a holistic GRN
// (union of per-celltype top TFs), and a TF co-regulation network.
// Runs in scenicplus.sif via its /opt/miniforge3/envs/gt python env.

process PLOT_TF_GENE_NETWORK {

    tag "tf_gene_network_plots"
    label 'process_medium'
    publishDir "${params.outdir}/enhancer_footprinting/network/plots", mode: 'copy'

    input:
    path adjacency_tsv          // tf_gene_adjacency.tsv from BUILD_TF_GENE_NETWORK

    output:
    path "grn_holistic.pdf",             emit: holistic_pdf,    optional: true
    path "grn_holistic.png",             emit: holistic_png,    optional: true
    path "grn_holistic_legend.png",      emit: holistic_legend, optional: true
    path "grn_tf_coregulation.pdf",      emit: coreg_pdf,       optional: true
    path "grn_tf_coregulation.png",      emit: coreg_png,       optional: true
    path "per_celltype/*.pdf",           emit: per_ct_pdfs,     optional: true
    path "per_celltype/*.png",           emit: per_ct_pngs,     optional: true

    script:
    def top_n        = params.enhancer_footprinting.network_plot_top_n     ?: 5
    def max_targets  = params.enhancer_footprinting.network_plot_max_targets ?: 15
    def min_shared   = params.enhancer_footprinting.network_plot_min_shared  ?: 5
    """
    export HDF5_USE_FILE_LOCKING=FALSE

    GT_PYTHON="/opt/miniforge3/envs/gt/bin/python"
    if [ ! -x "\${GT_PYTHON}" ]; then
        echo "WARNING: graph-tool conda env not found at \${GT_PYTHON}, falling back to system python"
        GT_PYTHON="python"
    fi

    \${GT_PYTHON} ${projectDir}/bin/plot_tf_gene_network.py \\
        --adjacency '${adjacency_tsv}' \\
        --top-n ${top_n} \\
        --max-targets-per-tf ${max_targets} \\
        --min-shared-coreg ${min_shared} \\
        --outdir .
    """
}
