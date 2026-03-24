// modules/multiome/scprinter_dorc.nf

nextflow.enable.dsl=2

process SCPRINTER_DORC {

    tag { "scprinter_dorc" }

    publishDir "${params.outdir}/dorc", mode: 'copy', overwrite: true

    // Use scprinter-capable env (should contain scprinter, scenic+, etc.)
    // Conda handled in config

    input:
    // ATAC AnnData with peaks x cells (e.g., ATAC_FINAL.out.peak_matrix)
    path atac_h5ad
    // RNA AnnData with counts for same cells (e.g., PLOT_POST_SCANVI.out.annotated_updated)
    path rna_h5ad
    // Optional label / covariate name (string); can be "none"
    val covariate_label

    output:
    path "dorc_all.tsv.gz",          emit: dorc_all
    path "dorc_significant.tsv.gz",  emit: dorc_sig
    path "dorc_scores.h5ad",         emit: dorc_scores
    path "dorc_jplot_genes.txt",     emit: dorc_genes, optional: true
    path "dorc_jplot.pdf",           emit: dorc_jplot_pdf, optional: true
    path "dorc_jplot.png",           emit: dorc_jplot_png, optional: true
    path "dorc_umap_plots/*.pdf",    emit: dorc_umap_pdfs, optional: true
    path "dorc_umap_plots/*.png",    emit: dorc_umap_pngs, optional: true
    path "dorc_logs.txt",            emit: dorc_log,   optional: true

    script:
    """
    python ${projectDir}/bin/run_scprinter_dorc.py \\
        --atac '${atac_h5ad}' \\
        --rna '${rna_h5ad}' \\
        --genome '${params.scprinter.genome}' \\
        --covariate '${covariate_label}' \\
        --cache-dir '${params.scprinter.cache_dir}' \\
        --outdir '.'
    """
}
