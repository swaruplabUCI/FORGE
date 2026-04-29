// modules/atac/post_qc_report.nf
//
// Post-hoc QC report — re-reads the merged peak_matrix_annotated.h5ad and
// renders QC visualizations that the in-line pipeline does not currently
// produce (e.g. FORGE-style TSS-vs-log10(nFragment) marginal density scatter).
// Designed to be invokable from both the main pipeline and the viz-only entry
// workflow.

process POST_QC_REPORT {
    tag "post_qc"
    label 'process_medium'
    publishDir "${params.outdir}/atac/post_qc", mode: 'copy'

    input:
    path "samples/*"                // per-sample h5ad files from ATAC_INITIAL_QC.individual_samples
    val  tsse_min
    val  nfrag_min
    val  nfrag_max

    output:
    path "*.pdf",                  emit: plots
    path "*.csv",                  emit: tables,   optional: true
    path "post_qc_manifest.json",  emit: manifest

    script:
    """
    python ${projectDir}/bin/post_qc_report.py \\
        --h5ad samples \\
        --tsse-min ${tsse_min} \\
        --nfrag-min ${nfrag_min} \\
        --nfrag-max ${nfrag_max} \\
        --outdir .
    """
}
