// modules/scprint/extract_ccan_enhancers.nf
//
// EXTRACT_CCAN_ENHANCERS — Recipe Step A4
// Extract distal enhancer peaks from Cicero CCANs by intersecting with
// TSS annotations.  Non-promoter peaks within a CCAN are labelled as
// putative enhancers with gene linkage via the CCAN's promoter anchor.
//
// Supports two modes:
//   'ccan'         — CCAN membership (permissive, original Recipe A)
//   'pairwise_95'  — Shi et al. 95th percentile co-accessibility threshold

process EXTRACT_CCAN_ENHANCERS {

    tag "extract_ccan_enhancers${condition_label ? '_' + condition_label : ''}"
    label 'process_medium'
    publishDir "${params.outdir}/enhancer_footprinting/ccan_enhancers", mode: 'copy'

    input:
    path cicero_connections   // Cicero co-accessibility connections TSV
    path ccan_assignments     // CCAN assignment file (can be NO_CCAN for pairwise_95)
    val gtf_path              // GTF path for TSS extraction
    val condition_label       // Optional: condition label for disease-stratified mode

    output:
    path "ccan_enhancer_peaks*.bed.gz",       emit: enhancer_peaks
    path "ccan_enhancer_gene_links*.tsv",     emit: gene_links
    path "ccan_enhancer_summary*.txt",        emit: summary

    script:
    def mode = params.enhancer_footprinting.enhancer_mode ?: 'ccan'
    def pct = params.enhancer_footprinting.percentile_threshold ?: 95
    def ccan_arg = mode == 'ccan' ? "--ccan '${ccan_assignments}'" : ""
    def cond_arg = condition_label ? "--condition-label '${condition_label}'" : ""
    """
    python ${projectDir}/bin/extract_ccan_enhancers.py \\
        --connections '${cicero_connections}' \\
        ${ccan_arg} \\
        --gtf '${gtf_path}' \\
        --promoter-upstream ${params.scprinter.promoter_upstream} \\
        --promoter-downstream ${params.scprinter.promoter_downstream} \\
        --mode ${mode} \\
        --percentile-threshold ${pct} \\
        ${cond_arg} \\
        --outdir .
    """
}
