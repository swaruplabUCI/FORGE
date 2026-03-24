// modules/scprint/extract_ccan_enhancers.nf
//
// EXTRACT_CCAN_ENHANCERS — Recipe Step A4
// Extract distal enhancer peaks from Cicero CCANs by intersecting with
// TSS annotations.  Non-promoter peaks within a CCAN are labelled as
// putative enhancers with gene linkage via the CCAN's promoter anchor.

process EXTRACT_CCAN_ENHANCERS {

    tag "extract_ccan_enhancers"
    label 'process_medium'
    publishDir "${params.outdir}/enhancer_footprinting/ccan_enhancers", mode: 'copy'

    input:
    path cicero_connections   // Cicero co-accessibility connections TSV
    path ccan_assignments     // CCAN assignment file from CICERO_FULL
    val gtf_path              // GTF path for TSS extraction

    output:
    path "ccan_enhancer_peaks.bed.gz",       emit: enhancer_peaks
    path "ccan_enhancer_gene_links.tsv",     emit: gene_links
    path "ccan_enhancer_summary.txt",        emit: summary

    script:
    """
    python ${projectDir}/bin/extract_ccan_enhancers.py \\
        --connections '${cicero_connections}' \\
        --ccan '${ccan_assignments}' \\
        --gtf '${gtf_path}' \\
        --promoter-upstream ${params.scprinter.promoter_upstream} \\
        --promoter-downstream ${params.scprinter.promoter_downstream} \\
        --outdir .
    """
}
