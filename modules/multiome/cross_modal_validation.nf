// modules/multiome/cross_modal_validation.nf
//
// CROSS_MODAL_VALIDATION — Recipe Step B6
// Correlate footprint depth with TF expression and target gene expression
// across cell states.  Classifies each TF as activator/repressor/indirect.

process CROSS_MODAL_VALIDATION {

    tag "cross_modal_validation"
    label 'process_medium'
    publishDir "${params.outdir}/multiome/cross_modal_validation", mode: 'copy'

    input:
    path footprint_h5ads      // Collected enhancer footprint h5ads
    path rna_h5ad             // Integrated RNA h5ad
    path target_genes_json    // eregulon_target_genes.json or DORC linkage

    output:
    path "cross_modal_validation.tsv",     emit: validation_table
    path "validation_plots/",              emit: plots, optional: true
    path "validation_summary.txt",         emit: summary

    script:
    """
    python ${projectDir}/bin/cross_modal_validation.py \\
        --footprint-dir '.' \\
        --footprint-h5ads ${footprint_h5ads} \\
        --rna-h5ad '${rna_h5ad}' \\
        --target-genes '${target_genes_json}' \\
        --outdir .
    """
}
