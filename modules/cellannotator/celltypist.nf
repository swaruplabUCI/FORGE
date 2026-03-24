// modules/cellannotator/celltypist.nf
//
// CellTypist-based cell type annotation for RNA data.
// Used as fallback when no reference atlas is provided for scANVI.

process RUN_CELLTYPIST {
    tag "celltypist"
    label 'process_medium'
    publishDir "${params.outdir}/cell_annotation", mode: 'copy'

    input:
    path rna_h5ad

    output:
    path "celltypist_annotated.h5ad", emit: annotated_h5ad
    path "*_celltypist_report.txt",   emit: report

    script:
    def model_arg = params.celltypist?.model ? "--model ${params.celltypist.model}" : ""
    """
    export HOME=/tmp

    python ${projectDir}/bin/run_cell_typist.py \\
        --input ${rna_h5ad} \\
        --output celltypist_annotated.h5ad \\
        --species ${params.species} \\
        ${model_arg} \\
        --majority_voting
    """
}
