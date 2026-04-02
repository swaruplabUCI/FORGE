// modules/atac/atac_celltypist.nf
//
// CellTypist-based ATAC cell type annotation using gene activity scores.
// Default mode when no custom marker file is provided.
// Reuses the same run_cell_typist.py script used for RNA annotation.

process ATAC_CELLTYPIST {
    tag "atac_celltypist"
    label 'process_medium'
    publishDir "${params.outdir}/atac/final", mode: 'copy'

    input:
    path gene_matrix   // gene_matrix.h5ad (log-normalized gene activity scores)

    output:
    path "atac_celltypist_annotations.h5ad", emit: annotated_gene_matrix

    script:
    def model_arg = params.celltypist?.model ? "--model ${params.celltypist.model}" : ""
    """
    export HOME=/tmp

    python ${projectDir}/bin/run_cell_typist.py \\
        --input ${gene_matrix} \\
        --output atac_celltypist_annotations.h5ad \\
        --species ${params.species} \\
        ${model_arg} \\
        --majority_voting
    """
}
