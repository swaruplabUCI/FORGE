// modules/atac/celltype_annotation.nf
//
// Tissue-aware ATAC cell type annotation using marker gene activity scores.
// Supports pbmc, brain, and general tissue types via configs/marker_genes.json.

process ATAC_CELLTYPE_ANNOTATION {
    label 'process_low'
    publishDir "${params.outdir}/atac/final", mode: 'copy'

    input:
    path cluster_scores  // cluster_avg_gene_scores.csv

    output:
    path "celltype_annotations.json", emit: annotations

    script:
    def tissue_arg = params.atac?.tissue_type ? "--tissue_type ${params.atac.tissue_type}" : "--tissue_type pbmc"
    def marker_arg = params.atac?.marker_file ? "--marker_file ${params.atac.marker_file}" : ""
    """
    python ${projectDir}/bin/atac_celltype_annotation.py \\
      --cluster_scores ${cluster_scores} \\
      --species ${params.species} \\
      --resolution 0_5 \\
      ${tissue_arg} \\
      ${marker_arg} \\
      --out celltype_annotations.json
    """
}
