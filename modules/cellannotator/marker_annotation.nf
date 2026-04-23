// modules/cellannotator/marker_annotation.nf
//
// Marker-gene-based cell type annotation for RNA data.
// Used when no CellTypist model fits the tissue (e.g., mouse kidney).
// Scores each cell against curated per-cell-type marker lists with
// sc.tl.score_genes and assigns labels by argmax under min-score and
// top-vs-second-best margin gates. Emits the same (annotated_h5ad, report)
// tuple as RUN_CELLTYPIST so downstream wiring is interchangeable.

process RUN_MARKER_ANNOTATION {
    tag "marker_annotation"
    label 'process_medium'
    publishDir "${params.outdir}/cell_annotation", mode: 'copy'

    input:
    path rna_h5ad
    path marker_file

    output:
    path "marker_annotated.h5ad",  emit: annotated_h5ad
    path "*_marker_report.txt",    emit: report

    script:
    def min_score  = params.rna?.marker_min_score    ?: 0.0
    def min_margin = params.rna?.marker_score_margin ?: 0.1
    """
    export HOME=/tmp

    python ${projectDir}/bin/run_marker_annotation.py \\
        --input ${rna_h5ad} \\
        --markers ${marker_file} \\
        --output marker_annotated.h5ad \\
        --min_score ${min_score} \\
        --min_margin ${min_margin}
    """
}
