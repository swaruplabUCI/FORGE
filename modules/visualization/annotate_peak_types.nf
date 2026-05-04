// modules/visualization/annotate_peak_types.nf
// Classify peaks as Promoter / Exonic / Intronic / Distal from a GTF.
// Feeds the 2D and 2E figure-equivalent plots.

nextflow.enable.dsl=2

process ANNOTATE_PEAK_TYPES {

    tag "annotate_peak_types"
    label 'process_medium'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent/peak_annotation", mode: 'copy'

    input:
    path peaks            // BED.gz or h5ad with peak ids
    path gtf

    output:
    path "peaks_annotated.tsv",     emit: tsv
    path "peaks_annotated.summary", emit: summary

    script:
    def window = params.shi_figures?.promoter_window ?: 2000
    """
    python ${projectDir}/bin/annotate_peak_types.py \\
        --peaks '${peaks}' \\
        --gtf '${gtf}' \\
        --promoter-window ${window} \\
        --out peaks_annotated.tsv \\
        --summary peaks_annotated.summary
    """
}
