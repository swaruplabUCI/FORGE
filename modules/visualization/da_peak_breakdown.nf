// modules/visualization/da_peak_breakdown.nf
// Stacked-bar peak-type breakdown for differential peaks (Shi Fig 2D).

nextflow.enable.dsl=2

process DA_PEAK_BREAKDOWN {

    tag "da_peak_breakdown"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent/2D_peak_type_breakdown", mode: 'copy'

    input:
    path "diff/*"           // DA_peaks_<ct>__<trt>_vs_<ctrl>.csv (collected)
    path peak_annotation    // peaks_annotated.tsv

    output:
    path "da_peak_breakdown.{pdf,png}", emit: plot
    path "da_peak_breakdown.tsv",       emit: table

    script:
    def trt = params.shi_figures?.treatment ?: params.differential?.treatment_condition
    def ctrl = params.shi_figures?.control ?: params.differential?.control_condition
    def fdr = params.shi_figures?.fdr_cutoff ?: 0.05
    def lfc = params.shi_figures?.lfc_cutoff ?: 0.5
    """
    python ${projectDir}/bin/plot_da_peak_breakdown.py \\
        --diff-dir diff \\
        --peak-annotation '${peak_annotation}' \\
        --treatment '${trt}' \\
        --control '${ctrl}' \\
        --fdr-cutoff ${fdr} \\
        --lfc-cutoff ${lfc} \\
        --outdir .
    """
}
