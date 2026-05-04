// modules/visualization/coacc_correlation_matrix.nf
// Co-accessibility correlation summary (Shi Fig 2C cleanup).

nextflow.enable.dsl=2

process COACC_CORRELATION_MATRIX {

    tag "coacc_correlation_matrix"
    label 'process_medium'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent/2C_coacc_correlation", mode: 'copy'

    input:
    path connections_ctrl, stageAs: 'cicero_connections_ctrl.tsv.gz'
    path connections_trt,  stageAs: 'cicero_connections_trt.tsv.gz'
    path peak_annotation

    output:
    path "coacc_correlation_matrix.{pdf,png}", emit: figure
    path "coacc_correlation_summary.tsv",      emit: summary

    script:
    def trt = params.shi_figures?.treatment ?: params.differential?.treatment_condition
    def ctrl = params.shi_figures?.control ?: params.differential?.control_condition
    """
    python ${projectDir}/bin/plot_coaccessibility_matrix.py \\
        --connections-ctrl '${connections_ctrl}' \\
        --connections-trt  '${connections_trt}' \\
        --peak-annotation  '${peak_annotation}' \\
        --control-label    '${ctrl}' \\
        --treatment-label  '${trt}' \\
        --outdir .
    """
}
