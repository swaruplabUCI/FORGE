// modules/visualization/tf_differential_volcano.nf
// Volcano plots of differential TF activity per cell type (Shi Fig 4B/5A/5C/5E).

nextflow.enable.dsl=2

process TF_DIFFERENTIAL_VOLCANO {

    tag "tf_differential_volcano"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent/4B5A5C5E_tf_volcano", mode: 'copy'

    input:
    path "diff_tf/*"          // tf_differential_<ct>_<trt>_vs_<ctrl>.csv

    output:
    path "*_volcano.{pdf,png}",         emit: volcanoes
    path "combined_volcano.{pdf,png}",  emit: combined
    path "volcano_manifest.json",       emit: manifest

    script:
    def trt = params.shi_figures?.treatment ?: params.differential?.treatment_condition
    def ctrl = params.shi_figures?.control ?: params.differential?.control_condition
    def fdr = params.shi_figures?.fdr_cutoff ?: 0.05
    def lfc = params.shi_figures?.lfc_cutoff ?: 0.5
    def labels = params.shi_figures?.volcano_top_labels ?: 10
    """
    python ${projectDir}/bin/plot_tf_differential_volcano.py \\
        --diff-tf-dir diff_tf \\
        --treatment '${trt}' \\
        --control '${ctrl}' \\
        --fdr-cutoff ${fdr} \\
        --lfc-cutoff ${lfc} \\
        --top-labels ${labels} \\
        --outdir .
    """
}
