// modules/visualization/locus_tf_binding.nf
// Locus-restricted differential TF binding bars (Shi Fig 4E approximation).

nextflow.enable.dsl=2

process LOCUS_TF_BINDING {

    tag "locus_tf_binding"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent/4E_locus_tf_binding", mode: 'copy'

    input:
    path adjacency
    path "diff_tf/*"
    path candidates

    output:
    path "*_locus_tf_binding.{pdf,png}",  emit: plots
    path "locus_tf_binding_summary.tsv",  emit: summary

    script:
    def trt = params.shi_figures?.treatment ?: params.differential?.treatment_condition
    def ctrl = params.shi_figures?.control ?: params.differential?.control_condition
    def top_tfs = params.shi_figures?.locus_top_tfs ?: 15
    """
    python ${projectDir}/bin/plot_locus_tf_binding.py \\
        --adjacency '${adjacency}' \\
        --diff-tf-dir diff_tf \\
        --candidates '${candidates}' \\
        --treatment '${trt}' \\
        --control '${ctrl}' \\
        --top-tfs-per-gene ${top_tfs} \\
        --outdir .
    """
}
