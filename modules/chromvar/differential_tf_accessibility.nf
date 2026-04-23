// modules/chromvar/differential_tf_accessibility.nf
//
// Per-celltype TF motif accessibility testing on the ChromVAR deviation matrix.
// Two modes (selected by params.differential_tf.mode):
//   'differential' — Shi et al. style Wilcoxon of TF devs between two
//                    conditions within a cell type.
//                    Input tuple:  (cell_type, treatment, control)
//   'descriptive'  — Cell type vs rest-of-cells; single-condition friendly.
//                    Input tuple:  (cell_type, condition_filter_or_"all", "rest")
//                    where the 3rd value is an unused placeholder for interface
//                    symmetry with differential mode.

process DIFFERENTIAL_TF_ACCESSIBILITY {
    tag "${cell_type}_${trt}_vs_${ctrl}"
    label 'process_medium'
    publishDir "${params.outdir}/differential_tf", mode: 'copy'

    input:
    path chromvar_dev                            // chromvar_matrix_deviations.h5ad
    path peak_matrix_annotated                   // for cell_type + condition metadata
    tuple val(cell_type), val(trt), val(ctrl)

    output:
    path "tf_differential_*.csv",     emit: tf_diff,      optional: true
    path "tf_differential_*.summary", emit: summaries,    optional: true
    path "tf_descriptive_*.csv",      emit: tf_desc,      optional: true
    path "tf_descriptive_*.summary",  emit: desc_summary, optional: true

    script:
    def mode      = params.differential_tf?.mode ?: 'differential'
    def cond_key  = params.differential_tf?.condition_key ?:
                    (params.differential?.condition_key ?: 'condition')
    def min_cells = params.differential_tf?.min_cells ?:
                    (params.differential?.min_cells ?: 50)
    def fdr_cut   = params.differential_tf?.fdr_cutoff ?:
                    (params.differential?.fdr_cutoff ?: 0.05)
    def trt_arg   = (mode == 'descriptive' && (!trt || trt == 'all')) ?
                    '' : "--treatment \"${trt}\""
    def ctrl_arg  = (mode == 'descriptive') ? '' : "--control \"${ctrl}\""
    """
    python ${projectDir}/bin/differential_tf_accessibility.py \\
        --chromvar-dev ${chromvar_dev} \\
        --annotated-peaks ${peak_matrix_annotated} \\
        --cell-type "${cell_type}" \\
        --mode ${mode} \\
        --condition-key ${cond_key} \\
        ${trt_arg} \\
        ${ctrl_arg} \\
        --min-cells ${min_cells} \\
        --fdr-cutoff ${fdr_cut}
    """
}
