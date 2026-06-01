// modules/visualization/render_cicero_lollipop.nf
//
// Cicero lollipop chart: TF-motif enhancer → promoter CCAN arcs per gene,
// comparing ctrl vs trt for a single cell type.
//
// For each condition the module reads:
//   {ccan_base}/{condition}/cicero_connections.tsv.gz
// These are constructed in the script from the two gz path inputs by symlinking
// them into a ccan_base/ directory with the expected subdirectory layout.
//
// The TF-motif BED (motif_peaks_bed) is produced upstream by
// MOTIF_SCAN_ENHANCERS or MOTIF_SCAN_ENHANCERS_UNION for the target TF; it
// lists every enhancer peak that contains a binding motif for that TF.
//
// Output per (cell_type, tf):
//   cicero_lollipop_{tf}__{ct}.pdf/.png
//   cicero_lollipop_{tf}__{ct}_genes.tsv
//
// Gate: only meaningful in differential runs; caller must ensure both gz files
// are real (not NO_FILE sentinels) before emitting into this channel.

nextflow.enable.dsl=2

process RENDER_CICERO_LOLLIPOP {

    tag "${cell_type}__${tf}"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/cicero_lollipop/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}",
               mode: 'copy'

    input:
    tuple val(cell_type), val(tf), path(motif_peaks_bed),
          path(ccan_ctrl_gz, stageAs: 'ccan_ctrl.tsv.gz'),
          path(ccan_trt_gz,  stageAs: 'ccan_trt.tsv.gz')
    path gtf
    val  ctrl_condition
    val  trt_condition

    output:
    tuple val(cell_type), val(tf),
          path("cicero_lollipop_*.pdf"),  emit: pdf,  optional: true
    tuple val(cell_type), val(tf),
          path("cicero_lollipop_*.png"),  emit: png,  optional: true
    tuple val(cell_type), val(tf),
          path("cicero_lollipop_*_genes.tsv"), emit: tsv, optional: true

    script:
    def safe_ct   = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_tf   = tf.replaceAll(/[\/\s\(\)]+/, '_')
    def out_tag   = "cicero_lollipop_${safe_tf}__${safe_ct}"
    def min_coacc = params.cis_rewiring?.min_coacc ?: 0.15
    def top_n     = params.cis_rewiring?.top_n_genes ?: 40
    """
    # Build the directory structure expected by render_cicero_lollipop.py:
    #   ccan_base/{ctrl_condition}/cicero_connections.tsv.gz
    #   ccan_base/{trt_condition}/cicero_connections.tsv.gz
    mkdir -p "ccan_base/${ctrl_condition}" "ccan_base/${trt_condition}"
    cp -L "ccan_ctrl.tsv.gz" "ccan_base/${ctrl_condition}/cicero_connections.tsv.gz"
    cp -L "ccan_trt.tsv.gz"  "ccan_base/${trt_condition}/cicero_connections.tsv.gz"

    python ${projectDir}/bin/render_cicero_lollipop.py \\
        --motif-peaks     '${motif_peaks_bed}' \\
        --ccan-base       ccan_base \\
        --cell-type       '${cell_type}' \\
        --gtf             '${gtf}' \\
        --ctrl-condition  '${ctrl_condition}' \\
        --trt-condition   '${trt_condition}' \\
        --min-coacc       ${min_coacc} \\
        --top-n           ${top_n} \\
        --out-tag         '${out_tag}' \\
        --outdir          .
    """

    stub:
    def safe_ct = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_tf = tf.replaceAll(/[\/\s\(\)]+/, '_')
    def tag     = "cicero_lollipop_${safe_tf}__${safe_ct}"
    """
    touch "${tag}.pdf" "${tag}.png" "${tag}_genes.tsv"
    """
}
