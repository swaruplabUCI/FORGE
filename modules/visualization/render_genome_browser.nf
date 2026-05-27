// modules/visualization/render_genome_browser.nf
//
// Multi-track ATAC genome browser rendered with matplotlib + pyBigWig.
// No pyGenomeTracks, no CCAN arc overlay.
//
// Modes:
//   absolute    — per-CT fill tracks, colored by cell class.
//   differential— per-CT overlay tracks:
//                   grey   = shared accessibility (min of ctrl and trt)
//                   red    = trt-unique
//                   blue   = ctrl-unique
//                 Requires condition-split BigWigs in manifest.by_condition.
//                 These are produced by EXPORT_ATAC_BIGWIGS when --condition-col
//                 is provided (params.enhancer_viz.condition_col set).
//
// BigWig files are staged into bigwigs/ so manifest relative paths resolve.
// An optional DA peaks BED adds a colored annotation strip at the top.
//
// Output: browser_{gene}_{mode}.png/.pdf

nextflow.enable.dsl=2

process RENDER_GENOME_BROWSER {

    tag "${gene}__${mode}"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/genome_browser", mode: 'copy'

    input:
    tuple val(gene), path(bw_manifest), path("bigwigs/*"), val(cell_types), val(mode)
    path gtf
    val  ctrl_condition
    val  trt_condition
    path da_peaks    // pass file('NO_FILE') when not available

    output:
    path "browser_*.png", emit: png
    path "browser_*.pdf", emit: pdf, optional: true

    script:
    def mode_flag     = mode ?: 'absolute'
    def safe_gene     = gene.replaceAll(/[\/\s\(\)]+/, '_')
    def window        = params.browser_viz?.window ?: 15000
    def n_bins        = params.browser_viz?.n_bins  ?: 700
    def ctrl_arg      = (ctrl_condition && ctrl_condition != 'none') ?
                        "--ctrl-condition '${ctrl_condition}'" : ''
    def trt_arg       = (trt_condition  && trt_condition  != 'none') ?
                        "--trt-condition '${trt_condition}'"  : ''
    def manifest_arg  = "--bw-manifest '${bw_manifest}'"
    def da_arg        = (da_peaks && da_peaks.name != 'NO_FILE') ?
                        "--da-peaks '${da_peaks}'" : ''
    def max_val_arg   = (params.browser_viz?.max_value) ?
                        "--max-value ${params.browser_viz.max_value}" : ''
    """
    python ${projectDir}/bin/render_genome_browser.py \\
        --gene        '${gene}' \\
        --gtf         '${gtf}' \\
        ${manifest_arg} \\
        --bw-dir      bigwigs \\
        --cell-types  '${cell_types}' \\
        --mode        '${mode_flag}' \\
        --window      ${window} \\
        --n-bins      ${n_bins} \\
        ${ctrl_arg} \\
        ${trt_arg} \\
        ${da_arg} \\
        ${max_val_arg} \\
        --out-png     'browser_${safe_gene}_${mode_flag}.png'
    """

    stub:
    def safe_gene = gene.replaceAll(/[\/\s\(\)]+/, '_')
    def m         = mode ?: 'absolute'
    """
    touch "browser_${safe_gene}_${m}.png" "browser_${safe_gene}_${m}.pdf"
    """
}
