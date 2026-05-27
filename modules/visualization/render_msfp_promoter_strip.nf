// modules/visualization/render_msfp_promoter_strip.nf
//
// Promoter MSFP strip: 150 bp zoom window centred on the best TF-motif hit,
// with motif logo inset, reference sequence, and genomic coordinates.
// One figure per (cell_type, TF) grouping; all target genes stacked vertically.
//
// Modes:
//   absolute    — arr[0] heatmap only
//   differential— Δ = arr[1] − arr[0]
//   all_three   — ctrl | trt | Δ stacked (requires ≥2 condition slices in h5ad)
//
// Input h5ads come from PROMOTER_MSFP_PER_CT (run_promoter_msfp_per_condition.py)
// which writes one h5ad per gene into params.outdir/promoter_msfp_overlay/per_cond_h5ad/.
// They are staged into scan_dir/ so the script can glob-find them.
//
// Output: {gene_tag}__promoter_strip_{mode}.png/.pdf

nextflow.enable.dsl=2

process RENDER_MSFP_PROMOTER_STRIP {

    tag "${cell_type}__${tfs}__${mode}"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/msfp_strips/promoter/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}",
               mode: 'copy'

    input:
    tuple val(cell_type), val(tfs), val(genes), path(h5ads, stageAs: 'scan_dir/*'), val(mode)
    val ctrl_condition
    val trt_condition

    output:
    path "promoter_strip_*.png", emit: png
    path "promoter_strip_*.pdf", emit: pdf, optional: true

    script:
    def safe_ct   = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_tfs  = tfs.replaceAll(/[\/\s\(\)]+/, '_')
    def mode_flag = mode ?: 'absolute'
    def ctrl_arg  = (ctrl_condition && ctrl_condition != 'none') ?
                    "--control-condition '${ctrl_condition}'"   : ''
    def trt_arg   = (trt_condition  && trt_condition  != 'none') ?
                    "--treatment-condition '${trt_condition}'"  : ''
    def genome    = params.scprinter?.genome ?: 'mm10'
    def max_scale = params.promoter_overlay?.max_scale ?: 30
    """
    python ${projectDir}/bin/render_msfp_promoter_strip.py \\
        --scan-dir      scan_dir \\
        --genes         '${genes}' \\
        --tfs           '${tfs}' \\
        --cell-type     '${cell_type}' \\
        --pfm           '${params.scprinter.pfms}' \\
        --cache-dir     '${params.scprinter.cache_dir}' \\
        --genome        '${genome}' \\
        --max-scale     ${max_scale} \\
        --mode          '${mode_flag}' \\
        ${ctrl_arg} \\
        ${trt_arg} \\
        --out-png       'promoter_strip_${safe_tfs}__${safe_ct}__${mode_flag}.png'
    """

    stub:
    def safe_ct  = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_tfs = tfs.replaceAll(/[\/\s\(\)]+/, '_')
    def m        = mode ?: 'absolute'
    """
    touch "promoter_strip_${safe_tfs}__${safe_ct}__${m}.png" \
          "promoter_strip_${safe_tfs}__${safe_ct}__${m}.pdf"
    """
}
