// modules/visualization/render_msfp_enhancer_strip.nf
//
// Enhancer MSFP strip: per-peak multi-scale footprint heatmaps for a target
// gene's Cicero-linked CCAN enhancers, with motif logos and reference
// coordinates.
//
// The h5ad comes from ENHANCER_FOOTPRINTING_PER_CT (one file per (CT, TF)):
//   results/scprinter/footprints/enhancer_footprints_{CT}_{TF}.h5ad
// Cicero connections and GTF are optional; when absent, TSS proximity
// filtering and region linking are skipped gracefully.
//
// Modes:
//   absolute    — arr[0] heatmap (single condition)
//   differential— Δ = arr[1] − arr[0] (requires ≥2 condition slices)
//
// Output: enhancer_strip_{TF}__{CT}__{mode}.png/.pdf

nextflow.enable.dsl=2

process RENDER_MSFP_ENHANCER_STRIP {

    tag "${cell_type}__${tf}__${target_gene}__${mode}"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/msfp_strips/enhancer/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}",
               mode: 'copy'

    input:
    tuple val(cell_type), val(tf), val(target_gene), path(enhancer_h5ad), path(cicero_gz), val(mode)
    path gtf
    val ctrl_condition
    val trt_condition

    output:
    path "enhancer_strip_*.png", emit: png
    path "enhancer_strip_*.pdf", emit: pdf, optional: true

    script:
    def safe_ct   = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_tf   = tf.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_gene = target_gene.replaceAll(/[\/\s\(\)]+/, '_')
    def mode_flag = mode ?: 'absolute'
    def genome    = params.scprinter?.genome ?: 'mm10'
    def ctx_bp    = params.msfp_strip?.context_bp ?: 500000
    def gtf_arg   = (gtf && gtf.name != 'NO_FILE') ? "--gtf '${gtf}'" : ''
    def cic_arg   = (cicero_gz && cicero_gz.name != 'NO_FILE') ?
                    "--cicero-connections '${cicero_gz}'" : ''
    def gene_arg  = target_gene ? "--target-gene '${target_gene}'" : ''
    """
    python ${projectDir}/bin/render_msfp_enhancer_strip.py \\
        --enhancer-h5ad '${enhancer_h5ad}' \\
        --tfs           '${tf}' \\
        --pfm           '${params.scprinter.pfms}' \\
        --cache-dir     '${params.scprinter.cache_dir}' \\
        --genome        '${genome}' \\
        --context-bp    ${ctx_bp} \\
        --mode          '${mode_flag}' \\
        ${gene_arg} \\
        ${gtf_arg} \\
        ${cic_arg} \\
        --out-png 'enhancer_strip_${safe_tf}__${safe_ct}__${safe_gene}__${mode_flag}.png'
    """

    stub:
    def safe_ct   = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_tf   = tf.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_gene = target_gene.replaceAll(/[\/\s\(\)]+/, '_')
    def m         = mode ?: 'absolute'
    """
    touch "enhancer_strip_${safe_tf}__${safe_ct}__${safe_gene}__${m}.png" \
          "enhancer_strip_${safe_tf}__${safe_ct}__${safe_gene}__${m}.pdf"
    """
}
