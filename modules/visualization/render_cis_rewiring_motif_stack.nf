// modules/visualization/render_cis_rewiring_motif_stack.nf
//
// RENDER_CIS_REWIRING_MOTIF_STACK — companion to MOTIF_IN_GAINED_CCANS.
// Per-TF stacked bar plot keyed on Δ-link ranking (NOT motif-match ranking,
// per user's explicit ordering preference). The treatment bar is split into
// shared / gained-no-motif / gained-with-motif so the figure preserves the
// finding "newly accessible enhancers don't host the TF directly" without
// sorting by it. Same upstream inputs as motif_in_gained_ccans (per-condition
// stratified Cicero CCAN→gene-link TSVs + union-peakset motif BED).

process RENDER_CIS_REWIRING_MOTIF_STACK {

    tag "${tf}"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/cis_rewiring/${tf}", mode: 'copy'

    input:
    tuple val(tf), path(motif_beds, stageAs: 'in_bed_*.bed'), path(gene_links_ctrl), path(gene_links_trt)
    val control_label
    val treatment_label

    output:
    tuple val(tf), path("cis_rewiring_motif_stack_${tf}.tsv"), emit: tsv
    tuple val(tf), path("cis_rewiring_motif_stack_${tf}.png"), emit: png

    script:
    def panels = params.cis_rewiring?.panels
    def panel_strs = (panels && !panels.isEmpty())
        ? panels.collect { "--panel '${it}'" }.join(' ')
        : ''
    def panel_file_arg = (params.cis_rewiring?.panel_file)
        ? "--panel-file '${params.cis_rewiring.panel_file}'"
        : ''
    def use_top_n = !panel_strs && !panel_file_arg
    def top_n_arg = use_top_n
        ? "--top-n-by-delta ${(params.cis_rewiring?.top_n_genes ?: 25) as Integer}"
        : ''
    def min_delta = (params.cis_rewiring?.min_delta ?: 1) as Integer
    def excl_pg = (params.cis_rewiring?.exclude_pseudogenes ?: true) ? '--exclude-pseudogenes' : ''
    """
    set -eo pipefail
    # Per-CT motif BEDs all describe the same union peakset; concat + dedup
    # to a single TF-wide motif BED that the renderer expects.
    cat in_bed_*.bed \\
        | awk 'BEGIN{OFS="\\t"} \$1!~/^#/ && NF>=3 {print \$1,\$2,\$3}' \\
        | sort -k1,1 -k2,2n -u \\
        > tf_motif_${tf}.bed

    python ${projectDir}/bin/render_cis_rewiring_motif_stack.py \\
        --links-control   '${gene_links_ctrl}' \\
        --links-treatment '${gene_links_trt}' \\
        --motif-bed       'tf_motif_${tf}.bed' \\
        --motif-name      '${tf}' \\
        --control-label   '${control_label}' \\
        --treatment-label '${treatment_label}' \\
        --min-delta       ${min_delta} \\
        ${excl_pg} \\
        ${top_n_arg} \\
        ${panel_strs} \\
        ${panel_file_arg} \\
        --out-tsv 'cis_rewiring_motif_stack_${tf}.tsv' \\
        --out-png 'cis_rewiring_motif_stack_${tf}.png'
    """

    stub:
    """
    touch "cis_rewiring_motif_stack_${tf}.tsv" "cis_rewiring_motif_stack_${tf}.png"
    """
}
