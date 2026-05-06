// modules/visualization/render_promoter_msfp_overlay.nf
//
// Per-(cell_type, gene) PNG overlay: stacked per-condition multi-scale
// footprint heatmaps (clipped to TF-binding band, default <=30 bp) above
// a JASPAR motif-hit track for the chosen TF list.
//
// TFs are ranked upstream from AGGREGATE_FP_STATS.out.triple_csv by
// (n_sites_in_gene_window, bind_dip_depth_mean) and capped at
// params.promoter_overlay.top_n_tfs. Palette is cycled from a shared
// default; per-TF colors can be overridden via --tf-colors.

nextflow.enable.dsl=2

process RENDER_PROMOTER_MSFP_OVERLAY {

    tag "${cell_type}__${gene}"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/promoter_msfp_overlay/png/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}",
               mode: 'copy'

    input:
    tuple val(cell_type), val(gene), path(fp_h5ad), val(tfs_csv), val(tf_colors_csv)

    output:
    path "${gene}__${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}__msfp_motif_overlay.png", emit: png

    script:
    def safe_ct = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    def score_threshold_arg = (params.promoter_overlay?.score_threshold != null) ?
        "--score-threshold ${params.promoter_overlay.score_threshold}" : ''
    """
    python ${projectDir}/bin/render_promoter_msfp_overlay.py \\
        --gene '${gene}' \\
        --cell-type '${cell_type}' \\
        --fp-h5ad '${fp_h5ad}' \\
        --out-png '${gene}__${safe_ct}__msfp_motif_overlay.png' \\
        --pfm '${params.scprinter.pfms}' \\
        --cache-dir '${params.scprinter.cache_dir}' \\
        --genome '${params.scprinter.genome}' \\
        --max-scale ${params.promoter_overlay.max_scale} \\
        --tfs '${tfs_csv}' \\
        --tf-colors '${tf_colors_csv}' \\
        ${score_threshold_arg}
    """

    stub:
    def safe_ct = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    """
    touch "${gene}__${safe_ct}__msfp_motif_overlay.png"
    """
}
