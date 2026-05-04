// modules/visualization/da_log2fc_heatmaps.nf
// Per-celltype log2FC heatmaps split into Promoter and Distal blocks (Shi Fig 2E).

nextflow.enable.dsl=2

process DA_LOG2FC_HEATMAPS {

    tag "da_log2fc_heatmaps"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent/2E_log2fc_heatmaps", mode: 'copy'

    input:
    path "diff/*"           // DA_peaks_<ct>__<trt>_vs_<ctrl>.csv
    path peak_annotation    // peaks_annotated.tsv
    path gene_links         // ccan_enhancer_gene_links.tsv

    output:
    path "*_log2fc.{pdf,png}",            emit: heatmaps
    path "log2fc_heatmaps_summary.tsv",   emit: summary

    script:
    def trt = params.shi_figures?.treatment ?: params.differential?.treatment_condition
    def ctrl = params.shi_figures?.control ?: params.differential?.control_condition
    def fdr = params.shi_figures?.fdr_cutoff ?: 0.05
    def lfc = params.shi_figures?.lfc_cutoff ?: 0.5
    def top_genes = params.shi_figures?.heatmap_top_genes ?: 8
    def max_rows = params.shi_figures?.heatmap_max_rows ?: 80
    """
    python ${projectDir}/bin/plot_da_log2fc_heatmaps.py \\
        --diff-dir diff \\
        --peak-annotation '${peak_annotation}' \\
        --gene-links '${gene_links}' \\
        --treatment '${trt}' \\
        --control '${ctrl}' \\
        --fdr-cutoff ${fdr} \\
        --lfc-cutoff ${lfc} \\
        --top-genes-per-block ${top_genes} \\
        --max-rows-per-block ${max_rows} \\
        --outdir .
    """
}
