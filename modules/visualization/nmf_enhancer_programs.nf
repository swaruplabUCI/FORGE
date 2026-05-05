// modules/visualization/nmf_enhancer_programs.nf
// Decompose distal enhancer accessibility into NMF programs (Shi Fig 2B).

nextflow.enable.dsl=2

process NMF_ENHANCER_PROGRAMS {

    tag "nmf_enhancer_programs"
    label 'process_high'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent/2B_nmf_enhancer", mode: 'copy'

    input:
    path peak_matrix          // consolidated_qc/peak_matrix_annotated.h5ad
    path enhancer_bed         // ccan_enhancer_peaks.bed.gz
    val  cell_type_col        // broad ATAC cell-type obs column on peak_matrix
                              // (wired from atac_broad_cell_type_key via SHI_FIGURES take:)
                              // TODO (QOL refactor): replace with sidecar emit from
                              // MERGE_ANNOTATIONS so the producer dictates the column name.

    output:
    path "nmf_program_heatmap.{pdf,png}", emit: heatmap
    path "nmf_elbow.{pdf,png}",           emit: elbow
    path "nmf_W.tsv",                     emit: W
    path "nmf_H.tsv",                     emit: H
    path "nmf_summary.json",              emit: summary

    script:
    def kgrid = params.shi_figures?.nmf_k_grid ?: '6,8,10,12'
    def pick_k = params.shi_figures?.nmf_pick_k ?: 0
    def min_cells = params.shi_figures?.nmf_min_cells ?: 200
    def top_cts = params.shi_figures?.nmf_top_celltypes ?: 12
    def cond_key = params.shi_figures?.condition_key ?: 'condition'
    def sample_key = params.shi_figures?.sample_key ?: 'sample'
    """
    python ${projectDir}/bin/nmf_enhancer_programs.py \\
        --peak-matrix '${peak_matrix}' \\
        --enhancers '${enhancer_bed}' \\
        --cell-type-key ${cell_type_col} \\
        --condition-key ${cond_key} \\
        --sample-key ${sample_key} \\
        --k-grid '${kgrid}' \\
        --pick-k ${pick_k} \\
        --min-cells ${min_cells} \\
        --top-celltypes ${top_cts} \\
        --outdir .
    """
}
