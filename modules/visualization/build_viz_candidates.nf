// modules/visualization/build_viz_candidates.nf
//
// Produces a CSV of (gene, tf, cell_type) triples for COMPOSITE_ENHANCER_VIZ
// to fan out over. Mirrors the Swarup PiD/AD chain:
//   gene-window enhancer ∩ TF motif (in ct)
//   ∧ TF differential FDR < cutoff (in ct)
//   ∧ |log2FC| > cutoff
//   → top-N TFs per (gene, ct) by |log2FC|
//
// Replaces the 17,490-panel cartesian fan-out (genes × TFs) at
// main.nf:958-968. Consumer downstream uses splitCsv to fan out.

nextflow.enable.dsl=2

process BUILD_VIZ_CANDIDATES {

    tag "build_viz_candidates"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/enhancer_viz/filter", mode: 'copy'

    input:
    path track_manifest          // track_manifest.json
    path motif_scan_manifest     // region_set_manifest.json
    path motif_bed_dir           // tf_enhancer_region_sets/ directory from MOTIF_SCAN_ENHANCERS
    path "diff_csvs/*"           // tf_differential_<ct>_<trt>_vs_<ctrl>.csv (collected)

    output:
    path "candidates.csv", emit: candidates
    path "candidates.log", emit: log

    script:
    def f = params.enhancer_viz?.composite_filter ?: [:]
    def top_n      = f.top_n_per_gene_per_ct ?: 3
    def fdr        = f.fdr_cutoff            ?: 0.05
    def lfc        = f.lfc_cutoff            ?: 0.5
    def min_sites  = f.min_motif_sites       ?: 1
    def trt        = params.differential?.treatment_condition ?: ''
    def ctrl       = params.differential?.control_condition   ?: ''
    """
    python ${projectDir}/bin/build_viz_candidates.py \\
        --track-manifest '${track_manifest}' \\
        --motif-scan-manifest '${motif_scan_manifest}' \\
        --motif-bed-dir '${motif_bed_dir}' \\
        --diff-csv-dir diff_csvs \\
        --treatment '${trt}' \\
        --control '${ctrl}' \\
        --top-n-per-gene-per-ct ${top_n} \\
        --fdr-cutoff ${fdr} \\
        --lfc-cutoff ${lfc} \\
        --min-motif-sites ${min_sites} \\
        --out-csv candidates.csv \\
        --debug-log candidates.log
    """
}
