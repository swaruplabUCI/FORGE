// modules/visualization/aggregate_fp_stats.nf
//
// Per-(cell_type, TF) scalar metrics rolled up from ENHANCER_FOOTPRINTING
// h5ads, joined with diff CSVs / tf_gene_adjacency / motif counts. Surfaces
// fp_depth_*, bind_dip_depth_*, n_sites_high_quality, chromvar_log2fc/p_adj,
// n_genes_in_network, score_chromvar/footprint/combined, evidence_tier.
//
// Two outputs:
//   fp_stats_per_pair.csv   — one row per (ct, TF), 657 rows
//   fp_stats_per_triple.csv — joined to gene windows, thousands of rows
//
// Used as a shortlist for users to find interesting MSFPs outside the
// top-3 per-(gene, ct) viz selection from BUILD_VIZ_CANDIDATES.

nextflow.enable.dsl=2

process AGGREGATE_FP_STATS {

    tag "aggregate_fp_stats"
    label 'process_medium'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/enhancer_footprinting/summary", mode: 'copy'

    input:
    path "fp_h5ads/*"            // enhancer_footprints_*.h5ad (collected)
    path "tfbs_h5ads/*"          // enhancer_tfbs_*.h5ad (collected)
    path "diff_csvs/*"           // tf_differential_<ct>_<trt>_vs_<ctrl>.csv (collected)
    path tf_gene_adjacency       // tf_gene_adjacency.tsv from BUILD_TF_GENE_NETWORK
    path motif_bed_dir           // tf_enhancer_region_sets/ from MOTIF_SCAN_ENHANCERS
    path region_set_manifest     // region_set_manifest.json
    path track_manifest          // track_manifest.json (for triple-form gene joining)
    // Tier-1 promoter evidence: per-(ct, gene) scPRINTER outputs from
    // SCPRINTER_FOOTPRINTING. Pass the published root directory as a single
    // symlink (same trick as COMPOSITE staging fix) — 39k h5ads in there.
    path promoter_fp_dir

    output:
    path "fp_stats_per_pair.csv",   emit: pair_csv
    path "fp_stats_per_triple.csv", emit: triple_csv, optional: true
    path "fp_stats_README.md",      emit: readme
    path "fp_stats.log",            emit: log

    script:
    def trt        = params.differential?.treatment_condition ?: ''
    def ctrl       = params.differential?.control_condition   ?: ''
    def dip_thr    = params.enhancer_viz?.fp_stats?.dip_threshold ?: 0.05
    """
    python ${projectDir}/bin/aggregate_fp_stats.py \\
        --fp-h5ads fp_h5ads/*.h5ad \\
        --tfbs-h5ads tfbs_h5ads/*.h5ad \\
        --diff-csv-dir diff_csvs \\
        --treatment '${trt}' \\
        --control '${ctrl}' \\
        --tf-gene-adjacency '${tf_gene_adjacency}' \\
        --motif-bed-dir '${motif_bed_dir}' \\
        --region-set-manifest '${region_set_manifest}' \\
        --track-manifest '${track_manifest}' \\
        --promoter-fp-dir '${promoter_fp_dir}' \\
        --dip-threshold ${dip_thr} \\
        --out-pair-csv fp_stats_per_pair.csv \\
        --out-triple-csv fp_stats_per_triple.csv \\
        --out-readme fp_stats_README.md \\
        --debug-log fp_stats.log
    """
}
