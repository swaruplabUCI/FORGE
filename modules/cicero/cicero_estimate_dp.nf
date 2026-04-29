// modules/cicero/cicero_estimate_dp.nf
//
// Cicero Phase 3 — step 1 of 3: build shared cicero CDS + estimate global
// distance_parameter. Replaces the single-threaded run_cicero() bottleneck
// in the old CICERO_FULL. Downstream fan-out is CICERO_FULL_CHROM.
// Ported from src_FORGE 2026-04-23; used verbatim for plain + stratified legs
// via include-as aliasing in main.nf.

process CICERO_ESTIMATE_DP {
    label 'process_medium'
    errorStrategy 'terminate'

    input:
    val triplets_path

    output:
    path "distance_parameter.txt",  emit: dp
    path "cicero_cds_shared.rds",   emit: cicero_cds
    path "gene_annotation.rds",     emit: gene_ann
    path "input_cds_ordered.rds",   emit: ordered_cds
    path "cds_summary.txt",         emit: summary, optional: true
    path "umap_*.pdf",              emit: umap_plots, optional: true

    script:
    """
    export HOME=/tmp/container_home
    mkdir -p /tmp/container_home
    export R_LIBS_USER=""
    cicero_estimate_dp.R \\
      --triplets "${triplets_path}" \\
      --outdir "." \\
      --num_dim ${params.cicero.num_dim} \\
      --sample_num ${params.cicero.sample_num} \\
      --window_bp 500000 \\
      --distance_constraint 250000 \\
      ${params.cicero.use_partition ? "--use_partition" : ""}
    """
}
