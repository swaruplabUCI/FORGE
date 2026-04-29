// modules/cicero/cicero_join.nf
//
// Cicero Phase 3 — step 3 of 3: rbind per-chromosome connections, filter,
// generate CCANs, produce pseudotime/regional plots. Emits the same downstream
// contract as the old CICERO_FULL (cicero_connections.tsv.gz,
// CCAN_assignments.tsv.gz, input_cds_ordered.rds).
// Adapted from src_FORGE 2026-04-23: added `val variant` to parameterize
// publishDir so Shi-stratified aliases land under stratified/${condition}/.

process CICERO_JOIN {
    tag "${variant ?: 'plain'}"
    label 'process_medium'
    errorStrategy 'terminate'
    publishDir "${params.cicero.outdir}/${variant ?: 'full_ConnsCoAccPseudoTAcc'}", mode: 'copy'

    input:
    path chrom_conns  // collected list of conns_*.tsv.gz from CICERO_FULL_CHROM
    path ordered_cds  // input_cds_ordered.rds from CICERO_ESTIMATE_DP
    val  gtf_path
    val  variant      // "" for plain, "stratified/<cond>" for Shi legs

    output:
    path "cicero_connections.tsv.gz", emit: connections
    path "CCAN_assignments.tsv.gz",   emit: ccan
    path "input_cds_ordered.rds",     emit: cds
    path "*", emit: all_files

    script:
    """
    export HOME=/tmp/container_home
    mkdir -p /tmp/container_home
    export R_LIBS_USER=""
    # Ensure the published CDS name matches the downstream contract
    if [ "${ordered_cds}" != "input_cds_ordered.rds" ]; then
        cp -L "${ordered_cds}" "input_cds_ordered.rds"
    fi
    cicero_join.R \\
      --conns_glob "conns_*.tsv.gz" \\
      --cds "input_cds_ordered.rds" \\
      --gtf "${gtf_path}" \\
      --connections_cutoff ${params.cicero.connections_cutoff} \\
      --ccan_min_coaccess ${params.cicero.ccan_min_coaccess} \\
      --outdir "."
    """
}
