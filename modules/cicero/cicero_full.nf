// modules/cicero/cicero_full.nf
//
// Run full Cicero pipeline: connections, CCAN, co-accessibility, pseudotime.

process CICERO_FULL {
    label 'process_high'
    errorStrategy 'ignore'
    publishDir "${params.cicero.outdir}/full_ConnsCoAccPseudoTAcc", mode: 'copy'

    input:
    val triplets_path
    val gtf_path
    val sample_num

    output:
    path "cicero_connections.tsv.gz", emit: connections
    path "CCAN_assignments.tsv.gz",  emit: ccan
    path "input_cds_ordered.rds",    emit: cds
    path "*",                        emit: all_files

script:
    """
    export HOME=/tmp/container_home
    mkdir -p /tmp/container_home
    export R_LIBS_USER=""
    run_cicero_full.R \\
      --triplets "${triplets_path}" \\
      --outdir "." \\
      --gtf "${gtf_path}" \\
      --num_dim ${params.cicero.num_dim} \\
      --sample_num ${sample_num} \\
      --connections_cutoff ${params.cicero.connections_cutoff} \\
      --ccan_min_coaccess ${params.cicero.ccan_min_coaccess} \\
      ${params.cicero.use_partition ? "--use_partition" : ""}
    """
}
