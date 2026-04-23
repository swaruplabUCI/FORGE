// modules/scprint/build_tf_gene_network.nf
//
// BUILD_TF_GENE_NETWORK — Phase 3
// Construct global TF-gene adjacency matrix from scPrinter binding scores
// and Cicero co-accessibility (Shi et al. TF_Net equivalent).

process BUILD_TF_GENE_NETWORK {

    tag "tf_gene_network"
    label 'process_medium'
    publishDir "${params.outdir}/enhancer_footprinting/network", mode: 'copy'

    input:
    path tf_targets_json       // tf_target_genes.json from MAP_TF_TO_TARGET_GENES
    path enhancer_gene_links   // ccan_enhancer_gene_links.tsv
    path cicero_connections    // Cicero connections TSV
    path "summaries/fp_summary_?.csv"  // Collected summaries; ? resolves name collisions
    val gtf_path               // GTF annotation

    output:
    path "tf_gene_adjacency.tsv",         emit: adjacency
    path "tf_gene_adjacency_matrix.tsv",  emit: matrix
    path "tf_gene_network.graphml",       emit: graphml, optional: true
    path "network_summary.txt",           emit: summary

    script:
    def threshold = params.enhancer_footprinting.binding_threshold ?: 'otsu'
    """
    python ${projectDir}/bin/build_tf_gene_network.py \\
        --tf-targets '${tf_targets_json}' \\
        --enhancer-gene-links '${enhancer_gene_links}' \\
        --cicero-connections '${cicero_connections}' \\
        --summary-dir summaries \\
        --gtf '${gtf_path}' \\
        --binding-threshold '${threshold}' \\
        --outdir .
    """
}
