// RANK_ENHANCER_STRIP_GENES
// Ranks per-(CT, TF) target genes for MSFP enhancer strip visualization using
// scPrinter binding scores (enhancer_tfbs obsm) × Cicero co-accessibility.
//
// Gating: only fires when msfp_enabled and msfp_strip.enabled. Receives
// binding_scores h5ads that already passed all upstream gates (min cell count,
// ChromVAR z-score floor) — no additional gating needed here. CTs that produce
// zero Cicero-ranked genes are logged and skipped by the downstream RENDER step.
//
// One task processes all (CT, TF) h5ads for all CTs together (shared Cicero
// and GTF indexes built once). Output: per_ct_genes.csv, per_ct_tf_genes.json,
// strip_gene_ranking.json.

process RANK_ENHANCER_STRIP_GENES {
    tag "rank_strip_genes"
    label 'process_medium'

    publishDir "${params.outdir}/enhancer_footprinting_per_ct/strip_gene_ranking",
               mode: 'copy', overwrite: true

    input:
    path tfbs_h5ads
    path cicero_connections
    path gtf
    val  top_n_regions
    val  top_k_genes

    output:
    path "per_ct_genes.csv",      emit: per_ct_genes
    path "per_ct_tf_genes.json",  emit: per_ct_tf_genes
    path "strip_gene_ranking.json", emit: ranking

    script:
    def h5ad_arg = (tfbs_h5ads instanceof List ? tfbs_h5ads : [tfbs_h5ads]).collect { "${it}" }.join(' ')
    """
    python ${projectDir}/bin/rank_enhancer_strip_genes.py \\
        --tfbs-h5ads ${h5ad_arg} \\
        --cicero-connections '${cicero_connections}' \\
        --gtf '${gtf}' \\
        --top-n-regions ${top_n_regions} \\
        --top-k-genes   ${top_k_genes} \\
        --min-coaccess  0.05 \\
        --promoter-window 2000
    """
}
