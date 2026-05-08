// modules/visualization/motif_in_gained_ccans.nf
//
// MOTIF_IN_GAINED_CCANS — directionality proxy for cis-rewiring.
// For one TF at a time: of the enhancers that gained CCAN links to a gene
// under treatment (gained = treatment-only enhancer→gene tuples), how many
// host the TF's motif? Compares the per-gene gained-with-motif fraction
// against the genome-wide background fraction across the union enhancer set.
// Motif BED comes from MOTIF_SCAN_ENHANCERS_UNION (region_sets/<TF>.bed) so
// treatment-only enhancers are covered by the scan (the methodological fix).

process MOTIF_IN_GAINED_CCANS {

    tag "${tf}"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/cis_rewiring/${tf}", mode: 'copy'

    input:
    tuple val(tf), path(motif_beds, stageAs: 'in_bed_*.bed'), path(gene_links_ctrl), path(gene_links_trt)
    val control_label
    val treatment_label

    output:
    tuple val(tf), path("cis_rewiring_${tf}.tsv"), emit: tsv
    tuple val(tf), path("cis_rewiring_${tf}.png"), emit: png

    script:
    def gene_set_arg  = (params.cis_rewiring?.gene_set)  ? "--gene-set '${params.cis_rewiring.gene_set.join(',')}'" : ''
    def gene_file_arg = (params.cis_rewiring?.gene_file) ? "--gene-file '${params.cis_rewiring.gene_file}'" : ''
    def top_n         = (params.cis_rewiring?.top_n_genes ?: 25) as Integer
    def min_delta     = (params.cis_rewiring?.min_delta ?: 1) as Integer
    """
    set -eo pipefail
    # Per-CT motif BEDs all describe the same union peakset; concat + dedup
    # to a single TF-wide motif BED that motif_in_gained_ccans.py expects.
    cat in_bed_*.bed \\
        | awk 'BEGIN{OFS="\\t"} \$1!~/^#/ && NF>=3 {print \$1,\$2,\$3}' \\
        | sort -k1,1 -k2,2n -u \\
        > tf_motif_${tf}.bed

    python ${projectDir}/bin/motif_in_gained_ccans.py \\
        --links-control   '${gene_links_ctrl}' \\
        --links-treatment '${gene_links_trt}' \\
        --motif-bed       'tf_motif_${tf}.bed' \\
        --control-label   '${control_label}' \\
        --treatment-label '${treatment_label}' \\
        --motif-name      '${tf}' \\
        --top-n           ${top_n} \\
        --min-delta       ${min_delta} \\
        ${gene_set_arg} \\
        ${gene_file_arg} \\
        --out-tsv 'cis_rewiring_${tf}.tsv' \\
        --out-png 'cis_rewiring_${tf}.png'
    """

    stub:
    """
    touch "cis_rewiring_${tf}.tsv" "cis_rewiring_${tf}.png"
    """
}
