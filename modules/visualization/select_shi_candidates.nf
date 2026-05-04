// modules/visualization/select_shi_candidates.nf
// Pick top TFs per cell type and curated gene loci that drive the
// Shi 4C/5B/5D/5F (curated networks) and 4E (locus TF binding) panels.

nextflow.enable.dsl=2

process SELECT_SHI_CANDIDATES {

    tag "select_shi_candidates"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent", mode: 'copy'

    input:
    path "diff_tf/*"          // tf_differential_<ct>_<trt>_vs_<ctrl>.csv (collected)
    path connections_ctrl, stageAs: 'cicero_connections_ctrl.tsv.gz'
    path connections_trt,  stageAs: 'cicero_connections_trt.tsv.gz'
    path gene_links_ctrl
    path gene_links_trt
    path adjacency            // tf_gene_adjacency.tsv (optional but recommended)

    output:
    path "candidates.json", emit: candidates
    path "tf_allowlist.tsv", emit: tf_allowlist

    script:
    def trt = params.shi_figures?.treatment ?: params.differential?.treatment_condition
    def ctrl = params.shi_figures?.control ?: params.differential?.control_condition
    def top_tfs = params.shi_figures?.top_tfs_per_celltype ?: 3
    def top_loci = params.shi_figures?.top_loci ?: 5
    def fdr = params.shi_figures?.fdr_cutoff ?: 0.05
    def seed = params.shi_figures?.seed_genes ?: 'Adam10,App,Apoe,Mapt,Trem2,Cx3cr1,Bin1,Picalm,Spi1,Srebf1,Hes1,Maf'
    """
    awk -F'\\t' 'NR>1 {print \$1}' '${adjacency}' | sort -u > tf_allowlist.tsv

    python ${projectDir}/bin/select_shi_candidates.py \\
        --diff-tf-dir diff_tf \\
        --connections-ctrl '${connections_ctrl}' \\
        --connections-trt  '${connections_trt}' \\
        --gene-links       '${gene_links_ctrl}' \\
        --gene-links-trt   '${gene_links_trt}' \\
        --treatment '${trt}' \\
        --control   '${ctrl}' \\
        --top-tfs ${top_tfs} \\
        --top-loci ${top_loci} \\
        --fdr-cutoff ${fdr} \\
        --seed-genes '${seed}' \\
        --tf-allowlist tf_allowlist.tsv \\
        --adjacency '${adjacency}' \\
        --out candidates.json
    """
}
