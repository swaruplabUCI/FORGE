// modules/visualization/curated_tf_networks.nf
// Curated top-3-TF subgraphs per cell type (Shi Fig 4C/5B/5D/5F).

nextflow.enable.dsl=2

process CURATED_TF_NETWORKS {

    tag "curated_tf_networks"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent/4C5B5D5F_curated_networks", mode: 'copy'

    input:
    path adjacency           // tf_gene_adjacency.tsv
    path candidates          // candidates.json
    path "diff_tf/*"         // tf_differential_<ct>_*.csv (for direction)

    output:
    path "*_curated_network.{pdf,png}",     emit: networks
    path "curated_networks_manifest.json",  emit: manifest

    script:
    def trt = params.shi_figures?.treatment ?: params.differential?.treatment_condition
    def ctrl = params.shi_figures?.control ?: params.differential?.control_condition
    def max_targets = params.shi_figures?.network_max_targets ?: 12
    """
    python ${projectDir}/bin/plot_tf_network_curated.py \\
        --adjacency '${adjacency}' \\
        --candidates '${candidates}' \\
        --diff-tf-dir diff_tf \\
        --treatment '${trt}' \\
        --control '${ctrl}' \\
        --max-targets-per-tf ${max_targets} \\
        --outdir .
    """
}
