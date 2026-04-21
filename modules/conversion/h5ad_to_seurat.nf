#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CONVERT_H5AD_TO_SEURAT {
    label 'process_medium'
    tag "h5ad_conversion"

    publishDir "${params.outdir}/converted", mode: 'copy', enabled: params.save_converted

    // Persist zellkonverter's basilisk Miniforge+anndata env to DFS.
    // Default XDG_CACHE_HOME=/tmp/cache (node-local tmpfs) forces a ~79 MB
    // Miniforge download from GitHub on every task; flaky worker egress caused
    // truncated downloads → exit 1. Pinning to ${launchDir}/cache/basilisk
    // downloads once per instance and reuses forever.
    containerOptions "--env XDG_CACHE_HOME=${launchDir}/cache/basilisk --env BASILISK_EXTERNAL_DIR=${launchDir}/cache/basilisk/R/basilisk"

    input:
    path h5ad
    val cell_type_key

    output:
    path "integrated_seurat.rds", emit: seurat_rds
    path "integrated_seurat_celltypes.txt", emit: cell_types

    script:
    """
    export HOME=\$(mktemp -d /tmp/nf_home_XXXXXX)

    convert_h5ad_to_seurat.R \\
        --h5ad ${h5ad} \\
        --output integrated_seurat.rds \\
        --cell_type_key ${cell_type_key} \\
        --assay_name RNA
    """
}
