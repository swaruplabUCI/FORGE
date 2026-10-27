#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// BETA (2026-10-23): entry point for params.onramp.seurat_rds. Converts a
// user-supplied Seurat object into the integrated-RNA h5ad that main.nf
// injects at ch_integrated_rna.
//
// Unlike CONVERT_H5AD_TO_SEURAT this process needs no containerOptions
// cache pinning: bin/convert_seurat_to_h5ad.R writes HDF5 through rhdf5
// rather than zellkonverter, so there is no Basilisk Miniforge environment
// to download or persist.
// process_medium (64 GB) matches the sibling CONVERT_H5AD_TO_SEURAT and fits
// on a standard node. The conversion streams the sparse slots straight to
// HDF5 rather than densifying, so peak memory is roughly the object itself
// plus one dgCMatrix coercion. An atlas large enough to exceed this wants a
// `withName: CONVERT_SEURAT_TO_H5AD { memory = ... }` override, not a bigger
// default that makes every small run unschedulable.
process CONVERT_SEURAT_TO_H5AD {
    label 'process_medium'
    tag "seurat_rds_onramp"

    publishDir "${params.outdir}/converted", mode: 'copy', enabled: params.save_converted

    input:
    path seurat_rds

    output:
    path "onramp_rna_integrated.h5ad",           emit: h5ad
    path "onramp_rna_integrated_celltypes.txt",  emit: cell_types

    script:
    def assay      = params.onramp_seurat?.assay      ?: 'RNA'
    def sample_key = params.onramp_seurat?.sample_key ?: 'orig.ident'
    // Omit the flag entirely when unset so the converter runs its own
    // detection precedence rather than being handed the string "null".
    def ct_key     = params.onramp_seurat?.cell_type_key
                     ? "--cell_type_key '${params.onramp_seurat.cell_type_key}'" : ''
    // 2026-10-27: hdWGCNA and CellChat now run on on-ramped objects
    // (RNA_NETWORKS), so the columns their differential tiers are configured
    // against must actually exist. Assert that during conversion -- seconds
    // here instead of a failure inside HDWGCNA_DIFFERENTIAL after the whole
    // network build. Only keys the user actually configured are checked.
    def req = [params.hdwgcna?.run  ? params.hdwgcna?.condition_key  : null,
               params.cellchat?.run ? params.cellchat?.condition_key : null]
              .findAll { it && it != 'none' }.unique()
    def req_keys = req ? "--require_keys '${req.join(',')}'" : ''
    """
    export HOME=\$(mktemp -d /tmp/nf_home_XXXXXX)

    convert_seurat_to_h5ad.R \\
        --rds ${seurat_rds} \\
        --output onramp_rna_integrated.h5ad \\
        --assay '${assay}' \\
        --sample_key '${sample_key}' \\
        ${ct_key} \\
        ${req_keys}
    """
}
