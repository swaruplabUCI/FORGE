process MULTIVI_VALIDATE {
    tag "MultiVI_validation"
    label 'process_high'
    publishDir "${params.outdir}/multiome/multivi/validation", mode: 'copy'

    input:
    path gap_filled_mudata   // Gap-filled MuData from MULTIVI_GAP_FILL
    path rna_h5ad            // Post-integration RNA h5ad (with UMAP)
    path atac_h5ad           // Annotated ATAC h5ad (with UMAP)

    output:
    path "validation/*.pdf", emit: plots
    path "validation/*.csv", emit: metrics
    path "validation/validation_summary.json", emit: summary

    script:
    """
    python ${projectDir}/bin/multivi_validate.py \
        --gap_filled_mudata ${gap_filled_mudata} \
        --rna_h5ad ${rna_h5ad} \
        --atac_h5ad ${atac_h5ad} \
        --output_dir validation \
        --cell_type_key ${params.multivi.cell_type_key} \
        --k_neighbors 50
    """
}
