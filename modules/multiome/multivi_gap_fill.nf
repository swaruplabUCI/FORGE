process MULTIVI_GAP_FILL {
    tag "MultiVI_gap_fill"
    label 'process_high_memory'
    publishDir "${params.outdir}/multiome/multivi/gap_fill", mode: 'copy'

    input:
    path mudata_h5mu         // Paired-only integrated MuData (from BUILD_MUDATA)
    path rna_h5ad            // Post-QC RNA h5ad (all cells, pre-intersection)
    path atac_h5ad           // Post-QC ATAC h5ad (all cells, pre-intersection)

    output:
    path "gap_fill/multivi_gap_filled.h5mu", emit: gap_filled
    path "gap_fill/gap_fill_stats.json", emit: stats
    path "gap_fill/multivi_gap_fill_model", emit: model, optional: true
    path "gap_fill/*.pdf", emit: plots

    script:
    """
    python ${projectDir}/bin/multivi_gap_fill.py \
        --mudata ${mudata_h5mu} \
        --rna_h5ad ${rna_h5ad} \
        --atac_h5ad ${atac_h5ad} \
        --output_dir gap_fill \
        --n_epochs ${params.multivi.n_epochs} \
        --batch_key ${params.multivi.batch_key} \
        --n_latent ${params.multivi.n_latent} \
        --modality_weights ${params.multivi.modality_weights} \
        --modality_penalty ${params.multivi.modality_penalty} \
        --n_posterior_samples ${params.multivi.gap_fill.n_posterior_samples} \
        --min_confidence ${params.multivi.gap_fill.min_confidence} \
        --cell_type_key ${params.multivi.cell_type_key}
    """
}
