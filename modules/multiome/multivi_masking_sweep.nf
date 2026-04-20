process MULTIVI_MASKING_SWEEP_ONE {
    tag "MultiVI_masking_sweep_frac${fraction}_seed${seed}"
    label 'process_high_memory'

    input:
    tuple val(fraction), val(seed), path(mudata_h5mu)

    output:
    path "fit_*_rna.csv",          emit: rna
    path "fit_*_atac.csv",         emit: atac
    path "fit_*_integration.json", emit: integration

    script:
    """
    python ${projectDir}/bin/multivi_masking_sweep_one.py \
        --mudata ${mudata_h5mu} \
        --output_dir . \
        --fraction ${fraction} \
        --seed ${seed} \
        --n_epochs ${params.multivi.n_epochs} \
        --batch_key ${params.multivi.batch_key} \
        --cell_type_key ${params.multivi.cell_type_key} \
        --n_latent ${params.multivi.n_latent} \
        --modality_weights ${params.multivi.modality_weights} \
        --modality_penalty ${params.multivi.modality_penalty}
    """
}

process MULTIVI_MASKING_SWEEP_AGGREGATE {
    tag "MultiVI_masking_sweep_aggregate"
    label 'process_high_memory'
    publishDir "${params.outdir}/multiome/multivi/masking_sweep", mode: 'copy'

    input:
    path fit_outputs, stageAs: 'fits/*'

    output:
    path "masking_sweep/masking_sweep_results.csv",       emit: results
    path "masking_sweep/masking_sweep_summary.json",      emit: summary
    path "masking_sweep/masking_degradation_curves.pdf",  emit: curves
    path "masking_sweep/*.pdf",                           emit: plots

    script:
    """
    python ${projectDir}/bin/multivi_masking_sweep_aggregate.py \
        --inputs fits \
        --output_dir masking_sweep
    """
}
