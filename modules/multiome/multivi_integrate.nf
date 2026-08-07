process MULTIVI_INTEGRATE {
    tag "MultiVI_integration"
    label 'process_high'
    publishDir "${params.outdir}/multiome/multivi", mode: 'copy'

    input:
    path mudata_h5mu

    output:
    path "multivi_model", emit: model
    path "multivi_integrated.h5mu", emit: integrated

    script:
    """
    echo "=== MULTIVI_INTEGRATE Debugging ==="
    echo "Working directory: \$(pwd)"
    echo "Input MuData: ${mudata_h5mu}"
    ls -lh

    python ${projectDir}/bin/run_multivi_integration.py \
        --mudata ${mudata_h5mu} \
        --out_model multivi_model \
        --out_mudata multivi_integrated.h5mu \
        --n_epochs ${params.multivi.n_epochs} \
        --batch_key ${params.multivi.batch_key} \
        --modality_weights ${params.multivi.modality_weights} \
        --modality_penalty ${params.multivi.modality_penalty} \
        --latent_dim ${params.multivi.n_latent} \
        --accelerator ${params.scvi_accelerator ?: 'auto'}
    """
}
