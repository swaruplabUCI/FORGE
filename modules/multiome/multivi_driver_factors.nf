process MULTIVI_DRIVER_FACTORS {
    tag "MultiVI_driver_factors"
    label 'process_high_memory'
    publishDir "${params.outdir}/multiome/multivi/driver_factors", mode: 'copy'

    input:
    path mudata_integrated   // MultiVI-integrated MuData
    path multivi_model       // Trained MultiVI model directory
    path mofa_mudata         // MOFA-integrated MuData (optional, for comparison)

    output:
    path "driver_factors/driver_factor_weights_rna.csv", emit: rna_weights
    path "driver_factors/driver_factor_weights_atac.csv", emit: atac_weights
    path "driver_factors/top_driver_genes_per_dim.csv", emit: top_features
    path "driver_factors/modality_specificity.csv", emit: specificity
    path "driver_factors/driver_factors_summary.json", emit: summary
    path "driver_factors/*.pdf", emit: plots

    script:
    def mofa_arg = mofa_mudata.name != 'NO_MOFA' ? "--mofa_mudata ${mofa_mudata}" : ""
    """
    python ${projectDir}/bin/multivi_driver_factors.py \
        --mudata ${mudata_integrated} \
        --model_dir ${multivi_model} \
        --output_dir driver_factors \
        --cell_type_key ${params.multivi.cell_type_key} \
        --method ${params.multivi.driver_factors.method} \
        --n_top_features ${params.multivi.driver_factors.n_top_features} \
        ${mofa_arg}
    """
}
