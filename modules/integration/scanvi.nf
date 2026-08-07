// modules/integration/scanvi.nf  v2.2.0
// FIX-33c: early_stopping_monitor fix is in bin/train_scanvi.py.

process TRAIN_SCANVI {
    tag "scANVI_${params.species}"
    label 'process_gpu'
    publishDir "${params.outdir}/scanvi", mode: 'copy'
    // Conda handled in config
    
    input:
    path prepared_ref       // Prepared reference h5ad
    path prepared_query     // Prepared query h5ad
    path label_key_file     // Label key text file
    path scvi_model_dir     // Directory containing trained SCVI model
    
    output:
    path "annotated_*.h5ad", emit: annotated
    path "scanvi_*", emit: model_dir, optional: true
    
    script:
    """
    python3 ${baseDir}/bin/train_scanvi.py \\
        --prepared_ref ${prepared_ref} \\
        --prepared_query ${prepared_query} \\
        --label_key_file ${label_key_file} \\
        --scvi_model_dir ${scvi_model_dir} \\
        --results_dir . \\
        --epochs ${params.scanvi_epochs} \\
        --accelerator ${params.scvi_accelerator ?: 'auto'}
    """
}
