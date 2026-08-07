// modules/integration/scvi.nf  v2.2.0
//
// FIX-33a: Changed --batch_key from 'source' to 'sample'.
//   'source' was a binary column (reference/query) providing almost no
//   batch correction.  'sample' maps to Donor ID in the reference and
//   sample_id in the query, giving scVI proper per-donor correction.
//   This is consistent with the scArches reference-mapping best practices.

process TRAIN_SCVI {
    label 'process_gpu'
    publishDir "${params.outdir}/rna/scvi", mode: 'copy'

    input:
    path prepared_ref

    output:
    path "scvi_*e", emit: scvi_model_dir
    path "reference_scvi_*.h5ad", emit: ref_with_embedding, optional: true

    script:
    """
    python3 ${baseDir}/bin/train_scvi.py \\
        --prepared_ref ${prepared_ref} \\
        --results_dir . \\
        --epochs ${params.n_epochs_scvi ?: 100} \\
        --batch_key sample \\
        --accelerator ${params.scvi_accelerator ?: 'auto'} \\
        --save_adata
    """
}
