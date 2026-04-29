// modules/atac/atac_qc_auto_thresh.nf

process ATAC_QC_AUTO_THRESH {
    label 'process_medium'
    publishDir "${params.outdir}/initial_qc", mode: 'copy'

    input:
    path dataset   // atac_initial_qc.h5ads
    path stats     // sample_statistics.csv

    output:
    path "thresholds.json", emit: thresholds

    script:
    """
    python ${projectDir}/bin/atac_qc_auto_thresh.py \
      --dataset ${dataset} \
      --stats ${stats} \
      --mode ${params.atac.threshold_mode ?: 'gmm'} \
      --out thresholds.json
    """
}
