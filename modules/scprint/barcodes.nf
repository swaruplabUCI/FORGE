process SCPRINTER_BARCODES {
    label 'process_small'
    publishDir "${params.scprinter.barcode_dir}", mode: 'copy'

    input:
    val qc_anndatas
    val sample_names

    output:
    path "*_barcodes.txt", emit: barcodes

    script:
    def qc_list
    if( qc_anndatas instanceof List ) {
        qc_list = qc_anndatas
    } else {
        qc_list = [ qc_anndatas ]
    }

    def sample_list
    if( sample_names instanceof List ) {
        sample_list = sample_names
    } else {
        sample_list = [ sample_names ]
    }

    def qc_args     = qc_list.collect     { "'${it}'" }.join(' ')
    def sample_args = sample_list.collect { "'${it}'" }.join(' ')

    """
    python ${projectDir}/bin/generate_scprinter_barcodes.py \\
      --qc-anndatas ${qc_args} \\
      --sample-names ${sample_args} \\
      --barcode-dir .
    """
}
