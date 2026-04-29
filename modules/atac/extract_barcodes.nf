process EXTRACT_BARCODES {
    label 'process_low'
    publishDir "${params.outdir}/barcodes", mode: 'copy'
    
    input:
    path h5ad_files
    
    output:
    path "*_barcodes.txt", emit: barcode_files
    
    script:
    """
    python ${baseDir}/bin/extract_barcodes.py \\
        --h5ad_files ${h5ad_files} \\
        --output_dir .
    """
}
