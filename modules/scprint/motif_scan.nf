process SCPRINTER_MOTIF_SCAN {
    label 'process_high'
    publishDir "${params.outdir}/scprinter/motifs/${cell_type}", mode: 'copy'

    input:
    path peak_matrix
    path da_peaks
    path printer
    val cell_type
    val footprints_done

    output:
    path "motif_enrichment_*.csv", emit: enrichment
    path "tfbs_*.h5ad", emit: tfbs
    path "motif_plots/*.png", emit: plots, optional: true

    script:
    def annotation_key = params.atac.marker_file ? 'cell_type' : (params.atac.annotation_method == 'scatanno' ? 'cell_type_prediction' : 'celltypist_prediction')
    """
    # FIX-35c: Persistent model cache (see footprinter.nf for details)
    export SCPRINTER_DATA='${params.scprinter.cache_dir}'

    run_scprinter_motif_scan.py \\
        --peak-matrix '${peak_matrix}' \\
        --da-peaks '${da_peaks}' \\
        --cell-type '${cell_type}' \\
        --cell-type-key ${annotation_key} \\
        --cache-dir '${params.scprinter.cache_dir}' \\
        --genome '${params.scprinter.genome}' \\
        --fdr ${params.scprinter.fdr_threshold} \\
        --printer-path '${printer}' \\
        --control-condition '${params.differential.control_condition}' \\
        --treatment-condition '${params.differential.treatment_condition}' \\
        --cpus ${task.cpus}
    """
}
