process SCPRINTER_MOTIF_SCAN {
    label 'process_high'
    // FIX-43: Sanitize cell_type for filesystem paths (spaces, slashes, parens)
    publishDir "${params.outdir}/scprinter/motifs/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}", mode: 'copy'

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
    // FIX-46: Resolve obs column for cell type (scATAnno: cell_type_prediction)
    def annotation_key = params.atac.marker_file ? 'cell_type' : (params.atac.annotation_method == 'scatanno' ? (params.atac.cell_type_col ?: 'cell_type_prediction') : 'celltypist_prediction')
    """
    # FIX-35c: Persistent model cache (see footprinter.nf for details)
    export SCPRINTER_DATA='${params.scprinter.cache_dir}'

    # FIX-40: Copy printer to local work dir to avoid concurrent write conflicts.
    # Multiple motif scan tasks share the same symlinked printer h5ad;
    # scprinter writes binding scores into the file, causing h5py collisions.
    PRINTER_LOCAL="printer_local_copy.h5ad"
    echo "FIX-40: Copying printer to local work dir..."
    cp -L '${printer}' "\${PRINTER_LOCAL}"
    if [ -d '${printer.baseName}_supp' ]; then
        cp -rL '${printer.baseName}_supp' "printer_local_copy_supp"
    fi

    run_scprinter_motif_scan.py \\
        --peak-matrix '${peak_matrix}' \\
        --da-peaks ${da_peaks} \\
        --cell-type '${cell_type}' \\
        --cell-type-key ${annotation_key} \\
        --cache-dir '${params.scprinter.cache_dir}' \\
        --genome '${params.scprinter.genome}' \\
        --fdr ${params.scprinter.fdr_threshold} \\
        --printer-path "\${PRINTER_LOCAL}" \\
        --control-condition '${params.differential.control_condition}' \\
        --treatment-condition '${params.differential.treatment_condition}' \\
        --cpus ${task.cpus}

    # FIX-40: Clean up local copy to save disk space
    rm -f "\${PRINTER_LOCAL}"
    rm -rf "printer_local_copy_supp"
    """
}
