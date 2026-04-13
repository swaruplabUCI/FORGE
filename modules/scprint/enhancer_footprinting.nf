// modules/scprint/enhancer_footprinting.nf
//
// ENHANCER_FOOTPRINTING — Recipe Steps A9+A10
// Run scPRINTER binding scores + multiscale footprints on enhancer region
// sets.  One task per (cell_type, TF) pair from the region set manifest.

process ENHANCER_FOOTPRINTING {

    tag "${cell_type}_${tf_name}"
    label 'process_high'
    maxForks 7
    // FIX-43: Sanitize cell_type/tf_name — slashes/spaces/parens in names break filesystem paths
    publishDir "${params.outdir}/enhancer_footprinting/footprints/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}/${tf_name.replaceAll(/[\/\s\(\)]+/, '_')}", mode: 'copy'

    input:
    path region_set_bed       // Per-TF BED file from MOTIF_SCAN_ENHANCERS or EXTRACT_EREGULON_REGIONS
    path printer              // scPrinter printer h5ad
    path peak_matrix          // Filtered peak matrix h5ad (for barcode groupings)
    val cell_type             // Cell type label
    val tf_name               // TF name label
    path cicero_connections   // Cicero connections TSV (gzipped)
    val cell_type_col         // FIX-43: obs column name for cell types

    output:
    // FIX-43: Sanitize cell_type/tf_name for filenames — '/' → '_'
    path "enhancer_footprints_*.h5ad",  emit: footprints, optional: true
    path "enhancer_tfbs_*.h5ad",        emit: binding_scores, optional: true
    path "enhancer_plots/*.png",                               emit: plots, optional: true
    path "enhancer_fp_summary.csv",                            emit: summary

    script:
    """
    # v2.9/FIX-82: paired MSFP+TFBS panels, cis-filter, CCAN colorbar
    # Disable HDF5 file locking — multiple tasks read the same printer h5ad
    export HDF5_USE_FILE_LOCKING=FALSE
    # Symlink scPRINTER cache to avoid network fetch on import
    _XDG="\${XDG_CACHE_HOME:-/tmp/cache}"
    mkdir -p "\${_XDG}"
    [ ! -e "\${_XDG}/scprinter" ] && ln -s '${params.scprinter.cache_dir}' "\${_XDG}/scprinter" || true

    # Local copy of printer — source is read-only to prevent corruption,
    # but Rust anndata backend requires read-write access
    PRINTER_LOCAL="printer_local_copy.h5ad"
    cp -L '${printer}' "\${PRINTER_LOCAL}"
    chmod u+w "\${PRINTER_LOCAL}"

    # FIX-75: Local copy of peak matrix — multiple tasks read the same file
    # concurrently on NFS, causing HDF5 locking errors (errno 11)
    PEAK_LOCAL="peak_matrix_local_copy.h5ad"
    cp -L '${peak_matrix}' "\${PEAK_LOCAL}"
    chmod u+w "\${PEAK_LOCAL}"

    python ${projectDir}/bin/run_enhancer_footprinting.py \\
        --region-set '${region_set_bed}' \\
        --printer-path "\${PRINTER_LOCAL}" \\
        --peak-matrix "\${PEAK_LOCAL}" \\
        --cell-type '${cell_type}' \\
        --tf-name '${tf_name}' \\
        --cache-dir '${params.scprinter.cache_dir}' \\
        --genome '${params.scprinter.genome}' \\
        --cpus ${task.cpus} \\
        --pfm-path '${params.scprinter.pfms}' \\
        --gtf '${params.scprinter.gtf_human}' \\
        --cicero-connections '${cicero_connections}' \\
        --cell-type-col '${cell_type_col}'
    """
}
