// modules/scprint/motif_scan_enhancers.nf
//
// MOTIF_SCAN_ENHANCERS — Recipe Steps A5+A6
// Scan enhancer regions for chromVAR-nominated TF motifs using scPRINTER's
// MOODS-based motif scanner.  Constructs per-(cell_type, TF) BED region sets
// for downstream footprinting.

process MOTIF_SCAN_ENHANCERS {

    tag "motif_scan_enhancers"
    label 'process_high'
    publishDir "${params.outdir}/enhancer_footprinting/motif_scan", mode: 'copy'

    input:
    path enhancer_peaks       // ccan_enhancer_peaks.bed.gz from EXTRACT_CCAN_ENHANCERS
    path motif_list_json      // per_celltype_motifs.json from EXTRACT_CHROMVAR_MOTIFS

    output:
    path "enhancer_motif_scan.tsv.gz",       emit: motif_scan
    path "tf_enhancer_region_sets/",         emit: region_sets
    path "region_set_manifest.json",         emit: manifest
    path "motif_scan_summary.txt",           emit: summary

    script:
    """
    export HDF5_USE_FILE_LOCKING=FALSE

    # Set up scPRINTER cache for genome FASTA access
    _XDG="\${XDG_CACHE_HOME:-/tmp/cache}"
    mkdir -p "\${_XDG}"
    [ ! -e "\${_XDG}/scprinter" ] && ln -s '${params.scprinter.cache_dir}' "\${_XDG}/scprinter" || true

    python ${projectDir}/bin/motif_scan_enhancers.py \\
        --enhancer-peaks '${enhancer_peaks}' \\
        --motif-list '${motif_list_json}' \\
        --pfms '${params.scprinter.pfms}' \\
        --cache-dir '${params.scprinter.cache_dir}' \\
        --genome '${params.scprinter.genome}' \\
        --cpus ${task.cpus} \\
        --min-regions 10 \\
        --outdir .
    """
}
