// modules/scprint/motif_scan_enhancers_union.nf
//
// MOTIF_SCAN_ENHANCERS_UNION — sibling of MOTIF_SCAN_ENHANCERS that scans the
// (ctrl ∪ trt) enhancer peakset. Identical to motif_scan_enhancers.nf except
// for tag + publishDir (motif_scan_union/) so outputs do not clobber the
// global motif scan. Produces per-(cell_type, TF) BED region sets used by the
// cis-rewiring motif-presence panels.

process MOTIF_SCAN_ENHANCERS_UNION {

    tag "motif_scan_enhancers_union"
    label 'process_high'
    publishDir "${params.outdir}/enhancer_footprinting/motif_scan_union", mode: 'copy'

    input:
    path enhancer_peaks       // ccan_enhancer_peaks_union.bed.gz from BUILD_UNION_ENHANCER_PEAKS
    path motif_list_json      // per_celltype_motifs.json from EXTRACT_CHROMVAR_MOTIFS

    output:
    path "enhancer_motif_scan.tsv.gz",       emit: motif_scan
    path "tf_enhancer_region_sets/",         emit: region_sets
    path "region_set_manifest.json",         emit: manifest
    path "motif_scan_summary.txt",           emit: summary

    script:
    """
    export HDF5_USE_FILE_LOCKING=FALSE

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

    stub:
    """
    mkdir -p tf_enhancer_region_sets
    touch enhancer_motif_scan.tsv.gz region_set_manifest.json motif_scan_summary.txt
    """
}
