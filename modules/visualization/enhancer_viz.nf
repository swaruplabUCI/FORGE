// modules/visualization/enhancer_viz.nf
//
// Recipe D: Composite Enhancer Visualization
// Process 1: Convert upstream outputs to pyGenomeTracks format
// Process 2: Per-locus composite figure rendering

nextflow.enable.dsl=2

process PREPARE_ENHANCER_VIZ_TRACKS {

    tag "enhancer_viz_tracks"
    label 'process_low'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/enhancer_viz/tracks", mode: 'copy'

    input:
    path cicero_connections    // Cicero connections CSV
    path enhancer_peaks        // CCAN enhancer peaks BED
    path "bigwigs/*"           // Per-cell-type bigWig files (staged into bigwigs/ subdir)
    path gtf                   // Gene annotation GTF
    val  target_genes          // List of target gene names (may be empty)
    path motif_scan_manifest   // region_set_manifest.json from MOTIF_SCAN_ENHANCERS

    output:
    path "track_manifest.json",        emit: track_manifest
    path "cicero_links.bedpe",         emit: links_bedpe
    path "enhancer_annotations.bed",   emit: enhancer_bed
    path "tracks_*.ini",               emit: track_inis

    script:
    def genes_arg = target_genes ? target_genes.join(',') : ''
    """
    python ${projectDir}/bin/prepare_enhancer_viz_tracks.py \\
        --cicero-connections '${cicero_connections}' \\
        --enhancer-peaks '${enhancer_peaks}' \\
        --bigwig-dir bigwigs \\
        --gtf '${gtf}' \\
        --target-genes '${genes_arg}' \\
        --motif-scan-manifest '${motif_scan_manifest}' \\
        --coaccess-threshold ${params.enhancer_viz.coaccess_threshold} \\
        --window-kb ${params.enhancer_viz.window_kb} \\
        --outdir .
    """
}

process COMPOSITE_ENHANCER_VIZ {

    tag "${gene}_${tf_name}"
    label 'process_high'
    maxForks 7
    errorStrategy 'ignore'
    publishDir "${params.outdir}/enhancer_viz/composites/${gene}/${tf_name}", mode: 'copy'

    input:
    path track_manifest        // track_manifest.json from PREPARE_ENHANCER_VIZ_TRACKS
    path "track_inis/*"        // FIX-95: Stage .ini files so pyGenomeTracks can find them
    path motif_scan            // Enhancer motif scan TSV
    val  gene                  // Gene name
    val  tf_name               // TF name
    path footprints_dir        // Directory with pre-computed footprint PNGs

    output:
    path "composite_${gene}_${tf_name}.png",  emit: composite_png, optional: true
    path "composite_${gene}_${tf_name}.pdf",  emit: composite_pdf, optional: true
    path "browser_${gene}.png",               emit: browser_png, optional: true
    path "summary_${gene}_${tf_name}.json",   emit: summary

    script:
    // FIX-97: Removed printer/peak_matrix inputs (caused 4h hang loading 1.9GB h5ad)
    """
    python ${projectDir}/bin/composite_enhancer_viz.py \\
        --track-manifest '${track_manifest}' \\
        --motif-scan '${motif_scan}' \\
        --gene '${gene}' \\
        --tf-name '${tf_name}' \\
        --footprints-dir '${footprints_dir}' \\
        --dpi ${params.enhancer_viz.dpi} \\
        --outdir .
    """
}
