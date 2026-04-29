// modules/visualization/enhancer_viz.nf
//
// Recipe D (ATAC-only port from FORGE): composite enhancer visualization via
// pyGenomeTracks. SDas-specific wiring sources per-cell-type bigWigs from
// EXPORT_ATAC_BIGWIGS (snap.ex.export_coverage) instead of FORGE's
// PYCISTOPIC_PREPARE. The two Python entry points are unchanged from
// src_FORGE/bin/ — they operate on a generic --bigwig-dir.

nextflow.enable.dsl=2

process PREPARE_ENHANCER_VIZ_TRACKS {

    tag "enhancer_viz_tracks"
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/enhancer_viz/tracks", mode: 'copy'

    input:
    path cicero_connections    // Cicero connections CSV (plain or stratified)
    path enhancer_peaks        // CCAN enhancer peaks BED
    path "bigwigs/*"           // Per-cell-type bigWig files (staged into bigwigs/)
    path gtf                   // Gene annotation GTF
    val  target_genes          // List of target gene names (may be empty)
    path motif_scan_manifest   // region_set_manifest.json from MOTIF_SCAN_ENHANCERS

    output:
    path "track_manifest.json",        emit: track_manifest
    path "cicero_links.bedpe",         emit: links_bedpe
    path "enhancer_annotations.bed",   emit: enhancer_bed
    path "tracks_*.ini",               emit: track_inis
    path "gene_windows/",              emit: gene_windows, optional: true

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
        --coaccess-threshold ${params.enhancer_viz.coaccess_threshold ?: 0.25} \\
        --window-kb ${params.enhancer_viz.window_kb ?: 100} \\
        --outdir .
    """
}

process COMPOSITE_ENHANCER_VIZ {

    // FIX-49: Sanitize gene + tf_name + cell_type for filesystem paths. JASPAR
    // carries dimer aliases like 'Tbx2/3' that have a slash; pass those through
    // raw and publishDir/output globs explode. Sanitize the same way as
    // enhancer_footprinting.nf:13 / scprinter/footprinter.nf:36. Raw names
    // still go to the Python script via the CLI args (composite logic does
    // its own case-insensitive lookups on raw labels).
    //
    // Per-ct split (Step 3): each task renders one (gene, TF, cell_type)
    // panel — slim ini (only that ct's bigwig) + that ct's footprint
    // composite. The legacy multi-ct stack remains available by passing an
    // empty cell_type (composite_enhancer_viz.py falls back).
    tag "${gene.replaceAll(/[\/\s\(\)]+/, '_')}_${tf_name.replaceAll(/[\/\s\(\)]+/, '_')}_${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}"
    label 'process_high'
    maxForks 7
    errorStrategy 'terminate'
    publishDir "${params.outdir}/enhancer_viz/composites/${gene.replaceAll(/[\/\s\(\)]+/, '_')}/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}/${tf_name.replaceAll(/[\/\s\(\)]+/, '_')}", mode: 'copy'

    input:
    path track_manifest        // track_manifest.json from PREPARE_ENHANCER_VIZ_TRACKS
    path "track_inis/*"        // .ini files staged so pyGenomeTracks can find them
    path motif_scan            // enhancer motif scan TSV
    val  gene                  // gene name
    val  tf_name               // TF name
    val  cell_type             // cell_type label (per-ct split)
    // Step 8: pass the published footprint root as a SINGLE directory symlink
    // instead of .collect()'ing 1,971 individual PNGs (1.4 GB) into every work
    // dir. Reduces .command.run staging time from ~14 min to <1 s per panel.
    // The script's find_footprint_pngs() already uses recursive glob.
    path "footprints"          // directory containing {ct}/{tf}/enhancer_global/*.png
    // Tier-1 promoter MSFP panel restoration (post-FIX-97). Same dir-symlink
    // trick — passed as scprinter/footprints publish root containing
    // {ct}/plots/{ct}_{gene}_composite.png. Empty NO_FILE sentinel when the
    // run lacks SCPRINTER_FOOTPRINTING outputs.
    path "scprinter_promoter_dir"

    output:
    path "composite_*.png",  emit: composite_png, optional: true
    path "composite_*.pdf",  emit: composite_pdf, optional: true
    path "browser_*.png",    emit: browser_png,   optional: true
    // Standalone Panel A (vector). Primary consumer-facing browser asset —
    // composite no longer embeds it (was compressed to ~3-6 in tall).
    path "browser_*.pdf",    emit: browser_pdf,   optional: true
    path "tracks_*.ini",     emit: slim_ini,      optional: true
    path "summary_*.json",   emit: summary

    script:
    def safe_gene = gene.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_tf   = tf_name.replaceAll(/[\/\s\(\)]+/, '_')
    def safe_ct   = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    def jaspar = params.scprinter?.pfms ?: ''
    """
    python ${projectDir}/bin/composite_enhancer_viz.py \\
        --track-manifest '${track_manifest}' \\
        --motif-scan '${motif_scan}' \\
        --gene '${gene}' \\
        --tf-name '${tf_name}' \\
        --cell-type '${cell_type}' \\
        --gene-safe '${safe_gene}' \\
        --tf-name-safe '${safe_tf}' \\
        --cell-type-safe '${safe_ct}' \\
        --jaspar-pfms '${jaspar}' \\
        --footprints-dir footprints \\
        --promoter-dir scprinter_promoter_dir \\
        --dpi ${params.enhancer_viz.dpi ?: 200} \\
        --outdir .
    """
}
