// ============================================================================
// TF_CHIP_CORROBORATION — empirical ChIP corroboration of FORGE TF calls
//
// Gap-B closer vs MAESTRO. MAESTRO ranks TFs against CistromeDB ChIP-seq via
// GIGGLE/LISA; FORGE ranks by motif (chromVAR) + footprint (scPrinter) + AUCell
// (SCENIC+). This process adds the *empirical ChIP* layer as an orthogonal
// corroboration input, producing a triangulated (motif x chip [x footprint])
// TF table — a strictly stronger claim than MAESTRO's ChIP-only ranking.
//
// GIGGLE/LOLA-style Fisher's-exact overlap test: for each FORGE-nominated
// (cell_type, TF), tests whether the TF's external ChIP peaks are enriched in
// that CT's characteristic ATAC peaks, over the consensus-peak universe.
//
// Nominations: EXTRACT_CHROMVAR_MOTIFS's per_celltype_motifs.json (+ optional
// SCENIC+ eRegulon TF list). Query peaks: one-vs-rest accessibility from the
// annotated peak matrix (no dependency on differential-peak calls).
//
// Container: inline `scgpu` (has anndata/scipy/pandas) — set inline rather than
// via a configs/*.config withName block, which would bust upstream process
// caches (see feedback_config_cache_stability).
// ============================================================================

process TF_CHIP_CORROBORATION {

    tag "tf_chip_corroboration"
    label 'process_medium'
    container params.containers.scgpu

    publishDir "${params.outdir}/validation/tf_chip_corroboration", mode: 'copy'

    input:
    path motif_json          // per_celltype_motifs.json (from EXTRACT_CHROMVAR_MOTIFS)
    path peak_matrix         // annotated peak matrix h5ad (query peak source)
    val  chip_ref_dir        // directory of per-TF ChIP BEDs (<TF>.bed[.gz])
    path eregulon_json       // optional SCENIC+ eRegulon TFs (or a NO_FILE stub)
    val  cell_type_col       // obs column for cell type

    output:
    path "tf_chip_corroboration.tsv",         emit: corroboration
    path "tf_triangulation.tsv",              emit: triangulation
    path "missing_chip_refs.txt",             emit: missing
    path "tf_chip_corroboration_summary.txt", emit: summary

    when:
    // Gate on both an enable flag AND a real reference dir being configured.
    (params.tf_chip?.run ?: false) && chip_ref_dir && chip_ref_dir != 'none'

    script:
    def ereg = (eregulon_json && eregulon_json.name != 'NO_FILE')
                ? "--eregulon-json '${eregulon_json}'" : ""
    def ctcol = (cell_type_col && cell_type_col != 'none')
                ? "--cell-type-col '${cell_type_col}'" : ""
    """
    tf_chip_corroboration.py \\
        --motif-json '${motif_json}' \\
        --peak-matrix '${peak_matrix}' \\
        --chip-ref-dir '${chip_ref_dir}' \\
        ${ereg} \\
        ${ctcol} \\
        --top-peaks-per-ct ${params.tf_chip?.top_peaks_per_ct ?: 2000} \\
        --min-cells ${params.qc?.cell_type_resolution?.min_cells ?: 50} \\
        --min-pct ${params.qc?.cell_type_resolution?.min_pct ?: 0.01} \\
        --fdr ${params.tf_chip?.fdr ?: 0.05} \\
        --min-odds-ratio ${params.tf_chip?.min_odds_ratio ?: 1.5} \\
        --min-overlap ${params.tf_chip?.min_overlap ?: 3} \\
        --outdir .
    """
}
