// modules/multiome/extract_eregulon_regions.nf
//
// EXTRACT_EREGULON_REGIONS — Recipe Steps B3/C4
// Parse SCENIC+ eRegulon TSVs into per-TF enhancer region sets and target
// gene lists.  Optionally intersect with DORC peaks.

process EXTRACT_EREGULON_REGIONS {

    tag "extract_eregulon_regions"
    label 'process_low'
    publishDir "${params.outdir}/enhancer_footprinting/eregulon_regions", mode: 'copy'

    input:
    path ereg_direct          // eRegulon_direct.tsv from SCENICPLUS_RUN
    path r2g_adj              // region_to_gene_adj.tsv from SCENICPLUS_RUN
    path dorc_sig             // dorc_significant.tsv.gz from SCPRINTER_DORC (optional, may be NO_FILE)

    output:
    path "eregulon_region_sets/",          emit: region_sets
    path "eregulon_target_genes.json",     emit: target_genes
    path "eregulon_manifest.json",         emit: manifest
    path "eregulon_dorc_intersection/",    emit: dorc_intersection, optional: true
    path "eregulon_summary.txt",           emit: summary

    script:
    def dorc_arg = !dorc_sig.name.startsWith('NO_FILE') ? "--dorc-sig '${dorc_sig}'" : ""

    """
    python ${projectDir}/bin/extract_eregulon_regions.py \\
        --ereg-direct '${ereg_direct}' \\
        --r2g '${r2g_adj}' \\
        ${dorc_arg} \\
        --min-regions 10 \\
        --outdir .
    """
}
