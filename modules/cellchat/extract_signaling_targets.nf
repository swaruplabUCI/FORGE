// modules/cellchat/extract_signaling_targets.nf
//
// EXTRACT_SIGNALING_TARGETS — Recipe Step C4
// For validated (receiver, TF) pairs, look up enhancer targets with
// 3-tier fallback: eRegulon → DORC-motif intersection → CCAN-motif fallback.

process EXTRACT_SIGNALING_TARGETS {

    tag "extract_signaling_targets"
    label 'process_low'
    publishDir "${params.outdir}/cellchat/signaling_targets", mode: 'copy'

    input:
    path validated_pairs      // validated_receiver_tf_pairs.tsv
    path eregulon_regions     // eregulon_region_sets/ directory (may be NO_FILE)
    path eregulon_manifest    // eregulon_manifest.json (may be NO_FILE)
    path ccan_region_sets     // tf_enhancer_region_sets/ directory (may be NO_FILE)
    path ccan_manifest        // region_set_manifest.json (may be NO_FILE)

    output:
    path "signaling_footprint_targets/",    emit: region_sets
    path "signaling_targets_manifest.json", emit: manifest
    path "signaling_metadata.tsv",          emit: metadata

    script:
    def ereg_arg = !eregulon_regions.name.startsWith('NO_FILE') ?
        "--eregulon-regions '${eregulon_regions}' --eregulon-manifest '${eregulon_manifest}'" : ""
    def ccan_arg = !ccan_region_sets.name.startsWith('NO_FILE') ?
        "--ccan-regions '${ccan_region_sets}' --ccan-manifest '${ccan_manifest}'" : ""

    """
    python ${projectDir}/bin/extract_signaling_targets.py \\
        --validated-pairs '${validated_pairs}' \\
        ${ereg_arg} \\
        ${ccan_arg} \\
        --min-regions 10 \\
        --outdir .
    """
}
