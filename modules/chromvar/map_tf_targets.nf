// ============================================================================
// MAP_TF_TO_TARGET_GENES
//
// Maps ChromVAR-enriched TFs to their putative target genes using:
//   1. motif-peak match matrix (varm["motif_match"] from chromvar_scan)
//   2. Cicero CCAN assignments (peak co-accessibility networks)
//   3. GTF annotation (gene promoter windows)
//
// For each (cell_type, TF): finds peaks containing the TF's motif, maps those
// peaks to CCANs, then identifies genes whose promoters share the same CCANs.
// This replaces the previous approach of treating TF names as gene names
// (footprinting Jun at Jun's own promoter) with actual TF→target mapping
// (footprinting Jun's binding sites at AGRN, ATAD3A, etc.).
// ============================================================================

process MAP_TF_TO_TARGET_GENES {

    label 'process_medium'

    publishDir "${params.outdir}/chromvar/tf_target_mapping", mode: 'copy'

    input:
    path chromvar_matrix     // chromvar_matrix.h5ad (with varm["motif_match"])
    path motif_json          // per_celltype_motifs.json
    path ccan_assignments    // CCAN_assignments.tsv.gz
    val  gtf_path            // GTF annotation file path

    output:
    path "tf_target_genes.json",         emit: tf_targets
    path "tf_target_mapping_report.txt", emit: report

    script:
    """
    map_tf_to_target_genes.py \
        --chromvar-matrix '${chromvar_matrix}' \
        --motif-json '${motif_json}' \
        --ccan '${ccan_assignments}' \
        --gtf '${gtf_path}' \
        --max-targets-per-tf 20 \
        --outdir .
    """
}
