// ============================================================================
// hdWGCNA Differential Module Eigengene (DME) Analysis
// FIX-29: Condition-comparative analysis — FindDMEs, module-trait correlation,
//         enrichment analysis per cell type
// ============================================================================

process HDWGCNA_DIFFERENTIAL {
    tag "${cell_type}"
    label 'process_medium'
    container params.containers.r_seurat

    // FIX: sanitize cell_type for filesystem paths
    publishDir "${params.outdir}/hdwgcna/differential/${cell_type_safe}", mode: 'copy'

    input:
    path seurat_rds
    val  cell_type
    val  cell_type_key
    val  condition_key
    val  control_condition
    val  treatment_condition
    val  traits              // comma-separated list of trait columns for module-trait correlation
    path group_mapping       // sample -> condition_group JSON (or 'NO_FILE'); recovers condition_key when absent from the object

    output:
    path "dme_results_*.csv",         emit: dme_results,         optional: true
    path "dme_plots/",                emit: dme_plots,           optional: true
    path "module_trait_*.pdf",        emit: module_trait_plots,   optional: true
    path "enrichment_results/",       emit: enrichment_results,  optional: true
    path "hdwgcna_diff_${cell_type_safe}.log", emit: log

    script:
    cell_type_safe = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    def traits_arg = traits ? "--traits ${traits}" : ""
    def mapping_arg = group_mapping.name != 'NO_FILE' ? "--group_mapping ${group_mapping}" : ""
    """
    run_hdwgcna_differential.R \\
        --seurat_rds ${seurat_rds} \\
        --cell_type "${cell_type}" \\
        --cell_type_key ${cell_type_key} \\
        --condition_key ${condition_key} \\
        --control "${control_condition}" \\
        --treatment "${treatment_condition}" \\
        ${traits_arg} \\
        ${mapping_arg} \\
        --output_prefix hdwgcna_diff 2>&1 | tee "hdwgcna_diff_${cell_type_safe}.log"
    """
}

process HDWGCNA_ENRICHMENT {
    tag "${cell_type}"
    label 'process_medium'
    container params.containers.r_seurat

    // FIX: sanitize cell_type for filesystem paths
    publishDir "${params.outdir}/hdwgcna/enrichment/${cell_type_safe}", mode: 'copy'

    input:
    path seurat_rds
    val  cell_type
    val  species

    output:
    path "enrichr_plots/",            emit: enrichr_plots,      optional: true
    path "enrichr_results.csv",       emit: enrichr_table,      optional: true
    path "network_plots/",            emit: network_plots,      optional: true
    path "module_umap.pdf",           emit: module_umap,        optional: true
    path "enrichment_${cell_type_safe}.log", emit: log

    script:
    cell_type_safe = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    """
    run_hdwgcna_enrichment.R \\
        --seurat_rds ${seurat_rds} \\
        --cell_type "${cell_type}" \\
        --species ${species} \\
        --output_prefix enrichment 2>&1 | tee "enrichment_${cell_type_safe}.log"
    """
}
