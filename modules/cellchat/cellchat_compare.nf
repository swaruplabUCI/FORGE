// ============================================================================
// CellChat Comparative Analysis Module
// FIX-29: Per-condition CellChat + mergeCellChat comparison
// ============================================================================

process CELLCHAT_PER_CONDITION {
    tag "${condition}"
    label 'process_medium'
    container params.containers.r_cellchat

    publishDir "${params.outdir}/cellchat/per_condition/${condition}", mode: 'copy'

    input:
    path h5ad
    val  condition
    val  cell_type_key
    val  condition_key
    val  species

    output:
    path "cellchat_${condition}.rds",    emit: cellchat_rds
    path "cellchat_${condition}_plots/", emit: plots
    val  condition,                      emit: condition_label

    script:
    """
    run_cellchat_per_condition.R \
        --h5ad ${h5ad} \
        --condition "${condition}" \
        --cell_type_key ${cell_type_key} \
        --condition_key ${condition_key} \
        --species ${species} \
        --output_rds cellchat_${condition}.rds \
        --output_dir cellchat_${condition}_plots
    """
}

process CELLCHAT_COMPARE {
    tag "compare"
    label 'process_medium'
    container params.containers.r_cellchat

    publishDir "${params.outdir}/cellchat/comparison", mode: 'copy'

    input:
    path cellchat_rds_files   // collected list of per-condition .rds files
    val  condition_labels     // ordered list of condition names
    val  species

    output:
    path "cellchat_comparison.rds",  emit: comparison_rds
    path "comparison_plots/",        emit: plots

    script:
    def labels_str = condition_labels.join(',')
    def rds_str = cellchat_rds_files.collect { it.name }.join(',')
    """
    run_cellchat_compare.R \
        --rds_files "${rds_str}" \
        --conditions "${labels_str}" \
        --species ${species} \
        --output_rds cellchat_comparison.rds \
        --output_dir comparison_plots
    """
}
