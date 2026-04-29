process CELLCHAT_PREP_H5AD {
    label 'process_low'
    publishDir "${params.outdir}/cellchat/prep", mode: 'copy'

    input:
    path gene_matrix
    path meta_h5ad

    output:
    path "cellchat_input.global.h5ad", emit: global
    path "cellchat_input.control.h5ad", emit: control
    path "cellchat_input.treatment.h5ad", emit: treatment

    script:
    def annotation_key = params.atac.marker_file ? 'cell_type' : (params.atac.annotation_method == 'scatanno' ? (params.atac.cell_type_col ?: 'cell_type_prediction') : 'celltypist_prediction')
    """
    python ${projectDir}/bin/cellchat_prep_h5ad.py \\
      --gene_matrix ${gene_matrix} \\
      --meta_h5ad ${meta_h5ad} \\
      --cell_type_key ${annotation_key} \\
      --condition_key ${params.cellchat.condition_key ?: 'condition'} \\
      --control '${params.cellchat.control}' \\
      --treatment '${params.cellchat.treatment}' \\
      --out_global cellchat_input.global.h5ad \\
      --out_control cellchat_input.control.h5ad \\
      --out_treatment cellchat_input.treatment.h5ad
    """
}

process CELLCHAT_INFER {
    label 'process_medium'
    publishDir "${params.outdir}/cellchat/infer", mode: 'copy'

    input:
    tuple val(tag), path(h5ad)

    output:
    tuple val(tag), path("${tag}_cellchat.rds"), emit: rds
    tuple val(tag), path("${tag}_cellchat_results.csv"), emit: csv

    script:
    def workers = params.cellchat.max_workers ?: task.cpus
    def annotation_key = params.atac.marker_file ? 'cell_type' : (params.atac.annotation_method == 'scatanno' ? (params.atac.cell_type_col ?: 'cell_type_prediction') : 'celltypist_prediction')
    """
    Rscript ${projectDir}/bin/run_cellchat_infer.R \\
      --input ${h5ad} \\
      --output_prefix ${tag} \\
      --cell_type_key ${annotation_key} \\
      --condition_key none \\
      --species ${params.species} \\
      --threads ${workers}
    """
}

process CELLCHAT_COMPARE {
    label 'process_medium'
    publishDir "${params.outdir}/cellchat/compare", mode: 'copy'

    input:
    path control_rds
    path treatment_rds

    output:
    path "cellchat_comparison.rds", emit: comparison_rds
    path "comparison_plots/",       emit: plots

    script:
    def rds_str    = "${control_rds.name},${treatment_rds.name}"
    def labels_str = "${params.cellchat.control},${params.cellchat.treatment}"
    """
    Rscript ${projectDir}/bin/run_cellchat_compare.R \\
      --rds_files "${rds_str}" \\
      --conditions "${labels_str}" \\
      --species ${params.species} \\
      --output_rds cellchat_comparison.rds \\
      --output_dir comparison_plots
    """
}
