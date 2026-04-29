process ATAC_DESCRIPTIVE_REPORT {
    publishDir "${params.outdir}/atac/descriptive", mode: 'copy'

    input:
    path peak_matrix
    path metadata

    output:
    path "summary.json", emit: summary
    path "cells_by_sample.csv", emit: by_sample
    path "cells_by_condition.csv", emit: by_condition
    path "cells_by_cell_type.csv", emit: by_cell_type
    path "cells_by_celltype_and_condition.csv", emit: by_celltype_condition, optional: true
    path "cell_types.txt", emit: cell_types

    script:
    def annotation_key = params.atac.marker_file ? 'cell_type' : (params.atac.annotation_method == 'scatanno' ? (params.atac.cell_type_col ?: 'cell_type_prediction') : 'celltypist_prediction')
    """
    python ${projectDir}/bin/atac_descriptive_report.py \\
      --peak_matrix ${peak_matrix} \\
      --metadata ${metadata} \\
      --condition_key "${params.differential.condition_key}" \\
      --cell_type_key "${annotation_key}" \\
      --sample_key "sample" \\
      --output_dir .
    """
}
