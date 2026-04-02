// modules/atac/merge_annotations.nf
//
// Merge cell type annotations + sample metadata into the ATAC peak matrix.
// Supports two annotation modes:
//   - 'marker':     cluster-level JSON from marker-based annotation → 'cell_type' column
//   - 'celltypist': per-cell CellTypist h5ad from gene activity scores → 'celltypist_prediction' column

process MERGE_ANNOTATIONS {
    label 'process_medium'
    publishDir "${params.outdir}/consolidated_qc", mode: 'copy'

    input:
    path peak_matrix
    path annotations       // JSON (marker mode) or CellTypist-annotated h5ad (celltypist mode)
    path metadata
    val  annotation_mode   // 'marker' or 'celltypist'

    output:
    path "peak_matrix_annotated.h5ad", emit: peak_matrix

    script:
    if (annotation_mode == 'celltypist') {
        """
        python ${projectDir}/bin/merge_annotations.py \\
            --peak-matrix ${peak_matrix} \\
            --celltypist-h5ad ${annotations} \\
            --metadata ${metadata} \\
            --output peak_matrix_annotated.h5ad
        """
    } else {
        """
        python ${projectDir}/bin/merge_annotations.py \\
            --peak-matrix ${peak_matrix} \\
            --annotations ${annotations} \\
            --metadata ${metadata} \\
            --resolution ${params.atac.annotation_resolution ?: 'leiden_0_5'} \\
            --output peak_matrix_annotated.h5ad
        """
    }
}
