process MERGE_ANNOTATIONS {
    publishDir "${params.outdir}/consolidated_qc", mode: 'copy'
    
    input:
    path peak_matrix
    path annotations_json
    path metadata  // Add metadata input
    
    output:
    path "peak_matrix_annotated.h5ad", emit: peak_matrix
    
    script:
    """
    python ${projectDir}/bin/merge_annotations.py \
        --peak-matrix ${peak_matrix} \
        --annotations ${annotations_json} \
        --metadata ${metadata} \
        --resolution ${params.atac.annotation_resolution ?: 'leiden_0_5'} \
        --output peak_matrix_annotated.h5ad
    """
}
