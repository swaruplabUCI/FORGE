// modules/atac/snapatac_diff.nf

// FIX-E7: Extract cell types from annotated peak matrix when params.differential.cell_types is empty.
// Reinforces wiring from celltypist/scanvi predictions into differential analysis.
process EXTRACT_ATAC_CELL_TYPES {
    label 'process_low'
    publishDir "${params.outdir}/differential", mode: 'copy'

    input:
    path peak_matrix

    output:
    path "atac_cell_types.txt", emit: cell_types

    script:
    def ct_key = params.differential.condition_key ?: 'cell_type'
    // Annotation key depends on mode: marker-based → 'cell_type', CellTypist → 'celltypist_prediction'
    def annotation_key = params.atac.marker_file ? 'cell_type' : 'celltypist_prediction'
    """
    python ${projectDir}/bin/extract_cell_types.py \\
        --h5ad ${peak_matrix} \\
        --cell_type_key ${annotation_key} \\
        --min_cells ${params.differential.min_cells} \\
        --output atac_cell_types.txt \\
        --counts_output atac_cell_type_counts.csv
    """
}

process SNAPATAC_DIFFERENTIAL {
    label 'process_medium'
    publishDir "${params.outdir}/differential", mode: 'copy'
    
    input:
    path peak_matrix
    path metadata
    tuple val(treatment), val(control)
    val cell_types    
    
    output:
    path "DA_peaks_${cell_type}.csv", emit: da_peaks
    path "plots_${cell_type}/*", emit: plots, optional: true
    
    script:
    """
    mkdir -p plots_${cell_type}
    
    python ${projectDir}/bin/run_snapatac_diff.py \\
        --peak_matrix ${peak_matrix} \\
        --metadata ${metadata} \\
        --cell_type "${cell_type}" \\
        --control "${params.differential.control_condition}" \\
        --treatment "${params.differential.treatment_condition}" \\
        --min_cells ${params.differential.min_cells} \\
        --fdr ${params.differential.fdr_threshold} \\
        --output_prefix "${cell_type}"
    """
}
