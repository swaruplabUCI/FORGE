// modules/atac/snapatac_diff.nf
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
