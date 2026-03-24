#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Run hdWGCNA for a single cell type
 */
process HDWGCNA_PER_CELLTYPE {
    tag "${cell_type}"
    label 'process_hdwgcna'
    
    publishDir "${params.outdir}/hdwgcna/${cell_type_safe}", mode: 'copy'
    
    input:
    tuple path(seurat_rds), val(cell_type)
    val cell_type_key
    val metadata_file
    
    output:
    tuple val("${cell_type}"), path("FiguresTables/**"), emit: figures, optional: true
    tuple val("${cell_type}"), path("*_seurat_with_wgcna.rds"), emit: seurat_obj, optional: true
    tuple val("${cell_type}"), path("*_complete_results.rds"), emit: results, optional: true
    tuple val("${cell_type}"), path("*_analysis_skipped.csv"), emit: skipped, optional: true
    tuple val("${cell_type}"), path("*.log"), emit: log
    
    script:
    // Sanitize cell type name for filenames
    cell_type_safe = cell_type.replaceAll(/[\/\s\(\)]+/, "_")
    def metadata_param = (metadata_file && metadata_file != 'NO_FILE') ? "--metadata ${metadata_file}" : ""
    
    """
    # Create output directory
    mkdir -p FiguresTables
    
    # Run hdWGCNA
    run_hdwgcna_celltype.R \\
        --seurat_rds ${seurat_rds} \\
        --cell_type "${cell_type}" \\
        --cell_type_key ${cell_type_key} \\
        --output_dir . \\
        --threads ${task.cpus} \\
        ${metadata_param} \\
        > "${cell_type_safe}.log" 2>&1 || {
            echo "Cell_Type,Status,Total_Cells" > "${cell_type_safe}_analysis_skipped.csv"
            echo "${cell_type},Failed,0" >> "${cell_type_safe}_analysis_skipped.csv"
            exit 0
        }
    """
}
