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
    set +e
    run_hdwgcna_celltype.R \\
        --seurat_rds ${seurat_rds} \\
        --cell_type "${cell_type}" \\
        --cell_type_key ${cell_type_key} \\
        --output_dir . \\
        --threads ${task.cpus} \\
        ${metadata_param} \\
        > "${cell_type_safe}.log" 2>&1
    rc=\$?
    set -e

    if [ \$rc -ne 0 ]; then
        # Capture the real error tail so the skipped CSV reflects what actually went wrong.
        # Sentinel skipped.csv files written *by the R script itself* (insufficient cells,
        # insufficient metacells) include a "Status" column already and use rc=0.
        last_err=\$(grep -i -m1 -E 'error|Execution halted' "${cell_type_safe}.log" | head -c 500 | tr '"|,' "'__")
        if [ -z "\$last_err" ]; then
            last_err=\$(tail -3 "${cell_type_safe}.log" | tr '\\n' ' ' | head -c 500 | tr '"|,' "'__")
        fi
        echo 'Cell_Type,Status,Reason,Exit_Code' > "${cell_type_safe}_analysis_skipped.csv"
        echo "\\"${cell_type}\\",Failed,\\"\$last_err\\",\$rc" >> "${cell_type_safe}_analysis_skipped.csv"
        # Per-CT failures are non-fatal: continue the pipeline so other cell types proceed.
        exit 0
    fi
    """
}
