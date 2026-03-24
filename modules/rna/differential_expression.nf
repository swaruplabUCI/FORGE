process ASSIGN_TEST_GROUPS {
    label 'process_medium'
    publishDir "${params.outdir}/rna_differential/group_assignment", mode: 'copy'
    
    input:
    path h5ad
    path group_mapping_json
    
    output:
    path "annotated_with_groups.h5ad", emit: h5ad
    path "group_assignment_summary.json", emit: summary
    
    script:
    """
    python ${projectDir}/bin/assign_test_groups.py \\
        --input ${h5ad} \\
        --output annotated_with_groups.h5ad \\
        --mapping ${group_mapping_json}
    """
}

process CONVERT_H5AD_FOR_MAST {
    tag 'convert_for_mast'

    input:
    path h5ad_with_groups   // from ASSIGN_TEST_GROUPS.out.h5ad

    output:
    path "seurat_for_mast.rds", emit: seurat_rds

    script:
    """
    # Use the same converter script used for hdWGCNA, but with a different output name
    convert_h5ad_to_seurat.R \\
        --h5ad ${h5ad_with_groups} \\
        --output seurat_for_mast.rds \\
        --cell_type_key condition_group \\
        --assay_name RNA
    """
}

process EXTRACT_CELL_TYPES_FOR_MAST {
    label 'process_low'
    publishDir "${params.outdir}/rna_differential/cell_types", mode: 'copy'

    input:
    // Use the annotated h5ad with groups from ASSIGN_TEST_GROUPS
    path h5ad_with_groups

    output:
    path "cell_types_for_mast.txt", emit: cell_types
    path "cell_type_counts_for_mast.csv", emit: cell_type_counts

    script:
    """
    python ${projectDir}/bin/extract_cell_types.py \\
        --h5ad ${h5ad_with_groups} \\
        --cell_type_key ${params.differential_rna.cell_type_key} \\
        --min_cells ${params.differential_rna.min_cells ?: 100} \\
        --output cell_types_for_mast.txt \\
        --counts_output cell_type_counts_for_mast.csv
    """
}

process RUN_MAST_DE {
    label 'process_high'
    publishDir "${params.outdir}/rna_differential/de_results", mode: 'copy'
    
    input:
    path seurat_rds
    tuple val(cell_type), val(group1), val(group2)
    
    output:
    path "${cell_type}_${group1}_vs_${group2}_DEResults.csv", emit: de_results
    
    script:
    """
    Rscript ${projectDir}/bin/run_mast_de.R \\
        --input ${seurat_rds} \\
        --output_dir . \\
        --cell_type "${cell_type}" \\
        --group1 "${group1}" \\
        --group2 "${group2}" \\
        --condition_key ${params.differential_rna.condition_key} \\
        --cell_type_key ${params.differential_rna.cell_type_key} \\
        --test_use MAST
    """
}

process CREATE_VOLCANO_PLOTS {
    label 'process_low'
    publishDir "${params.outdir}/rna_differential/volcano_plots", mode: 'copy'
    
    input:
    path de_results
    
    output:
    path "*.pdf", emit: plots
    
    script:
    def basename = de_results.baseName
    """
    Rscript ${projectDir}/bin/create_volcano_plots.R \\
        --input ${de_results} \\
        --output ${basename}_volcano.pdf \\
        --species ${params.species} \\
        --pval_cutoff 0.05 \\
        --fc_cutoff 0.25
    """
}

process RUN_GO_ENRICHMENT {
    label 'process_medium'
    publishDir "${params.outdir}/rna_differential/go_enrichment", mode: 'copy'
    
    input:
    path de_results
    
    output:
    path "${de_results.baseName}/*", emit: enrichment_results
    
    script:
    """
    mkdir -p ${de_results.baseName}
    
    Rscript ${projectDir}/bin/go_enrichment_single.R ${de_results}
    
    # Move results to named directory
    mv *.pdf *.csv *_summary.txt ${de_results.baseName}/ 2>/dev/null || true
    """
}
