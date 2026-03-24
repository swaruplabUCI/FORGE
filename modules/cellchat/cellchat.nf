#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process RUN_CELLCHAT {
    publishDir "${params.outdir}/cellchat", mode: 'copy'
    errorStrategy 'ignore'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(adata)
    val cell_type_key
    val condition_key
    val species

    output:
    path "${sample_id}_cellchat.rds",          emit: cellchat_rds
    path "${sample_id}_cellchat_results.csv", emit: csv,     optional: true
    path "${sample_id}_cellchat_plots/*",     emit: plots,   optional: true
    // (optional) add report later if you reintroduce HTML
    // path "${sample_id}_cellchat_report.html", emit: report, optional: true

    script:
    def condition_arg = condition_key != "none" ? "--condition_key ${condition_key}" : ""
    """
    export HOME=\$(mktemp -d /tmp/nf_home_XXXXXX)

    # 1) Inference: generate RDS + CSV
    Rscript ${projectDir}/bin/run_cellchat_infer.R \\
        --input ${adata} \\
        --output_prefix ${sample_id} \\
        --cell_type_key ${cell_type_key} \\
        ${condition_arg} \\
        --species ${species} \\
        --threads ${task.cpus}

    # 2) Plotting: load RDS and generate PDFs only if RDS exists
    Rscript ${projectDir}/bin/plot_cellchat_from_rds.R \\
        --input_rds ${sample_id}_cellchat.rds \\
        --output_prefix ${sample_id}
    """
}

process COMPARE_CELLCHAT {
    tag "comparison_${comparison_name}"
    publishDir "${params.outdir}/cellchat/comparison", mode: 'copy'
    
    input:
    path(cellchat_objects)
    val comparison_name
    
    output:
    path("${comparison_name}_comparison.rds"), emit: merged_object
    path("${comparison_name}_comparison_plots/*"), emit: plots
    path("${comparison_name}_comparison_report.html"), emit: report
    
    script:
    """
    export HOME=\$(mktemp -d /tmp/nf_home_XXXXXX)

    Rscript ${projectDir}/bin/compare_cellchat.R \\
        --input ${cellchat_objects} \\
        --output_prefix ${comparison_name} \\
        --threads ${task.cpus}
    """
}
