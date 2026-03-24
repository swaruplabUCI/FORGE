process MULTIVI_VISUALIZE {
    tag "MultiVI_visualization"
    label 'process_medium'
    publishDir "${params.outdir}/multiome/multivi/visualizations", mode: 'copy'

    input:
    path mudata_integrated
    path multivi_model

    output:
    path "multivi_plots/*", emit: plots
    path "multivi_imputation_metrics.csv", emit: metrics
    path "multivi_de_results.csv", emit: de_results, optional: true
    path "multivi_da_results.csv", emit: da_results, optional: true
    path "multivi_visualization_report.html", emit: report

    script:
    """
    python ${projectDir}/bin/multivi_visualize.py \
        --mudata ${mudata_integrated} \
        --model_dir ${multivi_model} \
        --output_dir multivi_plots \
        --cell_type_key ${params.multivi.cell_type_key} \
        --batch_key ${params.multivi.batch_key} \
        --run_imputation ${params.multivi.run_imputation} \
        --run_differential ${params.multivi.run_differential}
    """
}
