process MOFA_VISUALIZE {
    label 'hugemem'
    publishDir "${params.outdir}/mofa_visualization", mode: 'copy'
    
    input:
    path mofa_model
    path mudata_file
    path metadata_json
    
    output:
    path "integrated_*_mofa_complete.h5mu", emit: integrated_mudata
    path "figures/*", emit: plots
    path "mofa_*.csv", emit: tables
    path "mofa_integration_summary.json", emit: summary
    
    when:
    params.mofa.run
    
    script:
    """
    mkdir -p figures
    
    python ${projectDir}/bin/mofa_vis.py \\
        --mofa_model ${mofa_model} \\
        --mudata_file ${mudata_file} \\
        --metadata_file ${metadata_json} \\
        --output_dir .
    """
}
