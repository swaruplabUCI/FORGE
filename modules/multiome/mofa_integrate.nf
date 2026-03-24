process MOFA_INTEGRATE {
    label 'process_high'
    publishDir "${params.outdir}/multiome/mofa", mode: 'copy'
    
    input:
    path mudata
    
    output:
    path "mofa_integrated.h5mu", emit: integrated
    path "mofa_factors.csv", emit: factors
    path "model.hdf5", emit: model
    path "mofa_stats.json", emit: stats
    path "mofa_stats.json", emit: metadata 
    
    script:
    if (params.mofa.mode == 'high_memory') {
        """
        python ${projectDir}/bin/mofa_int.py \\
            --mudata ${mudata} \\
            --output_dir . \\
            --n_factors ${params.mofa.n_factors} \\
            --convergence_mode ${params.mofa.convergence_mode}
        """
    } else if (params.mofa.mode == 'bootstrap') {
        """
        python ${projectDir}/bin/bootstrap_mofa_integration.py \\
            --mudata_file ${mudata} \\
            --output_dir . \\
            --n_iterations ${params.mofa.bootstrap.n_iterations} \\
            --sample_fraction ${params.mofa.bootstrap.sample_fraction} \\
            --n_factors ${params.mofa.n_factors}
        """
    } else {
        error "Unknown MOFA mode: ${params.mofa.mode}"
    }
}
