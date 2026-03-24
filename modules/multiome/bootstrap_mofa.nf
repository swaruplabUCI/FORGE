// modules/multiome/bootstrap_mofa.nf

process BOOTSTRAP_MOFA_INTEGRATION {
    tag "${mudata.baseName}"
    label 'process_high_memory'
    publishDir "${params.outdir}/bootstrap_mofa", mode: 'copy'
    
    // Container assigned in nextflow.config via withName

    input:
    path mudata
    val n_iterations
    val sample_fraction
    val n_factors
    
    output:
    path "bootstrap_results/", emit: results_dir, type: 'dir'
    path "bootstrap_results/bootstrap_summary.json", emit: summary
    path "memory_log.jsonl", emit: memory_log
    path "bootstrap_results/mofa_model_iter_*.hdf5", emit: models
    path "bootstrap_results/iteration_*/factors.csv", emit: factors, optional: true
    
    script:
    """
    #!/usr/bin/env python3
    
    import sys
    import os
    
    # Add the bin directory to Python path
    sys.path.insert(0, '${projectDir}/bin')
    
    from bootstrap_mofa_integration import main
    
    # Set up arguments
    class Args:
        mudata_file = '${mudata}'
        n_iterations = ${n_iterations}
        sample_fraction = ${sample_fraction}
        n_factors = ${n_factors}
        output_dir = 'bootstrap_results'
        seed = 42
    
    sys.argv = ['bootstrap_mofa_integration.py',
                '--mudata_file', Args.mudata_file,
                '--n_iterations', str(Args.n_iterations),
                '--sample_fraction', str(Args.sample_fraction),
                '--n_factors', str(Args.n_factors),
                '--output_dir', Args.output_dir,
                '--seed', str(Args.seed)]
    
    exit_code = main()
    sys.exit(exit_code)
    """
    
    stub:
    """
    mkdir -p bootstrap_results/iteration_00
    touch bootstrap_results/bootstrap_summary.csv
    touch bootstrap_results/memory_log.jsonl
    touch bootstrap_results/iteration_00/mofa_model.hdf5
    touch bootstrap_results/iteration_00/factors.csv
    """
}

process CONSENSUS_ANALYSIS {
    tag "consensus"
    label 'process_medium'
    publishDir "${params.outdir}/bootstrap_mofa/consensus", mode: 'copy'
    
    // Container assigned in nextflow.config via withName

    input:
    path bootstrap_dir
    
    output:
    path "consensus_results/", emit: results_dir
    path "consensus_results/stability_analysis.png", emit: stability_plot
    path "consensus_results/consensus_summary.json", emit: summary
    
    script:
    """
    #!/usr/bin/env python3
    
    import sys
    sys.path.insert(0, '${projectDir}/bin')
    
    from consensus_analysis import main
    
    sys.argv = ['consensus_analysis.py',
                '--bootstrap_dir', '${bootstrap_dir}',
                '--output_dir', 'consensus_results']
    
    exit_code = main()
    sys.exit(exit_code)
    """
    
    stub:
    """
    mkdir -p consensus_results
    touch consensus_results/stability_analysis.png
    touch consensus_results/consensus_summary.json
    """
}

process ANALYZE_MEMORY_LOG {
    tag "memory_analysis"
    label 'process_low'
    publishDir "${params.outdir}/bootstrap_mofa/memory_analysis", mode: 'copy'
    
    input:
    path memory_log
    
    output:
    path "memory_analysis.txt", emit: report
    path "memory_plot.png", emit: plot, optional: true
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend
    import matplotlib.pyplot as plt
    
    # Load memory log
    records = []
    with open('${memory_log}', 'r') as f:
        for line in f:
            records.append(json.loads(line))
    
    df = pd.DataFrame(records)
    
    # Print columns for debugging
    print(f"Available columns: {df.columns.tolist()}")
    
    # Generate report
    with open('memory_analysis.txt', 'w') as f:
        f.write("=" * 80 + "\\n")
        f.write("MEMORY USAGE ANALYSIS\\n")
        f.write("=" * 80 + "\\n\\n")
        
        # Use correct column names from bootstrap script
        f.write(f"Peak memory usage: {df['process_memory_gb'].max():.2f} GB\\n")
        f.write(f"Average memory usage: {df['process_memory_gb'].mean():.2f} GB\\n")
        
        # Calculate available memory (total - used)
        df['system_available_gb'] = df['system_total_gb'] - df['system_used_gb']
        f.write(f"System memory available: {df['system_available_gb'].min():.2f} GB (minimum)\\n")
        f.write(f"System memory total: {df['system_total_gb'].iloc[0]:.2f} GB\\n")
        
        # Memory at key steps
        f.write("\\nMemory at key steps:\\n")
        for stage in df['stage'].unique():
            stage_mem = df[df['stage'] == stage]['process_memory_gb'].iloc[0]
            f.write(f"  {stage}: {stage_mem:.2f} GB\\n")
    
    # Generate plot
    try:
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(range(len(df)), df['process_memory_gb'], marker='o', label='Process Memory')
        ax.axhline(y=df['system_available_gb'].min(), color='r', linestyle='--', 
                   label='Min Available System Memory')
        ax.set_xlabel('Step Index')
        ax.set_ylabel('Memory (GB)')
        ax.set_title('Memory Usage During Bootstrap MOFA Integration')
        ax.legend()
        plt.tight_layout()
        plt.savefig('memory_plot.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Plot saved successfully")
    except Exception as e:
        print(f"Could not generate plot: {e}")
        import traceback
        traceback.print_exc()
    """
    
    stub:
    """
    echo "Memory analysis stub" > memory_analysis.txt
    """
}
