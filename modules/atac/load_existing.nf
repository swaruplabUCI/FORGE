// modules/atac/load_existing.nf
process LOAD_EXISTING_DATA {
    label 'process_low'
    publishDir "${params.outdir}/loaded_data", mode: 'copy'
    
    input:
    path peak_matrix
    path metadata
    
    output:
    path "filtered_peak_matrix.h5ad", emit: peak_matrix
    path "filtered_metadata.csv", emit: metadata
    path "data_summary.txt", emit: summary
    
    script:
    // Convert Nextflow boolean to Python boolean string
    def test_mode_py = params.test_mode ? "True" : "False"
    def test_samples_py = params.test_samples.collect{ "'${it}'" }.join(", ")
    
    """
    #!/usr/bin/env python3
    import anndata as ad
    import pandas as pd
    
    # Load data
    print("Loading peak matrix...")
    adata = ad.read_h5ad("${peak_matrix}")
    print(f"Loaded matrix shape: {adata.shape}")
    
    meta = pd.read_csv("${metadata}")
    print(f"Loaded metadata shape: {meta.shape}")
    
    # Add condition information to adata.obs if not present
    if 'condition' not in adata.obs.columns:
        print("Adding condition information to adata.obs...")
        # Map samples to conditions
        sample_to_condition = {}
        for _, row in meta.iterrows():
            sample_to_condition[row['sample_id']] = row['condition']
        
        # Apply mapping
        adata.obs['condition'] = adata.obs['sample'].map(sample_to_condition)
        print(f"Conditions added: {adata.obs['condition'].value_counts()}")
    
    # Filter for test samples if in test mode
    test_mode = ${test_mode_py}
    if test_mode:
        test_samples = [${test_samples_py}]
        print(f"Test mode: filtering for samples {test_samples}")
        adata = adata[adata.obs['sample'].isin(test_samples)].copy()
        meta = meta[meta['sample_id'].isin(test_samples)]
        print(f"Filtered to {adata.n_obs} cells")
    
    # Save filtered data
    adata.write_h5ad("filtered_peak_matrix.h5ad")
    meta.to_csv("filtered_metadata.csv", index=False)
    
    # Summary
    with open("data_summary.txt", "w") as f:
        f.write(f"Cells: {adata.n_obs}\\n")
        f.write(f"Peaks: {adata.n_vars}\\n")
        f.write(f"Samples: {meta.shape[0]}\\n")
        if 'cell_type' in adata.obs.columns:
            f.write("\\nCell types:\\n")
            for ct, count in adata.obs['cell_type'].value_counts().items():
                f.write(f"  {ct}: {count}\\n")
    """
}
