process BUILD_MUDATA {
    tag "Build_MuData"
    label 'process_high'
    publishDir "${params.outdir}/multiome/mudata", mode: 'copy'
    
    input:
    path rna_h5ad_files
    path atac_h5ad_files
    path scanvi_file
    path metadata_csv
    
    output:
    path "integrated.h5mu", emit: mudata
    path "mudata_stats.json", emit: stats
    
    script:
    """
    echo "=== BUILD_MUDATA Debugging ==="
    echo "Working directory: \$(pwd)"
    echo "RNA h5ad files:"
    ls -lh *.h5ad | grep -v peak_matrix | grep -v atac_complete | head -20
    echo ""
    echo "Total h5ad files: \$(ls *.h5ad 2>/dev/null | wc -l)"
    echo ""
    
    # FIXED: Distinguish RNA vs ATAC by filename patterns
    # RNA files: contain "_filtered_" 
    # ATAC files: do NOT contain "_filtered_" and are NOT peak_matrix.h5ad
    
    python ${projectDir}/bin/build_mudata_batched.py \\
        --rna_files \$(ls *_filtered_*.h5ad 2>/dev/null | xargs) \\
        --atac_files \$(ls *.h5ad 2>/dev/null | grep -v "_filtered_" | grep -v "peak_matrix" | grep -v "atac_complete" | grep -v "annotated" | xargs) \\
        --scanvi_file ${scanvi_file.name} \\
        --metadata_file ${metadata_csv.name} \\
        --output_dir . \\
        --batch_size ${params.mudata.batch_size}
    """
}
