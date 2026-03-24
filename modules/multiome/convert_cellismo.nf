process CONVERT_CELLISMO {
    tag "$cellismo_file.baseName"
    publishDir "${params.outdir}/cellismo_converted", mode: 'copy', enabled: params.save_intermediate
    
    input:
    path cellismo_file
    
    output:
    path "*.h5mu", emit: h5mu
    
    when:
    params.cellismo.detect
    
    script:
    def output_name = cellismo_file.baseName.replaceAll('\\.cellismo$', '') + '.h5mu'
    """
    #!/usr/bin/env python3
    import shutil
    from pathlib import Path
    
    input_path = Path("${cellismo_file}")
    output_path = Path("${output_name}")
    
    # Cellismo files are just h5mu with different extension
    # Validate structure and convert
    try:
        import mudata as md
        mdata = md.read_h5mu(input_path)
        print(f"✓ Valid cellismo/h5mu file: {mdata.n_obs} cells, {list(mdata.mod.keys())} modalities")
        mdata.write_h5mu(output_path)
        print(f"✓ Converted to: {output_path}")
    except Exception as e:
        print(f"✗ Error: {e}")
        exit(1)
    """
}

process CONCAT_CELLISMO {
    label 'hugemem'
    publishDir "${params.outdir}/cellismo_combined", mode: 'copy'
    
    input:
    path h5mu_files  // List of converted h5mu files
    
    output:
    path "combined_cellismo.h5mu", emit: combined_mudata
    path "integration_report.json", emit: report
    
    when:
    params.cellismo.detect && h5mu_files.size() > 1
    
    script:
    """
    python ${projectDir}/bin/concat_h5mu.py \\
        --input_files ${h5mu_files.join(' ')} \\
        --output_file combined_cellismo.h5mu \\
        --output_dir .
    """
}
