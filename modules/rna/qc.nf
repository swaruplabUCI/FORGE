process RNA_QC {
    label 'process_medium'
    tag "$sample"
    publishDir "${params.outdir}/rna/qc", mode: 'copy'
    
    input:
    tuple val(sample), path(mex_zip), path(demux_metadata_file), path(demux_souporcell_dir)

    output:
    tuple val(sample), path("${sample}_filtered_*.h5ad"), emit: filtered_h5ad  // Changed to glob pattern
    path "${sample}_qc_metrics.csv", emit: metrics
    path "*.png", emit: plots, optional: true

    script:
    // Generic demux argument wiring:
    //   - demux_metadata_file: CSV with barcode + label columns (e.g., June BD batches)
    //   - demux_souporcell_dir: Souporcell output directory (e.g., July/Nov BD batches)
    // Both resolve to NO_FILE for datasets without demultiplexing (e.g., 10x PBMC).
    def metadata_arg = ""
    def souporcell_arg = ""

    if (demux_metadata_file.name != 'NO_FILE') {
        metadata_arg = "--metadata_file ${demux_metadata_file}"
    } else if (demux_souporcell_dir.name != 'NO_FILE') {
        // Extract lane number from sample name (pattern: L<N>_<batch>)
        def lane_match = sample =~ /L(\d+)_/
        if (lane_match) {
            def lane_num = lane_match[0][1]
            souporcell_arg = "--souporcell_dir ${demux_souporcell_dir}/L${lane_num}_Bioproduct"
        }
    }
        
    """
    python ${baseDir}/bin/rna_qc.py \
        --input ${mex_zip} \
        --sample ${sample} \
        --min_genes ${params.min_genes} \
        --min_cells ${params.min_cells} \
        --mt_threshold ${params.mt_threshold} \
        --output ${sample}_filtered.h5ad \
        --metrics ${sample}_qc_metrics.csv \
        --plots_dir . \
        --species ${params.species} \
        ${metadata_arg} \
        ${souporcell_arg}
    """
}
