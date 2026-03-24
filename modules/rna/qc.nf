process RNA_QC {
    tag "$sample"
    publishDir "${params.outdir}/rna/qc", mode: 'copy'
    
    input:
    tuple val(sample), path(mex_zip)
    path metadata_file
    path july_souporcell_dir    
    path nov_souporcell_dir 
    
    output:
    tuple val(sample), path("${sample}_filtered_*.h5ad"), emit: filtered_h5ad  // Changed to glob pattern
    path "${sample}_qc_metrics.csv", emit: metrics
    path "*.png", emit: plots, optional: true
    
    script:
    def is_june = sample.contains('june')
    def is_july = sample.contains('july')
    def is_nov = sample.contains('nov') 

    def metadata_arg = ""
    def souporcell_arg = ""
    
    if (is_june && metadata_file.name != 'NO_FILE') {
        metadata_arg = "--metadata_file ${metadata_file}"
    } else if (is_july && july_souporcell_dir.name != 'NO_FILE') {
        // Extract lane number from sample name
        def lane_match = sample =~ /L(\d+)_july/
        if (lane_match) {
            def lane_num = lane_match[0][1]
            souporcell_arg = "--souporcell_dir ${july_souporcell_dir}/L${lane_num}_Bioproduct"
        }
    } else if (is_nov && nov_souporcell_dir.name != 'NO_FILE') {
        // Extract lane number from November sample name
        def lane_match = sample =~ /L(\d+)_nov/
        if (lane_match) {
            def lane_num = lane_match[0][1]
            // Check if it's ATAC or Bioproduct directory
            def bioproduct_dir = "${nov_souporcell_dir}/L${lane_num}_Bioproduct"
            def atac_dir = "${nov_souporcell_dir}/L${lane_num}_ATAC"
            // Use Bioproduct if it exists, otherwise ATAC
            souporcell_arg = "--souporcell_dir ${bioproduct_dir}"
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
