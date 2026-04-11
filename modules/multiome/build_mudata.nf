process BUILD_MUDATA {
    tag "Build_MuData"
    label 'process_high'
    publishDir "${params.outdir}/multiome/mudata", mode: 'copy'

    input:
    path rna_h5ad_files
    path atac_h5ad_files
    path scanvi_file
    path metadata_csv
    path sample_map
    path atac_annotations

    output:
    path "integrated.h5mu", emit: mudata
    path "mudata_stats.json", emit: stats

    script:
    def atac_annot_arg = atac_annotations.name != 'NO_FILE' ? "--atac_annotations ${atac_annotations}" : ""
    """
    echo "=== BUILD_MUDATA ==="
    echo "Working directory: \$(pwd)"
    echo "Sample map:"
    cat ${sample_map}
    echo ""
    echo "Total h5ad files: \$(ls *.h5ad 2>/dev/null | wc -l)"
    echo ""

    python ${projectDir}/bin/build_mudata_batched.py \\
        --sample_map ${sample_map} \\
        --scanvi_file ${scanvi_file.name} \\
        --metadata_file ${metadata_csv.name} \\
        --output_dir . \\
        --batch_size ${params.mudata.batch_size} \\
        ${atac_annot_arg}
    """
}
