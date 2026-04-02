process RUN_CELL_ANNOTATOR {
    label 'process_medium'
    tag "cell_annotator"
    publishDir "${params.outdir}/cell_annotation", mode: 'copy'

    input:
    path rna_h5ad

    output:
    path "annotated_rna.h5ad", emit: annotated_h5ad
    path "cell_type_annotation_report.txt", emit: report

    script:
    """
    export PYTHONPATH="${projectDir}/tools/cell-annotator-main/src:\${PYTHONPATH:-}"
    
    python ${projectDir}/bin/run_cell_annotator.py \\
        --input ${rna_h5ad} \\
        --output annotated_rna.h5ad \\
        --species ${params.species}
    """
}
