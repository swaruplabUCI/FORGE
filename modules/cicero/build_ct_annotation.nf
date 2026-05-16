// modules/cicero/build_ct_annotation.nf
//
// Generates ct_annotation_v2.csv from peak_matrix_annotated.h5ad.
// Wraps bin/build_ct_annotation_v2.py. Runs once per pipeline execution;
// output is broadcast to all CICERO_TRIPLETS_PER_CT strata.

process BUILD_CT_ANNOTATION {
    label 'process_low'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/cicero/per_ct", mode: 'copy'

    input:
    path  peak_matrix     // peak_matrix_annotated.h5ad
    val   cell_type_col   // obs column → CT labels (e.g. 'cell_type_broad')
    val   condition_col   // obs column → condition labels (e.g. 'condition')
    val   min_cells       // global CT abundance floor (absolute)
    val   min_pct         // global CT abundance floor (fraction of total)

    output:
    path "ct_annotation_v2.csv", emit: csv

    script:
    """
    build_ct_annotation_v2.py \\
        --input         "${peak_matrix}" \\
        --output        "ct_annotation_v2.csv" \\
        --cell-type-col "${cell_type_col}" \\
        --condition-col "${condition_col}" \\
        --min-cells     ${min_cells} \\
        --min-pct       ${min_pct}
    """
}
