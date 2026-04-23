// modules/cicero/cicero_stratified.nf
//
// Disease-stratified Cicero pipeline (Shi et al. methodology).
// Runs Cicero separately per condition to capture condition-specific
// co-accessibility patterns. Mirrors:
//   ADnControl_cicero_parallel.R / PiDnControl_cicero_parallel.R

process CICERO_TRIPLETS_STRATIFIED {
    tag "cicero_triplets_${condition}"
    label 'process_medium'
    publishDir "${params.cicero.outdir}/stratified/${condition}", mode: 'copy'

    input:
    path peak_matrix
    val condition_key
    val condition

    output:
    path "cicero_triplets_${condition}.tsv.gz", emit: triplets
    val condition,                               emit: condition_label

    script:
    """
    make_cicero_triplets_stratified.py \\
        --peak-matrix "${peak_matrix}" \\
        --condition-key "${condition_key}" \\
        --condition-value "${condition}" \\
        --out "cicero_triplets_${condition}.tsv.gz"
    """
}


process CICERO_FULL_STRATIFIED {
    tag "cicero_full_${condition}"
    label 'process_high'
    publishDir "${params.cicero.outdir}/stratified/${condition}", mode: 'copy'

    input:
    val triplets_path
    val gtf_path
    val sample_num
    val condition

    output:
    path "cicero_connections.tsv.gz",  emit: connections
    path "CCAN_assignments.tsv.gz",    emit: ccan
    path "input_cds_ordered.rds",      emit: cds
    val condition,                     emit: condition_label
    path "*",                          emit: all_files

    script:
    """
    run_cicero_full.R \\
        --triplets "${triplets_path}" \\
        --outdir "." \\
        --gtf "${gtf_path}" \\
        --num_dim ${params.cicero.num_dim} \\
        --sample_num ${sample_num} \\
        --connections_cutoff ${params.cicero.connections_cutoff} \\
        --ccan_min_coaccess ${params.cicero.ccan_min_coaccess} \\
        ${params.cicero.use_partition ? "--use_partition" : ""}
    """
}


process COMPARE_COACCESSIBILITY {
    tag "compare_coaccess"
    label 'process_medium'
    publishDir "${params.outdir}/enhancer_footprinting/coaccessibility_comparison", mode: 'copy'

    input:
    path "ctrl/*"              // staged into ctrl/ to avoid name collision
    path "trt/*"               // staged into trt/  to avoid name collision
    val control_label
    val treatment_label

    output:
    path "coaccessibility_correlation.tsv",     emit: correlation
    path "condition_specific_enhancers.bed",     emit: condition_enhancers
    path "coaccessibility_comparison.png",       emit: plot
    path "coaccessibility_summary.txt",          emit: summary

    script:
    """
    python ${projectDir}/bin/compare_coaccessibility.py \\
        --connections-ctrl ctrl/cicero_connections.tsv.gz \\
        --connections-trt trt/cicero_connections.tsv.gz \\
        --control-label '${control_label}' \\
        --treatment-label '${treatment_label}' \\
        --outdir .
    """
}
