process PREPARE_REFERENCE {
    label 'process_high'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path query_h5ad
    val species

    output:
    path "Refmapping_*HVG_*.h5ad", emit: prepared_ref
    path "Refmapping_query.h5ad", emit: prepared_query
    path "label_key.txt", emit: label_key

    script:
    // FIX-29b: Removed single quotes that prevented Nextflow interpolation.
    //          '${params.ref_dir_human_integrated}' was passed literally to bash
    //          as the string "${params.ref_dir_human_integrated}" which bash can't resolve.
    def ref_path = species == 'human' ?
        params.ref_dir_human_integrated :
        params.ref_dir_mouse_integrated

    """
    python3 ${baseDir}/bin/prepare_reference.py \\
        --query ${query_h5ad} \\
        --species ${species} \\
        --ref_path ${ref_path} \\
        --results_dir . \\
        --n_top_genes ${params.n_top_genes ?: 12000}
    """
}
