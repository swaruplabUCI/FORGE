// modules/scprint/resolve_coordinates.nf

process RESOLVE_GENE_COORDINATES {
    label 'process_low'
    publishDir "${params.outdir}/scprinter/coordinates", mode: 'copy'
    
    input:
    val species
    val target_genes
    path manual_coords_json  // Optional manual overrides
    
    output:
    path "gene_coordinates.json", emit: coordinates
    path "gene_coordinates.tsv", emit: coordinates_tsv
    path "coordinate_report.txt", emit: report
    
    script:
    def genes_str = target_genes.join(' ')
    def manual_arg = manual_coords_json.name != 'NO_FILE' ? "--manual-coords ${manual_coords_json}" : ""
    
    """
    # Resolve coordinates
    gene_coordinate_resolver.py \\
        --genes ${genes_str} \\
        --species ${species} \\
        --gtf-human '${params.scprinter.gtf_human}' \\
        --gtf-mouse '${params.scprinter.gtf_mouse}' \\
        --promoter-upstream ${task.ext.promoter_upstream ?: params.scprinter.promoter_upstream} \\
        --promoter-downstream ${task.ext.promoter_downstream ?: params.scprinter.promoter_downstream} \\
        ${manual_arg} \\
        --output gene_coordinates.json \\
        --format both \\
        > coordinate_report.txt
    
    # Also create a simplified format for the footprinting script
    cat coordinate_report.txt
    """
}
