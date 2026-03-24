process CICERO_TARGET_PLOTS {
    tag "cicero_plots"
    label 'process_medium'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/cicero/plots", mode: 'copy'

    input:
    path connections
    path ccan
    path cds
    val gtf_path
    val target_genes

    output:
    path "*.pdf", emit: plots

    // Container assigned in nextflow.config via withName

    script:
    def gene_list = target_genes.join(',')
    """
    run_cicero_target_plots.R \\
        --connections ${connections} \\
        --ccan ${ccan} \\
        --cds ${cds} \\
        --gtf ${gtf_path} \\
        --genes "${gene_list}" \\
        --outdir . \\
        --windows "100000,250000" \\
        --thresholds "0.1,0.05" \\
        --labels "zoom_in,zoom_out"
    """
}
