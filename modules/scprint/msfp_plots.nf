// modules/scprint/msfp_plots.nf
process SCPRINTER_MSFP_PLOTS {
    label 'process_medium'
    publishDir "${params.outdir}/scprinter/msfp_celltype_plots", mode: 'copy'

    input:
    path printer_h5ad
    path peak_matrix
    path coordinates_json
    val  genes

    output:
    path "*.png", emit: plots

    script:
    // Ensure we pass a comma-separated string to match your argparse expectation
    def genes_str = (genes instanceof List) ? genes.join(',') : genes.toString()

    """
    python ${projectDir}/bin/scprinter_msfp_celltype_plots.py \\
      --printer '${printer_h5ad}' \\
      --peak-matrix '${peak_matrix}' \\
      --coords '${coordinates_json}' \\
      --genes '${genes_str}' \\
      --cell-type-key '${params.cellchat.cell_type_key ?: "cell_type"}' \\
      --outdir . \\
      --cpus ${task.cpus} \\
      --scale ${params.scprinter.msfp_scale ?: 50}
    """
}
