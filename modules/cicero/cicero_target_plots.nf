// modules/cicero/cicero_target_plots.nf
//
// Generate gene-specific Cicero connection plots for arbitrary target genes.

process CICERO_TARGET_PLOTS {
    tag "cicero_plots"
    label 'process_medium'
    // Viz overhaul §1: split outputs by purpose. PDFs go under plots/gene_specific,
    // debug TSVs under plots/debug, summary CSV at the top of plots/.
    publishDir "${params.cicero.outdir}/plots/gene_specific", mode: 'copy', pattern: "*.pdf"
    publishDir "${params.cicero.outdir}/plots/debug",         mode: 'copy', pattern: "*_debug_*.tsv"
    publishDir "${params.cicero.outdir}/plots",               mode: 'copy', pattern: "cicero_target_gene_summary.csv"

    input:
    path connections
    path ccan
    path cds
    val gtf_path
    val target_genes

    output:
    path "*.pdf",                              emit: gene_plots
    path "*_debug_*.tsv",                      emit: debug_tsvs, optional: true
    path "cicero_target_gene_summary.csv",     emit: summary,    optional: true

    script:
    def genes_arg = target_genes.join(',')
    def windows    = params.cicero.plot_windows    ?: "100000,250000,500000"
    def thresholds = params.cicero.plot_thresholds ?: "0.10,0.05,0.025"
    def labels     = params.cicero.plot_labels     ?: "strict,permissive,exploratory"
    """
    run_cicero_target_plots.R \\
      --connections "${connections}" \\
      --ccan "${ccan}" \\
      --cds "${cds}" \\
      --gtf "${gtf_path}" \\
      --genes "${genes_arg}" \\
      --outdir "." \\
      --windows "${windows}" \\
      --thresholds "${thresholds}" \\
      --labels "${labels}"
    """
}
