// modules/cicero/cicero_triplets.nf
//
// Create cicero_triplets.tsv.gz from a cellxpeak_matrix.h5ad

process CICERO_TRIPLETS {
    label 'process_medium'
    errorStrategy 'ignore'
    publishDir "${params.cicero.outdir}", mode: 'copy'

    input:
    path peak_matrix

    output:
    path "cicero_triplets.tsv.gz", emit: triplets

    script:
    """
    make_cicero_triplets.py \\
      --peak-matrix "${peak_matrix}" \\
      --out "cicero_triplets.tsv.gz"
    """
}
