// modules/visualization/export_atac_bigwigs.nf
//
// Produce per-cell-type pseudobulk bigWigs from the ATAC fragment AnnDataSet
// using snapatac2's snap.ex.export_coverage. The ATAC-only analog of FORGE's
// PYCISTOPIC_PREPARE — no RNA dependency.
//
// Feeds PREPARE_ENHANCER_VIZ_TRACKS via --bigwig-dir.

nextflow.enable.dsl=2

process EXPORT_ATAC_BIGWIGS {

    tag "export_atac_bigwigs"
    label 'process_high'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/enhancer_viz/bigwigs", mode: 'copy'

    input:
    path anndataset          // atac_complete.h5ads
    path annotated_peaks     // peak_matrix_annotated.h5ad
    val  cell_type_col       // e.g. 'cell_type_prediction' (scATAnno) / 'cell_type_broad'
    val  min_cells           // skip groups below this cell count

    output:
    path "bigwigs/*.bw",             emit: bigwigs
    path "bigwigs/manifest.json",    emit: manifest

    script:
    def blacklist_arg = params.enhancer_viz.blacklist ? "--blacklist ${params.enhancer_viz.blacklist}" : ''
    """
    python ${projectDir}/bin/export_atac_bigwigs.py \\
        --anndataset '${anndataset}' \\
        --annotated-peaks '${annotated_peaks}' \\
        --cell-type-col '${cell_type_col}' \\
        --min-cells ${min_cells} \\
        --bin-size ${params.enhancer_viz.bin_size ?: 10} \\
        --normalization ${params.enhancer_viz.normalization ?: 'RPKM'} \\
        ${blacklist_arg} \\
        --n-jobs ${task.cpus} \\
        --outdir bigwigs
    """
}
