// modules/atac/atac_scatanno.nf
//
// scATAnno-based ATAC cell type annotation using reference peak atlases.
// Replaces CellTypist when params.atac.annotation_method == 'scatanno'.
// Requires a SnapATAC2 h5ad with fragment access for peak re-counting.

process ATAC_SCATANNO {
    tag "atac_scatanno"
    label 'process_high'
    publishDir "${params.outdir}/atac/final", mode: 'copy'

    input:
    path snap_h5ad         // SnapATAC2 AnnDataSet or individual h5ad with fragment_paired
    path reference_atlas   // scATAnno reference atlas h5ad (e.g., PBMC_reference_atlas_final.h5ad)

    output:
    path "scatanno_annotations.h5ad", emit: annotations
    path "scatanno_annotations.csv", emit: annotations_csv

    script:
    def dist_thresh = params.scatanno?.distance_threshold ?: 95
    def unc_thresh  = params.scatanno?.uncertainty_threshold ?: 0.5
    def atlas_name  = params.scatanno?.atlas_name ?: 'general'
    def n_dims      = params.scatanno?.n_dims ?: 30
    def knn         = params.scatanno?.knn_neighbors ?: 30
    """
    python ${projectDir}/bin/run_scatanno.py \\
        --query ${snap_h5ad} \\
        --reference ${reference_atlas} \\
        --output scatanno_annotations.h5ad \\
        --csv-output scatanno_annotations.csv \\
        --atlas ${atlas_name} \\
        --distance-threshold ${dist_thresh} \\
        --uncertainty-threshold ${unc_thresh} \\
        --n-dims ${n_dims} \\
        --knn-neighbors ${knn}
    """
}
