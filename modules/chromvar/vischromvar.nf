// modules/chromvar/vischromvar.nf

process VIS_CHROMVAR {
    label 'process_high'

    tag "vis_chromvar"

    publishDir "${params.outdir}/chromvar/plots", mode: 'copy', overwrite: true

    input:
        path chromvar_dev  // chromvar_matrix.Chromvar.h5ad from GPU_CHROMVAR

    output:
        path "*.pdf", emit: pdfs
        path "*.tsv", emit: tables

    script:
    def annotation_key = params.atac.marker_file ? 'cell_type' : (params.atac.annotation_method == 'scatanno' ? (params.atac.cell_type_col ?: 'cell_type_prediction') : 'celltypist_prediction')
    """
    vis_chromvar_nf.py \
      --dev-h5ad ${chromvar_dev} \
      --out-dir chromvar_plots \
      --cell-type-key ${annotation_key}

    # Move generated files into CWD so publishDir picks them up
    if [ -d chromvar_plots ]; then
      mv chromvar_plots/* .
    fi
    """
}
