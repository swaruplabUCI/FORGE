// modules/chromvar/vischromvar.nf

process VIS_CHROMVAR {

    tag "vis_chromvar"

    publishDir "${params.outdir}/chromvar/plots", mode: 'copy', overwrite: true

    input:
        path chromvar_dev  // chromvar_matrix.Chromvar.h5ad from GPU_CHROMVAR

    output:
        path "*.pdf", emit: pdfs
        path "*.tsv", emit: tables

    script:
    """
    vis_chromvar_nf.py \
      --dev-h5ad ${chromvar_dev} \
      --out-dir chromvar_plots

    # Move generated files into CWD so publishDir picks them up
    if [ -d chromvar_plots ]; then
      mv chromvar_plots/* .
    fi
    """
}
