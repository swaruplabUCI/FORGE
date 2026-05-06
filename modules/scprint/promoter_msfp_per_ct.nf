// modules/scprint/promoter_msfp_per_ct.nf
//
// Per-(cell_type, gene) promoter MSFP compute, stratified by condition.
// One task per cell type loads the printer + peak matrix once and iterates
// every gene in params.scprinter.target_genes, writing one standalone h5ad
// per gene shaped for render_promoter_msfp_overlay.py (per-condition rows
// in obs["name"], single-locus obsm[<chr:start-end>], uns["scales"]).
//
// Gated upstream by params.promoter_overlay.enabled AND
// params.differential.{control,treatment}_condition both set AND
// non-empty params.scprinter.target_genes — single-condition runs and
// untargeted runs auto-skip.

nextflow.enable.dsl=2

process PROMOTER_MSFP_PER_CT {

    tag "${cell_type}"
    label 'process_high'
    maxForks 14
    errorStrategy 'terminate'
    publishDir "${params.outdir}/promoter_msfp_overlay/per_cond_h5ad",
               mode: 'copy'

    input:
    tuple val(cell_type), val(genes_csv)
    path printer
    path peak_matrix
    path coords_json
    val cell_type_col
    val condition_col
    val control_condition
    val treatment_condition

    output:
    tuple val(cell_type), path("promoter_fp_*.h5ad"), emit: fp, optional: true

    script:
    """
    export HDF5_USE_FILE_LOCKING=FALSE
    _XDG="\${XDG_CACHE_HOME:-/tmp/cache}"
    mkdir -p "\${_XDG}"
    [ ! -e "\${_XDG}/scprinter" ] && ln -s '${params.scprinter.cache_dir}' "\${_XDG}/scprinter" || true

    SCRATCH="\${PWD}/scratch_promoter_msfp"
    mkdir -p "\${SCRATCH}"
    trap 'rm -rf "\${SCRATCH}"' EXIT

    PRINTER_LOCAL="\${SCRATCH}/printer_local_copy.h5ad"
    cp -L '${printer}' "\${PRINTER_LOCAL}"
    chmod u+w "\${PRINTER_LOCAL}"
    PEAK_LOCAL="\${SCRATCH}/peak_matrix_local_copy.h5ad"
    cp -L '${peak_matrix}' "\${PEAK_LOCAL}"
    chmod u+w "\${PEAK_LOCAL}"

    python ${projectDir}/bin/run_promoter_msfp_per_condition.py \\
        --printer "\${PRINTER_LOCAL}" \\
        --peak-matrix "\${PEAK_LOCAL}" \\
        --coords '${coords_json}' \\
        --genes '${genes_csv}' \\
        --cell-type '${cell_type}' \\
        --cell-type-col '${cell_type_col}' \\
        --condition-col '${condition_col}' \\
        --control-condition '${control_condition}' \\
        --treatment-condition '${treatment_condition}' \\
        --genome '${params.scprinter.genome}' \\
        --cpus ${task.cpus} \\
        --outdir .
    """

    stub:
    def safe_ct = cell_type.replaceAll(/[\/\s\(\)]+/, '_')
    """
    for g in \$(echo '${genes_csv}' | tr ',' ' '); do
        touch "promoter_fp_\${g}__${safe_ct}.h5ad"
    done
    """
}
