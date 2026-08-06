// modules/scprint/enhancer_footprinting_per_ct_strip.nf
//
// Second-pass ENHANCER_FOOTPRINTING_PER_CT run that populates obsm in the
// output h5ads with full 3D MSFP tensors for per-CT strip target genes.
//
// Identical to ENHANCER_FOOTPRINTING_PER_CT except:
//   - Takes val(strip_genes) as an additional input (CT-specific gene list
//     from RANK_ENHANCER_STRIP_GENES instead of the config-level static list).
//   - publishDir points to enhancer_footprinting_per_ct_strip/ to keep the
//     first-pass h5ads (used by BUILD_TF_GENE_NETWORK / AGGREGATE_FP_STATS)
//     intact and cache-stable.
//   - Separate process name avoids Nextflow DAG collision with first pass.
//
// Gating: only fires for CTs with non-empty strip_genes (filter in main.nf
// before invoking). CTs already passed upstream cell-count and ChromVAR z-score
// gates via RANK_ENHANCER_STRIP_GENES input (binding_scores only exist when
// ENHANCER_FOOTPRINTING_PER_CT succeeded for that CT).

nextflow.enable.dsl=2

process ENHANCER_FOOTPRINTING_PER_CT_STRIP {

    tag "${cell_type}"
    label 'process_high'
    maxForks 35
    errorStrategy 'terminate'
    publishDir "${params.outdir}/enhancer_footprinting_per_ct_strip/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}",
               mode: 'copy'

    input:
    // strip_genes is part of the per-CT tuple to guarantee it stays synchronized
    // with the manifest and beds for that CT (no separate val channel needed).
    tuple val(cell_type), path(tf_bed_manifest_json), path("tf_bed_files/*"), val(strip_genes)
    path printer
    path peak_matrix
    val cell_type_col
    val control_condition
    val treatment_condition

    output:
    path "enhancer_footprints_*.h5ad",  emit: footprints,     optional: true
    path "enhancer_tfbs_*.h5ad",        emit: binding_scores, optional: true
    path "enhancer_fp_summary_*.csv",   emit: summary,        optional: true

    script:
    """
    export HDF5_USE_FILE_LOCKING=FALSE
    _XDG="\${XDG_CACHE_HOME:-/tmp/cache}"
    mkdir -p "\${_XDG}"
    [ ! -e "\${_XDG}/scprinter" ] && ln -s '${params.scprinter.cache_dir}' "\${_XDG}/scprinter" || true

    SCRATCH="\${PWD}/scratch_fp_per_ct"
    mkdir -p "\${SCRATCH}"
    trap 'rm -rf "\${SCRATCH}"' EXIT

    PRINTER_LOCAL="\${SCRATCH}/printer_local_copy.h5ad"
    cp -L '${printer}' "\${PRINTER_LOCAL}"
    chmod u+w "\${PRINTER_LOCAL}"
    PEAK_LOCAL="\${SCRATCH}/peak_matrix_local_copy.h5ad"
    cp -L '${peak_matrix}' "\${PEAK_LOCAL}"
    chmod u+w "\${PEAK_LOCAL}"

    python ${projectDir}/bin/run_enhancer_footprinting_per_ct.py \\
        --printer-path "\${PRINTER_LOCAL}" \\
        --peak-matrix "\${PEAK_LOCAL}" \\
        --cell-type '${cell_type}' \\
        --tf-bed-manifest '${tf_bed_manifest_json}' \\
        --cache-dir '${params.scprinter.cache_dir}' \\
        --genome '${params.scprinter.genome}' \\
        --cell-type-col '${cell_type_col}' \\
        --cpus ${task.cpus} \\
        --pfm-path '${params.scprinter.pfms}' \\
        --gtf '${params.species == "mouse" ? params.scprinter.gtf_mouse : params.scprinter.gtf_human}' \\
        --control-condition '${control_condition}' \\
        --treatment-condition '${treatment_condition}' \\
        --strip-target-genes '${strip_genes}' \\
        --strip-context-bp ${params.msfp_strip?.context_bp ?: 500000} \\
        --outdir .
    """
}
