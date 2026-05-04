// modules/scprint/enhancer_footprinting_per_ct.nf
//
// V5 PER-CT ARCHITECTURE — toggleable via params.enhancer_footprinting.use_per_ct.
//
// Drop-in replacement for ENHANCER_FOOTPRINTING that processes one cell type
// per task with all its TFs in a single Python session. Loads printer + models
// + peak matrix + barcode grouping ONCE per ct (not 20× as in the sharded
// version) and amortizes the ~3-5 min setup across all (ct, TF) pairs in that
// cell type. Cuts SLURM task count from 657 → 33 for typical mouse-brain runs.
//
// Wired in main.nf via the ENHANCER_FOOTPRINTING_RECIPES workflow; toggled by
// params.enhancer_footprinting.use_per_ct (default false → legacy sharded path).
//
// Architecture trade-offs and validation status: see the docstring of
// `bin/run_enhancer_footprinting_per_ct.py` (companion script).
//
// HOW TO WIRE IN (when porting to another project):
//   1. include { ENHANCER_FOOTPRINTING_PER_CT } from './modules/scprint/enhancer_footprinting_per_ct'
//   2. Group region_set_manifest by cell_type and emit one tuple per ct:
//      def ch_per_ct = MOTIF_SCAN_ENHANCERS.out.manifest
//          .splitJson(path: 'region_sets')
//          .map { entry -> tuple(entry.cell_type, [tf: entry.tf, bed_file: entry.bed_file]) }
//          .groupTuple()  // (ct, [{tf, bed_file}, ...])
//          .map { ct, entries ->
//              def manifest = new groovy.json.JsonOutput().toJson(entries)
//              tuple(ct, manifest)
//          }
//   3. Call: ENHANCER_FOOTPRINTING_PER_CT(ch_per_ct, printer_ch, peak_matrix_ch, ...)
//   4. Combine outputs same as ENHANCER_FOOTPRINTING.out.{footprints, binding_scores, summary}
//      then feed to BUILD_TF_GENE_NETWORK / AGGREGATE_FP_STATS unchanged.
//
// WHY NOT WIRE INTO SDas PRODUCTION:
//   - 657 sharded ENHANCER_FOOTPRINTING tasks are already cached for SDas v4
//   - Switching architecture invalidates that cache → ~9 cluster-hours recompute
//   - Per-ct savings (~3.5×) only matters for fresh runs
//   - Keep this module accessible for porting to other projects / fresh starts

nextflow.enable.dsl=2

process ENHANCER_FOOTPRINTING_PER_CT {

    tag "${cell_type}"
    label 'process_high'
    // 2026-04-28: bumped 7 → 35 because per-ct fan-out is bounded by ct count
    // (10 active broad cts; up to 33 in fine mode), not the 657 sharded tasks
    // it replaces. 35 lets the entire ct cohort launch in one wave; SLURM
    // queue handles the resource cap if cluster can't absorb all at once.
    maxForks 35
    errorStrategy 'terminate'
    publishDir "${params.outdir}/enhancer_footprinting_per_ct/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}",
               mode: 'copy'

    input:
    // Single tuple per task: (ct label, manifest JSON, list of bed files for
    // every TF in this ct). Channels in main.nf compose all three into one
    // synchronized stream via groupTuple()/map; staging the bed list under
    // tf_bed_files/ keeps the Python script's relative-path resolver
    // unchanged from the smoke-test build.
    tuple val(cell_type), path(tf_bed_manifest_json), path("tf_bed_files/*")
    path printer
    path peak_matrix
    val cell_type_col          // peak matrix obs column (broad scATAnno)
    val control_condition      // optional: enables per-condition outputs
    val treatment_condition    // optional: enables differential outputs

    output:
    // Same emission shape as ENHANCER_FOOTPRINTING (sharded) so downstream
    // BUILD_TF_GENE_NETWORK / AGGREGATE_FP_STATS work unchanged.
    path "enhancer_footprints_*.h5ad",     emit: footprints,     optional: true
    path "enhancer_tfbs_*.h5ad",           emit: binding_scores, optional: true
    path "enhancer_global/*",              emit: global_plots,        optional: true
    path "enhancer_per_condition/*",       emit: per_condition_plots, optional: true
    path "enhancer_differential/*",        emit: differential_plots,  optional: true
    // NEW vs sharded: per-ct summary roll-up (drops AGGREGATE_FP_STATS to a
    // simple .collect() of 33 CSVs instead of walking 657 h5ads)
    path "enhancer_fp_summary_*.csv",      emit: summary

    script:
    """
    # Same scratch-copy pattern as ENHANCER_FOOTPRINTING (sharded) — printer
    # h5ad needs RW for the rust-AnnData backend.
    # cache-bust 2026-05-03: backed=True restored at 3 footprint sites + module-level
    # monkey-patch on scp.tl.get_footprint_score with retry-3x then re-raise on
    # rust-anndata snap.read panic. Per-TF except now catches BaseException with
    # KI/SE/GE allowlist so PanicException is properly isolated.
    export HDF5_USE_FILE_LOCKING=FALSE
    _XDG="\${XDG_CACHE_HOME:-/tmp/cache}"
    mkdir -p "\${_XDG}"
    [ ! -e "\${_XDG}/scprinter" ] && ln -s '${params.scprinter.cache_dir}' "\${_XDG}/scprinter" || true

    SCRATCH="\${TMPDIR:-/tmp}/fp_per_ct_\${SLURM_JOB_ID:-\$\$}"
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
        --outdir .
    """
}
