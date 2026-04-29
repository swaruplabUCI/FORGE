// ============================================================================
// SCPRINTER_FOOTPRINTING  v2.2.0
// TF footprinting analysis with optional differential conditions.
//
// FIX-23:  Added control_condition and treatment_condition as explicit inputs.
// FIX-23b: cell_type per-cell-type labels in discovery mode.
// FIX-35:  HDF5 file lock contention — 7 concurrent tasks all read the same
//          printer h5ad symlink, causing POSIX lock conflicts on parallel FS
//          (DFS7/NFS/Lustre).  HDF5_USE_FILE_LOCKING=FALSE disables locking
//          for read-only access.  Safe because all tasks only read the printer.
// FIX-38:  Rewrote global mode for focal cell-type + TFBS binding scores.
//          New outputs: tfbs_*.h5ad, composite plots (A+B+C).
// FIX-40:  Copy printer h5ad to local work dir before opening.  scprinter's
//          Rust anndata backend opens h5ad read-write; if a task panics
//          mid-write the shared file is permanently corrupted.  Copying
//          isolates each task so crashes can't destroy the original.
// FIX-45:  PYTHONUNBUFFERED=1 so print() statements stream to Nextflow logs
//          in real-time instead of being buffered until process exit.
// ============================================================================

process SCPRINTER_FOOTPRINTING {
    label 'process_high'

    // FIX-36: Each task loads the printer h5ad (read-only) + model files +
    // peak matrix independently.  Each is a separate SLURM job on a separate
    // node, so 7 parallel tasks = 7 × 300 GB across 7 nodes (no shared mem).
    // Legacy workflow used hugemem/666 GB/maxForks 1; Option B parallelizes
    // across standard partition nodes instead.
    // FIX-38: Each task now also loads the TFBS binding score model (~200MB).
    // FIX-40: Each task copies the printer h5ad (~10-30 GB) to local scratch.
    // Memory budget: printer copy + peak matrix + disp model + TFBS model
    // + footprint computation buffers.  250 GB per node is sufficient.
    maxForks 7

    // FIX-43: Sanitize cell_type — slashes/spaces/parens in names break filesystem paths
    // Viz overhaul §3a: plots are routed by mode at write time into the
    // global/, per_condition/, differential/ subdirs — publishDir copies
    // recursively, preserving the structure under the per-cell-type leaf.
    publishDir "${params.outdir}/scprinter/footprints/${cell_type.replaceAll(/[\/\s\(\)]+/, '_')}", mode: 'copy'

    input:
    path peak_matrix
    path metadata
    path printer
    val cell_type
    val target_genes
    path coordinates_json
    val control_condition      // FIX-23: empty string → global mode
    val treatment_condition    // FIX-23: empty string → global mode
    path tf_mapping            // Optional: tf_target_genes.json for TF annotation on plots
    path cicero_connections    // Optional: Cicero connections TSV.gz for CCAN arc panel
    path pfm_path              // Optional: JASPAR PFM file for TF motif logo panel
    val cell_type_col          // FIX-46: obs column name for cell types

    output:
    // FIX-43: Output globs already use wildcards — safe for cell types with '/'
    path "footprints_*.h5ad", emit: footprints
    path "tfbs_*.h5ad", emit: tfbs, optional: true                           // FIX-38
    path "enhancer_footprints_*.h5ad", emit: enhancer_footprints, optional: true
    // Viz overhaul §3a: split plots by mode. Each subdir is optional because
    // global mode produces only global/, differential mode produces all three.
    path "global/*",         emit: global_plots,        optional: true
    path "per_condition/*",  emit: per_condition_plots, optional: true
    path "differential/*",   emit: differential_plots,  optional: true
    path "footprint_summary.csv", emit: summary
    path "enhancer_footprint_summary.csv", emit: enhancer_summary, optional: true

    script:
    def genes_str = target_genes.join(',')

    // FIX-23: Condition arguments are optional — omitted in global mode
    def condition_args = (control_condition && treatment_condition) ?
        "--control-condition '${control_condition}' --treatment-condition '${treatment_condition}'" :
        ""

    def mode_label = condition_args ? "differential" : "global"
    def tf_mapping_arg = (tf_mapping.name != 'NO_FILE') ? "--tf-mapping '${tf_mapping}'" : ""
    def cicero_arg = (cicero_connections.name != 'NO_FILE') ? "--cicero-connections '${cicero_connections}'" : ""
    def pfm_arg = (pfm_path.name != 'NO_FILE') ? "--pfm-path '${pfm_path}'" : ""

    """
    # FIX-45: Ensure Python print() streams to Nextflow logs in real-time
    export PYTHONUNBUFFERED=1

    # FIX-35: Disable HDF5 POSIX file locking to prevent contention when
    # multiple concurrent tasks read the same printer h5ad on a parallel FS.
    export HDF5_USE_FILE_LOCKING=FALSE

    # FIX-35c: scprinter's model loader (load_disp_model, TFBS models, etc.)
    # uses pooch and caches to SCPRINTER_DATA.  Under Singularity --contain,
    # this defaults to /tmp/cache/scprinter (ephemeral), causing all 7 tasks
    # to re-download ~500 MB of models every run.  Setting SCPRINTER_DATA to
    # the pre-downloaded ref dir on DFS7 skips downloads entirely.
    export SCPRINTER_DATA='${params.scprinter.cache_dir}'

    echo "FIX-23b: scPRINTER footprinting — cell_type=${cell_type}, mode=${mode_label}"

    # FIX-40: Copy the printer h5ad to a local file in the work directory.
    # Nextflow symlinks the printer into each task's work dir, but all 7 tasks
    # point to the same physical file on DFS7.  scprinter's Rust anndata backend
    # opens h5ad for read-write and modifies printer.uns["footprints"] during
    # get_footprint_score().  If any task panics mid-write (e.g. the FIX-39
    # HDF5 cache bug), the shared file is permanently corrupted with
    # "address of object past end of allocation".  Copying gives each task its
    # own disposable h5ad that can't corrupt the BUILD_PRINTER output.
    #
    # Also copy the _supp directory if it exists (scprinter stores backed
    # footprint/binding-score data in {printer_name}_supp/).
    PRINTER_LOCAL="printer_local_copy.h5ad"
    echo "FIX-40: Copying printer to local work dir..."
    cp -L '${printer}' "\${PRINTER_LOCAL}"
    chmod u+w "\${PRINTER_LOCAL}"
    if [ -d '${printer.baseName}_supp' ]; then
        cp -rL '${printer.baseName}_supp' "printer_local_copy_supp"
    fi
    echo "FIX-40: Copy complete (\$(du -sh \${PRINTER_LOCAL} | cut -f1))"

    run_scprinter_footprinting.py \\
        --peak-matrix '${peak_matrix}' \\
        --cell-type '${cell_type}' \\
        --target-genes '${genes_str}' \\
        --cache-dir '${params.scprinter.cache_dir}' \\
        --genome '${params.scprinter.genome}' \\
        --printer-path "\${PRINTER_LOCAL}" \\
        --coordinates-file '${coordinates_json}' \\
        --cell-type-label '${cell_type}' \\
        ${condition_args} \\
        ${tf_mapping_arg} \\
        ${cicero_arg} \\
        ${pfm_arg} \\
        --cell-type-col '${cell_type_col}' \
        --cpus ${task.cpus}

    # FIX-40: Clean up local copy to save disk space in work dir
    rm -f "\${PRINTER_LOCAL}"
    rm -rf "printer_local_copy_supp"
    """
}
