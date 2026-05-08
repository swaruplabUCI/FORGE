// modules/cicero/cicero_per_ct.nf
//
// Per-CT × condition Cicero pipeline.
// Runs Cicero once per (cell_type_v2, condition) stratum, producing independent
// co-accessibility models for each cell population under each condition.
//
// Strata are defined by ct_annotation.csv (bin/build_ct_annotation_v2.py):
//   any obs column (--cell-type-col) × condition. CTs below the abundance
//   floor or in --exclude-labels are marked 'EXCLUDED' in the CSV and never
//   reach this module (filtered in the main workflow strata_ch).
//
// Per-stratum floor: params.cicero_per_ct.min_cells_per_stratum (default 250).
// Strata with fewer than min_cells_per_stratum cells per (CT, condition) pair
// exit 77 → errorStrategy ignore → cleanly skipped without pipeline failure.
//
// Processes:
//   CICERO_TRIPLETS_PER_CT  — subset h5ad → triplets TSV.gz
//   CICERO_FULL_PER_CT      — triplets → connections, CCANs, plots


// ---------------------------------------------------------------------------
// Phase 1 — Build cell × peak triplets for one (cell_type, condition) stratum
// ---------------------------------------------------------------------------

process CICERO_TRIPLETS_PER_CT {
    tag { "${cell_type.replaceAll(' ', '_')}_${condition}" }
    label 'process_medium'
    errorStrategy { task.exitStatus == 77 ? 'ignore' : 'terminate' }
    publishDir "${params.outdir}/cicero/per_ct/${cell_type.replaceAll(/\s+/, '_').replaceAll('/', '-')}/${condition}", mode: 'copy'

    input:
    path  peak_matrix      // peak_matrix_annotated.h5ad (full dataset; broadcast value)
    path  annotation_csv   // ct_annotation_v2.csv  (barcode, condition, cell_type_v2)
    tuple val(cell_type),
          val(condition)   // one stratum pair emitted by strata_ch

    output:
    tuple val(cell_type), val(condition),
          path("cicero_triplets_*.tsv.gz"), emit: triplets

    script:
    def ct_safe  = cell_type.replaceAll(/\s+/, '_').replaceAll('/', '-')
    def min_cells = params.cicero_per_ct?.min_cells_per_stratum ?: 250
    """
    make_cicero_triplets_per_ct.py \\
        --peak-matrix    "${peak_matrix}" \\
        --annotation-csv "${annotation_csv}" \\
        --cell-type      "${cell_type}" \\
        --condition      "${condition}" \\
        --min-cells      ${min_cells} \\
        --out            "cicero_triplets_${ct_safe}_${condition}.tsv.gz"
    """
}


// ---------------------------------------------------------------------------
// Phase 2 — Run monolithic Cicero on the stratum triplets
//
// Using run_cicero_full.R (monolithic) rather than ESTIMATE_DP + FULL_CHROM
// because per-CT strata are small (200–50k cells); the chromosome fan-out
// speedup is worth its overhead only for the full pooled dataset.
// ---------------------------------------------------------------------------

process CICERO_FULL_PER_CT {
    tag { "${cell_type.replaceAll(' ', '_')}_${condition}" }
    label 'process_high'
    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 1
    publishDir "${params.outdir}/cicero/per_ct/${cell_type.replaceAll(/\s+/, '_').replaceAll('/', '-')}/${condition}", mode: 'copy'

    input:
    tuple val(cell_type), val(condition), path(triplets)
    val   gtf_path      // absolute path to GTF (params.cicero.gtf_full)
    val   sample_num    // Cicero k-NN (params.cicero.sample_num)

    output:
    tuple val(cell_type), val(condition),
          path("cicero_connections.tsv.gz"),  emit: connections
    path  "CCAN_assignments.tsv.gz",          emit: ccan
    path  "input_cds_ordered.rds",            emit: cds
    path  "*.pdf",                            emit: plots,   optional: true

    script:
    """
    export HOME=/tmp/container_home
    mkdir -p /tmp/container_home
    export R_LIBS_USER=""
    run_cicero_full.R \\
        --triplets "${triplets}" \\
        --outdir   "." \\
        --gtf      "${gtf_path}" \\
        --num_dim            ${params.cicero.num_dim} \\
        --sample_num         ${sample_num} \\
        --connections_cutoff ${params.cicero.connections_cutoff} \\
        --ccan_min_coaccess  ${params.cicero.ccan_min_coaccess} \\
        ${params.cicero.use_partition ? "--use_partition" : ""}
    """
}
