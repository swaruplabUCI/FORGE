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
//   CICERO_TRIPLETS_PER_CT       — subset h5ad → triplets TSV.gz
//   CICERO_ESTIMATE_DP_PER_CT    — build shared CDS + estimate distance parameter
//   CICERO_FULL_CHROM_PER_CT     — per-chromosome Cicero models (fan-out)
//   CICERO_JOIN_PER_CT           — rbind chrom connections → final outputs
//   CICERO_FULL_PER_CT           — legacy monolithic runner (retained for reference)


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
// Phase 2a — Build shared CDS + estimate global distance_parameter
//            (one per (cell_type, condition) stratum)
// ---------------------------------------------------------------------------

process CICERO_ESTIMATE_DP_PER_CT {
    tag { "${cell_type.replaceAll(' ', '_')}_${condition}" }
    label 'process_medium'
    errorStrategy 'terminate'

    input:
    tuple val(cell_type), val(condition), path(triplets)

    output:
    tuple val(cell_type), val(condition),
          path("distance_parameter.txt"),
          path("cicero_cds_shared.rds"),
          path("gene_annotation.rds"),
          path("input_cds_ordered.rds"),
          emit: all_out
    path "cds_summary.txt", optional: true
    path "umap_*.pdf",      optional: true

    script:
    """
    export HOME=/tmp/container_home
    mkdir -p /tmp/container_home
    export R_LIBS_USER=""
    cicero_estimate_dp.R \\
      --triplets "${triplets}" \\
      --outdir "." \\
      --num_dim ${params.cicero.num_dim} \\
      --sample_num ${params.cicero.sample_num} \\
      --window_bp 500000 \\
      --distance_constraint 250000 \\
      ${params.cicero.use_partition ? "--use_partition" : ""}
    """
}


// ---------------------------------------------------------------------------
// Phase 2b — Per-chromosome Cicero model, fanned out across all chromosomes
//            (one task per (cell_type, condition, chrom))
// ---------------------------------------------------------------------------

process CICERO_FULL_CHROM_PER_CT {
    tag { "${cell_type.replaceAll(' ', '_')}_${condition}_${chrom}" }
    label 'process_medium'
    errorStrategy 'terminate'

    input:
    tuple val(cell_type), val(condition), val(chrom),
          path(cicero_cds_rds), path(gene_ann_rds), path(dp_file)

    output:
    tuple val(cell_type), val(condition),
          val(chrom), path("conns_${chrom}.tsv.gz"),
          emit: chrom_conns

    script:
    """
    export HOME=/tmp/container_home
    mkdir -p /tmp/container_home
    export R_LIBS_USER=""
    DP=\$(cat ${dp_file})
    cicero_full_chrom.R \\
      --chrom "${chrom}" \\
      --cicero_cds "${cicero_cds_rds}" \\
      --gene_annotation "${gene_ann_rds}" \\
      --dp "\${DP}" \\
      --window_bp 500000 \\
      --outdir "."
    """
}


// ---------------------------------------------------------------------------
// Phase 2c — Collect per-chromosome connections, run CCAN, write final outputs
//            publishDir mirrors the monolithic CICERO_FULL_PER_CT layout
// ---------------------------------------------------------------------------

process CICERO_JOIN_PER_CT {
    tag { "${cell_type.replaceAll(' ', '_')}_${condition}" }
    label 'process_medium'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/cicero/per_ct/${cell_type.replaceAll(/\s+/, '_').replaceAll('/', '-')}/${condition}", mode: 'copy'

    input:
    tuple val(cell_type), val(condition),
          path(chrom_conns),   // staged list of conns_chr*.tsv.gz
          path(ordered_cds)
    val   gtf_path

    output:
    tuple val(cell_type), val(condition),
          path("cicero_connections.tsv.gz"),  emit: connections
    path "CCAN_assignments.tsv.gz",           emit: ccan
    path "input_cds_ordered.rds",             emit: cds
    path "*.pdf",                             emit: plots, optional: true

    script:
    """
    export HOME=/tmp/container_home
    mkdir -p /tmp/container_home
    export R_LIBS_USER=""
    if [ "${ordered_cds}" != "input_cds_ordered.rds" ]; then
        cp -L "${ordered_cds}" "input_cds_ordered.rds"
    fi
    cicero_join.R \\
      --conns_glob "conns_*.tsv.gz" \\
      --cds "input_cds_ordered.rds" \\
      --gtf "${gtf_path}" \\
      --connections_cutoff ${params.cicero.connections_cutoff} \\
      --ccan_min_coaccess ${params.cicero.ccan_min_coaccess} \\
      --outdir "."
    """
}


// ---------------------------------------------------------------------------
// Phase 2 (legacy) — Monolithic Cicero runner; superseded by the three-step
//                    fan-out above; retained here for reference only.
// ---------------------------------------------------------------------------

process CICERO_FULL_PER_CT {
    tag { "${cell_type.replaceAll(' ', '_')}_${condition}" }
    label 'process_high'
    errorStrategy { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    maxRetries 1
    publishDir "${params.outdir}/cicero/per_ct/${cell_type.replaceAll(/\s+/, '_').replaceAll('/', '-')}/${condition}", mode: 'copy'

    input:
    tuple val(cell_type), val(condition), path(triplets)
    val   gtf_path
    val   sample_num

    output:
    tuple val(cell_type), val(condition),
          path("cicero_connections.tsv.gz"),  emit: connections
    path  "CCAN_assignments.tsv.gz",          emit: ccan
    path  "input_cds_ordered.rds",            emit: cds
    path  "*.pdf",                            emit: plots, optional: true

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
