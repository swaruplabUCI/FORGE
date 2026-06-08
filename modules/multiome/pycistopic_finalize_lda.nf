// modules/multiome/pycistopic_finalize_lda.nf
//
// PYCISTOPIC_FINALIZE_LDA — Phase 3c of the pyCisTopic fan-out pipeline.
//
// Collects all Topic*.pkl files from the parallel LDA fan-out, selects the
// best model, binarizes topics, computes DARs and region_sets/, and
// optionally computes gene activity.  Produces the same outputs as the
// former monolithic PYCISTOPIC_MERGE_LDA (cistopic_obj.pkl, region_sets/,
// pycistopic_gene_activity.h5ad).

nextflow.enable.dsl=2

process PYCISTOPIC_FINALIZE_LDA {
    tag "pycistopic_finalize"

    publishDir "${params.outdir}/pycistopic_atac",
        mode: 'copy', overwrite: true,
        pattern: "cistopic_obj.pkl"
    publishDir "${params.outdir}/pycistopic_atac",
        mode: 'copy', overwrite: true,
        pattern: "region_sets"
    publishDir "${params.outdir}/pycistopic_atac",
        mode: 'copy', overwrite: true,
        pattern: "pycistopic_gene_activity.h5ad"

    input:
    path "topic_pkls/*"   // all Topic*.pkl files collected from PYCISTOPIC_RUN_LDA
    path merged_pkl
    path cell_metadata
    path qc_dir
    path blacklist
    val  species

    output:
    path "cistopic_obj.pkl",                              emit: cistopic_obj
    path "region_sets",                                   emit: region_sets
    path "pycistopic_gene_activity.h5ad", optional: true, emit: gene_activity

    script:
    def sel_topics  = params.pycistopic?.selected_topics ?: 40
    def n_cpu       = params.pycistopic?.n_cpu           ?: task.cpus
    def do_gene_act = params.pycistopic?.do_gene_activity ? '--do-gene-activity' : ''
    """
    set -euo pipefail
    echo "=== PYCISTOPIC_FINALIZE_LDA ==="
    echo "  Topic pkls: \$(ls topic_pkls/Topic*.pkl 2>/dev/null | wc -l)"

    python ${projectDir}/bin/run_pycistopic_finalize_lda.py \\
        --merged-pkl      '${merged_pkl}' \\
        --topic-pkl-dir   topic_pkls \\
        --cell-metadata   '${cell_metadata}' \\
        --tss-bed         '${qc_dir}/tss.bed' \\
        --blacklist       '${blacklist}' \\
        --cell-type-col   cell_type_safe \\
        --condition-col   ${params.pycistopic?.condition_col ?: 'condition'} \\
        --selected-topics ${sel_topics} \\
        --n-cpu           ${n_cpu} \\
        --species         ${species} \\
        ${do_gene_act} \\
        --outdir          .

    echo "=== done ==="
    ls -lh cistopic_obj.pkl region_sets/ || true
    """
}
