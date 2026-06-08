// modules/multiome/pycistopic_merge_lda.nf
nextflow.enable.dsl=2

process PYCISTOPIC_MERGE_LDA {
    label 'process_high_memory'
    tag "pycistopic_merge_lda"

    publishDir "${params.outdir}/pycistopic", mode: 'copy', overwrite: true, pattern: "cistopic_obj.pkl"
    publishDir "${params.outdir}/pycistopic", mode: 'copy', overwrite: true, pattern: "region_sets"
    publishDir "${params.outdir}/pycistopic", mode: 'copy', overwrite: true, pattern: "pycistopic_gene_activity.h5ad"

    input:
    path pkls            // collected list of per-group *.pkl files from PYCISTOPIC_PER_GROUP
    path phase1_dir      // Phase 1 output dir (provides cell_metadata, tss.bed, blacklist)
    val  species
    val  topics
    val  selected_topics

    output:
    path "cistopic_obj.pkl",                          emit: cistopic_obj
    path "region_sets",                                emit: region_sets
    path "pycistopic_gene_activity.h5ad", optional: true, emit: gene_activity

    script:
    """
    set -euo pipefail

    echo "=== PYCISTOPIC_MERGE_LDA ==="
    echo "PKL files staged:"
    ls -lh *.pkl 2>/dev/null || true
    N_PKLS=\$(ls *.pkl 2>/dev/null | wc -l)
    if [[ "\$N_PKLS" -eq 0 ]]; then
        echo "ERROR: no PKL files staged"
        exit 1
    fi

    # Stage PKLs into a subdirectory so --pkl-dir is unambiguous
    mkdir -p pkl_input
    for f in *.pkl; do mv "\$f" pkl_input/; done

    python ${projectDir}/bin/run_pycistopic_merge_lda.py \\
        --pkl-dir         pkl_input \\
        --cell-metadata   ${phase1_dir}/cell_metadata_for_pycistopic.safe.tsv \\
        --tss-bed         ${phase1_dir}/qc/tss.bed \\
        --blacklist       ${phase1_dir}/blacklist.bed \\
        --cell-type-col   cell_type_safe \\
        --condition-col   condition \\
        --topics          ${topics} \\
        --selected-topics ${selected_topics} \\
        --n-cpu           ${task.cpus} \\
        --species         ${species} \\
        --do-gene-activity \\
        --outdir          .

    echo "=== PYCISTOPIC_MERGE_LDA done ==="
    ls -lh cistopic_obj.pkl region_sets/ pycistopic_gene_activity.h5ad 2>/dev/null || true

    N_DARS=\$(ls region_sets/DARs_cell_type/*.bed 2>/dev/null | wc -l)
    echo "[merge_lda] \${N_DARS} DARs_cell_type BEDs produced"
    if [[ "\$N_DARS" -eq 0 ]]; then
        echo "WARNING: DARs_cell_type/ is empty — barcode reconciliation may have failed"
    fi
    """
}
