// modules/multiome/pycistopic_merge_objects.nf
//
// PYCISTOPIC_MERGE_OBJECTS — Phase 3a of the pyCisTopic fan-out pipeline.
//
// Merges all per-group CistopicObject pkl files from Phase 2 into a single
// atlas object with cell-type/condition annotations attached, then saves
// merged_cistopic.pkl for parallel LDA fan-out (one job per topic count).

nextflow.enable.dsl=2

process PYCISTOPIC_MERGE_OBJECTS {
    tag "pycistopic_merge_objects"

    publishDir "${params.outdir}/pycistopic_atac",
        mode: 'copy', overwrite: true,
        pattern: "merged_cistopic.pkl"

    input:
    path "pkls/*"       // all per-group *_cistopic.pkl from Phase 2
    path cell_metadata

    output:
    path "merged_cistopic.pkl", emit: merged_pkl

    script:
    """
    set -euo pipefail
    echo "=== PYCISTOPIC_MERGE_OBJECTS ==="
    echo "  pkl count: \$(ls pkls/*.pkl 2>/dev/null | wc -l)"

    python ${projectDir}/bin/run_pycistopic_merge_objects.py \\
        --pkl-dir       pkls \\
        --cell-metadata '${cell_metadata}' \\
        --outdir        .

    echo "=== done ==="
    ls -lh merged_cistopic.pkl
    """
}
