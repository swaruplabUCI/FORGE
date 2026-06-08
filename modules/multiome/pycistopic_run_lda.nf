// modules/multiome/pycistopic_run_lda.nf
//
// PYCISTOPIC_RUN_LDA — Phase 3b of the pyCisTopic fan-out pipeline.
//
// One task per topic count: loads merged_cistopic.pkl and runs MALLET LDA
// for exactly n_topics, saving Topic{n}.pkl.  All topic-count tasks run in
// parallel on separate SLURM nodes; wall-clock cost is that of the slowest
// single model rather than the sequential sum.

nextflow.enable.dsl=2

process PYCISTOPIC_RUN_LDA {
    tag "pycistopic_lda_${n_topics}"

    input:
    val  n_topics
    path merged_pkl

    output:
    path "Topic${n_topics}.pkl", emit: topic_pkl

    script:
    def n_cpu   = params.pycistopic?.n_cpu  ?: task.cpus
    def n_iter  = params.pycistopic?.n_iter ?: 500
    """
    set -euo pipefail
    echo "=== PYCISTOPIC_RUN_LDA n_topics=${n_topics} ==="

    python ${projectDir}/bin/run_pycistopic_run_lda_single.py \\
        --merged-pkl  '${merged_pkl}' \\
        --n-topics    ${n_topics} \\
        --n-cpu       ${n_cpu} \\
        --n-iter      ${n_iter} \\
        --outdir      .

    echo "=== done ==="
    ls -lh Topic${n_topics}.pkl
    """
}
