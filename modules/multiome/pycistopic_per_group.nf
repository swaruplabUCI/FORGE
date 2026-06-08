// modules/multiome/pycistopic_per_group.nf
nextflow.enable.dsl=2

process PYCISTOPIC_PER_GROUP {
    label 'process_medium'
    tag "${row.cell_type_safe}__${row.condition}"

    input:
    tuple val(row), path(phase1_dir)

    output:
    path "*.pkl", emit: pkl

    script:
    """
    set -euo pipefail

    echo "=== PYCISTOPIC_PER_GROUP: ${row.cell_type_safe}__${row.condition} (n_cells=${row.n_cells}) ==="

    python ${projectDir}/bin/run_pycistopic_create_group.py \\
        --cell-type      "${row.cell_type_safe}" \\
        --condition      "${row.condition}" \\
        --cell-metadata  ${phase1_dir}/cell_metadata_for_pycistopic.safe.tsv \\
        --fragments-map  ${phase1_dir}/fragments_map.tsv \\
        --consensus-bed  ${phase1_dir}/consensus_peak_calling/consensus_regions.bed \\
        --qc-dir         ${phase1_dir}/qc \\
        --blacklist      ${phase1_dir}/blacklist.bed \\
        --sample-id-col  sample_id \\
        --cell-type-col  cell_type_safe \\
        --condition-col  condition \\
        --n-cpu          ${task.cpus} \\
        --outdir         .

    echo "=== PKL files produced ==="
    ls -lh *.pkl 2>/dev/null || { echo "ERROR: no PKL file written"; exit 1; }
    """
}
