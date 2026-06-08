// modules/multiome/pycistopic_prepare_atac.nf
//
// PYCISTOPIC_ATAC_PREPARE — ATAC-only Phase 1 entrypoint.
//
// Identical to PYCISTOPIC_PHASE1 except cell metadata is built directly
// from the ATAC peak matrix h5ad (via build_atac_cell_metadata.py) rather
// than from the RNA h5ad. Used when running pyCisTopic in ATAC-only mode
// (params.pycistopic.atac_only = true), e.g. for datasets without paired RNA
// or for post-hoc ATAC re-analysis under new taxonomy labels.
//
// Outputs the same phase1_out/ directory structure as PYCISTOPIC_PHASE1 so
// that PYCISTOPIC_PER_GROUP and the Phase 3 sub-modules can be shared
// across both modes without modification.

nextflow.enable.dsl=2

process PYCISTOPIC_ATAC_PREPARE {
    label 'process_high_memory'
    tag "pycistopic_atac_prepare"

    publishDir "${params.outdir}/pycistopic/phase1", mode: 'copy', overwrite: true, pattern: "phase1_out/group_list.tsv"
    publishDir "${params.outdir}/pycistopic/phase1", mode: 'copy', overwrite: true, pattern: "phase1_out/consensus_peak_calling/consensus_regions.bed"
    publishDir "${params.outdir}/pycistopic/bigwigs", mode: 'copy', overwrite: true, pattern: "phase1_out/consensus_peak_calling/pseudobulk_bw_files/*.bw"

    input:
    path peak_matrix_h5ad   // ATAC peak matrix (cell types read from obs directly)
    path sample_metadata
    val  species
    path blacklist_bed
    path gtf_file
    val  cell_type_key
    val  condition_col
    val  min_cells

    output:
    path "phase1_out",                                                                    emit: phase1_dir
    path "phase1_out/group_list.tsv",                                                     emit: group_list
    path "phase1_out/cell_metadata_for_pycistopic.safe.tsv",                             emit: cell_metadata
    path "phase1_out/qc",                                                                 emit: qc_dir
    path "phase1_out/blacklist.bed",                                                      emit: blacklist
    path "phase1_out/consensus_peak_calling/pseudobulk_bw_files/*.bw", optional: true,   emit: pseudobulk_bigwigs

    script:
    """
    set -euo pipefail

    echo "=== PYCISTOPIC_ATAC_PREPARE ==="
    echo "Peak matrix: ${peak_matrix_h5ad}"
    echo "Species:     ${species}"
    echo "Condition:   ${condition_col}"
    echo "Min cells:   ${min_cells}"

    # 1) Build cell-level metadata from ATAC peak matrix (no RNA h5ad required)
    python ${projectDir}/bin/build_atac_cell_metadata.py \\
        --atac-h5ad       ${peak_matrix_h5ad} \\
        --sample-metadata ${sample_metadata} \\
        --cell-type-key   ${cell_type_key} \\
        --condition-key   ${condition_col} \\
        --out-tsv         cell_metadata_for_pycistopic.tsv

    echo "Cell metadata (head):"
    head cell_metadata_for_pycistopic.tsv || true

    # 2) Add cell_type_safe column
python - <<'PY'
import pandas as pd, re

df = pd.read_csv("cell_metadata_for_pycistopic.tsv", sep="\\t")

def safe(s):
    s = str(s)
    s = re.sub(r"\\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s

df["cell_type_safe"] = df["cell_type"].map(safe)
df.to_csv("cell_metadata_for_pycistopic.safe.tsv", sep="\\t", index=False)
print(df[["cell_type", "cell_type_safe"]].drop_duplicates().head(20))
PY

    # 3) Build fragments_map.tsv from sample_metadata
python - << EOF
import os
import pandas as pd

df = pd.read_csv("${sample_metadata}")
df = df.drop_duplicates(subset=["sample_id"]).copy()

missing = [c for c in ["sample_id", "fragment_file"] if c not in df.columns]
if missing:
    raise ValueError(f"sample_metadata missing columns: {missing}. Found: {list(df.columns)}")

paths = []
for _, row in df.iterrows():
    sample_id = str(row["sample_id"])
    frag_file = str(row["fragment_file"])
    data_dir  = str(row.get("data_dir", "."))
    full_path = frag_file if (os.path.isabs(frag_file) and os.path.exists(frag_file)) \\
                          else os.path.join(data_dir, frag_file)
    if not os.path.exists(full_path):
        raise FileNotFoundError(f"Fragments file not found for '{sample_id}': {full_path}")
    paths.append((sample_id, full_path))

pd.DataFrame(paths, columns=["sample_id", "fragments_path"]).to_csv(
    "fragments_map.tsv", sep="\\t", index=False)
print(f"Wrote fragments_map.tsv for {len(paths)} samples")
EOF

    echo "fragments_map.tsv (head):"
    head fragments_map.tsv || true

    # 4) Normalise blacklist (handle .gz)
    if [[ "${blacklist_bed}" == *.gz ]]; then
        gunzip -c "${blacklist_bed}" > blacklist.bed
        BLACKLIST="blacklist.bed"
    else
        BLACKLIST="${blacklist_bed}"
    fi

    # 5) Phase 1: consensus peak calling + QC + group_list
    mkdir -p phase1_out
    python ${projectDir}/bin/run_pycistopic_phase1.py \\
        --fragments-map   fragments_map.tsv \\
        --cell-metadata   cell_metadata_for_pycistopic.safe.tsv \\
        --species         ${species} \\
        --sample-id-col   sample_id \\
        --cell-type-col   cell_type_safe \\
        --variable        cell_type_safe \\
        --blacklist-bed   "\$BLACKLIST" \\
        --gtf             ${gtf_file} \\
        --n-cpu           ${task.cpus} \\
        --condition-col   ${condition_col} \\
        --min-cells       ${min_cells} \\
        --outdir          phase1_out

    # 6) Stage helper files into phase1_out so PYCISTOPIC_PER_GROUP can find them
    cp cell_metadata_for_pycistopic.safe.tsv phase1_out/
    cp fragments_map.tsv                     phase1_out/
    cp "\$BLACKLIST"                         phase1_out/blacklist.bed

    echo "=== group_list.tsv ==="
    head -11 phase1_out/group_list.tsv || true
    N_GROUPS=\$(tail -n +2 phase1_out/group_list.tsv | wc -l)
    echo "[PYCISTOPIC_ATAC_PREPARE] \${N_GROUPS} CT×condition groups will be processed in Phase 2"
    if [[ "\$N_GROUPS" -eq 0 ]]; then
        echo "ERROR: group_list.tsv is empty — no groups passed min_cells=${min_cells}"
        exit 1
    fi
    """
}
