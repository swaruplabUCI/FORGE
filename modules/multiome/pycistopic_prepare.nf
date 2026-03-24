// modules/multiome/pycistopic_prepare.nf
nextflow.enable.dsl=2

process PYCISTOPIC_PREPARE {
    tag "pycistopic_prepare"

    publishDir "${params.outdir}/pycistopic", mode: 'copy', overwrite: true, pattern: "cistopic_obj.pkl"
    publishDir "${params.outdir}/pycistopic", mode: 'copy', overwrite: true, pattern: "region_sets"
    publishDir "${params.outdir}/pycistopic", mode: 'copy', overwrite: true, pattern: "pycistopic_gene_activity.h5ad"
    publishDir "${params.outdir}/pycistopic/bigwigs", mode: 'copy', overwrite: true, pattern: "consensus_peak_calling/pseudobulk_bw_files/*.bw"

    input:
    // Sample-level metadata CSV (sample_id, rna_file, etc.)
    path sample_metadata
    // Integrated RNA h5ad
    path rna_h5ad
    // Species code for pycisTopic TSS (e.g., "hsapiens")
    val species
    // MuData stats JSON from BUILD_MUDATA (mudata_stats.json)
    path mudata_stats
    path blacklist_bed
    // Gencode GTF for local TSS generation (avoids BioMart)
    path gtf_file

    output:
    path "cistopic_obj.pkl", emit: cistopic_obj
    path "region_sets",        emit: region_sets
    path "pycistopic_gene_activity.h5ad", optional: true, emit: gene_activity
    path "consensus_peak_calling/pseudobulk_bw_files/*.bw", optional: true, emit: pseudobulk_bigwigs

    // Conda handled in nextflow.config

    script:
    """
    set -euo pipefail

    echo "=== PYCISTOPIC_PREPARE ==="
    echo "Sample metadata: ${sample_metadata}"
    echo "RNA h5ad:        ${rna_h5ad}"
    echo "Species:         ${species}"
    echo "CONDA_PREFIX:    ${CONDA_PREFIX}"

    # 1) Build cell-level metadata (barcode -> sample_id, cell_type)
    python ${projectDir}/bin/build_cell_metadata_for_pycistopic.py \\
        --rna-h5ad ${rna_h5ad} \\
        --sample-metadata ${sample_metadata} \\
        --out-tsv cell_metadata_for_pycistopic.tsv

    echo "Cell metadata (head):"
    head cell_metadata_for_pycistopic.tsv || true

    # 2) Build fragments_map.tsv from sample_metadata_unified.csv
    #    Uses columns: sample_id, fragment_file, batch
    echo "Building fragments_map.tsv from sample metadata..."

python - <<'PY'
import pandas as pd, re

df = pd.read_csv("cell_metadata_for_pycistopic.tsv", sep="\t")

def safe(s):
    s = str(s)
    s = re.sub(r"\s+", "_", s)          # spaces -> _
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)  # everything else -> _
    s = re.sub(r"_+", "_", s).strip("_")
    return s

df["cell_type_safe"] = df["cell_type"].map(safe)
df.to_csv("cell_metadata_for_pycistopic.safe.tsv", sep="\t", index=False)

print(df[["cell_type","cell_type_safe"]].drop_duplicates().head(20))
PY

python - << EOF
import os
import pandas as pd

sample_metadata_path = "${sample_metadata}"
out_map_path = "fragments_map.tsv"

df = pd.read_csv(sample_metadata_path)

# Restrict to demux-level ATAC samples only
df = df[df["sample_type"] == "demux"].copy()

required_cols = ["sample_id", "fragment_file", "batch"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(
        f"sample_metadata is missing required columns for pycisTopic coord fragments map: {missing}. "
        f"Found columns: {list(df.columns)}"
    )

batch_to_dir = {
    "june": "${params.atac.demux_coord_dir_june}",
    "july": "${params.atac.demux_coord_dir_july}",
    "nov":  "${params.atac.demux_coord_dir_nov}",
}

def is_lane_donor_id(x: str) -> bool:
    parts = x.split("_")
    if len(parts) != 3:
        return False
    lane, donor_word, donor_num = parts
    return (lane.startswith("L") and lane[1:].isdigit()
            and donor_word == "Donor" and donor_num.isdigit())

paths = []
for _, row in df.iterrows():
    sample_id = str(row["sample_id"])
    frag_root = str(row["fragment_file"])
    batch = str(row["batch"])

    sample_id_unique = sample_id

    if batch in batch_to_dir:
        # BD Rhapsody multi-batch: use batch-specific coord-sorted dir
        if batch == "june" and sample_id.startswith("S") and sample_id[1:].isdigit():
            sample_id_unique = f"{sample_id}_june"
        if batch in ("july", "nov") and is_lane_donor_id(sample_id):
            sample_id_unique = f"{sample_id}_{batch}"

        base_dir = batch_to_dir[batch]
        filename = f"{frag_root}_coord_sorted.bed.gz"
        full_path = os.path.join(base_dir, filename)
        if not os.path.exists(full_path):
            raise FileNotFoundError(
                f"Coord-sorted fragments file not found for sample_id '{sample_id_unique}': {full_path}"
            )
    else:
        # Generic (e.g. 10x): fragment_file is the filename, data_dir has the path
        data_dir = str(row.get("data_dir", "."))
        # Try the fragment file directly (may already be a full path)
        if os.path.isabs(frag_root) and os.path.exists(frag_root):
            full_path = frag_root
        else:
            full_path = os.path.join(data_dir, frag_root)
        if not os.path.exists(full_path):
            raise FileNotFoundError(
                f"Fragments file not found for sample_id '{sample_id}': {full_path}"
            )

    paths.append((sample_id_unique, full_path))

map_df = pd.DataFrame(paths, columns=["sample_id", "fragments_path"])
map_df.to_csv(out_map_path, sep="\\t", index=False)

print(f"[PYCISTOPIC_PREPARE] Wrote coord-sorted fragments map for {len(map_df)} samples to {out_map_path}")
print(map_df.head())
EOF
    echo "fragments_map.tsv (head):"
    head fragments_map.tsv || true

    # Make blacklist portable: stage from params, and ensure it's an uncompressed .bed
    BLACKLIST_SRC="${blacklist_bed}"

    if [[ "\$BLACKLIST_SRC" == *.gz ]]; then
      gunzip -c "\$BLACKLIST_SRC" > blacklist.bed
      BLACKLIST="blacklist.bed"
    else
      BLACKLIST="\$BLACKLIST_SRC"
    fi

python - <<'PY'
import json, pandas as pd
cm = pd.read_csv("cell_metadata_for_pycistopic.tsv", sep="\t")
stats = json.load(open("mudata_stats.json"))
mudata_samples = set([d["sample_id"] for d in stats["samples"]])
cm_samples = set(cm["sample_id"].unique())
print("cell_metadata samples:", len(cm_samples))
print("mudata samples:", len(mudata_samples))
print("intersection:", len(cm_samples & mudata_samples))
assert len(cm_samples & mudata_samples) > 0, "NO sample_id overlap between cell_metadata and MuData stats!"
PY

    # 3) Run pycisTopic wrapper with multi-sample fragments map
    echo "MuData stats file: ${mudata_stats}"
    python ${projectDir}/bin/run_pycistopic_prepare.py \
      --fragments-map fragments_map.tsv \
      --cell-metadata cell_metadata_for_pycistopic.safe.tsv \
      --species hsapiens \
      --sample-id-col sample_id \
      --cell-type-col cell_type_safe \
      --variable cell_type_safe \
      --mudata-stats "mudata_stats.json" \
      --blacklist-bed "\$BLACKLIST" \
      --gtf ${gtf_file} \
      --outdir .

    echo "PYCISTOPIC_PREPARE done."
    ls -lh
    """
}
