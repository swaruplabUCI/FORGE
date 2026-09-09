# Training a custom CellTypist model

Used by: `RUN_CELLTYPIST`, `ATAC_CELLTYPIST` — for any tissue the
[CellTypist model zoo](https://www.celltypist.org/models) does not already cover.

!!! info "This is not part of the FORGE pipeline"
    Nothing on this page runs inside Nextflow. A CellTypist model is an input to FORGE,
    built once, out of band, and then pointed at from your dataset config. The scripts
    below are provided so you can build one for your own tissue; they are reproduced in
    full rather than shipped in `bin/`, precisely because they are not pipeline code.

The authoritative instructions are upstream, in the CellTypist repository's
[**Supplemental guidance: generate a custom model**](https://github.com/Teichlab/celltypist#supplemental-guidance-generate-a-custom-model)
section. This page does not replace it — it records the one recipe FORGE has actually
used, with the parameter choices and the failure mode that cost us a job.

| Upstream resource | Link |
|---|---|
| CellTypist source and README | [github.com/Teichlab/celltypist](https://github.com/Teichlab/celltypist) |
| Custom-model guidance | [Supplemental guidance: generate a custom model](https://github.com/Teichlab/celltypist#supplemental-guidance-generate-a-custom-model) |
| `celltypist.train` API reference | [celltypist.readthedocs.io](https://celltypist.readthedocs.io/en/latest/celltypist.train.html) |
| Pre-trained model zoo | [celltypist.org/models](https://www.celltypist.org/models) |
| Citation | Domínguez Conde *et al.*, *Science* 376:eabl5197 (2022) — [doi:10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197) |

## When you need one

Check the model zoo first. Build your own when either is true:

- **No zoo model covers your tissue.** Zoo coverage is concentrated in immune, brain, and
  lung; many other tissues have no entry. Check the current list before assuming either
  way — models are added over time.
- **The zoo model's label resolution is wrong for your question.** A model that returns
  "Epithelial cell" is useless if your analysis turns on distinguishing PTS1 from PTS3.

A custom model needs an annotated **reference** dataset — someone else's atlas, with cell
type labels you trust. You are transferring their labels onto your cells, so the model is
only ever as good as that atlas.

---

## Prerequisites

- An annotated reference `.h5ad` with **raw integer counts** available (in `X` or `.raw`)
  and a label column in `.obs`
- A container or environment with `celltypist`, `scanpy`, and `anndata`. In FORGE this is
  `scgpu_extended.sif`, which already carries CellTypist.
- Memory proportional to the atlas: ~128 GB was comfortable for 112k cells × 16.5k genes.
  Training itself is short; **loading and normalizing the atlas dominates**.

Set a working directory once:

```bash
export CT=/path/to/ref/my_tissue_atlas
export SIF=/path/to/singularity_cache/scgpu_extended.sif
mkdir -p "$CT"
```

---

## Step 1 — download the reference atlas

Any annotated atlas will do. [CellxGene Discover](https://cellxgene.cziscience.com/) is the
most convenient source because its collections are already `.h5ad`, already have a stable
per-dataset download URL, and carry a curated `.obs` schema — including `is_primary_data`,
which matters in step 2.

```bash title="download.sh"
#!/bin/bash
#SBATCH --job-name=dl_atlas
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -euo pipefail
cd "${CT:?set CT to your atlas directory}"

echo "Starting download at $(date)"
# CellxGene Discover: copy the dataset's .h5ad URL from the collection page.
wget -c "https://datasets.cellxgene.cziscience.com/<DATASET-UUID>.h5ad" \
     -O atlas.h5ad \
     --progress=dot:giga
echo "Download complete at $(date)"
ls -lh atlas.h5ad
```

Record the **dataset UUID and the download date**. CellxGene collections are versioned and
the per-dataset URL points at a specific version; you cannot reconstruct which one you used
from the `.h5ad` alone.

## Step 2 — inspect before you train

Do this interactively, before writing the training script. Three things decide every
parameter downstream:

```bash
singularity exec --bind /path/to/data "$SIF" python3 - <<'PY'
import anndata
a = anndata.read_h5ad("atlas.h5ad", backed="r")
print(a)
print("\nlabel columns:")
for c in a.obs.columns:
    if a.obs[c].dtype.name in ("category", "object"):
        print(f"  {c:<30s} {a.obs[c].nunique()} values")
print("\nis_primary_data:", a.obs.get("is_primary_data", "ABSENT"))
PY
```

**Which label column?** Curated atlases usually carry both an ontology-normalized column
(`cell_type`, coarse, safe) and the authors' own (`author_cell_type`, finer, closer to the
biology you probably care about). Prefer the author labels when the fine distinctions are
the point, and accept that some will be rare.

**Is the atlas a meta-atlas?** Many CellxGene entries pool several studies and include cells
that also appear in another dataset. `is_primary_data == False` marks those. Training on
them double-counts donors and inflates the apparent size of your reference.

**Where are the raw counts?** This is the one that bites — see step 3.

## Step 3 — the counts trap

CellTypist requires expression **log1p-normalized to 10,000 counts per cell**, and enforces
it with `check_expression=True`. Curated atlases very often store *some* normalization in
`X` that is not to 10k, while keeping integer counts in `.raw`. Handing over `X` fails:

```text
🍳 Preparing data before training
ValueError: 🛑 Invalid expression matrix, expect log1p normalized expression
to 10000 counts per cell
```

That is CellTypist working correctly. The fix is to rebuild the AnnData from `.raw` and
normalize it yourself — not to set `check_expression=False`, which silently trains a model
on the wrong scale and produces a `.pkl` that misbehaves only at prediction time.

CellxGene `.h5ad` files also index `var` by ENSEMBL ID while carrying symbols in
`var['feature_name']`. FORGE's query data is symbol-indexed, so the swap is mandatory —
without it the model and your query share zero features and every prediction is garbage.

## Step 4 — train

```python title="train_celltypist.py"
#!/usr/bin/env python3
"""
Train a custom CellTypist model from an annotated reference atlas.

Strategy:
- Primary data only (drops cells duplicated from other studies in a meta-atlas)
- author_cell_type labels; drop 'Unknown'
- Use a.raw (integer counts); swap ENSEMBL var.index -> gene symbols (feature_name)
- Normalize raw counts to 10k + log1p (X is often pre-normalized, but not to 10k)
- SGD + mini-batch + balance_cell_type for >100k-cell scale
- Two-pass feature_selection (top_genes=300 per class)
- with_mean=False keeps the sparse matrix sparse during scaling
"""
import warnings; warnings.filterwarnings('ignore')
import anndata, scanpy as sc, celltypist
from pathlib import Path
import os

ATLAS   = Path(os.environ["CT"]) / "atlas.h5ad"
OUTDIR  = Path(os.environ["CT"])
LABELS  = "author_cell_type"        # <- the .obs column chosen in step 2
N_JOBS  = int(os.environ.get("SLURM_CPUS_PER_TASK", 16))

print("=" * 60)
print("CellTypist custom model training")
print("=" * 60)

# -- 1. Load ----------------------------------------------------------------
print("\n[1] Loading atlas...", flush=True)
a = anndata.read_h5ad(ATLAS)
print(f"    Full: {a.n_obs} cells x {a.n_vars} genes")

# -- 2. Primary data only ---------------------------------------------------
# Meta-atlases include cells re-used from constituent studies; training on them
# double-counts donors. Drop this block if your atlas has no such column.
if "is_primary_data" in a.obs:
    a = a[a.obs["is_primary_data"]].copy()
    print(f"    Primary: {a.n_obs} cells")

# -- 3. Drop Unknown --------------------------------------------------------
a = a[a.obs[LABELS] != "Unknown"].copy()
print(f"    After dropping Unknown: {a.n_obs} cells")
print("\n    Cell type distribution:")
for ct, n in a.obs[LABELS].value_counts().items():
    print(f"      {ct:<35s} {n:>6d}")

# -- 4. Build normalized AnnData from raw counts ----------------------------
# See step 3: X may be pre-normalized but not to 10k; a.raw holds the integers.
print("\n[2] Building raw-count AnnData from a.raw...", flush=True)
a_raw = a.raw.to_adata()
a_raw.obs = a.obs.copy()

# Swap ENSEMBL var.index -> gene symbols, so the model's features match query data
a_raw.var_names = a_raw.var["feature_name"].values
a_raw.var_names_make_unique()
print(f"    Example var names: {list(a_raw.var_names[:5])}")

print("\n[3] Normalizing raw counts to 10k + log1p...", flush=True)
sc.pp.normalize_total(a_raw, target_sum=1e4)
sc.pp.log1p(a_raw)
x_check = a_raw.X[:50, :].toarray().max() if hasattr(a_raw.X, "toarray") else a_raw.X[:50, :].max()
print(f"    X max after norm (first 50 cells): {x_check:.2f}  (expect <=9.21)")

# -- 5. Train ---------------------------------------------------------------
print("\n[4] Training CellTypist model...", flush=True)
print(f"    SGD mini-batch | feature_selection=True | n_jobs={N_JOBS} | balance_cell_type=True")

model = celltypist.train(
    X                 = a_raw,
    labels            = LABELS,        # obs column name
    with_mean         = False,         # keep sparse; saves several GB of RAM
    check_expression  = True,          # leave ON -- see step 3
    use_SGD           = True,
    mini_batch        = True,
    batch_number      = 100,
    batch_size        = 1000,
    epochs            = 10,
    balance_cell_type = True,
    feature_selection = True,
    top_genes         = 300,
    n_jobs            = N_JOBS,
    details           = "Describe the atlas, its constituent studies, and the label set.",
    source            = "CellxGene <DATASET-UUID> / <citation>",
    version           = "1.0",
)
print("    Training complete.")

# -- 6. Save ----------------------------------------------------------------
out_pkl = OUTDIR / "atlas_celltypist.pkl"
model.write(str(out_pkl))
print(f"\n[5] Saved model: {out_pkl}")
print(f"    Classes ({len(model.cell_types)}): {sorted(model.cell_types)}")
```

```bash title="train_celltypist.sh"
#!/bin/bash
#SBATCH --job-name=train_celltypist
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -euo pipefail

echo "Job started: $(date)"
echo "Node: $(hostname)"

singularity exec --bind /dfs7 \
  "${SIF:?set SIF to scgpu_extended.sif}" \
  python3 "${CT:?set CT to your atlas directory}/train_celltypist.py"

echo "Job finished: $(date)"
```

### Why these parameters

| Parameter | Value | Reason |
|---|---|---|
| `use_SGD`, `mini_batch` | `True` | The default logistic-regression solver does not scale past ~50k cells in reasonable memory. Above that, SGD with mini-batches is the documented path. |
| `batch_number` / `batch_size` / `epochs` | 100 / 1000 / 10 | Upstream's recommended starting point for >100k cells: 100 × 1000 = 100k cells seen per epoch. |
| `balance_cell_type` | `True` | Rare classes are otherwise never sampled into a mini-batch. Essential when class sizes span three orders of magnitude. |
| `feature_selection` | `True` | Two-pass: train, rank features per class, retrain on the top ones. Costs a second pass, materially improves the model. |
| `top_genes` | `300` | Per class, not total. Upstream default. |
| `with_mean` | `False` | Centering densifies the matrix. Keeping it sparse is the difference between 128 GB working and not. |
| `check_expression` | `True` | Never turn this off (step 3). |
| `details` / `source` / `version` | set them | These are written into the `.pkl` and are the only provenance that travels with the model. |

The walltime above is generous on purpose. Actual training was ~3.5 minutes; the request
covers atlas load and normalization on a slow filesystem.

## Step 5 — verify

```bash
singularity exec --bind /dfs7 "$SIF" python3 - <<'PY'
import celltypist, os
m = celltypist.models.Model.load(os.environ["CT"] + "/atlas_celltypist.pkl")
print("classes:", len(m.cell_types))
print(sorted(m.cell_types))
print("features:", len(m.features))
print("example features:", list(m.features[:10]))
PY
```

Check three things:

1. **Feature names are symbols, not ENSEMBL IDs.** If you see `ENSMUSG…`, step 3's `var`
   swap did not happen and the model will match nothing in your query.
2. **The class list is the one you expect**, and `Unknown` is absent.
3. **Class count matches labels minus dropped.** A silently missing class usually means it
   had too few cells to survive feature selection.

Then run it on a small query and look at the distribution — a model that assigns 95% of
cells to one label is telling you the reference and query are not comparable.

---

## Wiring it into FORGE

Point the dataset config at the `.pkl` by absolute path:

```groovy
params {
    celltypist {
        enabled = true
        model   = "/path/to/ref/my_tissue_atlas/atlas_celltypist.pkl"
    }
}
```

Any `params.celltypist.model` value that is not a bare zoo filename is passed straight to
`celltypist.models.Model.load()` as a path, so no other change is needed.

!!! warning "Custom models are not in the container"
    Zoo models are baked into `scgpu_extended.sif`; a custom `.pkl` is not. It lives on
    your filesystem, must be inside a bind-mounted path, and **must be archived alongside
    your results** — it is not reconstructible from the pipeline alone.

---

## Provenance: the model FORGE actually uses

One custom model has been built and used in this project.

| | |
|---|---|
| **Model file** | `kidney_atlas_celltypist.pkl` (1.1 MB, 35 classes) |
| **HPC3 path** | `/dfs7/swaruplab/lesolano/ref/mouse_kidney_snRNA_atlas/` |
| **Status** | `[custom-built]` |
| **Used by** | `Kidney_Mm_BD_r5` — `configs/datasets/kidney_mm_bd.config` |
| **Build date** | June 8, 2026 |
| **Container** | `scgpu_extended.sif` |

**Source atlas.** Mouse kidney snRNA-seq reference atlas, CellxGene Discover dataset
`945ee50d-d14e-4e86-baaa-febc7eb00409`
([PMC10238935](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10238935/)), downloaded
2026-06-08 as `kidney_atlas_PMC10238935.h5ad` (650,446,815 bytes). Pools eight studies:
Hinze20, Kirita20, Park18, Miao21, Wu19, Dumas20, Janosevic21, Conway20.

**Filtering.**

| Stage | Cells |
|---|---|
| Full atlas | 141,401 (× 16,575 genes) |
| `is_primary_data == True` (drops Kirita20 duplicates) | 115,551 |
| After dropping `Unknown` | **111,904** |

**Labels.** `author_cell_type` — 36 values in the filtered atlas, 35 after dropping
`Unknown`, all 35 retained in the trained model. Chosen over the ontology `cell_type`
column because the kidney analysis depends on proximal-tubule segment identity (PTS1 /
PTS2 / PTS3 / PTS3T2) and on vascular subtypes that the ontology column collapses.

**Expression.** Rebuilt from `.raw` integer counts, `var` re-indexed ENSEMBL →
`feature_name`, then `normalize_total(target_sum=1e4)` + `log1p`. Post-normalization max
over the first 50 cells: 8.12 (bound 9.21).

**Training.** `celltypist.train(use_SGD=True, mini_batch=True, batch_number=100,
batch_size=1000, epochs=10, balance_cell_type=True, feature_selection=True, top_genes=300,
with_mean=False, check_expression=True, n_jobs=16)`.

**Run record.** Two SLURM jobs on `standard`, 16 CPU / 128 GB:

| Job | Outcome |
|---|---|
| `52944727` | **Failed** — `ValueError: Invalid expression matrix, expect log1p normalized expression to 10000 counts per cell`. `X` was pre-normalized, but not to 10k. |
| `52944736` | **Succeeded** — 2026-06-08 01:30:42 → 01:34:08 on `hpc3-14-28` (~3.5 min). |

Job `52944727` is the origin of step 3 above, and is recorded here rather than discarded:
the fix changed which matrix the model was trained on, so it is a substantive part of how
this model came to exist.

!!! caution "Known limitation of this model"
    Class sizes span PTS2 at 28,380 cells down to MD at 7. `balance_cell_type=True` keeps
    the rare classes in training, but a 7-cell class does not support a trustworthy
    classifier. Treat predictions for the smallest classes as hypotheses to check against
    marker expression, not as calls. The per-class counts are in the job log at
    `train_celltypist_52944736.log`.

---

## Reporting this in a manuscript

Record: the atlas and its accession or CellxGene dataset UUID, the download date, the label
column used, every filtering step with resulting cell counts, the normalization applied,
the full `celltypist.train` parameter set, the CellTypist version, and the final class
count. Cite the CellTypist paper. The `.pkl` is a research resource in its own right —
deposit it, or cite this page as the build recipe.
