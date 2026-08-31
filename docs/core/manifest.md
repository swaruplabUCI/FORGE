# The manifest CSV

The manifest is the **only** file that describes your data to FORGE. Everything
else — which samples exist, which files belong to them, which experimental group
each belongs to, how many parallel branches the pipeline builds — is derived from
it. Get this file right and the architecture does the rest.

FORGE points at it with a single parameter:

```groovy
params.metadata_file = '/path/to/my_manifest.csv'
```

---

## Minimal working example

This is a complete, valid manifest for a single-sample human 10x multiome run:

```csv
sample_id,batch,sample_type,original_lane_id,rna_file,fragment_file,condition_group,data_dir
10k_PBMC,PBMC,lane,L1,10k_PBMC_raw_feature_bc_matrix.h5,10k_PBMC_atac_fragments.tsv.gz,ConditionA,/data/pbmc
```

And this is a 4 sample, two-condition mouse experiment. Note these columns for critical highlighted changes:
- sample_id
- condition

<pre><code>sample_id,batch,sample_type,original_lane_id,rna_file,fragment_file,condition_group,data_dir
<mark style="background-color: #d1fae5; color: #065f46;">AD_17p9_rep4</mark>,AD_multiome,lane,L1,AD_17p9_rep4_raw_feature_bc_matrix.h5,AD_17p9_rep4_atac_fragments.tsv.gz,<mark style="background-color: #fef3c7; color: #92400e; font-weight: bold;">TG</mark>,/data/ad
<mark style="background-color: #d1fae5; color: #065f46;">AD_2p5_rep2</mark>,AD_multiome,lane,L1,AD_2p5_rep2_raw_feature_bc_matrix.h5,AD_2p5_rep2_atac_fragments.tsv.gz,<mark style="background-color: #fef3c7; color: #92400e; font-weight: bold;">TG</mark>,/data/ad
<mark style="background-color: #d1fae5; color: #065f46;">WT_13p4_rep2</mark>,AD_multiome,lane,L1,WT_13p4_rep2_raw_feature_bc_matrix.h5,WT_13p4_rep2_atac_fragments.tsv.gz,<mark style="background-color: #dbeafe; color: #1e40af; font-weight: bold;">WT</mark>,/data/ad
<mark style="background-color: #d1fae5; color: #065f46;">WT_2p5_rep2</mark>,AD_multiome,lane,L1,WT_2p5_rep2_raw_feature_bc_matrix.h5,WT_2p5_rep2_atac_fragments.tsv.gz,<mark style="background-color: #dbeafe; color: #1e40af; font-weight: bold;">WT</mark>,/data/ad</code></pre>

!!! tip "Scaling is a manifest edit, not a code edit"
    FORGE has been run on 1 sample and on 11 samples with the same `main.nf`.
    Adding samples means adding rows. Per-sample QC, CellBender, and annotation
    fan out automatically from the row count.

---

## Column reference

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | **yes** | Unique identifier for the sample. Becomes the key that joins RNA and ATAC across the whole pipeline, and the prefix on nearly every output file. Rows with a blank `sample_id` are silently skipped. Duplicates are a pre-flight error. |
| `sample_type` | **yes** | Always `lane`, lowercase. Vestigial — see [below](#sample_type). |
| `batch` | **yes** in practice | Batch/group label. Used for batch correction, and as the lookup key into `params.batch_dirs` when `data_dir` is not given. |
| `rna_file` | for RNA runs | Filename (not a path) of the RNA count matrix, resolved relative to the row's data directory. 10x: `*_raw_feature_bc_matrix.h5`. BD Rhapsody: `*_RSEC_MolsPerCell_MEX.zip`. May also name a MEX **directory**. |
| `fragment_file` | for ATAC runs | Filename of the ATAC fragments file. 10x: `*_atac_fragments.tsv.gz`. BD: `*_ATAC_Fragments.bed.gz`. If you give a value with no extension, FORGE appends `.bed.gz`. |
| `condition_group` | for differential | Experimental group (e.g. `WT` / `TG`). This is the axis every differential comparison is built on. Single-condition datasets still need a value — use one label for all rows. Empty values warn and default to `Control`. |
| `data_dir` | see below | Absolute directory holding this row's files. If present it wins outright; if absent, the directory is looked up from `params.batch_dirs[batch]`. |
| `original_lane_id` | optional | Sequencing lane (e.g. `L1`). Only consulted when the batch is listed in `params.batch_dirs_use_lane_subdir`, in which case it becomes a subdirectory under the batch directory. |
| `coord_data_dir` | optional | Directory of coordinate-sorted ATAC fragments, when they live apart from the barcode-sorted ones. Falls back to `params.atac_coord_batch_dirs[batch]`. |

### `sample_type`

Set this to exactly **`lane`**, lowercase. The value is validated
case-sensitively, so capitalized variants fail the pre-flight check.

!!! warning "This column is vestigial, and `lane` is not the best name for it"
    `lane` does not mean a sequencing lane. The name is a holdover from an earlier 
    workflow and will eventually be removed in a pending update.

    **Write `lane` today.** It is the only value the pre-flight check accepts, and
    a manifest using anything else will be rejected before any work starts.

### Paths: `data_dir` vs `batch_dirs`

FORGE gives you two ways to say where files live:

=== "Per-row (simple)"

    Put an absolute full file path `data_dir` on every row. Best for one or two directories.

    ```csv
    sample_id,batch,sample_type,rna_file,data_dir
    S1,b1,lane,S1_raw_feature_bc_matrix.h5,/data/run1
    S2,b1,lane,S2_raw_feature_bc_matrix.h5,/data/run1
    ```

=== "Per-batch (scales better)"

    Omit `data_dir` and map batches to directories in your config. Best when
    many samples share a few locations.

    ```groovy
    params.batch_dirs = [
        june: '/data/june_run',
        july: '/data/july_run',
    ]
    ```

    ```groovy
    // Optional: append original_lane_id as a subdirectory for these batches
    params.batch_dirs_use_lane_subdir = ['july']
    // → /data/july_run/L1/<rna_file>
    ```

If a row has neither `data_dir` nor a matching `batch_dirs` entry, the run stops
with `No directory configured for batch '<batch>'`.

---

## What FORGE checks before it runs

FORGE validates the manifest in its **pre-flight checklist**, before a single
compute task is submitted. A typo in a manifest will get flagged in seconds,
not crash several hours into a run.

The checklist covers:

- **Required columns present** — `sample_id`, `sample_type`, and `rna_file`
- **`sample_id` uniqueness**
- **`rna_file` actually exists** on disk, per row
- **MEX directories are complete** — if `rna_file` names a directory, it must
  contain `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz`
- **`.h5` files are not truncated** — a size probe warns on suspiciously small files
- **`condition_group` is present** whenever a condition-aware workflow is enabled
- **Config agrees with the manifest** — e.g. `rna.run = true` against a manifest
  with no `rna_file` values is an error

A failure looks like this, with every problem reported at once rather than one
per re-run:

```text
================================================================================
PRE-FLIGHT CHECKLIST FAILED (3 error(s)):
================================================================================
  1. Manifest CSV not found: /path/to/datasets/pbmc_10x/10k_pbmc_manifest.csv
  2. atac.annotation_method='scatanno' requires params.scatanno.reference_atlas
     (path to a .h5ad reference). Set this explicitly in your dataset config.
  3. rna.run=true but the manifest contains no rows with a non-null rna_file.
     Either populate rna_file in the manifest, or set rna.run=false for
     ATAC-only runs.
================================================================================
```

### Misspelled column names are flagged, not accepted

If a required column is missing, FORGE looks for near-misses in your header and names them in the error:

```text
Manifest CSV missing required columns: [sample_id (did you mean: sample_ID?)].
Found: [sample_ID, batch, sample_type, rna_file]
```

This is **diagnostic only**. FORGE does not remap `sample_ID` to `sample_id`. Column matching is exact, the run will fail.

!!! note "`fragment_file` paths are not pre-flight checked"
    Per-row file-existence validation covers `rna_file` only. A missing or
    misnamed `fragment_file` is not caught by the checklist; **but it surfaces later** as ATAC channels are built:

    ```text
    Manifest validation failed: ATAC fragment file not found -> /data/pbmc/...
    ```

    That still happens before heavy compute, but it is not part of the initial
    all-errors-at-once report.

!!! tip "Run the checklist on its own"
    You can exercise all of this with no containers, no references, and no GPU:

    ```bash
    nextflow run main.nf -preview -c my_dataset.config
    ```

    It completes in seconds. See [Verifying FORGE works](../verification.md).

---

## Common mistakes

!!! failure "Full paths in `rna_file` / `fragment_file`"
    These columns take **filenames**, not paths. The directory comes from
    `data_dir` or `batch_dirs`. A full path here produces a doubled path.

!!! failure "Omitting `condition_group` on a single-condition dataset"
    Several downstream stages key off this column even when there is nothing to
    contrast. Give every row the same label (the PBMC example uses
    `ConditionA`) rather than leaving it blank.

!!! failure "`sample_id` values that differ between RNA and ATAC"
    `sample_id` is the join key across modalities. If the RNA and ATAC halves of
    one biological sample carry different IDs, they will be treated as two
    samples and never integrated.

!!! failure "Assuming barcodes match across layers"
    Raw barcodes differ in format between 10x and BD, and between RNA and ATAC
    (`sample:barcode` vs `barcode-sample`). FORGE normalizes these internally. 
    Do not adjust barcodes to try to help it.

---

## Where to go next

- [nextflow.config](config.md) — the parameters that decide what runs on this manifest
- [main.nf architecture](architecture.md) — how manifest rows become parallel channels
- [Quickstart](../quickstart.md) — a manifest you can actually run today
