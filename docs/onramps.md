# On-ramps & resuming

There are two different ways to avoid recomputing work, and they solve different
problems.

| Mechanism | Use it when | Scope |
|---|---|---|
| **`-resume`** | The same run failed, was interrupted, or you changed something downstream | Nextflow's task cache, within one project |
| **On-ramps** | You have intermediate files from *elsewhere* and want to start mid-pipeline | Explicit file injection via `params.onramp` |

---

## `-resume`

Always pass it:

```bash
nextflow run main.nf -profile cluster,singularity -c my_study.config -resume
```

Nextflow hashes each task's inputs, code, and container. Unchanged tasks are
restored from cache; only what actually changed re-runs. A run that dies at hour
six of eight continues from hour six.

### What silently busts the cache

Understanding this saves a lot of wasted compute:

- **Touching an upstream file's mtime.** Reading an h5ad in write mode (`r+`) is
  enough — it updates mtime, changes the hash, and re-runs everything downstream.
  Open intermediates read-only.
- **Adding a `withName:` block to a resource-tier config.** This shifts the
  configuration hash of processes below it. For one-off adjustments, prefer an
  inline directive in the module.
- **Editing a `bin/` script**, even in a comment — the script is a hashed input.

!!! warning "`-preview` and `-resume` interact badly"
    A `-preview` invocation is written to Nextflow's run history. A subsequent
    bare `-resume` can select that empty preview session and report nothing
    cached. Run previews from a separate directory, or resume an explicit session:

    ```bash
    nextflow log                        # list sessions
    nextflow run main.nf -resume <session-id> ...
    ```

---

## On-ramps

An on-ramp injects a pre-computed artifact so FORGE skips the stage that would
have produced it. This is how you reuse an expensive result across projects, or
demonstrate a downstream stage without paying for its upstream.

```groovy
params {
    onramp {
        rna_integrated_h5ad   = '/prior/results/integration/rna_integrated.h5ad'
        atac_peak_matrix_h5ad = '/prior/results/atac/final/peak_matrix.h5ad'
    }
}
```

### Available on-ramps

| Key | Skips | Notes |
|---|---|---|
| `rna_integrated_h5ad` | The whole RNA arm through integration | See the pairing rule below |
| `seurat_rds` | The whole RNA arm through integration | **BETA.** A Seurat `.rds` instead of an h5ad. Mutually exclusive with `rna_integrated_h5ad`; see below |
| `atac_peak_matrix_h5ad` | ATAC QC, peak calling, clustering | See the side-key warnings below |
| `mudata_h5mu` | Multiome integration | |
| `printer_h5ad` | `SCPRINTER_BUILD_PRINTER` | |
| `cicero_connections` + `cicero_ccan` + `cicero_cds` | Cicero | **All-or-none triple** |
| `chromvar_deviations` + `chromvar_raw` | ChromVAR | **All-or-none pair** |
| `atac_individual_samples_dir` | — | Side key: per-sample h5ads, needed by scPRINTER |
| `atac_anndataset` | — | Side key: SnapATAC2 binding object, needed by enhancer footprinting |
| `rna_per_sample_h5ads_dir` | — | Side key: per-sample RNA h5ads, needed by multiome integration |
| `rna_cellchat_csv` | — | Side key: CellChat input for footprinting recipes |

### `seurat_rds` — the Seurat on-ramp (BETA)

Hand FORGE an existing Seurat object and it converts to the integrated-RNA h5ad
before anything downstream runs. The converted object *is*
`rna_integrated_h5ad` from that point on, so every rule in this document that
applies to one applies to the other — including the multiome pairing rule.
Setting both is a hard error, because they are two routes to the same artifact
and there is no defensible winner:

```text
params.onramp.seurat_rds and params.onramp.rna_integrated_h5ad both set. They
are two routes to the same artifact (the integrated RNA h5ad) -- set exactly
one.
```

```groovy
params {
    onramp {
        seurat_rds = '/data/my_study/integrated_seurat.rds'
    }
    onramp_seurat {
        cell_type_key = 'subclass'      // null = auto-detect; required if hdWGCNA/CellChat run
        assay         = 'RNA'
        sample_key    = 'orig.ident'
    }
}
```

| `params.onramp_seurat` key | Default | Meaning |
|---|---|---|
| `cell_type_key` | `null` | `meta.data` column holding cell type labels. `null` auto-detects with the same precedence as `plot_postscanvi.py`: `cell_type` > `scanvi_prediction` > `celltypist_prediction` > `cell_type_prediction` > `cell_type_marker` > `seurat_annotations`. Set it explicitly for anything non-obvious — an Allen-style object carries `class`/`subclass`/`cluster` and only you know which level the analysis is about. **Required** (auto-detect refused) when `hdwgcna.run` or `cellchat.run` is true |
| `assay` | `'RNA'` | Assay to export. Must be gene x cell; a Signac `ChromatinAssay` is rejected with a pointer to `atac_peak_matrix_h5ad` |
| `sample_key` | `'orig.ident'` | `meta.data` column mapped onto `obs['sample']`. Falls back to a single synthetic sample when absent |

Objects older than Seurat 5 are run through `UpdateSeuratObject()` first, and v5
split layers (`counts.1`, `counts.2`, ...) are `JoinLayers()`'d — without that
you would silently export only the first sample.

The converted h5ad is published to `<outdir>/converted/onramp_rna_integrated.h5ad`
when `params.save_converted` is on, so you can inspect exactly what the rest of
the pipeline saw.

**What BETA means here.** The conversion is verified bit-exact for counts,
per-cell QC, per-gene totals, cell-type composition and cell/gene ordering;
embeddings agree to ~1e-14 (float64 round-trip, not data movement). Beyond the
matrix, the analyses that consume it are gated too: running hdWGCNA natively
and through a full `.rds -> h5ad -> Seurat` round trip gives **identical**
soft-power selection, module counts, module assignments, kME, module
eigengenes and hub genes across 12 cell types on one reference and 3 on
another, and CellChat agrees on all 14 scored quantities including every
inferred ligand-receptor interaction and its probability. Two runs of the same
arm on the same input also agree exactly, so "identical" here means identical,
not "within tolerance".

What is still *not* gated: MAST differential expression goes through the same
converter but has no equivalence test of its own, and neither do the multiome
and GRN consumers (MuData assembly, DORC, SCENIC+, pycisTopic). The on-ramp is
also unreleased. Check those paths against a native run before publishing off
them.

### hdWGCNA and CellChat on an on-ramped object

These used to live inside the RNA subworkflow, which **every** RNA on-ramp
skips — so they silently never ran on an on-ramped object, `rna_integrated_h5ad`
included. As of 2026-10-27 both are hoisted into `RNA_NETWORKS`, which hangs off
the integrated RNA h5ad regardless of what produced it. They now run normally.

Because the cell-type labels drive network construction, the column stops being
a plumbing detail and becomes a scientific choice — it *is* the level your
co-expression networks get built at. Auto-detection is therefore **not allowed**
when either analysis is enabled:

```text
params.onramp.seurat_rds with hdWGCNA/CellChat enabled requires
params.onramp_seurat.cell_type_key to be set explicitly — auto-detection is
not allowed when the labels drive network construction. Name the meta.data
column holding the annotation level you want networks built at
(e.g. 'subclass', not 'class').
```

Set `cell_type_key` and it proceeds. Leave networks disabled and auto-detection
still applies as before.

Every other `meta.data` column is carried into `obs` verbatim, so a
`condition_key` you configured for the hdWGCNA or CellChat differential tiers
either survived the conversion or was never in the object. The converter checks
those columns up front and dies in seconds rather than failing inside
`HDWGCNA_DIFFERENTIAL` after the whole network build:

```text
ERROR: required meta.data column(s) absent from the Seurat object: diagnosis.
These are configured as condition/group keys for hdWGCNA or CellChat, which
run on the on-ramped object.
```

### What X holds, and why your object should already be normalised

The on-ramp and the native RNA path write the integrated h5ad with **opposite
conventions**, which matters if you ever inspect one by hand:

| | `X` | `layers/` |
|---|---|---|
| native path (`plot_postscanvi.py`) | log-normalised | `counts` = raw counts |
| `seurat_rds` on-ramp | **raw counts** | `logcounts` = log-normalised |

`X` is raw on the on-ramp deliberately: DORC and SCENIC+ model count data, and
promoting normalised values into `X` would corrupt them silently. Note also
that a scanpy `/raw/X` group is *not* raw counts — it is `adata.raw`, a
pre-HVG-subset snapshot taken **after** `log1p`. The actual counts are in
`layers/counts`.

Consumers that need normalised values read `layers/logcounts` and fall back to
`X` when it is absent, so both conventions land in the same place. If you are
writing a new consumer, do the same — do not assume `X` is normalised.

**Normalise before you on-ramp.** If your object's `data` layer is identical to
its `counts` layer — i.e. you never ran `NormalizeData()` — the converter has no
normalised values to carry across and writes no `logcounts` layer. FORGE then
normalises internally (CPM + `log1p`) on both return paths — the Seurat
conversion that feeds hdWGCNA and MAST, and the readers that feed CellChat —
so nothing downstream is handed raw counts. That backstop is bit-identical to
`NormalizeData(normalization.method = "LogNormalize", scale.factor = 1e4)`,
but it is still a backstop rather than a substitute: it applies FORGE's
normalisation rather than yours, and it cannot reproduce a custom approach such
as SCTransform. Run
`NormalizeData()` (or your own equivalent) on the object first, and the values
you chose are what the analyses use.

### Bundles are all-or-none

Some artifacts are only meaningful together, and the pre-flight checklist enforces
that rather than letting you discover it mid-run:

```text
Cicero onramp is an all-or-none triple: cicero_connections + cicero_ccan +
cicero_cds. Got 2/3 set.
```

```text
ChromVAR onramp is a pair: chromvar_deviations + chromvar_raw. Got 1/2 set.
```

### Side keys travel with their main artifact

A peak matrix on its own is not enough for every downstream consumer, because
some of them need per-sample objects that the peak matrix does not contain. You
get warned:

```text
WARN: atac_peak_matrix_h5ad onramp set without atac_individual_samples_dir —
scPRINTer cannot build per-sample fragment binding.
```

```text
WARN: atac_peak_matrix_h5ad onramp set without atac_anndataset —
ENHANCER_FOOTPRINTING_RECIPES will be unable to bind snapatac2 anndataset.
```

One combination is a hard error rather than a warning, because there is no valid
way to proceed:

```text
rna_integrated_h5ad onramp + run_multiome_integration=true requires either
rna_per_sample_h5ads_dir (for fresh multiome) or mudata_h5mu (skip multiome).
Set one, or set run_multiome_integration=false.
```

### Forward-declared keys are rejected

`params.onramp` declares several keys that no consumer reads yet. Setting one is
an error, not a silent no-op — it would otherwise look like you had skipped a
stage that in fact still ran:

```text
params.onramp.cistopic_obj_pkl is set but is forward-declared only — no consumer
exists in any FORGE main.nf. Unset this key; the producer step will run
regardless.
```

Currently forward-declared only: `cistopic_obj_pkl`, `da_peaks_dir`, and the
four stratified Cicero keys (`cicero_connections_ctrl/_trt`,
`cicero_ccan_ctrl/_trt`).

`seurat_rds` used to be on this list and no longer is — it became a real
on-ramp, documented above.

---

## Re-rendering figures without recomputing

For plots alone, `VIZ_ONLY` is cheaper than an on-ramped full run:

```bash
nextflow run main.nf -entry VIZ_ONLY -c my_study.config
```

```groovy
params.viz_only {
    peak_matrix_h5ad   = 'results/atac/final/peak_matrix.h5ad'
    cicero_connections = 'results/cicero/connections.csv'
    cicero_ccan        = 'results/cicero/CCAN_assignments.tsv.gz'
    cicero_cds         = 'results/cicero/input_cds_ordered.rds'
    target_genes       = 'APOE,TREM2,CST7'
}
```

The three Cicero keys are required together to render Cicero target plots.

---

## Validate on-ramps before running

On-ramp mistakes are exactly what the pre-flight checklist is best at — every
rule above is checked in about ten seconds:

```bash
nextflow run main.nf -preview -c my_study.config
```

See [Verifying FORGE works](verification.md).
