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
| `atac_peak_matrix_h5ad` | ATAC QC, peak calling, clustering | See the side-key warnings below |
| `mudata_h5mu` | Multiome integration | |
| `printer_h5ad` | `SCPRINTER_BUILD_PRINTER` | |
| `cicero_connections` + `cicero_ccan` + `cicero_cds` | Cicero | **All-or-none triple** |
| `chromvar_deviations` + `chromvar_raw` | ChromVAR | **All-or-none pair** |
| `atac_individual_samples_dir` | — | Side key: per-sample h5ads, needed by scPRINTER |
| `atac_anndataset` | — | Side key: SnapATAC2 binding object, needed by enhancer footprinting |
| `rna_per_sample_h5ads_dir` | — | Side key: per-sample RNA h5ads, needed by multiome integration |
| `rna_cellchat_csv` | — | Side key: CellChat input for footprinting recipes |

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
params.onramp.seurat_rds is set but is forward-declared only — no consumer
exists in any FORGE main.nf. Unset this key; the producer step will run
regardless.
```

Currently forward-declared only: `cistopic_obj_pkl`, `seurat_rds`,
`da_peaks_dir`, and the four stratified Cicero keys
(`cicero_connections_ctrl/_trt`, `cicero_ccan_ctrl/_trt`).

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
