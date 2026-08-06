# Troubleshooting

Failures in FORGE fall into a few recurring shapes. This page is organized by what
you actually see.

**Before anything else**, re-run pre-flight — it costs ten seconds and catches
most configuration faults:

```bash
nextflow run main.nf -preview -c my_study.config
```

---

## Pre-flight checklist failed

Good news: nothing has run yet, and every problem is reported at once.

| Error | Fix |
|---|---|
| `Manifest CSV not found` | Check `params.metadata_file` — a relative path resolves against the launch directory, not the config file. |
| `missing required columns: [sample_id (did you mean: sample_ID?)]` | Fix the header to match **exactly**. The suggestion tells you which column to correct; FORGE does not accept near-misses. |
| `atac.annotation_method='scatanno' requires params.scatanno.reference_atlas` | Set the atlas path, or switch to `annotation_method = 'celltypist'`. |
| `rna.run=true but the manifest contains no rows with a non-null rna_file` | Populate `rna_file`, or set `rna.run = false`. |
| `Duplicate sample_id in lane rows` | `sample_id` must be unique. |
| `missing 'condition_group' column but a condition-aware workflow is enabled` | Add the column, or disable the differential/stratified workflows. |
| `RNA file not found for sample 'X'` | `rna_file` is a **filename**; the directory comes from `data_dir` or `batch_dirs`. A full path here yields a doubled path. |
| `MEX directory for sample 'X' missing files` | A MEX directory needs `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz`. |
| `Cicero onramp is an all-or-none triple` | Supply all of `cicero_connections`, `cicero_ccan`, `cicero_cds` — or none. |
| `params.onramp.X is forward-declared only` | Unset it. That stage runs regardless. |
| `resource_tier` rejected | Must be lowercase `small`, `medium`, `large`, or `auto`. |

---

## Nothing happened, and there was no error

Almost always a gate. Three causes, in order of likelihood:

1. **An inner gate without its outer gate.** `msfp_strip.enabled = true` does
   nothing unless `enhancer_footprinting.msfp_enabled = true` also holds, because
   strips render from footprints that were never computed. Same pattern for
   `promoter_overlay` and the `shi_figures` tiers.
2. **A misspelled parameter.** `-c` merges rather than validates, so an unknown key
   is silently ignored. Confirm what actually resolved:

    ```bash
    nextflow config -c my_study.config | grep -A5 msfp_strip
    ```

3. **A cell type filtered by the resolution floor.** Any per-cell-type fan-out
   skips cell types below `max(qc.cell_type_resolution.min_cells, min_pct × total)`.
   Per-condition Cicero additionally requires 250 cells *per stratum*. Skips are
   logged — search the log for the cell-type name.

---

## Fewer samples processed than the manifest has rows

Check `sample_type`. Rows whose value does not match a recognized branch are
dropped from the channels. Values are validated case-sensitively, so a
capitalized variant fails the checklist — but a genuinely unrecognized string can
pass validation and still match nothing.

---

## Out of memory / exit code 137

Exit 137 is a SLURM or cgroup kill, not an application bug.

```bash
grep -E "MULTIVI|OOM|137" logs/nextflow/trace.txt
```

Common cases:

- **MultiVI** needs far more *host* memory than its GPU footprint suggests. Raise
  host memory, not GPU memory. `run_imputation = true` increases this further.
- **SCENIC+** wants ≥ 256 GB and is the most memory-hungry stage in FORGE.
- **Footprinting on a tiny cell type.** Memory use is driven by region and mode
  counts, not cell count, so a 7-cell cell type can cost as much as a large one.
  Raise `qc.cell_type_resolution.min_cells` rather than raising memory.
- **A GPU partition memory ceiling.** Some partitions cap total host memory per
  job regardless of what you request; a large ask can be silently clipped. Check
  the partition's limits before assuming the request was honored.

Right-size from evidence — `logs/nextflow/report.html` gives peak RSS per process.

---

## Everything re-ran despite `-resume`

Nextflow hashes task inputs, code, and containers. Things that bust the cache
without looking like changes:

- **An upstream file's mtime changed.** Opening an h5ad in write mode (`r+`) is
  enough. Open intermediates read-only.
- **A new `withName:` block in a resource-tier config**, which shifts the
  configuration hash of processes below it. Prefer an inline directive for one-off
  changes.
- **Any edit to a `bin/` script**, including comments.
- **A previous `-preview` in the same directory.** A bare `-resume` can select the
  preview's empty session. Resume an explicit session:

    ```bash
    nextflow log
    nextflow run main.nf -resume <session-id> ...
    ```

Diagnose before rebuilding: enumerate the hashed inputs and compare upstream
mtimes rather than concluding Nextflow is nondeterministic.

---

## The run died with no useful error

If the log truncates mid-flush with a bare exit 1, a blank reason, and no
`Caused by:` block, **check your filesystem quota first**. Quota exhaustion kills
the Nextflow head JVM while it is writing, which produces exactly this signature.

```bash
df -h /path/to/workdir
dd if=/dev/zero of=/path/to/workdir/.probe bs=1M count=10 && rm /path/to/workdir/.probe
```

The `dd` probe is instant and definitive where quota tools are unavailable.

---

## Container and job-submission failures

| Symptom | Cause |
|---|---|
| `FATAL: could not open image` | `params.containers` path is wrong, or the `.sif` is not readable from compute nodes. |
| Container cannot see your data | Add bind mounts: `singularity.runOptions = '--nv -B /data -B /refs'`. |
| CUDA not available inside the container | `--nv` is missing, or the job landed on a non-GPU node. |
| `QOSMaxSubmitJobPerUserLimit` | Your QOS caps concurrent submissions. Add `maxForks` to the process **and** an error-strategy retry — `maxForks` alone does not prevent submit-time rejections. |
| `sbatch: error: ...` with no `.command.log` | A submit-time failure. There is no task working directory, so inspect the scheduler error rather than looking for process logs. |
| `Invalid partition` / `Invalid account` | Site-specific settings. See [Adapting to your cluster](setup/cluster.md). |

---

## Wrong-looking biology

Not a crash, but worth checking before you trust results:

- **CellTypist model mismatch.** `Immune_All_Low.pkl` is the default and is wrong
  for most tissues. A mismatched model gives confident, plausible, wrong labels.
- **Species/genome mismatch.** Pre-flight checks build consistency, but confirm
  `params.species` matches your GTF, blacklist, and motif databases.
- **Differential direction.** A positive `log2FC` means up in the **treatment**
  condition. Verify `control_condition` and `treatment_condition` are not swapped.
- **RNA/ATAC disagreement.** Expected to some degree — they are independent
  estimates by design. Large disagreement usually points at an annotation
  reference, not at the integration.

---

## Getting help

Include the following, which together identify almost any failure:

1. The pre-flight output: `nextflow run main.nf -preview -c my_study.config`
2. Resolved config: `nextflow config -c my_study.config`
3. The failing task's `.command.log` and `.command.err` from its work directory
4. The relevant `logs/nextflow/trace.txt` row

Open an issue at
[github.com/swaruplabUCI/FORGE/issues](https://github.com/swaruplabUCI/FORGE/issues).
