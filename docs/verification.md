# Verifying FORGE works

FORGE's published results come from four datasets that took roughly **1,000
compute-hours** in total with the goal of demonstrating the pipeline at scale. 

So we've split verification into three tiers. Each answers a different question, at
increasing costs.

| Tier | Question it answers | Time | Needs |
|---|---|---|---|
| **1. Pre-flight + DAG** | Is my configuration coherent? Does the whole graph wire up? | **~15 seconds** | Nextflow only |
| **2. Tiny dataset** | Does a real analysis arm produce real numbers? | **~3h** on 8 CPUs | 78.9 MB data, real containers, CPU only, ~15 GB free disk |
| **3. Full datasets** | Does it scale, and reproduce the published figures? | hours–days | Full references, GPU, HPC |

---

## Tier 1 — pre-flight and DAG construction

This requires **no containers, no reference files, no GPU, and no downloads**. It
runs on a laptop.

### Run it against the shipped fixture

FORGE ships a self-contained fixture — `test_data/` (a few hundred KB of
placeholder files) plus a config that enables every optional block:

```bash
nextflow run main.nf -profile test -preview \
    -c configs/datasets/test_preview.config
```

Expected result:

```text
WARN: Missing container files: [...]. Not required for -preview, but a real
      run needs them in singularity_cache/.

PRE-FLIGHT CHECKLIST PASSED (8 checks):
    [OK] Manifest schema (2 rows)
    [OK] Species/genome consistency
    [OK] GTF files (5 paths validated)
    [OK] scATAnno reference atlas (stub_atlas.h5ad)
    [OK] MOFA mode (high_memory)
    [OK] scPRINTER genome (hg38)
    [OK] CellTypist model: Immune_All_Low.pkl
    [OK] Resource tier (test)
  Warnings: 1 (see above)
```

The container warning is expected on a fresh clone. A
preview builds the process graph but launches no task, so no image is ever
entered. The check converts to an error for a real run. If you have already
built the containers you will see nine checks and no warnings instead.

`-profile test` strips every site-specific scheduler assumption (SLURM
partitions, accounts, QOS, `--gres` GPU strings) so this works on any machine,
and redirects `outdir` to `results_test/` so nothing can touch real results.

### Then run it against your own config

```bash
nextflow run main.nf -preview -c configs/datasets/my_study.config
```

`-preview` builds the entire workflow graph and runs FORGE's pre-flight checklist
without submitting a single task. You get one of two outcomes.

**Configuration problems, reported all at once:**

```text
================================================================================
PRE-FLIGHT CHECKLIST FAILED (3 error(s)):
================================================================================
  1. Manifest CSV not found: /path/to/10k_pbmc_manifest.csv
  2. atac.annotation_method='scatanno' requires params.scatanno.reference_atlas
     (path to a .h5ad reference). Set this explicitly in your dataset config.
  3. rna.run=true but the manifest contains no rows with a non-null rna_file.
     Either populate rna_file in the manifest, or set rna.run=false for
     ATAC-only runs.
================================================================================
```

Every error is actionable and names the parameter to fix. Fix them and re-run.

**Or a clean graph**, meaning your manifest parses, your species and genome
builds agree, every reference your enabled modules need exists, your on-ramp
bundles are complete, and the process graph resolves.

### Tier 1 covers:

- `-preview` builds the entire workflow graph
- It surfaces every construction-time defect:
 a. an unresolved module include
 b. a channel referenced outside the scope where it was declared
 c. a `val` input that evaluates to `null`, or a missing nested config block.
In practice these are specific actionable errors that are quickly fixed with minimal troubleshooting.

### Tier 1 missess:

Tier 1 will miss anything that only triggers once tasks run:

| Missed in Tier 1 | Why | Covered by |
|---|---|---|
| Runtime channel joins | e.g. the per-sample multiome join keys on output *filenames* | Tier 2 |
| Tool behaviour and numerical output | No process body executes | Tier 2 |
| Resource sizing (OOM, walltime) | Run must crash for these reasons | Tier 3 |

So a clean tier-1 result means "this configuration is coherent and the pipeline
will start.

You can also confirm what your layered config actually resolved to, without
launching anything:

```bash
nextflow -c configs/datasets/my_study.config config -profile cluster,singularity
```

!!! tip "This is the check to run after any config edit"

!!! warning "Don't follow `-preview` with a bare `-resume`"
    A `-preview` invocation is recorded in Nextflow's run history. A subsequent
    bare `nextflow run ... -resume` can pick up the preview's empty session and
    conclude there is nothing cached. Either run `-preview` from a separate
    directory, or resume from an explicit session ID.

---

## Tier 2 — the tiny dataset

Tier 1 proves the wiring. Tier 2 proves it can move biological information between tools as expected.

The tiny dataset is a subset of the public **10x Genomics 10k PBMC multiome**
sample: roughly 1,000 cells with fragments restricted to `chr21` and `chr22`,
**78.9 MB** in total. It is drawn from the same source data as the full PBMC
example, so it exercises the real code paths rather than a separate toy path.

The subset allows for efficient testing, but its' premise disallows interpreting pipeline inputs or outputs as meaningful biology. A few examples:

- **ATAC QC thresholds are set explicitly.** The shipped auto-thresholds are
  derived from whole-genome fragment counts and reject nearly every cell when
  fragments are restricted to two chromosomes.
- **Downstream cell-count floors are lowered** (`hdwgcna.min_cells = 50`). At
  1,000 cells a PBMC cell type is ~30–400 cells, so the shipped default of 100
  would skip most cell types.

**What it covers:** manifest parsing and pre-flight, RNA QC → CellBender →
CellTypist annotation, ATAC QC → peak calling → clustering → ATAC cell-type
annotation, Cicero co-accessibility, hdWGCNA (whole-dataset and per-cell-type),
CellChat, and MOFA+/MultiVI integration. Measured: **94/94 tasks, exit 0,
1 h 42 m 45 s wall-clock, 6.6 CPU-hours** on 8 CPUs.

Note that scVI/scANVI training does **not** run here — annotation goes through
CellTypist (`rna.annotation_method = 'celltypist'`), so `TRAIN_SCVI` and
`TRAIN_SCANVI` are not part of the 94 tasks. scvi-tools CPU viability and
seeding were verified separately.

### Tier 2 missess:

| Stage | Config | Why it is excluded |
|---|---|---|
| SCENIC+ / pycisTopic | `scenicplus.run`, `pycistopic.run` = `false` | Expensive compute and needs cisTarget databases |
| scPRINTER enhancer footprinting | `scprinter.run`, `enhancer_footprinting.run` = `false` | Expensive compute |
| ChromVAR | `chromvar.run` = `false` | `bin/gpu_chromvar_nf.py` imports `cupy`/`rmm` at module scope — a hard GPU dependency, and this tier is CPU-only |
| Differential workflows | `differential*.run` = `false` | The subset is a single-condition design, so every condition-aware workflow has nothing to contrast |

**No precomputed on-ramp bundle ships with the tutorial.** The
[on-ramp mechanism](onramps.md) is documented on its own terms instead.

!!! note "Status"
    The tiny dataset is **built and verified** — it runs end to end, twice,
    94/94 tasks, exit 0. Tier 1 is fully
    available today against any dataset config, and tier 3 is reproducible from
    the published configs in `configs/datasets/`.

TODO: Confirm data is available for Tier 2.
---

## Tier 3 — the published datasets

Four datasets, spanning two species, two assay platforms, and single- and
two-condition designs:

| Dataset | Species | Platform | Design |
|---|---|---|---|
| PBMC 10k | Human (hg38) | 10x | Single condition |
| Alzheimer's mouse model | Mouse (mm10) | 10x | WT vs TG |
| Mouse brain | Mouse | BD Rhapsody | Single condition |
| Mouse kidney | Mouse | BD Rhapsody | Single condition |

Their configs are in `configs/datasets/`. Measured cost across all four:

| Metric | Value |
|---|---|
| Unique tasks | 4,984 |
| Compute-hours | 1,012 |
| Node-hours | 6,547 |

---

## Recommended order

1. **`nextflow config`** — confirm your layered parameters resolved as intended.
2. **`-preview`** — confirm the manifest, references, and DAG are coherent.
3. **Tiny dataset** — confirm real outputs appear, end to end.
4. **Your data, minimal config** — RNA + ATAC only, expensive blocks off.
5. **Enable one block at a time**, re-running `-preview` after each config change.

Steps 1 and 2 cost seconds and catch the large majority of problems. Do not skip
them.
