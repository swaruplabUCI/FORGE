# Verifying FORGE works

FORGE's published results come from four datasets that took roughly **1,000
compute-hours** in total. That demonstrates the pipeline at scale, but it is a
poor way to convince yourself — or a reviewer — that the pipeline *functions*.

So verification is split into three tiers. Each answers a different question, at
a very different cost.

| Tier | Question it answers | Time | Needs |
|---|---|---|---|
| **1. Pre-flight** | Does my configuration make sense? Does the DAG wire up? | **~10 seconds** | Nextflow only |
| **2. Tiny dataset** | Does a real analysis arm produce real numbers? | **~15–60 minutes** | ~200 MB data, 2 containers, CPU only |
| **3. Full datasets** | Does it scale, and reproduce the published figures? | hours–days | Full references, GPU, HPC |

Start at tier 1. It costs nothing and catches most mistakes.

---

## Tier 1 — pre-flight validation

This requires **no containers, no reference files, and no GPU**. It runs on a
laptop.

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

Every error is actionable and names the parameter to fix. Fix them and re-run —
each cycle is seconds.

**Or a clean graph**, meaning your manifest parses, your species and genome
builds agree, every reference your enabled modules need exists, your on-ramp
bundles are complete, and the process graph resolves.

You can also confirm what your layered config actually resolved to, without
launching anything:

```bash
nextflow config -profile cluster,singularity -c configs/datasets/my_study.config
```

!!! tip "This is the check to run after any config edit"
    Ten seconds here regularly saves a multi-hour failure at hour three of a run.

!!! warning "Don't follow `-preview` with a bare `-resume`"
    A `-preview` invocation is recorded in Nextflow's run history. A subsequent
    bare `nextflow run ... -resume` can pick up the preview's empty session and
    conclude there is nothing cached. Either run `-preview` from a separate
    directory, or resume from an explicit session ID.

---

## Tier 2 — the tiny dataset

Tier 1 proves the wiring. Tier 2 proves the science runs.

The tiny dataset is a subset of the public **10x Genomics 10k PBMC multiome**
sample: roughly 1,000 cells with fragments restricted to two chromosomes, about
200 MB in total. It is the same dataset as the full PBMC example, so the tiny run
and the published run differ *only in scale* — which makes it a genuine test of
the pipeline rather than a separate toy path.

**What it covers:** manifest parsing and pre-flight, RNA QC → scVI → clustering →
CellTypist annotation, ATAC QC → peak calling → clustering, Cicero
co-accessibility, ChromVAR motif enrichment, and MOFA+/MultiVI integration.

**What it deliberately does not cover, and why:**

| Stage | Why it is excluded | How it is demonstrated instead |
|---|---|---|
| SCENIC+ / pycisTopic | Needs cisTarget databases that cannot be meaningfully subset to two chromosomes | Precomputed outputs injected via [on-ramps](onramps.md) |
| scPRINTER enhancer footprinting | Single most expensive process in FORGE — 54% of all compute across the four published datasets | Precomputed outputs injected via on-ramps |

Using on-ramps for those two is not a workaround bolted on for the tutorial —
it is the same mechanism FORGE uses in production to resume from checkpoints, so
exercising it here is itself part of the verification.

!!! note "Status"
    The tiny dataset is being packaged for release with a citable DOI. Until it
    lands, tier 1 is fully available today against any dataset config, and tier 3
    is reproducible from the published configs in `configs/datasets/`.

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

The distribution is extremely uneven, which is worth knowing before you enable
everything: **`ENHANCER_FOOTPRINTING_PER_CT` alone is 54% of all compute-hours.**
Per-cell-type cost ranged from 1.2 h (AD) to 12.2 h (Brain), tracking cells and
fragments per cell type rather than any configuration difference — all four runs
used identical selection caps.

This is why `enhancer_footprinting.msfp_enabled` ships as `false`. Enable it
knowingly.

---

## Recommended order

1. **`nextflow config`** — confirm your layered parameters resolved as intended.
2. **`-preview`** — confirm the manifest, references, and DAG are coherent.
3. **Tiny dataset** — confirm real outputs appear, end to end.
4. **Your data, minimal config** — RNA + ATAC only, expensive blocks off.
5. **Enable one block at a time**, re-running `-preview` after each config change.

Steps 1 and 2 cost seconds and catch the large majority of problems. Do not skip
them because the run "should" work.
