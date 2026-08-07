# Verifying FORGE works

FORGE's published results come from four datasets that took roughly **1,000
compute-hours** in total. That demonstrates the pipeline at scale, but it is a
poor way to convince yourself — or a reviewer — that the pipeline *functions*.

So verification is split into three tiers. Each answers a different question, at
a very different cost.

| Tier | Question it answers | Time | Needs |
|---|---|---|---|
| **1. Pre-flight + DAG** | Is my configuration coherent? Does the whole graph wire up? | **~15 seconds** | Nextflow only |
| **2. Tiny dataset** | Does a real analysis arm produce real numbers? | **~15–60 minutes** | ~200 MB data, 2 containers, CPU only |
| **3. Full datasets** | Does it scale, and reproduce the published figures? | hours–days | Full references, GPU, HPC |

Start at tier 1. It costs nothing and catches most mistakes.

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

Expected result, in about fifteen seconds:

```text
PRE-FLIGHT CHECKLIST PASSED (9 checks):
    [OK] Manifest schema (2 rows)
    [OK] Species/genome consistency
    [OK] GTF files (5 paths validated)
    [OK] scATAnno reference atlas (stub_atlas.h5ad)
    [OK] MOFA mode (high_memory)
    [OK] scPRINTER genome (hg38)
    [OK] CellTypist model: Immune_All_Low.pkl
    [OK] Containers (5 SIF files)
    [OK] Resource tier (test)
  No warnings.
```

`-profile test` strips every site-specific scheduler assumption (SLURM
partitions, accounts, QOS, `--gres` GPU strings) so this works on any machine,
and redirects `outdir` to `results_test/` so nothing can touch real results.

This is the artifact to hand a reviewer who wants to confirm the pipeline is real
without provisioning a cluster.

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

Every error is actionable and names the parameter to fix. Fix them and re-run —
each cycle is seconds.

**Or a clean graph**, meaning your manifest parses, your species and genome
builds agree, every reference your enabled modules need exists, your on-ramp
bundles are complete, and the process graph resolves.

### What tier 1 covers — and what it does not

`-preview` builds the entire workflow graph, so it catches more than parameter
validation. It surfaces every construction-time defect: an unresolved module
include, a channel referenced outside the scope where it was declared, a `val`
input that evaluates to `null`, or a missing nested config block. In practice
this is where most breakage lives, which is why it is worth running after any
edit.

What it cannot catch is anything that only happens once tasks run:

| Not covered by tier 1 | Why | Covered by |
|---|---|---|
| Runtime channel joins | e.g. the per-sample multiome join keys on output *filenames*, which only exist once processes have run | Tier 2 |
| Tool behaviour and numerical output | No process body executes | Tier 2 / 3 |
| Resource sizing (OOM, walltime) | Nothing is scheduled | Tier 3 |

So a clean tier-1 result means "this configuration is coherent and the pipeline
will start" — not "this run will finish."

You can also confirm what your layered config actually resolved to, without
launching anything:

```bash
nextflow -c configs/datasets/my_study.config config -profile cluster,singularity
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
