# Quickstart

This page is one linear path from an empty directory to a running pipeline. Work
through it top to bottom.

The early steps need almost nothing — you can validate a complete FORGE
configuration with only Nextflow installed, before committing to container builds
or reference downloads. Do those steps first.

---

## Step 1 — Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/     # or anywhere on your PATH
nextflow -version                    # need >= 23.04
```

On an HPC system Nextflow is often already available:

```bash
module avail nextflow
module load nextflow
```

## Step 2 — Clone FORGE

```bash
git clone https://github.com/swaruplabUCI/FORGE.git
cd FORGE
```

## Step 3 — Confirm the pipeline works, before installing anything

FORGE ships a self-contained fixture, so you can verify the whole pipeline
graph right now — no containers, no reference downloads, no GPU, no cluster:

```bash
nextflow run main.nf -profile test -preview \
    -c configs/datasets/test_preview.config
```

In about fifteen seconds you should see:

```text
PRE-FLIGHT CHECKLIST PASSED (9 checks):
    [OK] Manifest schema (2 rows)
    [OK] Species/genome consistency
    ...
  No warnings.
```

That confirms your Nextflow install, the repository, and the full process graph
are all sound. If you only do one thing from this page, do this.

You can also print the fully merged configuration at any point:

```bash
nextflow config -profile cluster,singularity
```

---

## Step 4 — Write a manifest

Create `my_manifest.csv`. One row per sample; filenames go in `rna_file` and
`fragment_file`, and the directory goes in `data_dir`:

```csv
sample_id,batch,sample_type,original_lane_id,rna_file,fragment_file,condition_group,data_dir
my_sample,batch1,lane,L1,my_sample_raw_feature_bc_matrix.h5,my_sample_atac_fragments.tsv.gz,ConditionA,/data/my_study
```

Even a single-condition dataset needs a `condition_group` value — use one label
for every row. Full column semantics: [The manifest CSV](core/manifest.md).

## Step 5 — Write a dataset config

Create `my_study.config`. Start deliberately small, with the expensive stages
off:

```groovy
params {
    species       = 'human'                  // or 'mouse' — required, no default
    metadata_file = '/path/to/my_manifest.csv'
    outdir        = 'results_my_study'
    resource_tier = 'small'

    // References (see docs/setup/references.md)
    gtf_human_full = '/refs/gencode.v38.annotation.gtf'
    blacklist_bed  = '/refs/hg38-blacklist.v2.bed'

    // Annotation. RNA uses CellTypist; ATAC requires a scATAnno atlas (or your
    // own atac.marker_file) — there is no atlas-free ATAC option.
    celltypist { model = 'Immune_All_Low.pkl' }
    scatanno   { reference_atlas = '/refs/scatanno_pbmc_atlas.h5ad' }

    // Leave the expensive arms off for the first pass
    enhancer_footprinting { msfp_enabled = false }
    scenicplus            { run = false }
    pycistopic            { run = false }
}
```

!!! tip "Copy a worked example instead"
    `configs/datasets/` holds the configs behind the four published datasets, and
    `configs/datasets/example_template.config` is a commented starting point.

## Step 6 — Validate before you run

This is the highest-value step on the page. It needs **no containers, no
references, and no GPU**:

```bash
nextflow run main.nf -preview -c my_study.config
```

FORGE builds the whole workflow graph and runs its pre-flight checklist without
submitting any work. Either you get a clean graph, or you get every problem at
once:

```text
================================================================================
PRE-FLIGHT CHECKLIST FAILED (3 error(s)):
================================================================================
  1. Manifest CSV not found: /path/to/my_manifest.csv
  2. atac.annotation_method='scatanno' requires params.scatanno.reference_atlas
     (path to a .h5ad reference). Set this explicitly in your dataset config.
  3. rna.run=true but the manifest contains no rows with a non-null rna_file.
================================================================================
```

Fix, re-run, repeat. Each cycle is seconds. Do not proceed until this is clean.

---

## Step 7 — Obtain the containers

FORGE runs in five Singularity containers. The **recipes** ship in this
repository (`docs/defs/`, about 40 KB of text), so you can build the images
yourself; you do not need to download anything large.

```bash
# On a machine where you have root or --fakeroot:
bash hpc_defs/BUILD_ON_HPC.sh all
```

If your cluster does not permit builds, either ask your administrators to run
them, or build on a workstation and copy the `.sif` files across — full guidance,
including CPU-architecture requirements, is in
[Containers](setup/containers.md#if-you-lack-root-or-fakeroot).

Place the resulting `.sif` files in `singularity_cache/`, or point
`params.containers` at wherever they live.

## Step 8 — Download references

FORGE needs external GTFs, blacklists, motif databases, and (for annotation and
GRN inference) reference atlases — roughly 600 GB for the complete set, though a
minimal RNA + ATAC run needs far less. The manifest of files, sizes, and sources
is in [Reference files](setup/references.md).

## Step 9 — Run

```bash
nextflow run main.nf \
    -profile cluster,singularity \
    -c my_study.config \
    -resume
```

Always pass `-resume`. Nextflow caches completed tasks, so an interrupted or
failed run continues rather than restarting.

!!! warning "One caveat about `-resume` after `-preview`"
    A `-preview` invocation is recorded in Nextflow's run history, and a
    subsequent bare `-resume` can latch onto that empty session and report
    nothing cached. Run previews from a separate directory, or resume from an
    explicit session ID.

---

## Step 10 — Read the output

On completion FORGE prints a map of every output location. Broadly:

| Path under `outdir` | Contents |
|---|---|
| `rna_qc/`, `integration/` | Per-sample QC and the integrated, annotated RNA object |
| `atac/final/` | Peak matrix, clustering, cell-type annotations |
| `cicero/`, `chromvar/`, `scprinter/` | Co-accessibility, motif deviations, footprints |
| `multiome/` | MOFA+ factors, MultiVI latent space, MuData export |
| `cellchat/`, `hdwgcna/` | Communication and co-expression networks |
| `logs/nextflow/` | Trace, timeline, and execution report |

`logs/nextflow/report.html` is the fastest way to see per-process runtime and
peak memory, which is what you need before enabling anything expensive.

---

## Then enable one thing at a time

Turn on a single block, re-run `-preview` to confirm the config is still
coherent, then run. Reasonable order:

1. `pycistopic.run` + `scenicplus.run` — GRN inference (needs cisTarget references)
2. `differential.run` / `differential_rna.run` — needs ≥ 2 `condition_group` values
3. `enhancer_footprinting.msfp_enabled` — **the expensive one.** This single stage
   was 54% of all compute across the four published datasets. Enable it knowingly,
   and read [the cost breakdown](verification.md#tier-3-the-published-datasets) first.

---

## Where to go next

- [The three core files](core/index.md) — the manifest, the config, the architecture
- [Verifying FORGE works](verification.md) — the full verification ladder
- [On-ramps & resuming](onramps.md) — skip stages you have already computed
- [Troubleshooting](troubleshooting.md) — common failures and their causes
