# Quickstart

This page is one linear path from an empty directory to a running pipeline **on
your own data**. Work through it top to bottom.

!!! tip "Want to watch FORGE run before committing your own data?"
    Do the [Tutorial](tutorial.md) instead. It runs the complete pipeline
    end-to-end on a 79 MB public PBMC dataset in about two hours on 8 CPUs, with
    **no GPU and no external reference downloads** — everything it needs ships
    with the dataset. Steps 1-3 below are still worth doing first; they install
    Nextflow and prove the pipeline graph is sound. After that, the Tutorial is
    the faster way to see real output, and you can return here when you are
    ready to point FORGE at your own samples.

The early steps need almost nothing — you can validate a complete FORGE
configuration with only Nextflow installed, before committing to container builds
or reference downloads. Do those steps first.

---

## Step 1 — Install Nextflow

FORGE is developed and tested against **Nextflow 25.10.0**. Pin it. The supported
window is `>= 25.04.0, < 26.0.0`. Setting `NXF_VER` below is the only thing that
actually protects you: an out-of-window Nextflow aborts while parsing the config,
before the pipeline's own version check can produce a friendlier message.

No `sudo` and no system-wide install is needed. Install into a directory you own:

```bash
# NOTE ON `~`: throughout this page, `~` means the root of the directory you
# intend to house FORGE in -- not necessarily your literal home directory.
# On most HPC systems $HOME is a small, quota-limited volume, and a full run
# will not fit there. Pick a workspace with room (scratch, a lab share, a
# project volume) and treat that as `~` for every step below, e.g.:
#     cd /path/to/your/workspace
# If your home directory really is where you want FORGE, `~` works literally.

export NXF_VER=25.10.0                      # set BEFORE the install; add to ~/.bashrc too
mkdir -p ~/bin
cd ~/bin
curl -s https://get.nextflow.io | bash      # writes ./nextflow here
export PATH="$HOME/bin:$PATH"               # add to ~/.bashrc to persist
nextflow -version                           # should report 25.10.0
```

`NXF_VER` matters, and it is set *before* the install on purpose. The
`get.nextflow.io` launcher always fetches the *newest* release, which is outside
the supported window; exporting `NXF_VER` first makes it download the pinned
version directly. Set it afterwards instead and the installer still works, but it
downloads a 26.x release you never use and prints that version at you. If you skip
it entirely, Step 3 fails with `Config parsing failed`.

Nextflow needs Java 17 (`java -version` to check). Most HPC systems provide it by
default; if not, `module load java/17`.

> **HPC note:** many clusters have no `nextflow` module — UC Irvine's HPC3, for
> example, provides `singularity` and `java` but not `nextflow`. Check with
> `module avail nextflow`, and if there is no hit, use the user-directory install
> above. It needs no administrator involvement.

## Step 2 — Clone FORGE

Step 1 left you in `~/bin`. Go back to your workspace root first, or the clone
lands inside your Nextflow install directory:

```bash
cd ~                                        # your workspace root -- see the note in Step 1
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
PRE-FLIGHT CHECKLIST PASSED (8 checks):
    [OK] Manifest schema (2 rows)
    [OK] Species/genome consistency
    ...
  Warnings: 1 (see above)
```

The one warning says the container images are missing. On a fresh clone that is
expected — a preview launches no task, so nothing is ever run inside an image.

That confirms your Nextflow install, the repository, and the full process graph
are all sound. If you only do one thing from this page, do this.

You can also print the fully merged configuration at any point:

```bash
nextflow config -profile cluster,singularity
```

---

## Step 4 — Write a manifest

The next two boxes are **file contents, not commands.** Create each file with a
text editor (`nano my_manifest.csv`, `vim my_manifest.csv`, or whatever you use)
and paste the box into it. Do not paste them straight into the shell — the line
wrapping in a browser will corrupt them.

Work from wherever you want this study to live; the commands below assume you
are inside the cloned `FORGE` directory (`cd ~/FORGE`), but any directory works
as long as the paths in the config point at the right places.

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

    // These three MUST be set explicitly, even though they repeat
    // gtf_human_full above. nextflow.config interpolates them at parse time,
    // which happens BEFORE your dataset config merges, so they do not inherit
    // your value -- they resolve to 'null' and pre-flight rejects the run.
    cicero    { gtf_full  = '/refs/gencode.v38.annotation.gtf' }
    scprinter {
        gtf_human = '/refs/gencode.v38.annotation.gtf'
        gtf_mouse = '/refs/gencode.vM10.annotation.gtf'
    }

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
ERROR ~
================================================================================
PRE-FLIGHT CHECKLIST FAILED (5 error(s)):
================================================================================
  1. Manifest CSV not found: /path/to/my_manifest.csv
  2. scprinter.gtf_human (scPRINTER/enhancer/chromVAR) file does not exist:
     /refs/gencode.v38.annotation.gtf. Verify the path and ensure the GTF is
     accessible from the execution host.
  3. cicero.gtf_full file does not exist: /refs/gencode.v38.annotation.gtf.
     Verify the path and ensure the GTF is accessible from the execution host.
  4. scATAnno reference atlas not found: /refs/scatanno_pbmc_atlas.h5ad
  5. rna.run=true but the manifest contains no rows with a non-null rna_file.
     Either populate rna_file in the manifest, or set rna.run=false for
     ATAC-only runs.
================================================================================
```

That is the real output of the Step 5 config, copied verbatim — five errors, not
a trimmed illustration. They are all the placeholder paths, working as intended:
every path is checked for *existence*, not merely for being set. Point them at
your real files and the list empties out.

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
`params.containers` at wherever they live. **The directory does not exist on a
fresh clone** — it is gitignored, so create it first:

```bash
mkdir -p singularity_cache
```

!!! note "This is where the quickstart gets expensive"
    Steps 7 and 8 are the only heavy parts of setup: five container images
    (~13 GB) and reference files (up to ~600 GB for the complete set). If your
    goal right now is to *see FORGE work* rather than to process your own
    samples, stop here and do the [Tutorial](tutorial.md) — it needs the
    containers but **none** of the references on this page.

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

- [Tutorial](tutorial.md) — run the whole pipeline end to end on a small public
  dataset, CPU-only, with or without SLURM. The fastest way to see FORGE work.
- [The three core files](core/index.md) — the manifest, the config, the architecture
- [Verifying FORGE works](verification.md) — the full verification ladder
- [On-ramps & resuming](onramps.md) — skip stages you have already computed
- [Troubleshooting](troubleshooting.md) — common failures and their causes
