# Tutorial — run FORGE end to end in about two hours

This is the fastest honest way to see FORGE actually work. You run the real
pipeline — real containers, real tools, real numbers — on a deliberately small
slice of a public 10x PBMC multiome dataset.

!!! warning "Draft — measured values pending"
    Every value written as **⟨TBD⟩** below is waiting on a completed
    re-measurement run. They are marked in the source with `<!-- FILL:... -->`
    comments so they can be found with a single grep:

    ```bash
    grep -rn "FILL:" docs/tutorial.md
    ```

    Do not publish this page with any **⟨TBD⟩** remaining. See
    [Regenerating the numbers on this page](#regenerating-the-numbers-on-this-page).

---

## What this is, and what it is not

**It is** a wiring-and-plumbing demo. It proves the pipeline is correctly
connected end to end: that every process runs, every join lands, every container
resolves, and real files come out the other side.

**It is not** a biological result, and you should not read it as one. The ATAC
side is restricted to **chr21 + chr22 — about 2.6% of hg38** — and the RNA side
is subsampled to ~1,000 cells. Two consequences worth stating plainly:

- **QC distributions will not look like a real run.** TSS enrichment, fragments
  per cell and the fragment-size distribution are all computed from 2.6% of the
  genome. This is also why the tutorial config sets ATAC QC thresholds
  explicitly rather than relying on the defaults — see the caveat section below.
- **Cell-type labels are over-fit.** CellTypist `Immune_All_Low.pkl` assigns
  ⟨TBD⟩ distinct labels to ~1,000 cells<!-- FILL:n-celltypist-labels (run 1 gave 45) -->,
  including classes that cannot exist in peripheral blood (thymocyte and
  progenitor populations). That is expected at this scale. Judge the pipeline by
  whether the stages run and join, not by the labels.

If you want biology, run one of the four published datasets in
`configs/datasets/` — see [Verifying FORGE works](verification.md) for the full
three-tier picture.

---

## Before you start

| Requirement | Value |
|---|---|
| Disk — results | ⟨TBD⟩<!-- FILL:results-size (run 1 gave 3.85 GB) --> |
| Disk — including `work/` | ⟨TBD⟩<!-- FILL:total-disk; recommend a rounded "have N GB free" --> |
| RAM — single heaviest task | ⟨TBD⟩<!-- FILL:peak-rss --> |
| Wall-clock | ⟨TBD⟩<!-- FILL:wallclock; state the CPU count it was measured at --> |
| CPU-hours | ⟨TBD⟩<!-- FILL:cpu-hours --> |
| GPU | **Not required.** The tutorial is CPU-only. |
| Network | Required. See the note below. |

!!! note "The run reaches the network"
    `RUN_CELLTYPIST` downloads CellTypist models at runtime from
    `celltypist.cog.sanger.ac.uk`. This is not optional and it is not cached in
    the container. If you are on an air-gapped cluster this step will fail.

You also need Nextflow and Singularity/Apptainer, and the FORGE containers. See
[Installation](setup/install.md) and [Containers](setup/containers.md).

---

## 1. Get the tutorial dataset

The dataset is **78.9 MB** and ships as a release asset — it is not in the git
repository.

```bash
cd /path/to/FORGE
# ⟨TBD⟩ download command + checksum verification
```
<!-- FILL:download-instructions — GitHub Release URL, tarball name, sha256 line,
     and the `tar xzf` invocation. Keep it Zenodo-swappable (single versioned
     tarball + README + checksums.txt). -->

Unpacked, it must land at `tutorial_data/` in the repository root:

```text
tutorial_data/
├── manifest.csv
├── samples/
│   ├── TUTORIAL_PBMC_raw_feature_bc_matrix.h5      9.1 MB
│   ├── TUTORIAL_PBMC_atac_fragments.tsv.gz        25.1 MB
│   └── TUTORIAL_PBMC_atac_fragments.tsv.gz.tbi
└── refs/
    ├── gencode_chr21_22.gtf                       44.7 MB
    └── blacklist_chr21_22.bed
```

Put it somewhere else with `--tutorial_data /your/path`.

??? info "What is actually in the dataset"
    A subset of the public 10x Genomics 10k PBMC Multiome sample (hg38,
    `chr`-prefixed).

    | Quantity | Value |
    |---|---|
    | Barcodes | 20,000 = **1,000 cells + 19,000 background** |
    | Features | 40,545 = 36,601 GEX + 3,944 chr21/22 peaks |
    | Fragments | 2,805,164 on chr21+chr22 |
    | Cell UMI range | 5,291 – 34,244 |

    The 19,000 background barcodes are deliberate, not padding: CellBender needs
    the ambient distribution to estimate contamination, and a cells-only matrix
    would silently break correction. All GEX features are kept — only the ATAC
    side is restricted — because annotation markers are genome-wide.

---

## 2. Validate before you run

Fifteen seconds here beats a long failure. This builds the full DAG and runs the
pre-flight checklist without executing anything:

```bash
nextflow run main.nf -preview \
    -profile tutorial,singularity \
    -c configs/datasets/tutorial_pbmc.config
```

You should see `PRE-FLIGHT CHECKLIST PASSED` with no warnings. If not, stop and
fix it — see [Troubleshooting](troubleshooting.md).

!!! warning "Do not follow `-preview` with a bare `-resume`"
    The preview is written to Nextflow's run history, and a subsequent bare
    `-resume` can select that empty session. Either run the preview from a
    separate directory, or resume an explicit session id.

---

## 3. Run it

Both paths run the same pipeline with the same `tutorial` profile.

=== "Without SLURM"

    On a workstation or any single machine:

    ```bash
    ./launch_tutorial.sh
    ```

    Or call Nextflow directly if you prefer:

    ```bash
    nextflow run main.nf \
        -profile tutorial,singularity \
        -c configs/datasets/tutorial_pbmc.config \
        -resume
    ```

    Nextflow will use the whole machine's CPU and RAM for its task pool, which
    is what you want here.

=== "With SLURM"

    ```bash
    sbatch -A <your-account> -p <your-partition> launch_tutorial.sh
    ```

    Account and partition are deliberately not hardcoded — there is no portable
    default. Everything else has one; override the defaults on the sbatch line
    if your site needs different values:

    ```bash
    #SBATCH --cpus-per-task=8
    #SBATCH --mem=48G
    #SBATCH --time=08:00:00
    ```

    The whole pipeline runs inside **one** allocation using Nextflow's local
    executor. It does not submit a job per process — that would require
    site-specific account, partition and QOS settings in every `withName` block,
    which no portable tutorial can supply.

    !!! note "Why the SLURM path loads an extra config"
        `configs/tutorial_slurm.config` sizes Nextflow's task pool to your
        allocation. Nextflow's local executor otherwise reads the *machine's*
        total CPU and RAM, not the cgroup it is confined to — so inside an
        allocation it will start more concurrent tasks than the allocation can
        hold, and the job is OOM-killed with no useful error. The launcher adds
        this file only when `SLURM_JOB_ID` is set.

**Launcher options**

| Option | Effect |
|---|---|
| `--outdir DIR` | Publish results to `DIR` (default `results_tutorial`) |
| `--no-resume` | Ignore the cache and run cold |
| `--preview` | Build the DAG and exit (~15 s, no compute) |

---

## 4. What you should see

<!-- FILL:results-table — populate from the re-measure run's trace.tsv and
     outputs. Everything here is a STRUCTURAL claim (counts, shapes) that should
     be stable run to run now that scvi-tools is seeded. Values from run 1 are
     noted where known; confirm each before publishing. -->

| Stage | Quantity | Expected |
|---|---|---|
| Pre-flight | checks passed | ⟨TBD⟩ |
| CellBender | ambient-corrected counts + report | ⟨TBD⟩ |
| RNA QC | cells passing | ⟨TBD⟩ |
| CellTypist | distinct labels | ⟨TBD⟩ *(run 1: 45)* |
| ATAC initial QC | cells passing | ⟨TBD⟩ *(run 1: 944 of 1,000)* |
| ATAC | median TSS enrichment | ⟨TBD⟩ *(run 1: 16.30)* |
| ATAC | peaks called | ⟨TBD⟩ |
| Cicero | cells / peaks entering | ⟨TBD⟩ *(run 1: 817 / 12,085)* |
| Cicero | connections, CCANs | ⟨TBD⟩ |
| MOFA+ | factors, variance explained | ⟨TBD⟩ |
| MultiVI | joint latent shape | ⟨TBD⟩ |
| hdWGCNA | cell types with modules | ⟨TBD⟩ *(run 1: 45)* |
| CellChat | interactions inferred | ⟨TBD⟩ |
| **Total** | **tasks succeeded** | ⟨TBD⟩ *(run 1: 94/94)* |

### The single most informative check

**RNA vs ATAC concordance.** The two arms are annotated *independently* —
different tools, different matrices, different `obs` columns. The ATAC side never
sees the RNA labels. So agreement between them is a measured result rather than
something the pipeline constructed, and it is the best evidence that the run is
genuinely working.

⟨TBD⟩<!-- FILL:concordance — report the RNA-vs-ATAC agreement and say how it was
     computed. State the honest expectation: at 2.6% of the genome and ~1,000
     cells this will be modest. Do not oversell it. -->

### Reference outputs

⟨TBD⟩<!-- FILL:reference-outputs — link the expected_results.json + reference
     PNGs release asset, and say how to diff against it. Ship metrics and figures
     only, NOT the h5ad/h5mu objects: they are hundreds of MB and go stale on
     every pipeline change. -->

---

## 5. Where the outputs land

```text
results_tutorial/
├── cellbender/            ambient-corrected counts + HTML report
├── rna/                   QC'd and integrated RNA h5ads
├── cell_annotation/       CellTypist labels
├── atac/
│   ├── initial_qc/        thresholds, QC plots, per-sample h5ad
│   └── final/             peak matrix, cell-type annotations
├── cicero/                connections, CCANs
├── multiome/
│   ├── mudata/            joint MuData object
│   ├── mofa/              factors + model
│   └── multivi/           joint latent + UMAPs
├── cellchat/              interaction networks
└── hdwgcna/               per-cell-type co-expression modules
```

!!! note "Where the execution report actually goes"
    `trace.tsv`, `report.html` and `timeline.html` are written to
    `results/pipeline_info/` — **not** to `results_tutorial/`. The path is
    interpolated when `nextflow.config` is parsed, before the `tutorial` profile
    sets `outdir`. Look there when you want per-process runtimes and peak memory.

---

## 6. What the tutorial does not cover

Four stages are switched off, and it is worth being clear about why rather than
leaving you to discover it:

| Stage | Why it is off |
|---|---|
| **ChromVAR** | `bin/gpu_chromvar_nf.py` imports `cupy` and `rmm` at module scope and calls `cp.cuda.set_allocator`. Hard GPU dependency, no CPU fallback. |
| **scPRINTER footprinting** | The single most expensive stage in FORGE — 54% of all compute across the four published datasets. Wildly out of proportion here. |
| **SCENIC+ / pycisTopic** | Needs cisTarget reference databases far larger than the whole tutorial dataset. |
| **Differential analyses** | The tutorial is one sample with one `condition_group`. There is nothing to contrast. |

!!! note "No ChromVAR on-ramp is shipped, deliberately"
    An earlier plan was to ship a precomputed ChromVAR bundle so you could
    exercise the stages downstream of it without a GPU. It was built and it is
    barcode-compatible — but under the tutorial's own settings it unlocks
    **nothing**. Every consumer of ChromVAR output is gated behind something the
    tutorial disables: `VIS_CHROMVAR` needs differential conditions,
    `EXTRACT_CHROMVAR_MOTIFS` and `MAP_TF_TO_TARGET_GENES` need
    `scprinter.run`, `DIFFERENTIAL_TF_ACCESSIBILITY` needs
    `differential_tf.run`, and `ENHANCER_FOOTPRINTING_RECIPES` needs the
    scPRINTER printer object.

    Shipping it would have meant either turning on the most expensive stage in
    the pipeline or inventing conditions for a single-condition dataset. Neither
    is worth it for a wiring demo, so the bundle is not distributed.

On-ramps themselves are a real and useful mechanism — just not one this tutorial
needs. See [On-ramps & resuming](onramps.md) for what they do, the all-or-none
bundle rules, and the pre-flight checks that enforce them.

!!! danger "If you do use on-ramps: artifacts are keyed to the run that produced them"
    ChromVAR output is a cells × motifs matrix keyed by ATAC barcodes, and its
    consumers join it back against the same run — per-cell-type motif extraction
    needs that run's `cell_type` column, and TF-to-target mapping joins against
    its Cicero CCANs. An artifact from a *different* run overlaps on almost
    nothing and fails quietly: near-empty results that still look like success.
    Never substitute ChromVAR output from the published PBMC run into a
    tutorial-scale run, or vice versa.

---

## 7. Caveats worth stating plainly

**The ATAC QC thresholds are dataset-specific and deliberately set.**
`configs/datasets/tutorial_pbmc.config` sets `atac.min_counts`, `max_counts`,
`min_tsse` and `min_fragments` explicitly. Leaving them unset does not fall back
to something sensible — `modules/atac/atac_initial_qc.nf` only passes
`--min_counts` when the parameter is truthy, so a null goes to the Python
default of 5,000 fragments. Against chr21+chr22 that retains about 5% of cells
and the ATAC arm collapses. Do not copy this config's ATAC thresholds to a
whole-genome dataset; they are ~40× lower than appropriate.

**CellBender's `total_droplets` must be strictly below the barcode count.**
It is 15,000 here against 20,000 barcodes. Setting it equal crashes inside
CellBender's prior computation with an `IndexError`.

**Timings do not extrapolate.** Wall-clock is dominated by a few long serial
tasks and by how many of the 45-way hdWGCNA and 25-way Cicero fan-outs can run
at once. More cores help those substantially; they do not help the serial tasks
at all.

---

## Regenerating the numbers on this page

Every ⟨TBD⟩ above comes from one clean run:

```bash
sbatch -A <account> -p <partition> launch_tutorial.sh \
    --no-resume --outdir results_tutorial_remeasure
```

Then:

1. **Runtime and memory** — `results/pipeline_info/trace.tsv`. Aggregate by
   process name; report summed realtime, longest single task, and peak RSS.
2. **Structural counts** — the per-stage JSON summaries under the output tree
   (`atac/initial_qc/atac_initial_qc_summary.json`, `multiome/*/`*`_stats.json`,
   and so on).
3. **Disk** — sum `stat -c %s` rather than using `du -sh`; `du` under-reports
   freshly written files on some parallel filesystems.
4. **Confirm allocations landed** — check `trace.tsv`'s `cpus`/`memory` columns
   against the tier. This failure mode is silent, so do not assume.

Working notes and the previous measured baseline live in
`dev_notes/phase3/T2_RESOURCE_BASELINE.md`.
