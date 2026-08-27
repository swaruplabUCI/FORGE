# Tutorial — run FORGE end to end in about two hours

This is the fastest honest way to see FORGE actually work. You run the real
pipeline — real containers, real tools, real numbers — on a deliberately small
slice of a public 10x PBMC multiome dataset.

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
  **45** distinct labels to ~1,000 cells, including classes that cannot exist in
  peripheral blood — `Double-negative_thymocytes`, `ELP`, `CD8a_a`,
  `Age-associated_B_cells`. That is expected at this scale. Judge the pipeline by
  whether the stages run and join, not by the labels.

If you want biology, run one of the four published datasets in
`configs/datasets/` — see [Verifying FORGE works](verification.md) for the full
three-tier picture.

---

## Before you start

| Requirement | Value |
|---|---|
| Disk — results | **3.9 GB** |
| Disk — including `work/` | **~15 GB free** (results 3.9 GB + `work/` 4.8 GB, plus headroom) |
| RAM — single heaviest task | **8.8 GB** (`ATAC_FINAL_PIPELINE`) |
| Wall-clock | **1 h 43 min** on 8 CPUs / 48 GB |
| CPU-hours | **6.6** |
| GPU | **Not required.** The tutorial is CPU-only. |
| Network | Required. See the note below. |

Wall-clock is dominated by a handful of long serial tasks, so more cores help
less than you would expect: 10 CPUs finished in 1 h 33 min, only ten minutes
faster than 8. What extra cores do buy is the 45-way hdWGCNA and 25-way Cicero
fan-outs. Summed across all 94 tasks the run is 2 h 49 min of task time.

!!! note "The run reaches the network"
    `RUN_CELLTYPIST` fetches its model at runtime from
    `celltypist.cog.sanger.ac.uk` — `Immune_All_Low.pkl`, 2.7 MB, a few seconds.
    It is not bundled in the container, so on an air-gapped cluster this step
    fails. Pass an absolute path to a pre-staged `.pkl` via `celltypist.model`
    to skip the fetch entirely.

You also need Nextflow and Singularity/Apptainer, and the FORGE containers. See
[Installation](setup/install.md) and [Containers](setup/containers.md).

---

## 1. Get the tutorial dataset

The dataset is **78.9 MB** and ships as a release asset — it is not in the git
repository.

```bash
cd /path/to/FORGE

REL=https://github.com/swaruplabUCI/FORGE/releases/download/tutorial-data-v1
curl -LO $REL/forge_tutorial_pbmc_v1.tar.gz
curl -LO $REL/forge_tutorial_pbmc_v1.tar.gz.sha256

sha256sum -c forge_tutorial_pbmc_v1.tar.gz.sha256   # must print: OK

mkdir -p tutorial_data
tar -xzf forge_tutorial_pbmc_v1.tar.gz -C tutorial_data/
```

Verify the checksum before unpacking rather than after. A truncated download
still extracts, and the failure surfaces hours later as an unreadable fragment
file rather than as a bad download.

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

You should see `PRE-FLIGHT CHECKLIST PASSED`. If not, stop and
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

<!-- Measured from SLURM job 55119808 (results_tutorial_remeasure) and confirmed
     identical against run 1. Regenerate by re-running the tutorial and reading
     pipeline_info/trace.tsv. -->

| Stage | Quantity | Expected |
|---|---|---|
| ATAC initial QC | cells passing | **944** of 1,000 |
| ATAC | median fragments / cell | **1,419.5** |
| ATAC | median TSS enrichment | **16.30** |
| ATAC final | peak matrix | **817 cells × 12,085 peaks** |
| Cicero | cells / peaks entering | **817 / 12,085** |
| CellTypist | distinct labels | **45** |
| MuData | cells in RNA ∩ ATAC | **767** (21,014 genes) |
| MOFA+ | factors | **3** (seed 42) |
| hdWGCNA | cell types with modules | **45** |
| **Total** | **tasks succeeded** | **94 / 94** |

!!! success "These numbers are reproducible"
    Every value in that table came out **identical** across two independent
    cold runs — including the median TSS enrichment agreeing to 15 significant
    figures (`16.300101023624922`). That is the `params.random_seed = 42`
    seeding of scvi-tools working end to end. If your run differs on a
    *structural* count, something is genuinely different about your inputs or
    container, and it is worth investigating rather than shrugging off.

    Runtime and peak memory are **not** equally stable — see the caveats.

!!! warning "Do not trust `atac/final/atac_pipeline_summary.json` thresholds"
    That file reports `min_counts: 5000, min_tsse: 6` — the Python script's
    argparse defaults, not what was applied. The module never passes those
    flags; filtering uses the computed per-sample values in
    `atac/initial_qc/sample_thresholds.json` (here `min_counts 660.15`,
    `max_counts 14886.7`, `min_tsse 6.31`), which is what takes 944 cells down
    to 817. Read the latter file.

### RNA vs ATAC concordance — a wiring check, not a result

The two arms are annotated *independently* — different tools, different matrices,
different `obs` columns. The ATAC side never sees the RNA labels. So the fact that
both labels land on a common cell index is real evidence the cross-modal join
works.

The *agreement* between them, at this scale, is not evidence of anything. Here is
the measured outcome, so you can check your own run against it rather than wonder:

| Quantity | Value |
|---|---|
| Cells carrying both labels | 767 (766 scored) |
| RNA labels (CellTypist) | 35 raw → 10 broad classes |
| ATAC labels (marker path) | **1 raw → 1 broad class** |
| L1 broad-class agreement | **3.79%** |
| Chance floor | **3.79%** |

**Agreement equals the chance floor exactly, and that is the expected result.**
The ATAC arm assigns `Plasma_cells` to all 817 cells, so its labels carry no
information and "agreement" collapses to whatever share of RNA cells happen to
fall in the matching broad class (29/766). This is not a pipeline failure — it is
the direct consequence of restricting ATAC to `chr21`+`chr22`: the marker panel
driving ATAC annotation is genome-wide, and on 2.6% of the genome it has almost
nothing to score against.

Everything upstream of the labelling is healthy, which is the part worth checking:
ATAC clustering resolves **4 / 7 / 9 Leiden clusters** at resolutions 0.5 / 1.0 /
2.0, and the RNA side produces a sensible PBMC composition (448 monocytes, 105
CD4+ T, 91 DC, 34 CD8+ T, …). The structure is there; only the label assignment
degenerates.

!!! warning "Do not quote this number as a quality metric"
    It measures the subset's genome restriction, not FORGE's cross-modal
    agreement. The published datasets — run on whole genomes with the scATAnno
    atlas — are where concordance is meaningful; PBMC scores **0.917** there.
    Reproduce the number below to confirm your run matches, and nothing more.

Reproduce it with:

```bash
singularity exec --bind "$PWD" singularity_cache/snapatac_extended.sif \
    python3 bin/tutorial_concordance.py \
        --h5mu results_tutorial/multiome/multivi/multivi_integrated.h5mu \
        --outdir concordance
```

It reads only `obs` from the h5mu (never materialises a matrix), maps both
vocabularies onto the same 10 broad classes used for the published datasets, and
writes `summary.json`, `L1_confusion_matrix.csv` and `per_cell_results.csv`.

### Reference outputs

The same release carries two artifacts for checking your run against ours:

```bash
REL=https://github.com/swaruplabUCI/FORGE/releases/download/tutorial-data-v1
curl -LO $REL/expected_results.json
curl -LO $REL/figures.tar.gz
```

`expected_results.json` splits into two blocks, and the split is the point:

- **`structural`** — cell counts, peak counts, task count, factor counts. These
  are deterministic. Two cold runs with `params.random_seed = 42` reproduced them
  exactly. A difference here means something real changed about your inputs or
  your containers, and is worth chasing.
- **`informational`** — wall-clock, peak RSS, directory sizes. These are *not*
  stable between runs. They are recorded so you know the scale to expect. Do not
  assert on them.

Compare the structural block against the numbers your own run printed:

```bash
python3 -c "import json; d=json.load(open('expected_results.json')); print(json.dumps(d['structural'], indent=2))"
```

`figures.tar.gz` holds 12 reference figures spanning every arm — RNA QC and UMAPs,
ATAC QC, MultiVI, MOFA+ and Cicero — plus a `CHECKSUMS.txt` you can verify with
`sha256sum -c`. They are for eyeballing shape and sanity, not for pixel diffing:
figure rendering is not byte-reproducible across matplotlib versions.

Deliberately **not** shipped: the `.h5ad` and `.h5mu` objects. They are hundreds
of megabytes and go stale on every pipeline change, which would make them a
liability rather than a reference.

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
tasks — `HDWGCNA_ENRICHMENT` 53 min, `CELLBENDER` 23 min, `RUN_CELLCHAT`
22 min, `MULTIVI_INTEGRATE` 19 min. Going from 8 to 10 CPUs saved only ten
minutes. Extra cores buy the 45-way hdWGCNA and 25-way Cicero fan-outs and
nothing else.

**Peak memory varies between runs; the structural counts do not.** Across two
cold runs `CELLBENDER` peaked at 1.50 GB then 2.90 GB, and `ATAC_INITIAL_QC` at
5.80 GB then 7.10 GB — while every cell count and cluster count was identical.
If you are sizing a machine, do not trust a single observation: allow ~50%
headroom over the numbers above.

---

## Regenerating the numbers on this page

Every measured number above comes from one clean run:

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

The measured baseline for the reference run is summarised in the table at the
top of this page and, in full, in the `expected_results.json` release asset.
