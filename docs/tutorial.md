# Tutorial — run FORGE end to end in about three hours

This is the fastest honest way to see FORGE actually work. You run the real
pipeline with its' real containers which put real tools to work. You'll generate real outputs from a deliberately small slice of the same public 10x PBMC multiome dataset used in the manuscript validation.

---

## Goals and intentions

**It is** a wiring-and-plumbing demo. It demonstrates the pipeline is correctly
connected end to end. You'll inspect input formats, observe processes as they run, verify that containers resolve on your local hardware, and will get an idea for what kind of real files will come out the other side.

**It is not** intended to generate a meaningful biological result. The ATAC
side is restricted to **chr21 + chr22 — about 2.6% of hg38** — and the RNA side
is subsampled to ~1,000 cells.

- **QC distributions will not look like a real run.** (TSS enrichment, fragments
  per cell, fragment-size distribution, etc) Note this is also why the tutorial config sets ATAC QC thresholds explicitly rather than relying on the defaults. 
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
| Disk — results | **~320 MB** |
| Disk — including `work/` | **~15 GB free** (results ~320 MB + `work/` ~4.5 GB, plus headroom) |
| RAM — single heaviest task | **9.1 GB** (`ATAC_FINAL_PIPELINE`) |
| Wall-clock | **42 min** measured on 50 CPUs / 300 GB; **~1 h 10 min** estimated on 8 CPUs / 48 GB — see below |
| CPU-hours | **4.6** |
| GPU | **Not required.** The tutorial is CPU-only. |
| Network | Required — one 2.8 MB CellTypist model is fetched at runtime. |

### About the wall-clock figure

The intention here is to give you an idea of gained efficiency from increases in compute resources:

| Allocation | Concurrency | Wall-clock | |
|---|---|---|---|
| 50 CPUs / 300 GB | peak 25 concurrent tasks | **42 min** | measured |
| 8 CPUs / 48 GB | few tasks at a time | **~2 h 20 min** | measured |

Extra cores do not speed up the long serial stages; what they buy is throughput
during the hdWGCNA and 25-way Cicero fan-outs. Summed across every task the run
is just over **2 h 0 min** of task time. The two longest single tasks are `RUN_CELLCHAT` (22.7 min) and
`CELLBENDER` (16.5 min). Extra cores won't shorten these times significantly.

!!! note "The run reaches the network"
    `RUN_CELLTYPIST` fetches its model at runtime from
    `celltypist.cog.sanger.ac.uk` (one file, 2.8 MB). The model cache lives at
    `/tmp/.celltypist` **on the node executing the task**, because containers run
    with `--home /tmp`; it is therefore not shared between nodes and not
    persisted between runs.

You also need Nextflow and Singularity/Apptainer installed, and the FORGE containers built. See
[Installation](setup/install.md) and [Containers](setup/containers.md).

!!! tip "On a module-based cluster, load Singularity first"
    `launch_tutorial.sh` loads the module itself, so the pipeline run works
    without you doing anything. But every command on this page that you type
    **manually** (verification steps in sections 4&5) needs
    `singularity` on your `PATH`. Without it you get
    `bash: singularity: command not found`.

    ```bash
    module avail singularity          # see what your site provides
    module load singularity           # or: module load apptainer
                                      # or pin a version: module load singularity/3.11.3
    ```

    Skip this if Singularity is installed system-wide.

---

## 1. Get the tutorial dataset

The download is **36 MB**; unpacked it is **78.9 MB**. It ships as a release
asset.

Run these from inside your cloned `FORGE` directory. Replace the `cd` below with
your actual path.

```bash
cd /path/to/FORGE

REL=https://github.com/swaruplabUCI/FORGE/releases/download/tutorial-data-v1
curl -LO $REL/forge_tutorial_pbmc_v1.tar.gz
curl -LO $REL/forge_tutorial_pbmc_v1.tar.gz.sha256

sha256sum -c forge_tutorial_pbmc_v1.tar.gz.sha256   # must print: OK

mkdir -p tutorial_data
tar -xzf forge_tutorial_pbmc_v1.tar.gz -C tutorial_data/
```

Running the checksum before unpacking rather than after is recommended as it helps prevent errors due to truncated downloads.

Check what you got against what you should have:

```bash
tree tutorial_data/          # or, if tree is not installed:
find tutorial_data -type f | sort
```

Unpacked, it must land at `tutorial_data/` in the repository root, and your
output should match this tree exactly — six files, no more:

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

You should see `PRE-FLIGHT CHECKLIST PASSED (7 checks)` and `No warnings.`
If not, stop and fix it — see [Troubleshooting](troubleshooting.md).

The process list `-preview` prints is the workflow graph. To inspect it properly,
ask for the DAG — this still executes nothing:

```bash
nextflow run main.nf -preview -profile tutorial,singularity \
    -c configs/datasets/tutorial_pbmc.config \
    -with-dag tutorial_dag.mmd
```

`nextflow.config` sets no `dag.file`, so **no DAG is written unless you pass
`-with-dag`**. The tutorial's graph should contain **25 distinct processes**:

```text
ATAC_CELLTYPE_ANNOTATION  ATAC_DESCRIPTIVE_REPORT  ATAC_FINAL_PIPELINE
ATAC_INITIAL_QC           ATAC_MAKE_THRESHOLDS     BUILD_MUDATA
CELLBENDER                CICERO_ESTIMATE_DP       CICERO_FULL_CHROM
CICERO_JOIN               CICERO_TRIPLETS          CONCAT_BATCHES
CONVERT_H5AD_TO_SEURAT    EXPORT_MUDATA_RNA        HDWGCNA_ENRICHMENT
HDWGCNA_PER_CELLTYPE      MERGE_ANNOTATIONS        MOFA_INTEGRATE
MOFA_VISUALIZE            MULTIVI_INTEGRATE        MULTIVI_VISUALIZE
PLOT_POST_SCANVI          RNA_QC                   RUN_CELLCHAT
RUN_CELLTYPIST
```

List them from your own DAG with:

```bash
grep -oE '\["[A-Z][A-Z0-9_]+"\]' tutorial_dag.mmd | tr -d '["]' | sort -u
```

`.mmd` is Mermaid — GitHub renders it inline, as does <https://mermaid.live>.
Use `-with-dag tutorial_dag.html` instead for a self-contained interactive
graph.

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

### Values that must match

These are counts and order statistics over counts. They are invariant to
hardware, thread count and allocation size, and have been confirmed identical
across independent runs on different machines. **Start by checking here.**

| Stage | Quantity | Expected |
|---|---|---|
| ATAC initial QC | cells passing | **944** of 1,000 |
| ATAC | median fragments / cell | **1,419.5** |
| ATAC | median TSS enrichment | **16.300101** |
| ATAC final | peak matrix | **817 cells × 12,085 peaks** |
| ATAC | Leiden clusters @ 0.5 / 1.0 / 2.0 | **4 / 7 / 9** |
| Cicero | cells / peaks entering | **817 / 12,085** |
| Cicero | triplets | **531,088** |
| RNA | cells annotated | **924** |
| MuData | cells in RNA ∩ ATAC | **767** |
| MOFA+ | factors | **3** |
| Concordance | ATAC distinct labels | **1** (degenerate — see section 4) |
| Concordance | L1 agreement vs chance floor | **equal** |

If one of these differs, something is likely different about your inputs or
your containers.

### Values to inspect, not match

The pipeline contains numerically sensitive stages. These
quantities shift by small amounts between machines with different configurations such as core counts.
The known shifts are tiny and they do not
change any structural count above. However, they do propagate into derived values.

| Stage | Quantity | Reference run |
|---|---|---|
| RNA | genes after filtering | 21,014 |
| CellTypist | distinct raw labels | 45 |
| hdWGCNA | cell types with modules | 45 |
| Cicero | connections / CCAN assignments | 3,498 / 166 |
| CellChat | interactions | 2,029 |
| Concordance | cells scored | 766 |
| Total | tasks succeeded | 94 |

Read these and sanity check that they are the same order of magnitude. Small differences here are expected on different hardware configurations. `params.random_seed = 42`
controls every stage that exposes a seed; these stages do not expose one. Moving forward, if we are able to stabalize these into deterministic values, we will shift them to the values in the table above.

!!! warning "`atac/final/atac_pipeline_summary.json` thresholds do not reflect biology"
    That file reports `min_counts: 5000, min_tsse: 6` — the Python script's
    argparse defaults, not what was applied. The module never passes those
    flags; filtering uses the computed per-sample values in
    `atac/initial_qc/sample_thresholds.json` (here `min_counts 660.15`,
    `max_counts 14886.7`, `min_tsse 6.31`), which is what takes 944 cells down
    to 817. Read the latter file.

### RNA vs ATAC concordance. This result is a wiring check

A hallmark of FORGE is that the two modality arms are annotated *independently* — different tools, different matrices,
different `obs` columns. The ATAC side never sees the RNA labels. Ideally,
both labels land on a common cell index as real convergent evidence the cross-modal join
works.

The *agreement* between them, at this scale, is not meaningful, powered, nor interpretable evidence. Here is the measured outcome:

| Quantity | Value |
|---|---|
| Cells carrying both labels | 767 (766 scored) |
| RNA labels (CellTypist) | 35 raw → 10 broad classes |
| ATAC labels (marker path) | **1 raw → 1 broad class** |
| L1 broad-class agreement | **3.79%** |
| Chance floor | **3.79%** |

**Agreement will equal the chance floor exactly.**
The ATAC arm assigns `Plasma_cells` to all 817 cells, so its labels carry no
information and "agreement" collapses to whatever share of RNA cells happen to
fall in the matching broad class (29/766). This is
a direct consequence of restricting ATAC to `chr21`+`chr22`. 

ATAC clustering resolves **4 / 7 / 9 Leiden clusters** at resolutions 0.5 / 1.0 /
2.0, and the RNA side produces a sensible PBMC composition (448 monocytes, 105
CD4+ T, 91 DC, 34 CD8+ T, …). The structure is there; only the label assignment
degenerates.

!!! warning "Do not quote this number as a quality metric"
    It reflects the subset's restricted genome, not FORGE's cross-modal agreement.
    Reproducing the number below only confirms your run matches and is wired correctly up to this point.

Reproduce it with:

```bash
module load singularity        # or `module load apptainer` — see note below
singularity exec --bind "$PWD" singularity_cache/snapatac_extended.sif \
    python3 bin/tutorial_concordance.py \
        --h5mu results_tutorial/multiome/multivi/multivi_integrated.h5mu \
        --outdir concordance
```

!!! warning "`singularity: command not found`?"
    `launch_tutorial.sh` loads the module for you, so the pipeline itself runs
    without you ever doing it by hand. Commands you run **manually**, like the
    one above, need it on your `PATH` yourself. On a module-based cluster:
    `module load singularity` (or `apptainer`). Check what your site provides
    with `module avail singularity`. If Singularity is installed system-wide,
    no module is needed.

It reads only `obs` from the h5mu (never materialises a matrix), maps both
vocabularies onto the same 10 broad classes used for the published datasets, and
writes `summary.json`, `L1_confusion_matrix.csv` and `per_cell_results.csv`.

### Reference outputs

The same release additionally carries three artifacts for checking your run against ours.

1. A structural contract.

2. A per-file checksums.

3. Reference figures.

```bash
REL=https://github.com/swaruplabUCI/FORGE/releases/download/tutorial-data-v1
curl -LO $REL/expected_results.json
curl -LO $REL/checksums_data.txt
curl -LO $REL/figures.tar.gz
```

`expected_results.json` splits into two blocks:

- **`structural`** — cell counts, peak counts, task count, factor counts. These
  are deterministic. A difference here means something is wrong.
- **`informational`** — wall-clock, peak RSS, directory sizes. These are *not*
  stable between runs.  

Compare the structural block against the numbers your own run printed:

```bash
python3 -c "import json; d=json.load(open('expected_results.json')); print(json.dumps(d['structural'], indent=2))"
```

### Checksums — verify the numbers and inspect the figures

`figures.tar.gz` holds 12 reference figures spanning every arm — RNA QC and UMAPs,
ATAC QC, MultiVI, MOFA+ and Cicero. They are a **sanity check**. Do not checksum them:

!!! warning "Figures never hash-match, even when they are correct"
    PDF figures embed `/CreationDate` and `/ModDate`. Two runs that produce
    pixel-identical plots still hash differently. We confirmed that 7 of the 12
    reference figures are otherwise byte-identical. A `sha256sum -c` over the figure set confoundingly will report `FAILED`.

Checksum the **numeric** outputs instead. `checksums_data.txt` covers **168
files** — every `.json`, `.csv`, `.tsv` and `.tsv.gz` the run publishes,
including the ATAC thresholds and QC summaries, the CellTypist labels, all 45
hdWGCNA module tables, the MuData/MOFA+ stats and factors, the CellChat
interaction table, and the Cicero connections, CCANs and triplets:

```bash
REL=https://github.com/swaruplabUCI/FORGE/releases/download/tutorial-data-v1
curl -LO $REL/checksums_data.txt

python3 bin/verify_tutorial_outputs.py \
    --results results_tutorial \
    --checksums checksums_data.txt
```

Expected output is `168/168 matched`. 

#### Verify the rest by inspection

The outputs that carry a timestamp cannot be checksummed, so we recommend a human
verification rather than a machine one. There are only two places this applies:

| Output | How to verify |
|---|---|
| The 12 reference figures in `figures.tar.gz` | Open yours side by side with ours. You are looking for the same *shape* — cluster structure in the UMAPs, the same QC distributions, comparable MOFA+ variance bars. |
| `mofa_visualization/mofa_integration_summary.json` | Ignore the `timestamp` and `output_file` fields, which differ every run by design. Read the numbers beside them and check them against the `mofa` block of `expected_results.json`. |

If the 168 checksums pass and those two look right, your run reproduced ours. You are now ready for a real FORGE run!

!!! note "Why a script and not `sha256sum -c`"
    gzip embeds the source mtime in its member header, so
    `cicero_connections.tsv.gz` hashes differently between runs even when the TSV
    inside is byte-identical. The script hashes gzip members **decompressed**,
    which is stable.

    Regenerate the manifest after a legitimate pipeline change with
    `--write`.

We chose **not** to ship the `.h5ad` and `.h5mu` objects. They are large and go stale on every pipeline change. 

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

Four stages are intentionally switched off:

| Stage | Why it is off |
|---|---|
| **ChromVAR** | `bin/gpu_chromvar_nf.py` imports `cupy` and `rmm` at module scope and calls `cp.cuda.set_allocator`. Hard GPU dependency, no CPU fallback. |
| **scPRINTER footprinting** | The single most expensive stage in FORGE — 54% of all compute across the four published datasets. |
| **SCENIC+ / pycisTopic** | Needs cisTarget reference databases far larger than the whole tutorial dataset. |
| **Differential analyses** | The tutorial is one sample with one `condition_group`. There is nothing to contrast. |

On-ramps themselves are a real and useful mechanism. The tutorial is meant to quickly get you up and running. For on-ramps, move to full-fledged FORGE. See [On-ramps & resuming](onramps.md) for what they do, the all-or-none bundle rules, and the pre-flight checks that enforce them.
---

## 7. Caveats 

**The ATAC QC thresholds.**
`configs/datasets/tutorial_pbmc.config` sets `atac.min_counts`, `max_counts`,
`min_tsse` and `min_fragments` explicitly. Leaving them unset does not fall back
to something sensible. `modules/atac/atac_initial_qc.nf` only passes
`--min_counts` when the parameter is truthy, so a null goes to the Python
default of 5,000 fragments. Against chr21+chr22 that retains about 5% of cells
and the ATAC arm collapses. Do not copy this config's ATAC thresholds to a
whole-genome dataset.

**CellBender's `total_droplets` must be strictly below the barcode count.**
It is 15,000 here against 20,000 barcodes. Setting it equal crashes inside
CellBender's prior computation with an `IndexError`.

**Wall-clock is dominated by a few long serial
tasks.**  — `HDWGCNA_ENRICHMENT` 53 min, `CELLBENDER` 23 min, `RUN_CELLCHAT`
22 min, `MULTIVI_INTEGRATE` 19 min. 

**Peak memory may vary between runs; the structural counts do not.** During our own "tutorial" testing,  `CELLBENDER` peaked between 1.50 GB and 2.90 GB. `ATAC_INITIAL_QC` between
5.80 GB and 7.10 GB. That said, every cell count and cluster count was identical. On large clusters, we suspect runs landing on slightly different nodes with slightly different hardware configurations will lead to these variations without fundamentally impacting compute outputs. When sizing a run to a cluster, allow for reasonable headroom.

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
