# FORGE tutorial — Phase 3 handoff

**Written:** 2026-08-06. **Purpose:** resume the tutorial/verification work from
this document alone, with no prior conversation context.

---

## 1. Orientation

Two directories, easily confused:

| Path | What it is |
|---|---|
| `/dfs7/swaruplab/lesolano/src_FORGE` | **The git repo.** Public-facing. All work below happens here. |
| `/dfs7/swaruplab/lesolano/FORGE` | **Execution instances** — `AD_Mm_10X_r5`, `Brain_Mm_BD_r5`, `Kidney_Mm_BD_r5`, `PBMC_Hs_10X_r5`, plus `oneOff/`. Each has its OWN copy of `main.nf` + `modules/`. Nothing here reads from `src_FORGE`. |

**The overall mission:** a GitHub Pages tutorial for FORGE, plus a
quickly-verifiable small dataset proving FORGE works in minutes rather than the
~1,000 compute-hours the four published datasets took.

**Verification is a three-tier ladder:**

| Tier | Proves | Cost | Status |
|---|---|---|---|
| T1 `-preview` | Config coherent + full DAG constructs | ~15 s | **DONE, working** |
| T2 tiny PBMC | A real analysis arm produces real numbers | ~15–60 min | **THIS PHASE** |
| T3 four datasets | Scale + published figures | hours–days | already existed |

---

## 2. What is already done

### Committed on branch `dev` (8 commits from this session)

```
7d081c2 FEAT: dependency-free verification fixture + config-default fixes
2c96ffb DOCS: MkDocs Material documentation site
fd79e33 DOCS: fix README clone URL and correct the modules/ layout
4f1b6e1 FEAT: MSFP enhancer-strip discovery architecture + TF ChIP corroboration
bbf8330 SYNC: absorb Kidney-instance improvements (broad map + auto topic selection)
2bf5996 PERF: right-size MOFA_VISUALIZE; gitignore container build logs
010710a FIX: pycisTopic merge mangles cell_names, breaking SCENIC+ GEX/ACC join
af84fa2 FIX: hdWGCNA DME — condition recovery, active WGCNA name, FC direction
```

### T1 works today — this is the reviewer-facing artifact

```bash
cd /dfs7/swaruplab/lesolano/src_FORGE
export PATH=/dfs7/swaruplab/lesolano/tools:$PATH
nextflow run main.nf -profile test -preview -c configs/datasets/test_preview.config
# → PRE-FLIGHT CHECKLIST PASSED (9 checks)   No warnings.   ~15 s, exit 0
```

No containers, references, GPU, or downloads. Backed by:
`configs/datasets/test_preview.config`, `configs/resource_tiers/test.config`,
the `test` profile in `nextflow.config`, and the `test_data/` fixture.

### Docs site

MkDocs Material in `docs/`, deploys via `.github/workflows/docs.yml`
(`mkdocs build --strict`, so broken links fail CI). Build it locally with the
venv at `/tmp/lesolano/mkdocsvenv` **if it still exists** — otherwise recreate
per §6. Pages: `index`, `quickstart`, `verification`,
`core/{index,manifest,config,architecture}`,
`setup/{install,containers,references,cluster}`,
`guides/{rna,atac,regulatory,integration,visualization}`, `onramps`,
`troubleshooting`.

### Phase 2 was pivoted — do not redo it

An exhaustive `-stub-run` tier was prototyped (a generator wrote 121 `stub:`
blocks) and **deliberately reverted**. `modules/` is untouched. Reason: stub
bodies must reproduce per-sample output **filenames**, because FORGE's multiome
join is filename-keyed and fails with
`Fail-fast: Missing ATAC for sample X. Check ATAC h5ad basenames vs manifest sample_id`.
Not mechanical, and needs lockstep upkeep across 134 processes. `-preview`
already caught every defect without that cost. The generator is kept for
reference only at `dev_notes/phase3/gen_stubs_REVERTED.py`.

---

## 3. Phase 3 decisions (locked by the user 2026-08-06)

1. **ChromVAR is excluded from T2.** `bin/gpu_chromvar_nf.py` imports `cupy` and
   `rmm` at module top level and calls `cp.cuda.set_allocator` — a hard GPU
   dependency with no fallback. Set `chromvar.run = false` **and**
   `cicero.use_chromvar_targets = false` (it defaults `true` and would otherwise
   demand ChromVAR). Ship precomputed ChromVAR output for an on-ramp demo.
2. **Hosting: GitHub first, Zenodo later.** Start with GitHub Release assets;
   move the final version to Zenodo for a DOI. Structure the deposit so that
   move is a straight copy (single versioned tarball + `README` + `checksums.txt`).
3. **ATAC subset = chr21 + chr22** (~2.6% of hg38).
4. **Redistribution: proceed.** User will confirm 10x / cisBP / JASPAR terms
   before anything is made public.

---

## 4. Verified facts (measured, do not re-derive)

**Source data** — `/dfs7/swaruplab/lesolano/FORGE/PBMC_Hs_10X_r5/data/`:

| File | Size |
|---|---|
| `10k_PBMC_Multiome_nextgem_Chromium_X_atac_fragments.tsv.gz` | 3.3 GB |
| `…atac_fragments.tsv.gz.tbi` | present — **tabix works** |
| `10k_PBMC_Multiome_nextgem_Chromium_X_raw_feature_bc_matrix.h5` | 205 MB |

- hg38, **`chr`-prefixed**; `chr21` and `chr22` both present in the tabix index
  (index also lists scaffolds like `KI270728.1` — exclude those).
- `tabix` binary available at `/pub/lesolano/miniconda3/bin/tabix` (base conda —
  a binary, not pip; prefer a container's tabix if convenient).

**References are cheap — this was the big win:**

| Reference | Size | Plan |
|---|---|---|
| `/dfs7/swaruplab/lesolano/ref/snapatac2/cisBP_2.00_human.meme` | **0.79 MB** | ship as-is |
| `/dfs7/swaruplab/lesolano/ref/scprinter/JASPAR2022_core_nonredundant.jaspar` | **744 KB** | ship as-is |
| GTF (`ref/Gencode_GRCh38/gencode.v38.annotation.gtf`) | subset to chr21/22 → ~2 MB | ship subset |
| Blacklist | subset → KB | ship subset |
| **scATAnno PBMC atlas** | **2.76 GB** | **AVOIDED** — see below |

**The atlas is avoided via the marker path.** ATAC annotation supports only
`scatanno` (needs the 2.76 GB atlas) **or** `atac.marker_file`. There is **no
ATAC CellTypist mode** — it was removed; `bin/` has no `ATAC_CELLTYPIST` module,
and setting `annotation_method = 'celltypist'` aborts with
`ATAC_CELLTYPIST on gene activity has been removed — use scATAnno instead`.
`configs/marker_genes.json` already contains human PBMC markers → use
`atac.marker_file`.

**GPU dependence.** `train_scvi.py`, `train_scanvi.py`,
`run_multivi_integration.py` all call `.train(max_epochs=…)` with **no
`accelerator=` argument**, so scvi-tools defaults to `"auto"` (GPU if present,
else CPU). `cellbender_wrapper.py` delegates to `cellbender.base_cli.main` and
never passes `--cuda`. MOFA+ takes `gpu_mode`. **Theory says all of these are
CPU-capable; only ChromVAR is a hard blocker.** NOT yet proven empirically —
that is Step 0.

---

## 5. Phase 3 plan

### Step 0 — verify CPU fallback — ✅ **DONE 2026-08-06**

**Result: T2 can run CPU-only.** Verified on a GPU-less node (`nvidia-smi` absent)
inside `scgpu_extended.sif` **without `--nv`**, so CUDA was genuinely unavailable
(`torch.cuda.is_available() == False`, `GPU available: False, used: False`):

| Stage | Result |
|---|---|
| scVI | PASS — trained on CPU |
| scANVI | PASS |
| MultiVI | PASS (`latent shape: (300, 5)`) |
| MOFA+ | PASS with `gpu_mode=False` |
| CellBender | imports; CLI is CPU unless `--cuda`, which FORGE never passes |
| **ChromVAR** | **hard GPU dep — excluded from T2 by decision** |

Two gotchas found while doing this, both recorded so they are not repeated:

- The first run reported `MultiVI FAIL`. That was a **bug in the test**, not a CPU
  problem: it used the deprecated `MULTIVI.setup_anndata`. FORGE correctly uses
  `MULTIVI.setup_mudata` on a MuData with `rna`/`atac` modalities
  (`bin/run_multivi_integration.py:57`). `dev_notes/phase3/multivi_cpu_test.py`
  mirrors FORGE's real path and passes.
- Wrapping a test in a shell command ending in `echo "EXIT=$?"` masks the real
  exit code from the task notification. Check the script's own reported exit.

**The device is now EXPLICIT rather than inferred** (`accelerator='auto'` was
being relied on implicitly). New param:

```groovy
scvi_accelerator = 'auto'   // 'auto' | 'cpu' | 'gpu'
```

Threaded as `--accelerator` through `bin/train_scvi.py`, `bin/train_scanvi.py`,
`bin/run_multivi_integration.py` and their three modules
(`modules/integration/scvi.nf`, `modules/integration/scanvi.nf`,
`modules/multiome/multivi_integrate.nf`). Default `'auto'` preserves existing GPU
behaviour exactly; the tutorial config sets `'cpu'`. Both values verified against
scvi-tools 1.4.2 by `dev_notes/phase3/accelerator_arg_test.py`.

**Still to handle in Step 2:** `TRAIN_SCVI`, `TRAIN_SCANVI` and
`MULTIVI_INTEGRATE` carry `label 'process_gpu'`, and the production tiers attach
GPU `clusterOptions` to them. The `tutorial` tier must override those so a
CPU-only run does not request a GPU it will not use.

### Step 1 — build the subset (~2–3 h)

Work in a new dir, e.g. `/dfs7/swaruplab/lesolano/FORGE/oneOff/20260806_tutorial/`
(**not** the session scratchpad — it is wiped between sessions).

1. **Choose barcodes.** From the raw h5, rank barcodes by total UMI; take the top
   ~1,000 as "cells". Intersect with barcodes present in the chr21/22 fragments
   so both modalities cover the same cells (`sample_id` is the join key
   everywhere in FORGE).
2. **CRITICAL — keep empty droplets.** FORGE feeds CellBender a *raw*
   (unfiltered) matrix and it needs the ambient background distribution.
   Retain ~1,000 cells **plus ~19,000 low-count barcodes**. Then set
   `cellbender.expected_cells = 1000`, `total_droplets = 20000`.
   Subsetting to cells only would silently break ambient correction.
3. **RNA:** subset barcodes, **keep all genes** — CellTypist markers are
   genome-wide, so restricting genes to chr21/22 would wreck annotation.
   Write 10x-format `.h5` named `<sample_id>_raw_feature_bc_matrix.h5`.
4. **ATAC:** `tabix frags.tsv.gz chr21 chr22` → filter to the chosen barcodes →
   `bgzip` → `tabix -p bed` re-index. Estimated well under 100 MB.
5. **References:** `awk` the GTF and blacklist down to chr21/22; copy cisBP and
   JASPAR verbatim.

### Step 2 — config + resource tier (~1 h)

- `configs/datasets/tutorial_pbmc.config`
- `configs/resource_tiers/tutorial.config` — a NEW tier. Do **not** reuse the
  `test` tier: it is deliberately 1 CPU / 1 GB and cannot run real work. Add
  `tutorial` to the `allowedTiers` list in `main.nf` (~line 855, currently
  `['small','medium','large','auto','test']`) and to the tier `includeConfig`
  chain at the bottom of `nextflow.config`.

Required settings (each for a reason found the hard way):

```groovy
params {
    species       = 'human'
    metadata_file = "${projectDir}/tutorial_data/manifest.csv"
    resource_tier = 'tutorial'

    // ATAC annotation without the 2.76 GB atlas
    atac { marker_file     = "${projectDir}/configs/marker_genes.json"
           tissue_type     = 'pbmc'
           sample_metadata = "${projectDir}/tutorial_data/manifest.csv" }  // REQUIRED, see §7

    chromvar { run = false }                    // hard cupy dependency
    cicero   { use_chromvar_targets = false }   // defaults true → would demand ChromVAR

    pycistopic { run = false }                  // cisTarget DBs can't be subset
    scenicplus { run = false }
    enhancer_footprinting { run = false; msfp_enabled = false }   // 54% of all compute

    // scprinter.gtf_{human,mouse} MUST be set explicitly here — see §7
    scprinter { gtf_human = '<subset gtf>'; pfms = '<JASPAR>' }

    n_epochs_scvi = 20                          // tiny data; keep runtime low
    multivi { n_epochs = 20 }

    scvi_accelerator = 'cpu'                    // explicit CPU; verified in Step 0
}
```

**In scope for T2:** manifest + pre-flight, CellBender, RNA QC → scVI →
clustering → CellTypist, ATAC QC → peak calling → clustering → marker
annotation, Cicero, MOFA+, MultiVI, CellChat, hdWGCNA.
**Out of scope:** ChromVAR, SCENIC+/pycisTopic, scPRINTER footprinting — all
demoed via `params.onramp`.

### Step 3 — run and iterate (~half a day)

```bash
cd /dfs7/swaruplab/lesolano/src_FORGE
export PATH=/dfs7/swaruplab/lesolano/tools:$PATH
# always validate first — 15 s vs a long failure
nextflow run main.nf -preview -c configs/datasets/tutorial_pbmc.config
# then the real thing
nextflow run main.nf -profile standard,singularity \
    -c configs/datasets/tutorial_pbmc.config -resume
```

Expect several rounds of genuine failures. Record wall-time and peak RSS per
process from `logs/nextflow/trace.txt` — publish measured numbers, not guesses.

### Step 4 — on-ramp bundle (optional)

Trim precomputed ChromVAR / SCENIC+ / footprinting outputs from the real PBMC
run so users can exercise downstream visualization. Respect the all-or-none
bundle rules (Cicero is a triple, ChromVAR a pair) — `docs/onramps.md` documents
them.

### Step 5 — publish + document

- Package `tutorial_data/` + subset refs as a versioned tarball with
  `checksums.txt`; attach to a GitHub Release. Keep it Zenodo-ready.
- Rewrite the **Tier 2** section of `docs/verification.md` with measured numbers
  and remove the `!!! note "Status"` block saying it is "being packaged".
- Add `docs/tutorial.md` (add it to `nav:` in `mkdocs.yml`) and link from
  `docs/quickstart.md`.
- Document the honest caveat: **restricting to chr21+chr22 means TSS-enrichment
  and QC distributions come from ~2.6% of the genome**, so QC plots will not
  resemble a real run. State it plainly.

---

## 6. Environment gotchas

- **Nextflow:** `/dfs7/swaruplab/lesolano/tools/nextflow` (standalone launcher,
  v25.10.0). Not on `PATH` by default — `export PATH=/dfs7/swaruplab/lesolano/tools:$PATH`.
- **`-c` is a TOP-LEVEL Nextflow option** and must precede the subcommand:
  `nextflow -c f.config config -profile x`. `nextflow config -c f.config` fails
  with `Unknown option: -c`.
- **Never `pip install` against base conda.** Bare `python3` resolves to
  `/pub/lesolano/miniconda3/bin/python3`. Use a container
  (`singularity exec … scgpu_extended.sif python3 …`) or a dedicated env.
- **All exploratory Python goes through a container**, even one-liners.
  Containers are in `src_FORGE/singularity_cache/` (5 × `.sif`).
- **DFS7 is slow on cold metadata.** First `ls`/`git status` in a directory can
  take minutes. Use long timeouts and background long commands rather than
  assuming a hang.
- **Don't follow `-preview` with a bare `-resume`** — the preview's empty
  session gets picked from history. Run previews from a separate directory or
  resume an explicit session id.
- **Adding a `stub:` block does NOT bust the resume cache** (verified: identical
  task hash, `cached: 2`).
- **Docs venv:** `/tmp/lesolano/mkdocsvenv` may be gone after a reboot.
  Recreate outside base conda, then
  `pip install -r docs-requirements.txt`. Pygments is pinned to **2.18.0** on
  purpose: mkdocs-material 9.5.x needs `pymdown-extensions~=10.2`, and that pair
  crashes on Pygments ≥2.19 with
  `AttributeError: 'NoneType' object has no attribute 'replace'`.

---

## 7. Config traps that cost real time

These are unguarded by pre-flight and produce cryptic errors:

| Symptom | Cause |
|---|---|
| `Argument of \`file()\` function cannot be null` | `atac.run = true` without `atac.sample_metadata` (main.nf calls `file()` on it unguarded), or `differential_rna.run = true` without `differential_rna.group_mapping`. |
| `scprinter.gtf_human … resolves to 'null'` | `nextflow.config` interpolates it from `params.gtf_human_full` at parse time, **before** your `-c` file merges. Every dataset config must re-declare `scprinter.gtf_human/gtf_mouse` explicitly. |
| `Cannot get property 'run' on null object` | A nested block missing from base config. Fixed in 7d081c2 for `multivi.masking_sweep/driver_factors/gap_fill`; watch for others. |
| `A process input channel evaluates to null -- Invalid declaration 'val …'` | A Nextflow `val` cannot be null. The disabled sentinel is the **string** `"none"` (e.g. `scenicplus.bc_transform_func`). |
| Inner gate on, nothing happens | `msfp_strip.enabled` needs `enhancer_footprinting.msfp_enabled` too. Outer gate wins. |

---

## 8. Outstanding, NOT part of Phase 3

1. **`dev` has never been merged to `main`.** `main` is at `e3a3465`; `dev` is
   **17 commits ahead**, and **13 commits are unpushed** to `origin/dev`. The
   public GitHub face therefore lags everything since late May, including the
   docs site. User wants to review the diff before any push.
2. **`sample_type` cleanup (user is handling).** Docs currently say exactly
   `lane`. The user wants the accepted values to become `rna`/`atac`. The
   validator currently checks case-sensitively against `['lane','demux']` and
   `demux` is archival. When the code changes, update
   `docs/core/manifest.md` — the `sample_type` column row and its short
   subsection.
3. **The MSFP enhancer-strip absorb (4f1b6e1) has never run against real data.**
   Its graph construction is now verified by T1 (that is how the missing
   `ch_per_ct_input` hoist was found), but no real execution has exercised it.
4. **`example_template.config` says `sample_metadata = null` is fine**
   ("Usually same as params.metadata_file"). It is not — see §7. Worth fixing.

---

## 9. First actions on resume

Step 0 is **done** — CPU-only is confirmed viable and the accelerator is now an
explicit param. Resume at **Step 1**:

1. `cd /dfs7/swaruplab/lesolano/src_FORGE && git log --oneline -3` — confirm you
   are on `dev`.
2. Start **Step 1** (build the chr21+chr22 subset) in
   `/dfs7/swaruplab/lesolano/FORGE/oneOff/20260806_tutorial/`. Remember the two
   traps: keep ~19,000 empty droplets for CellBender, and keep all genes on the
   RNA side.
3. Re-run the tests in `dev_notes/phase3/` any time the container or scvi-tools
   version changes — they are cheap and they are what proved CPU viability.
