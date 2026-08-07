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

## 3b. INVARIANT — the tutorial must never break the working pipeline

The four published datasets are the product of ~1,000 compute-hours. **No
tutorial-enabling setting may change their behaviour.** Enforce it structurally,
not by care:

1. **Every tutorial setting lives in the tutorial config or the `tutorial` tier
   — never in `nextflow.config`'s defaults, never in a module.** Reduced epochs,
   `scvi_accelerator = 'cpu'`, disabled blocks: all scoped to
   `configs/datasets/tutorial_pbmc.config`.
2. **New base-config params must be behaviour-preserving by construction.**
   `scvi_accelerator` defaults to `'auto'`, which is exactly what scvi-tools did
   implicitly before. Modules read it as
   `${params.scvi_accelerator ?: 'auto'}` so an older config lacking the key
   still works.
3. **Never edit `small` / `medium` / `large` tiers for tutorial purposes.** Add a
   separate `tutorial` tier. Note also that adding a `withName:` block to an
   existing tier shifts the config hash of processes below it and can invalidate
   cached work on `-resume`.
4. **GPU is attached by two independent mechanisms** — get both, or a CPU run
   will queue for a GPU it never uses:
   - `label 'process_gpu'` on `TRAIN_SCVI` and `TRAIN_SCANVI`
   - tier `withName:` blocks with `containerOptions = '--nv'` and
     `accelerator = 1` (this is how `MULTIVI_INTEGRATE` gets one — it has **no**
     `process_gpu` label)

**Regression check to run after any change in this area** (verified clean at
commit `1b4027d`):

```bash
cd /dfs7/swaruplab/lesolano/src_FORGE
export PATH=/dfs7/swaruplab/lesolano/tools:$PATH
for c in ad_mm_10x pbmc_10x_10k kidney_mm_bd brain_mm_bd; do
  nextflow -c configs/datasets/$c.config config -profile cluster,gpu,singularity 2>/dev/null \
    | grep -E "scvi_accelerator|resource_tier ="
done
# Expect: scvi_accelerator = 'auto' for all four, and their original tiers
# (ad=medium, the rest=small). Anything else is a regression.
```

Also confirm `git diff` on the three scvi modules shows **only** the added
`--accelerator` line, and that `label 'process_gpu'` is still present in
`modules/integration/scvi.nf` and `modules/integration/scanvi.nf`.

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

### Step 1 — build the subset — ✅ **DONE & VERIFIED 2026-08-06**

Everything lives in `/dfs7/swaruplab/lesolano/FORGE/oneOff/20260806_tutorial/`:

```
build/                          # reusable, re-runnable build scripts
  01_extract_fragments.py       # tabix chr21+22, tally fragments/barcode
  02_subset_h5.py               # pick barcodes, write subset multiome .h5
  03_finalize.py                # filter frags, bgzip+index, subset refs, manifest
  04_verify.py                  # validate the result is a real FORGE input
  frags_chr21_22.tsv (499 MB)   # intermediate, do NOT ship
  frag_counts.tsv, selected_barcodes.tsv
out/                            # ← THIS is the shippable dataset (~79 MB)
  manifest.csv
  samples/TUTORIAL_PBMC_raw_feature_bc_matrix.h5        9.1 MB
  samples/TUTORIAL_PBMC_atac_fragments.tsv.gz          25.1 MB  (+ .tbi)
  refs/gencode_chr21_22.gtf                            44.7 MB
  refs/blacklist_chr21_22.bed                            ~KB
```

**Measured outcome:**

| Quantity | Value |
|---|---|
| Features kept | 40,545 = 36,601 GEX + 3,944 chr21/22 peaks |
| Barcodes kept | 20,000 = **1,000 cells + 19,000 background** |
| Cell UMI range | 5,291 – 34,244 |
| Background UMI range | 1 – 99 |
| Nonzeros | 3,716,357 (2.82% of the source's 131.5 M) |
| Fragments | 2,805,164 of 11,630,797 on chr21+22 (24.1%) |
| GTF records | 101,557 of 3,150,424 |

**Design decisions worth not re-litigating:**

- **The source `.h5` is a COMBINED multiome matrix** (36,601 Gene Expression +
  111,743 Peaks over 733,612 barcodes) — not RNA-only. The subset preserves that
  structure so it is a drop-in: `bin/rna_qc.py` calls `sc.read_10x_h5`, whose
  `gex_only=True` default filters to GEX itself, and CellBender consumes the raw
  `.h5` directly. Peaks are restricted to chr21/22 to stay consistent with the
  fragment subset; **all** GEX features are kept because annotation markers are
  genome-wide.
- **19,000 background barcodes are deliberate, not padding.** CellBender needs the
  ambient distribution; a cells-only matrix would silently break correction.
- **The selection validated itself:** 10,023 barcodes cleared "≥200 chr21/22
  fragments AND UMI>0", which matches the ~10k cells this dataset is named for.
- Containers: `pysam 0.23.3` is in **snapatac_extended.sif** (used for tabix
  reads and bgzf writing — neither container ships the `bgzip`/`tabix` CLIs).
  `h5py`/`scipy` work from **scgpu_extended.sif**. No base-conda tooling needed.

**Re-run the whole build** (idempotent) with the four scripts in order; see the
`singularity exec` invocations in this file's §6 for the container flags.

**Chromosome choice is a flag, not an edit.** All three build scripts take
`--chroms` (default `chr21,chr22`, matching what is built and shipped). Passing a
single chromosome was considered and rejected: the ATAC side is where the subset
bites, and dropping to one chromosome halves the gene content available to
gene-activity / marker annotation for only ~25 MB of savings. If you do change
it, pass the SAME `--chroms` to 01, 02 and 03 — they are independent processes,
not a shared config.

**Verification** — `build/04_verify.py`. First run reported one FAIL:

```
[FAIL] tabix range query returns rows — 0 in chr21:1–2Mb
```

**That was a bug in the test, not the data.** `chr21:1–2 Mb` is the acrocentric
short arm — N-rich and unmappable, so it legitimately holds zero fragments. The
real data is fine: chr21 spans 5,030,700–46,699,859 (1,102,772 fragments) and
chr22 spans 10,519,343–50,808,172 (1,702,392), both sorted, and a query at
chr21:30–31 Mb returns 19,316 rows.

Two checks after it (`5 columns per row`, `coordinates sorted`) had therefore
been passing **vacuously** on an empty list. The script now derives its probe
window from the data (midpoint of the contig's actual coordinate span), guards
those two checks against an empty result, and adds a whole-contig ordering check.

Re-run it any time the dataset is rebuilt:

```bash
W=/dfs7/swaruplab/lesolano/FORGE/oneOff/20260806_tutorial
S=/dfs7/swaruplab/lesolano/src_FORGE/singularity_cache
singularity exec --contain --home /tmp --bind /dfs7 --bind /tmp \
  --env PYTHONNOUSERSITE=1 --env HDF5_USE_FILE_LOCKING=FALSE \
  --env NUMBA_CACHE_DIR=/tmp/nb --env MPLCONFIGDIR=/tmp/mpl --env XDG_CACHE_HOME=/tmp/c \
  $S/snapatac_extended.sif python3 $W/build/04_verify.py
```

It checks: scanpy `read_10x_h5` works as `rna_qc.py` uses it; cells AND ambient
background both present; fragments are bgzf + tabix-queryable + sorted +
5-column; RNA/ATAC barcodes overlap (the cross-modality join key); GTF is
non-empty and confined to chr21/22. Exits non-zero listing failures.

**Current result: ALL 14 CHECKS PASSED, exit 0.** Notable values:

```
[PASS] GEX-only read works — 36,601 genes          (20,000 barcodes)
[PASS] has real cells (UMI>=1000) — 1,000
[PASS] has ambient background (0<UMI<100) — 19,000
[PASS] tabix range query — 23,330 rows in a populated 1 Mb window
[PASS] whole-contig ordering sorted — 1,102,772 fragments on chr21
[PASS] ATAC barcodes are a subset of RNA — 0 ATAC-only
TOTAL DATASET SIZE: 78.9 MB
```

**The dataset is a valid FORGE input. It is NOT yet a working pipeline run** —
that is Step 3, and it is where real failures are expected.

**Open size question:** the GTF is 44.7 MB of the ~79 MB total. Gzipping would cut
it substantially, but it is unclear which T2 consumers accept `.gz`, so it was
left plain rather than guessed at. Revisit if release size matters.

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
| **`resource_tier` in a `-c` file is IGNORED** | The `includeConfig` chain at the bottom of `nextflow.config` is evaluated while that file is parsed — **before** `-c` files merge. Verified: a `-c` file declaring `'tutorial'` still got small.config's 200 GB `MULTIVI_INTEGRATE` block. Set the tier in a **profile** (like `test` / `tutorial`) instead. The per-instance `launch.sh` scripts work around it by scraping the tier out of the config and re-passing `--resource_tier` on the CLI. |
| Clearing `containerOptions` globally | Wipes the `--env R_LIBS_USER=/dev/null` that `nextflow.config` sets for `RUN_CELLCHAT`/`COMPARE_CELLCHAT`, without which basilisk fails on a conda lookup. Clear GPU's `--nv` with a targeted `withName` block, never a wildcard. |

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

## 9. Progress and first actions on resume

| Step | Status |
|---|---|
| 0 — CPU-only viability + explicit `accelerator` | ✅ done, verified |
| 1 — build chr21+chr22 subset | ✅ done, **all 14 checks pass**, 78.9 MB |
| 2 — tutorial config + `tutorial` tier + profile | ✅ done, `-preview` passes (7 checks) |
| **3 — run end-to-end and iterate** | ⬜ **RESUME HERE** |
| 4 — on-ramp bundle (optional) | ⬜ |
| 5 — publish (GitHub Release) + docs | ⬜ |

**Step 2 shipped** (all verified):

- `configs/datasets/tutorial_pbmc.config` — the dataset config
- `configs/resource_tiers/tutorial.config` — CPU-only tier, real resources
- `tutorial` **profile** in `nextflow.config` — this is what selects the tier
- `'tutorial'` added to `allowedTiers` in `main.nf`

Verified: `-preview` passes 7 checks / no warnings; GPU fully stripped
(`gres=gpu` 0, non-empty `clusterOptions` 0, the GPU process set resolves to
`accelerator = null`); containers and CellChat's `R_LIBS_USER` guard preserved;
and the resolved **production** config is **byte-identical to HEAD** (1351 lines,
empty diff) with GPU intact (`gres=gpu` 10, `accelerator = 1` 6, `--nv` 10).

**Resume at Step 3 — the first real run:**

```bash
cd /dfs7/swaruplab/lesolano/src_FORGE
export PATH=/dfs7/swaruplab/lesolano/tools:$PATH
nextflow run main.nf -profile tutorial,singularity \
    -c configs/datasets/tutorial_pbmc.config -resume
```

Note `-profile tutorial,singularity` — **not** `test` (which disables containers)
and **not** relying on `resource_tier` in the `-c` file (which does not select a
tier; see §7). Expect several rounds of genuine failures; record wall-time and
peak RSS per process from `logs/nextflow/trace.txt` so the published numbers are
measured rather than guessed.

Before committing anything that touches shared config, **re-run the §3b
regression check** — ideally the byte-diff form:

```bash
timeout 400 nextflow -c configs/datasets/pbmc_10x_10k.config config \
  -profile cluster,gpu,singularity > /tmp/after.txt
git stash push -q -- nextflow.config main.nf      # tracked files only
timeout 400 nextflow -c configs/datasets/pbmc_10x_10k.config config \
  -profile cluster,gpu,singularity > /tmp/before.txt
git stash pop -q
diff /tmp/before.txt /tmp/after.txt && echo "no production impact"
```

(`git stash push` with a pathspec fails outright if any listed path is untracked —
list only tracked files, or the stash silently isn't created and the diff is
meaningless.)

Re-run the tests in `dev_notes/phase3/` whenever the containers or scvi-tools
version change — they are cheap, and they are what proved CPU viability.
