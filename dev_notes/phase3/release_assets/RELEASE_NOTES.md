# FORGE tutorial dataset — v1 (DRAFT)

> **Status: draft. Nothing has been uploaded to GitHub.** These assets are staged
> under `dev_notes/phase3/release_assets/`. Cutting the actual Release is a
> separate, deliberate step — see "When you cut the release" at the bottom.

Suggested tag: `tutorial-data-v1`
Suggested title: **FORGE tutorial dataset v1 (10k PBMC multiome, chr21+chr22)**

---

## What this release provides

The tiny dataset behind FORGE's [Tier 2 verification](../../../docs/verification.md)
and [tutorial](../../../docs/tutorial.md). It runs the real pipeline — real
containers, real tools — CPU-only, in about an hour and three quarters, with no
cluster, no GPU, and no 600 GB reference bundle.

| Asset | Size | Purpose |
|---|---|---|
| `forge_tutorial_pbmc_v1.tar.gz` | 36 MB (78.9 MB unpacked) | the dataset |
| `expected_results.json` | 2 KB | structural contract to diff your run against |
| `figures.tar.gz` | ~2 MB | 12 reference figures spanning every arm |

Source: a subset of the public
[10x Genomics 10k PBMC multiome (Chromium X)](https://www.10xgenomics.com/datasets)
sample — ~1,000 cells, ATAC fragments restricted to `chr21` and `chr22`, plus
chr21/22-subset GENCODE and blacklist references.

## Install

```bash
mkdir -p tutorial_data
tar -xzf forge_tutorial_pbmc_v1.tar.gz -C tutorial_data/
```

Yielding `tutorial_data/{manifest.csv,refs/,samples/}`. That is the layout
`configs/datasets/tutorial_pbmc.config` expects by default; keep it elsewhere and
pass `--tutorial_data /path/to/data`.

Verify the download:

```bash
sha256sum -c forge_tutorial_pbmc_v1.tar.gz.sha256
```

```
3ad17cd76ed301da1d345b815af70c729ce8b55b77f4e8fadfa4b76b5d12c883  forge_tutorial_pbmc_v1.tar.gz
```

## Run

```bash
nextflow run . -profile tutorial,singularity -c configs/datasets/tutorial_pbmc.config
```

`tutorial` selects the CPU-only resource tier; `singularity` keeps containers on.
Budget ~15 GB free disk (3.88 GB results + 4.81 GB work).

## Checking your run

`expected_results.json` splits into two blocks, and the distinction matters:

- **`structural`** — deterministic. Two cold runs reproduced these exactly, with
  scvi-tools seeded (`params.random_seed = 42`). If your run differs here,
  something is genuinely different about your inputs or container and is worth
  investigating.
- **`informational`** — wall-clock, peak RSS, directory sizes. **Not** stable run
  to run. Recorded for scale only; do not assert on them.

Headline structural values:

| | |
|---|---|
| Tasks completed | 94 / 94, exit 0 |
| RNA | 924 cells × 21,014 genes (CellTypist) |
| ATAC | 944 → 817 cells, 12,085 peaks |
| Paired multiome | 767 cells |
| MOFA+ | 3 factors |
| Cicero | 3,498 connections, 166 CCAN assignments |
| CellChat | 2,029 interactions |
| hdWGCNA | 45 cell types with output |

## Known limitations — read before interpreting anything

**This dataset demonstrates that the pipeline runs. It does not produce
biologically meaningful results, and is not meant to.**

Restricting ATAC to `chr21`+`chr22` means TSS enrichment, peak counts and QC
distributions are computed from ~2.6% of the genome. They will not resemble a
real run.

The clearest consequence: **ATAC cell-type annotation degenerates to a single
label** (`Plasma_cells`) across all 817 cells, because the marker panel driving it
is genome-wide. RNA-vs-ATAC concordance is therefore **3.79%, exactly equal to the
chance floor** — the join is real and worth verifying, but the agreement number
carries no information. Do not quote it as a quality metric. On the published
whole-genome datasets, PBMC concordance is 0.917.

Structure upstream of labelling is healthy: ATAC clustering resolves 4 / 7 / 9
Leiden clusters, and the RNA arm gives a sensible PBMC composition.

Several stages are off by design — ChromVAR (hard GPU dependency),
SCENIC+/pycisTopic (cisTarget databases cannot be subset), scPRINTER footprinting
(the single most expensive process in FORGE), and all differential workflows
(single-condition design). No on-ramp bundle ships to stand in for them.

## Reproducing these assets

The tarball is built deterministically — fixed mtime, ownership and member order —
so rebuilding from the same `tutorial_data/` reproduces the same sha256:

```bash
tar --sort=name --owner=0 --group=0 --numeric-owner \
    --mtime='2026-08-07 00:00:00 UTC' \
    -czf forge_tutorial_pbmc_v1.tar.gz \
    -C tutorial_data manifest.csv refs samples
```

`expected_results.json` and the reference figures are generated, never
hand-written, so they cannot drift from the run they describe:

```bash
singularity exec --bind /dfs7,/tmp singularity_cache/snapatac_extended.sif \
    python3 dev_notes/phase3/build_release_assets.py \
        --results results_tutorial_remeasure \
        --outdir dev_notes/phase3/release_assets
```

Per-figure sha256 values are in `figures/CHECKSUMS.txt`.

---

## When you cut the release

Two things are worth knowing before this goes out:

1. **On a private repo, release asset URLs require authentication** — an anonymous
   `wget` gets a 404. The download instructions above only become true for other
   people once the repo is public.
2. **GitHub Pages does not serve private repos** outside Enterprise Cloud, and
   `.github/workflows/docs.yml` triggers on pushes to `main` while all tutorial
   work currently sits on `dev`. The docs site needs both resolved before it is
   reachable.

The command, once you have decided:

```bash
gh release create tutorial-data-v1 \
    --repo swaruplabUCI/FORGE \
    --title "FORGE tutorial dataset v1 (10k PBMC multiome, chr21+chr22)" \
    --notes-file dev_notes/phase3/release_assets/RELEASE_NOTES.md \
    dev_notes/phase3/release_assets/forge_tutorial_pbmc_v1.tar.gz \
    dev_notes/phase3/release_assets/forge_tutorial_pbmc_v1.tar.gz.sha256 \
    dev_notes/phase3/release_assets/expected_results.json
```

Afterwards the two remaining `<!-- FILL: -->` slots in `docs/tutorial.md`
(`download-instructions`, `reference-outputs`) can be filled with the real asset
URLs.
