# Reference File Manifest

This document lists every external reference file required by FORGE, grouped by tool. For each entry it records the file size, canonical download source, the path used on UCI HPC3 (Swarup Lab), and any notes about custom build steps or modifications.

**Legend**
- `[standard]` — File downloaded as-is from the source; no modification.
- `[custom-built]` — File assembled, filtered, or reformatted from raw source data. Build scripts are provided.
- `[runtime-subsampled]` — Full file is downloaded; subsampling happens automatically inside the pipeline at runtime.

---

## Genome annotations (Gencode GTF)

Used by: Cicero, pycisTopic, scPRINTER, SCENIC+

| File | Assembly | Size | Source | HPC3 path | Note |
|------|----------|------|--------|-----------|------|
| `gencode.v38.annotation.gtf` | GRCh38 (hg38) | 1.46 GB | [Gencode release 38](https://www.gencodegenes.org/human/release_38.html) | `/dfs7/swaruplab/lesolano/ref/Gencode_GRCh38/` | `[standard]` |
| `gencode.vM10.annotation.gtf` | GRCm38 (mm10) | 802 MB | [Gencode mouse release M10](https://www.gencodegenes.org/mouse/release_M10.html) | `/dfs7/swaruplab/lesolano/ref/Gencode_GRCm38/` | `[standard]` |
| `gencode.vM37.annotation.gtf` | GRCm39 (mm39) | — | [Gencode mouse release M37](https://www.gencodegenes.org/mouse/release_M37.html) | `/dfs7/swaruplab/lesolano/ref/Gencode_GRCm39/` | `[standard]` |

### Download

```bash
# Human GRCh38 v38
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
gunzip gencode.v38.annotation.gtf.gz

# Mouse GRCm38 vM10
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz
gunzip gencode.vM10.annotation.gtf.gz
```

---

## Genome FASTA and indices

Used by: SCENIC+ (sequence context), scPRINTER (bias model)

| File | Assembly | Size | Source | HPC3 path | Note |
|------|----------|------|--------|-----------|------|
| `hg38.fa` + `.fai` | GRCh38 | ~3 GB | Gencode / UCSC | `/dfs7/swaruplab/lesolano/ref/Gencode_GRCh38/GRCh38_FaGFF3/` | `[standard]` |
| `GRCm38.primary_assembly.genome.fa` + `.fai` | GRCm38 | 2.77 GB | [Gencode](https://www.gencodegenes.org/mouse/release_M10.html) | `/dfs7/swaruplab/lesolano/ref/Gencode_GRCm38/` | `[standard]` |
| `GRCm39.primary_assembly.genome.fa` + `.fai` | GRCm39 | — | [Gencode](https://www.gencodegenes.org/mouse/release_M37.html) | `/dfs7/swaruplab/lesolano/ref/Gencode_GRCm39/` | `[standard]` |

---

## ATAC blacklist regions

Used by: SnapATAC2 peak calling (pycisTopic `blacklist_bed` parameter)

| File | Assembly | Source | HPC3 path | Note |
|------|----------|--------|-----------|------|
| `hg38-blacklist.v2.bed.gz` | hg38 | [Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist) | `/dfs7/swaruplab/lesolano/ref/blacklist_ATAC_multispecies/Blacklist/lists/` | `[standard]` |
| `mm10-blacklist.v2.bed.gz` | mm10 | [Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist) | `/dfs7/swaruplab/lesolano/ref/blacklist_ATAC_multispecies/Blacklist/lists/` | `[standard]` |

---

## CellTypist pre-trained models

Used by: `RUN_CELLTYPIST`, `ATAC_CELLTYPIST`

CellTypist models are downloaded from the [CellTypist model zoo](https://www.celltypist.org/models) and are bundled inside the `scgpu_extended.sif` container. **No modifications were made to any model.** Users extending the pipeline to new tissues should download the appropriate `.pkl` file from the model zoo and set `params.celltypist.model` accordingly.

| Model file | Tissue | Config parameter | Note |
|------------|--------|-----------------|------|
| `Immune_All_Low.pkl` | Human PBMC / immune cells | `params.celltypist.model` | `[standard]` — bundled in container |
| `Mouse_Whole_Brain.pkl` | Mouse brain | `params.celltypist.model` | `[standard]` — bundled in container |

Additional available models (not used in current examples):
- `Pan_Fetal_Human.pkl` — human fetal tissues
- `Developing_Mouse_Brain.pkl` — developmental mouse brain

To use a model not bundled in the container, set `params.celltypist.model` to the full path of a locally downloaded `.pkl` file.

---

## scATAnno reference atlases

Used by: `ATAC_SCATANNO` (reference-based ATAC cell-type annotation)

scATAnno is a SnapATAC2-based tool that projects query ATAC cells onto a reference atlas for annotation. The tool ships with limited built-in references. **Both atlases used by FORGE were custom-built** from public source data, as described below.

### PBMC atlas (Human)

| File | Size | HPC3 path |
|------|------|-----------|
| `PBMC_reference_atlas_final.h5ad` | 2.76 GB | `/dfs7/swaruplab/lesolano/ref/scATAnno_PBMC/` |

**Status:** `[custom-built]`  
**Source data:** Healthy adult PBMC multiome data (10x Genomics).  
**Build date:** April 8, 2025.  
**Why custom:** The scATAnno package provides a healthy adult reference atlas, but we rebuilt it from curated source data to ensure compatibility with our SnapATAC2 version and to control cell-type label resolution.

An alternative version (`Healthy_Adult_reference_atlas.h5ad`, 3.0 GB) is available in the same directory.

> Build scripts and notebook are available in `/dfs7/swaruplab/lesolano/ref/scATAnno_PBMC/`. A user-facing build tutorial will be added to `docs/containers.md` in a future update.

### Mouse brain atlas

| File | Size | HPC3 path |
|------|------|-----------|
| `mouse_brain_reference_atlas.h5ad` | 1.96 GB | `/dfs7/swaruplab/lesolano/ref/scATAnno_mouse_brain/` |

**Status:** `[custom-built]`  
**Source data:** GEO accession [GSE246791](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246791) — Whole Mouse Brain (WMB) SnapATAC2 anndata (126 GB compressed tar archive).  
**Build date:** April 10, 2025.  
**Why custom:** No mouse brain scATAnno reference atlas existed at the time of development. We extracted the WMB SnapATAC2 anndata, mapped subclass labels using a hand-curated name mapping (`sa2.subclass.names.map.csv`), and built a compact atlas suitable for reference-based annotation.

Build scripts: `build_reference.py`, `build_reference.sh` in `/dfs7/swaruplab/lesolano/ref/scATAnno_mouse_brain/`.

### Download source data for custom build

```bash
# Mouse brain — download WMB SnapATAC2 anndata from GEO
# GSE246791: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246791
# File: GSE246791_wmb_SnapATAC2_anndata.tar.gz (126 GB)
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246791/suppl/GSE246791_wmb_SnapATAC2_anndata.tar.gz"
```

---

## Allen Brain atlas (mouse)

Used by: `PREPARE_REFERENCE` → `TRAIN_SCVI` → `TRAIN_SCANVI` (Path A RNA annotation)

| File | Size | HPC3 path |
|------|------|-----------|
| `AllenRef_mouse10xv2.h5ad` | 145.7 GB | `/dfs7/swaruplab/lesolano/ref/ALLEN_mouse/` |

**Status:** `[runtime-subsampled]`  
**Source:** [Allen Brain Cell Atlas (ABCA)](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas) — mouse 10x v2 transcriptomics (1.17 million cells).  
**Build scripts:** `buildref.sh`, `buildref.py`, `BuildAllenMouseRef.ipynb` in `/dfs7/swaruplab/lesolano/ref/ALLEN_mouse/`.

**Important:** The full 145.7 GB `.h5ad` is required on disk. FORGE does **not** pre-downsample this file. Instead, the `PREPARE_REFERENCE` process automatically subsamples to a maximum of 50,000 cells (stratified by cell type) at runtime before training scVI/scANVI. This subsampling is necessary because loading the full 1.17M-cell reference into GPU memory for training is not feasible on current HPC nodes.

### Download

The Allen Brain Cell Atlas data is distributed through AWS S3. See the [ABC Atlas Data Access tutorial](https://alleninstitute.github.io/abc_atlas_access/intro.html) for programmatic download instructions.

---

## SEA-AD human brain atlas

Used by: `PREPARE_REFERENCE` → `TRAIN_SCVI` → `TRAIN_SCANVI` (Path A RNA annotation, human brain datasets)

| File | Size | HPC3 path |
|------|------|-----------|
| `SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad` | 36.3 GB | `/dfs7/swaruplab/lesolano/ref/SEA_AD/` |

**Status:** `[standard]` — downloaded from the SEA-AD portal, no modification.  
**Source:** [Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)](https://portal.brain-map.org/explore/seattle-alzheimers-disease) — Middle Temporal Gyrus (MTG) snRNAseq, final QC'd nuclei, 2024-02-13 release.

Like the Allen mouse atlas, this reference is `[runtime-subsampled]` to 50,000 cells during `PREPARE_REFERENCE`.

> **Note:** SEA-AD is not currently used in the 10x example configs (PBMC uses Path B / CellTypist-only; the Alzheimer's mouse example uses the Allen mouse atlas). It is the preferred reference for human brain multiome datasets.

---

## ChromVAR motif library (cisBP)

Used by: `GPU_CHROMVAR`

| File | Organism | Source | HPC3 path | Note |
|------|----------|--------|-----------|------|
| `cisBP_2.00_human.meme` | Human | [cisBP v2.00](http://cisbp.ccbr.utoronto.ca/) | `/dfs7/swaruplab/lesolano/ref/snapatac2/` | `[standard]` |
| `cisBP_2.00_mouse.meme` | Mouse | [cisBP v2.00](http://cisbp.ccbr.utoronto.ca/) | `/dfs7/swaruplab/lesolano/ref/snapatac2/` | `[standard]` |

These files are also distributed with the scATAnno package. Alternatively, download directly from the cisBP database in MEME format.

---

## JASPAR 2022 motif database

Used by: `SCPRINTER_FOOTPRINTING`, `SCPRINTER_MOTIF_SCAN`

| File | Size | Source | HPC3 path | Note |
|------|------|--------|-----------|------|
| `JASPAR2022_core_nonredundant.jaspar` | small | [jaspar.elixir.no](https://jaspar.elixir.no/) | `/dfs7/swaruplab/lesolano/ref/scprinter/` | `[standard]` |

```bash
# Download
wget https://jaspar.elixir.no/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.zip
unzip JASPAR2022_CORE_non-redundant_pfms_jaspar.zip
```

---

## SCENIC+ / cistarget rankings and scores

Used by: `SCENICPLUS_RUN`

These are large pre-computed region-vs-motif ranking and score matrices from the [aertslab cistarget resources](https://resources.aertslab.org/cistarget/). Download from the `v10nr_clust` collection.

### Human (hg38)

| File | Size | HPC3 path | Note |
|------|------|-----------|------|
| `hg38_screen_v10_clust.regions_vs_motifs.rankings.feather` | 35.2 GB | `/dfs7/swaruplab/lesolano/ref/scenic_plus_resources/human/` | `[standard]` (symlink → `hg38_region_based.rankings.feather`) |
| `hg38_screen_v10_clust.regions_vs_motifs.scores.feather` | 13.9 GB | `/dfs7/swaruplab/lesolano/ref/scenic_plus_resources/human/` | `[standard]` (symlink → `hg38_region_based.scores.feather`) |
| `hg38_screen_v10_clust.genes_vs_motifs.rankings.feather` | 311 MB | same | `[standard]` — gene-based alternative |
| `hg38_screen_v10_clust.genes_vs_motifs.scores.feather` | 342 MB | same | `[standard]` |
| `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` | 98.7 MB | `/dfs7/swaruplab/lesolano/ref/scenic_plus_resources/` | `[standard]` — HGNC motif annotations |

### Mouse (mm10)

| File | Size | HPC3 path | Note |
|------|------|-----------|------|
| `mm10_screen_v10_clust.regions_vs_motifs.rankings.feather` | 17.8 GB | `/dfs7/swaruplab/lesolano/ref/scenic_plus_resources/mouse/` | `[standard]` |
| `mm10_screen_v10_clust.regions_vs_motifs.scores.feather` | 8.2 GB | same | `[standard]` |
| `mm10_screen_v10_clust.genes_vs_motifs.rankings.feather` | 237 MB | same | `[standard]` |
| `mm10_screen_v10_clust.genes_vs_motifs.scores.feather` | 266 MB | same | `[standard]` |
| `motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl` | 113.1 MB | `/dfs7/swaruplab/lesolano/ref/scenic_plus_resources/` | `[standard]` — MGI motif annotations |

### Mouse (mm39)

| File | Size | HPC3 path | Note |
|------|------|-----------|------|
| `mm39_screen_v10_clust.regions_vs_motifs.rankings.feather` | 22.6 GB | `/dfs7/swaruplab/lesolano/ref/scenic_plus_resources/mouse/` | `[standard]` |
| `mm39_screen_v10_clust.regions_vs_motifs.scores.feather` | 10.9 GB | same | `[standard]` |

### Download

```bash
# Human rankings (region-based, hg38) — ~35 GB
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather

# Human scores
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather

# Motif annotations (HGNC)
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

# Mouse mm10 rankings
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather
```

For the full aertslab resource listing, see: https://resources.aertslab.org/cistarget/databases/

---

## scPRINTER cache

Used by: `SCPRINTER_BUILD_PRINTER`, `SCPRINTER_FOOTPRINTING`, `SCPRINTER_MOTIF_SCAN`

scPRINTER maintains a persistent local cache of genome files, dispersion models, bias correction tracks, and cisBP motif libraries. The cache is ~95 GB for a single species.

| Key file | Size | Note |
|----------|------|------|
| `hg38_bias_v2.bw` | 17.5 GB | Tn5 insertion bias bigwig for hg38 |
| `dispersion_model_py_v2.h5` | small | Per-cell read depth dispersion model |
| `gencode_v41_GRCh38.fa.gz` | 844 MB | Genome sequence (human) |
| `gencode_v25_GRCm38.fa.gz` | — | Genome sequence (mouse mm10) |
| `CisBP_Human_FigR/` | — | Pre-formatted human cisBP motifs |
| `CisBP_Mouse_FigR/` | — | Pre-formatted mouse cisBP motifs |
| `JASPAR2022_core_nonredundant.jaspar` | small | (also referenced directly — see above) |

Set `params.scprinter.cache_dir` to a directory with sufficient disk space. On first run, scPRINTER will populate the cache automatically by downloading from its CDN. For air-gapped or restricted-network environments, pre-populate the cache manually following the [scPRINTER documentation](https://github.com/pinellolab/scPRINTER).

**HPC3 path:** `/dfs7/swaruplab/lesolano/ref/scprinter/`

---

## Storage summary

| Category | Total size | Species |
|----------|-----------|---------|
| Allen Brain mouse atlas | 145.7 GB | Mouse |
| SEA-AD human brain atlas | 36.3 GB | Human |
| SCENIC+ feather indices | ~70 GB | Human + Mouse |
| scPRINTER cache | ~95 GB | Human + Mouse |
| scATAnno atlases | ~5 GB | Human + Mouse |
| Genome FASTA files | ~6 GB | Human + Mouse |
| Gencode GTF files | ~2.3 GB | Human + Mouse |
| Motif databases | ~500 MB | Human + Mouse |
| Blacklists | small | Human + Mouse |
| **Total (approximate)** | **~360 GB** | (without SEA-AD and with one species of scPRINTER) |
| **Total (full, both species)** | **~600 GB** | |

---

## Config parameter quick reference

The table below maps each reference file to the `params` key that points to it in `nextflow.config` / your dataset config.

| `params` key | Reference | Tool |
|--------------|-----------|------|
| `gtf_human_full` | Gencode v38 GTF | Cicero, SCENIC+, scPRINTER |
| `gtf_mouse_full` | Gencode vM10 GTF | Cicero, SCENIC+, scPRINTER |
| `blacklist_bed` | ENCODE blacklist BED | pycisTopic |
| `ref_dir_human_integrated` | SEA-AD directory | scANVI (Path A) |
| `ref_dir_mouse_integrated` | Allen mouse directory | scANVI (Path A) |
| `scatanno.reference_atlas` | scATAnno .h5ad | scATAnno |
| `atac.cisbp_human` | cisBP human.meme | ChromVAR |
| `atac.cisbp_mouse` | cisBP mouse.meme | ChromVAR |
| `scenicplus.ctx_rankings` | cistarget rankings .feather | SCENIC+ |
| `scenicplus.ctx_scores` | cistarget scores .feather | SCENIC+ |
| `scenicplus.motif_annotations` | motifs .tbl | SCENIC+ |
| `scenicplus.fai` | genome .fa.fai | SCENIC+ |
| `scprinter.cache_dir` | scPRINTER cache directory | scPRINTER |
| `scprinter.pfms` | JASPAR .jaspar | scPRINTER |
| `cicero.gtf_full` | Gencode GTF | Cicero |
| `pycistopic.gtf` | Gencode GTF | pycisTopic |
| `pycistopic.blacklist_bed` | ENCODE blacklist | pycisTopic |
