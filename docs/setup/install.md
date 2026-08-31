# Installation and Getting Started

## Prerequisites

| Requirement | Minimum version | Notes |
|-------------|-----------------|-------|
| [Nextflow](https://www.nextflow.io/) | **25.10.0** (supported: `>=25.04.0, <26.0.0`) | Install to a user directory, no `sudo` — see [Quickstart Step 1](../quickstart.md). Pin with `export NXF_VER=25.10.0`; the installer otherwise fetches a release too new to run FORGE. |
| [Singularity](https://sylabs.io/singularity/) or [Apptainer](https://apptainer.org/) | 3.8 | Must be on `$PATH` |
| Java | 17 | Required by Nextflow |
| SLURM | any | HPC job submission |
| NVIDIA GPU | **optional** — A30-class only for scVI/scANVI | FORGE runs CPU-only end to end (see the [Tutorial](../tutorial.md)). V100-class is ample for CellBender, ChromVAR, MOFA+ and MultiVI. See [GPU tiers](../index.md#gpu-tiers). |
| High-memory node | ≥ 256 GB RAM | SCENIC+, large reference atlas loading |
| Disk space | ~600 GB (references) + ~13 GB (containers) | See reference manifest |

> **HPC note:** FORGE is designed for SLURM-managed clusters with Singularity/Apptainer. Local execution works for development, dry-runs, and the tutorial, but is not practical for full datasets. Stages that need high memory will fail locally unless that memory is available; GPU stages fall back to CPU rather than failing.

---

## 1. Clone the repository

```bash
git clone https://github.com/swaruplabUCI/FORGE.git
cd FORGE
```

---

## 2. Obtain Singularity containers

FORGE uses five containers ranging from 1.1 GB to 4.7 GB each (~13 GB total). GitHub's file size limit makes hosting them in the repository impractical. Two options are available:

### Option A: Download pre-built containers (recommended)

Pre-built `.sif` files are available at: **[download link — to be provided]**

Place all five files in `singularity_cache/` inside the cloned repository:

```
FORGE/
└── singularity_cache/
    ├── scgpu_extended.sif      (~3.6 GB)
    ├── snapatac_extended.sif   (~4.7 GB)
    ├── scenicplus.sif          (~1.8 GB)
    ├── cicero.sif              (~1.9 GB)
    └── seurat_extended.sif     (~1.1 GB)
```

```bash
mkdir -p singularity_cache
# Download each .sif into singularity_cache/ using the link above
```

### Option B: Build from source

Complete Singularity `.def` files and step-by-step build instructions are in [docs/containers.md](containers.md). Building is also the right path if you need to extend a container with additional tools.

---

## 3. Download reference files

FORGE requires several large external references. The complete manifest with download sources, sizes, and notes on custom-built files is in [docs/references.md](references.md).

The tables below list the minimum references required for each 10x example.

### Human PBMC (hg38)

| Reference | Size | Source |
|-----------|------|--------|
| Gencode v38 GTF | 1.46 GB | [gencodegenes.org — human release 38](https://www.gencodegenes.org/human/release_38.html) |
| ENCODE hg38 blacklist v2 | small | [github.com/Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist/tree/master/lists) |
| cisBP 2.00 human.meme | small | bundled with scATAnno (see references.md) |
| PBMC scATAnno atlas | 2.76 GB | custom-built — see [references.md §scATAnno](references.md#scatanno-reference-atlases) |
| SCENIC+ hg38 rankings (.feather) | 35.2 GB | [aertslab cistarget resources](https://resources.aertslab.org/cistarget/) |
| SCENIC+ hg38 scores (.feather) | 13.9 GB | [aertslab cistarget resources](https://resources.aertslab.org/cistarget/) |
| SCENIC+ HGNC motif annotations | 98.7 MB | [aertslab cistarget resources](https://resources.aertslab.org/cistarget/) |
| JASPAR 2022 core non-redundant | small | [jaspar.elixir.no](https://jaspar.elixir.no/) |
| scPRINTER cache (hg38) | ~95 GB | auto-populated on first run or manual download via scPRINTER |

### Mouse Alzheimer's model (mm10)

| Reference | Size | Source |
|-----------|------|--------|
| Gencode vM10 GTF | 802 MB | [gencodegenes.org — mouse release M10](https://www.gencodegenes.org/mouse/release_M10.html) |
| ENCODE mm10 blacklist v2 | small | [github.com/Boyle-Lab/Blacklist](https://github.com/Boyle-Lab/Blacklist/tree/master/lists) |
| cisBP 2.00 mouse.meme | small | bundled with scATAnno (see references.md) |
| Mouse brain scATAnno atlas | 1.96 GB | custom-built from GEO GSE246791 — see [references.md §scATAnno](references.md#scatanno-reference-atlases) |
| Allen Brain mouse atlas | 145.7 GB | custom-built from Allen Brain Cell Atlas — see [references.md §Allen](references.md#allen-brain-atlas-mouse) |
| SCENIC+ mm10 rankings (.feather) | 17.8 GB | [aertslab cistarget resources](https://resources.aertslab.org/cistarget/) |
| SCENIC+ mm10 scores (.feather) | 8.2 GB | [aertslab cistarget resources](https://resources.aertslab.org/cistarget/) |
| SCENIC+ MGI motif annotations | 113.1 MB | [aertslab cistarget resources](https://resources.aertslab.org/cistarget/) |
| JASPAR 2022 core non-redundant | small | [jaspar.elixir.no](https://jaspar.elixir.no/) |
| scPRINTER cache (mm10) | ~95 GB | auto-populated on first run or manual download via scPRINTER |

---

## 4. Configure reference paths

Open your dataset config (e.g., `examples/nextflow_PBMC_Hs_10X.config`) and update all absolute paths to match your local filesystem. The key parameters to change:

```groovy
params {
    // Genome annotation
    gtf_human_full = '/your/path/to/gencode.v38.annotation.gtf'

    // Blacklist
    // (set in atac and pycistopic blocks)

    // scATAnno
    scatanno {
        reference_atlas = '/your/path/to/PBMC_reference_atlas_final.h5ad'
    }

    // ChromVAR motifs
    atac {
        cisbp_human = '/your/path/to/cisBP_2.00_human.meme'
    }

    // Cicero
    cicero {
        gtf_full = '/your/path/to/gencode.v38.annotation.gtf'
        gtf_plot = '/your/path/to/gencode.v38.annotation.gtf'
    }

    // SCENIC+
    scenicplus {
        ctx_rankings      = '/your/path/to/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather'
        ctx_scores        = '/your/path/to/hg38_screen_v10_clust.regions_vs_motifs.scores.feather'
        motif_annotations = '/your/path/to/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
        gtf               = '/your/path/to/gencode.v38.annotation.gtf'
        fai               = '/your/path/to/hg38.fa.fai'
    }

    // scPRINTER
    scprinter {
        cache_dir = '/your/path/to/scprinter_cache'
        pfms      = '/your/path/to/JASPAR2022_core_nonredundant.jaspar'
    }
}
```

---

## 5. Prepare input data and metadata

### Metadata CSV

Create a `metadata.csv` file with one row per sample:

```csv
sample_id,fragment_file,condition_group
pbmc_10k,/path/to/pbmc_10k_atac_fragments.tsv.gz,Control
```

Required columns: `sample_id`, `fragment_file`.  
Additional columns (e.g., `condition_group`, `batch`) are passed through to output objects.

Set `params.metadata_file = '/path/to/metadata.csv'` in your dataset config.

### RNA input

For 10x Multiome data, set `params.batch_dirs` to a map of sample names to their CellRanger-ARC output directories:

```groovy
params {
    batch_dirs = [
        'pbmc_10k': '/path/to/cellranger_arc_output/pbmc_10k'
    ]
}
```

---

## 6. Run the pipeline

```bash
# PBMC 10x example
nextflow run main.nf \
    -profile cluster,singularity \
    -c examples/nextflow_PBMC_Hs_10X.config \
    -resume

# AD mouse 10x example
nextflow run main.nf \
    -profile cluster,singularity \
    -c examples/nextflow_AD_Mm_10X.config \
    -resume
```

`-resume` enables Nextflow's task-level caching. Restarting after a failure or parameter change will skip already-completed steps.

### Resource tier

The default tier (`auto`, equivalent to `small`) is appropriate for PBMC-scale data (~10k cells, single sample). Override at the command line or in your config:

```bash
# Override on the command line
nextflow run main.nf -profile cluster,singularity \
    -c my_dataset.config \
    --resource_tier medium \
    -resume
```

| Tier | Intended scale | Notes |
|------|---------------|-------|
| `small` / `auto` | ≤ 20k cells, 1–5 samples | Default; safe for any dataset |
| `medium` | 20k–100k cells | Intermediate memory/CPU |
| `large` | > 100k cells, 50+ samples | BD-scale; hugemem nodes required |

---

## Example datasets

### 10x PBMC — Human (hg38)

- **Dataset:** 10x Genomics Multiome PBMC 10k
- **Data type:** Paired GEX + ATAC (10x Chromium Multiome)
- **Download:** [10x Genomics Datasets](https://www.10xgenomics.com/datasets) — *10k Human PBMCs, Multiome v1.0, Chromium X*
- **Config:** `examples/nextflow_PBMC_Hs_10X.config`
- **Entry script:** `examples/main_PBMC_Hs_10X.nf`
- **Annotation path:** CellTypist (`Immune_All_Low.pkl`) + scATAnno (PBMC atlas)

### 10x Alzheimer's Mouse Model (mm10)

- **Dataset:** [download link — to be provided]
- **Data type:** Paired snRNA-seq + snATAC-seq (10x Chromium Multiome)
- **Config:** `examples/nextflow_AD_Mm_10X.config`
- **Entry script:** `examples/main_AD_Mm_10X.nf`
- **Annotation path:** CellTypist (`Mouse_Whole_Brain.pkl`) + scATAnno (mouse brain atlas) + scANVI against Allen Brain mouse atlas

---

## Cluster configuration

FORGE targets [UCI HPC3](https://rcic.uci.edu/hpc3/) by default. The `cluster` profile loads `configs/profiles/hpc3_cluster.config`.

For a different SLURM cluster, copy `configs/profiles/hpc3_cluster.config`, adapt the partition/account names and GPU node settings, add a new profile entry to `nextflow.config`, and invoke with `-profile yourcluster,singularity`.

Key SLURM parameters (can be overridden in your dataset config or on the CLI):

```groovy
params {
    slurm_account               = 'your_lab_account'
    slurm_account_gpu           = 'your_gpu_account'
    slurm_partition_cpu         = 'standard'
    slurm_partition_gpu         = 'gpu'
    slurm_partition_gpu_hugemem = 'gpu-hugemem'
    slurm_gpu_type              = 'A30'
    slurm_gpu_count             = 1
}
```

---

## On-ramp: resuming from pre-computed intermediates

Any major output artifact can be injected at startup to skip upstream stages. This is useful when iterating on downstream analysis without re-running QC or integration.

```groovy
params {
    onramp {
        // Skip CellBender + RNA QC + scVI/scANVI — start from integrated RNA
        rna_integrated_h5ad   = '/path/to/rna_integrated.h5ad'

        // Skip all ATAC processing — start from peak matrix
        atac_peak_matrix_h5ad = '/path/to/atac_peak_matrix.h5ad'

        // Skip Cicero — inject pre-computed connections
        cicero_connections    = '/path/to/cicero_connections.csv'
        cicero_ccan           = '/path/to/CCAN_assignments.tsv.gz'
        cicero_cds            = '/path/to/input_cds_ordered.rds'

        // Skip ChromVAR
        chromvar_deviations   = '/path/to/chromvar_deviations.h5ad'
        chromvar_raw          = '/path/to/chromvar_raw.h5ad'

        // Skip scPRINTER printer build
        printer_h5ad          = '/path/to/printer.h5ad'
    }
}
```

All on-ramp keys default to `null` (pipeline runs from raw data). See `nextflow.config` for the complete list.

---

## Troubleshooting

**Exit code 137 (out of memory)**  
Increase resource tier (`--resource_tier large`) or add a process-level override in your dataset config:
```groovy
process {
    withName: 'SCENICPLUS_RUN' { memory = '512.GB' }
}
```

**Singularity bind/permission errors**  
Check that `singularity.runOptions` in `nextflow.config` includes `--bind /your/data/path`. The default binds `/dfs7` (UCI HPC3 scratch); adjust for your filesystem.

**GPU processes not dispatched to GPU queue**  
Verify your cluster profile sets `--gres=gpu:...` in `clusterOptions` for `process_gpu`-labeled processes. Check `configs/profiles/hpc3_cluster.config` as a reference.

**`Invalid resource_tier value`**  
Values are case-sensitive: `small`, `medium`, `large`, `auto`. `Medium` or `LARGE` will cause a preflight error.

**`sample_id_regex did not match barcode sample id`**  
The scPRINTER manifest builder expects barcode filenames to contain the sample ID. Set `params.scprinter.sample_id_regex` to a named-capture regex whose `id` group extracts your metadata `sample_id`. Example: `'^(?<id>.+?)_batch\\d+$'`.
