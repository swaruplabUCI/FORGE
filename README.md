# FORGE

**FORGE** is a Nextflow pipeline for end-to-end analysis of single-cell and single-nucleus multiome (RNA + ATAC) data on HPC infrastructure. It integrates ambient RNA correction, QC, unsupervised clustering, reference-based cell-type annotation, co-accessibility network inference, TF motif enrichment, footprinting, gene regulatory network reconstruction, and multi-modal visualization into a single reproducible workflow.

Developed at the [Swarup Lab](https://swaruplab.bio.uci.edu/), University of California, Irvine.

📖 **[Documentation](https://swaruplabUCI.github.io/FORGE/)** — quickstart, the three core files, and per-arm guides.

> **Verify FORGE works in ~15 seconds — no containers, references, GPU, or downloads:**
> ```bash
> nextflow run main.nf -profile test -preview -c configs/datasets/test_preview.config
> ```
> This runs against the self-contained `test_data/` fixture that ships with the
> repository. It validates the manifest, every reference and parameter check, and
> constructs the complete process graph:
> ```
> PRE-FLIGHT CHECKLIST PASSED (8 checks):   Warnings: 1 (containers absent)
> ```
> Swap in `-c your_dataset.config` to validate your own setup the same way. FORGE
> reports every configuration problem at once, before any compute is submitted.
> See [Verifying FORGE works](docs/verification.md).

---

## Pipeline overview

```
Raw 10x/BD data
  │
  ├─ RNA arm ──────────────────────────────────────────────────────────────────┐
  │   CellBender → RNA QC → scVI → scANVI + CellTypist → DE (MAST) → hdWGCNA  │
  │                                                                            │
  ├─ ATAC arm ─────────────────────────────────────────────────────────────────┤
  │   SnapATAC2 QC → peak calling → scATAnno / CellTypist annotation           │
  │   → Cicero (CCANs) → ChromVAR (GPU) → SCENIC+/pycisTopic (GRN)            │
  │   → scPRINTER (TF footprinting at enhancers)                               │
  │                                                                            │
  ├─ Multiome integration ─────────────────────────────────────────────────────┤
  │   MOFA+ → MultiVI → MuData export                                          │
  │                                                                            │
  └─ Communication & visualization ───────────────────────────────────────────┘
      CellChat → enhancer footprinting recipes (A/B/C) → genome browser tracks
```

---

## Key capabilities

| Stage | Tools |
|-------|-------|
| Ambient RNA correction | CellBender |
| RNA QC, integration, batch correction | scanpy, scVI |
| Reference-based cell-type annotation (RNA) | scANVI, CellTypist |
| ATAC QC, peak calling, clustering | SnapATAC2 |
| Reference-based cell-type annotation (ATAC) | scATAnno |
| Co-accessibility networks | Cicero (R) |
| TF motif enrichment | ChromVAR (GPU-accelerated, cisBP) |
| Gene regulatory networks | SCENIC+ / pycisTopic |
| TF footprinting at enhancers | scPRINTER (JASPAR 2022) |
| RNA co-expression networks | hdWGCNA |
| Cell-cell communication | CellChat |
| Multi-modal integration | MOFA+, MultiVI |
| Differential expression | MAST |
| Differential accessibility | SnapATAC2 |

---

## Supported genomes

| Assembly | Species | Annotation |
|----------|---------|------------|
| GRCh38 (hg38) | Human | Gencode v38 |
| GRCm38 (mm10) | Mouse | Gencode vM10 |
| GRCm39 (mm39) | Mouse | Gencode vM37 |

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [Singularity](https://sylabs.io/singularity/) or [Apptainer](https://apptainer.org/) ≥ 3.8
- SLURM-managed HPC cluster
- NVIDIA GPU (A30 or equivalent) — required for CellBender, scVI/scANVI, ChromVAR, MOFA+, MultiVI
- High-memory nodes (≥ 256 GB RAM) — required for SCENIC+ and large reference atlases
- ~600 GB disk space for reference files (see [docs/setup/references.md](docs/setup/references.md))

---

## Quick start

```bash
# 1. Clone
git clone https://github.com/swaruplabUCI/FORGE.git
cd FORGE

# 2. Obtain containers (see docs/setup/install.md)
#    Place all .sif files in singularity_cache/

# 3. Download reference files (see docs/setup/references.md)
#    Update paths in your dataset config

# 4. Run — PBMC 10x example
nextflow run main.nf \
    -profile cluster,singularity \
    -c examples/nextflow_PBMC_Hs_10X.config \
    -resume
```

Full instructions: [docs/setup/install.md](docs/setup/install.md)

---

## Example datasets

FORGE ships with two fully annotated 10x Genomics example configurations (BD configurations are available but require collaborator metadata not yet released):

| Dataset | Species | Config | Entry script |
|---------|---------|--------|--------------|
| PBMC — healthy adult | Human (hg38) | `examples/nextflow_PBMC_Hs_10X.config` | `examples/main_PBMC_Hs_10X.nf` |
| Alzheimer's mouse model | Mouse (mm10) | `examples/nextflow_AD_Mm_10X.config` | `examples/main_AD_Mm_10X.nf` |

Data download links and setup: [docs/setup/install.md#example-datasets](docs/setup/install.md#example-datasets)

---

## Containers

FORGE runs in five Singularity containers. The **build recipes** (`.def` files, ~40 KB total) ship in [docs/defs/](docs/defs/), so the containers are fully reproducible from this repository — you do not need to download anything large. The built `.sif` **images** are 1.8–4.7 GB each (~14 GB total), too large for git.

Building a `.sif` needs root or `--fakeroot`. If your cluster restricts that, either ask your administrators to run the build, or build on a machine where you do have privileges and copy the images across — running a `.sif` requires no elevated privileges. See [docs/setup/containers.md](docs/setup/containers.md).

| Container | Role |
|-----------|------|
| `scgpu_extended.sif` | Python/GPU: CellBender, scVI, scANVI, CellTypist, MOFA+, MultiVI |
| `snapatac_extended.sif` | SnapATAC2, scATAnno, scPRINTER, ChromVAR |
| `scenicplus.sif` | SCENIC+, pycisTopic |
| `cicero.sif` | R: Cicero (Monocle3), Signac |
| `seurat_extended.sif` | R: Seurat, MAST, CellChat, hdWGCNA |

---

## Reference files

FORGE depends on large external reference atlases and motif databases that cannot be bundled with the repository (~600 GB total). See [docs/setup/references.md](docs/setup/references.md) for:

- Complete file manifest with sizes and download sources
- Notes on custom-built references (scATAnno mouse brain atlas, Allen Brain subsampling)
- HPC paths for Swarup Lab (UCI HPC3)

---

## Repository layout

```
FORGE/
├── main.nf                    # Pipeline entry point
├── nextflow.config            # Base config (all parameter defaults)
├── modules/                   # Per-tool Nextflow process definitions
│   ├── rna/                   # CellBender, QC, scVI/scANVI, MAST
│   ├── atac/                  # SnapATAC2 QC, peak calling, differential
│   ├── cellannotator/         # CellTypist, scATAnno, annotation merging
│   ├── multiome/              # MOFA+, MultiVI, MuData, pycisTopic, SCENIC+
│   ├── integration/           # Cross-modal integration & validation
│   ├── cicero/                # Co-accessibility (CCANs)
│   ├── chromvar/              # TF motif enrichment
│   ├── scprint/               # scPRINTER footprinting & enhancer strips
│   ├── cellchat/              # Cell-cell communication
│   ├── hdwgcna/               # RNA co-expression networks
│   ├── conversion/            # Format conversion (h5ad/Seurat/MuData)
│   ├── validation/            # Orthogonal validation (e.g. TF ChIP corroboration)
│   └── visualization/         # Genome browser, composites, bigWig export
├── configs/
│   ├── datasets/              # Per-experiment parameter overrides
│   ├── profiles/              # HPC cluster SLURM profile
│   └── resource_tiers/        # small / medium / large resource presets
├── bin/                       # Standalone Python helper scripts
├── examples/                  # Worked 10x example configs & launch scripts
├── singularity_cache/         # Built .sif images (not tracked in git)
├── mkdocs.yml                 # Documentation site config
└── docs/                      # Documentation site source
    ├── index.md               # Landing page
    ├── quickstart.md          # Linear start-to-finish walkthrough
    ├── verification.md        # How to verify FORGE works
    ├── core/                  # Manifest CSV, nextflow.config, main.nf architecture
    ├── setup/                 # Install, containers, references, cluster adaptation
    ├── guides/                # Per-arm how-tos
    └── defs/                  # Singularity build recipes (.def)
```

---

## On-ramp: resuming from checkpoints

Any major intermediate can be injected directly, skipping upstream stages:

```groovy
params {
    onramp {
        rna_integrated_h5ad   = '/path/to/rna_integrated.h5ad'
        atac_peak_matrix_h5ad = '/path/to/atac_peak_matrix.h5ad'
        // cicero_connections, chromvar_deviations, printer_h5ad, ...
    }
}
```

---

## Citation

If you use FORGE in your research, please cite:

> Solano LE, Swarup V, et al. FORGE: a Nextflow pipeline for end-to-end single-cell multiomics analysis. *[manuscript in preparation]*, 2026.

---

## License

BSD 3-Clause — see [LICENSE](LICENSE).

---

## Contact

Luis Enrique Solano · lesolano@uci.edu  
Swarup Lab · University of California, Irvine · https://swaruplab.bio.uci.edu/
