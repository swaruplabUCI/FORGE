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

| Stage | Tool | Docs | Paper |
|-------|------|------|-------|
| Ambient RNA correction | CellBender | [docs](https://cellbender.readthedocs.io/) | [Fleming et al. 2023, *Nat Methods*](https://doi.org/10.1038/s41592-023-01943-7) |
| RNA QC, clustering | scanpy | [docs](https://scanpy.readthedocs.io/) | [Wolf et al. 2018, *Genome Biol*](https://doi.org/10.1186/s13059-017-1382-0) |
| RNA integration, batch correction | scVI | [scvi-tools](https://docs.scvi-tools.org/) | [Lopez et al. 2018, *Nat Methods*](https://doi.org/10.1038/s41592-018-0229-2) |
| Reference-based annotation (RNA) | scANVI | [scvi-tools](https://docs.scvi-tools.org/) | [Xu et al. 2021, *Mol Syst Biol*](https://doi.org/10.15252/msb.20209620) |
| Reference-based annotation (RNA) | CellTypist | [celltypist.org](https://www.celltypist.org/) | [Domínguez Conde et al. 2022, *Science*](https://doi.org/10.1126/science.abl5197) |
| ATAC QC, peak calling, clustering | SnapATAC2 | [docs](https://kzhang.org/SnapATAC2/) | [Zhang et al. 2024, *Nat Methods*](https://doi.org/10.1038/s41592-023-02139-8) |
| Reference-based annotation (ATAC) | scATAnno | [PyPI](https://pypi.org/project/scATAnno/) | [Jiang et al. 2025, *Genom Proteom Bioinform*](https://doi.org/10.1093/gpbjnl/qzaf108) |
| Co-accessibility networks | Cicero | [docs](https://cole-trapnell-lab.github.io/cicero-release/) | [Pliner et al. 2018, *Mol Cell*](https://doi.org/10.1016/j.molcel.2018.06.044) |
| TF motif enrichment | chromVAR | [docs](https://greenleaflab.github.io/chromVAR/) | [Schep et al. 2017, *Nat Methods*](https://doi.org/10.1038/nmeth.4401) |
| Gene regulatory networks | SCENIC+ | [docs](https://scenicplus.readthedocs.io/) | [Bravo González-Blas et al. 2023, *Nat Methods*](https://doi.org/10.1038/s41592-023-01938-4) |
| Topic modelling of ATAC | pycisTopic | [docs](https://pycistopic.readthedocs.io/) | [Bravo González-Blas et al. 2023, *Nat Methods*](https://doi.org/10.1038/s41592-023-01938-4) |
| TF footprinting at enhancers | scPRINTER | [GitHub](https://github.com/buenrostrolab/scPrinter) | [Hu et al. 2025, *Nature*](https://doi.org/10.1038/s41586-024-08443-4) |
| RNA co-expression networks | hdWGCNA | [docs](https://smorabit.github.io/hdWGCNA/) | [Morabito et al. 2023, *Cell Rep Methods*](https://doi.org/10.1016/j.crmeth.2023.100498) |
| Cell–cell communication | CellChat | [GitHub](https://github.com/jinworks/CellChat) | [Jin et al. 2021, *Nat Commun*](https://doi.org/10.1038/s41467-021-21246-9) |
| Multi-modal integration | MOFA+ | [docs](https://biofam.github.io/MOFA2/) | [Argelaguet et al. 2020, *Genome Biol*](https://doi.org/10.1186/s13059-020-02015-1) |
| Multi-modal integration | MultiVI | [scvi-tools](https://docs.scvi-tools.org/) | [Ashuach et al. 2023, *Nat Methods*](https://doi.org/10.1038/s41592-023-01909-9) |
| Differential expression | MAST | [GitHub](https://github.com/RGLab/MAST) | [Finak et al. 2015, *Genome Biol*](https://doi.org/10.1186/s13059-015-0844-5) |
| Differential accessibility | SnapATAC2 | [docs](https://kzhang.org/SnapATAC2/) | [Zhang et al. 2024, *Nat Methods*](https://doi.org/10.1038/s41592-023-02139-8) |

**Motif databases.** Footprinting scores against
[JASPAR 2022](https://jaspar.elixir.no/)
([Castro-Mondragon et al. 2022, *NAR*](https://doi.org/10.1093/nar/gkab1113));
motif enrichment uses [CIS-BP](https://cisbp.ccbr.utoronto.ca/)
([Weirauch et al. 2014, *Cell*](https://doi.org/10.1016/j.cell.2014.08.009)).

---

## Validated genomes

These are the assemblies FORGE has been run and validated against. They are not
a whitelist — the pipeline takes genome, annotation, and motif resources as
configurable paths, so other assemblies work once you supply the matching
reference set. See [docs/setup/references.md](docs/setup/references.md) for what a new
assembly needs.

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
- **GPU — optional.** FORGE runs CPU-only end to end. See the GPU tiers below.
- High-memory nodes (≥ 256 GB RAM) — required for SCENIC+ and large reference atlases
- ~600 GB disk space for reference files (see [docs/setup/references.md](docs/setup/references.md))

### GPU tiers

FORGE is CPU-only by default and no stage requires a GPU to produce results — the
[tutorial](docs/tutorial.md) completes the full 94-task pipeline without one. Where a
GPU is available, these are the stages that use it:

| Tier | Stages | Notes |
|---|---|---|
| **A30-class** (≥ 24 GB VRAM) | `TRAIN_SCVI`, `TRAIN_SCANVI` | The only genuine A30 requirement. Latent-space training on full atlases exhausts smaller cards. |
| **V100-class** (16 GB VRAM) | `CELLBENDER`, `GPU_CHROMVAR`, `MULTIVI_*`, `MOFA_INTEGRATE` | Ample at these sizes; the shipped resource tiers already pin V100 for most MultiVI processes. |
| **CPU-only** | everything else (~100 of ~108 processes) | No GPU code path at all. |

Set `params.slurm_gpu_type` to whatever your site provides — see
[docs/setup/cluster.md](docs/setup/cluster.md) for the `--gres` syntax.

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

> Solano LE, Swarup V, et al. Flow Orchestrated Regulatory Genomics Engine
> (FORGE): A Configurable Nextflow Pipeline for End-to-End snMultiome
> Analysis. *[manuscript in preparation]*, 2026.

---

## License

BSD 3-Clause — see [LICENSE](LICENSE).

---

## Contact

Luis Enrique Solano · lesolano@uci.edu  
Swarup Lab · University of California, Irvine · https://swaruplab.bio.uci.edu/
