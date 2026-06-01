# FORGE containers

This document describes the five Singularity/Apptainer containers used by the
FORGE single-cell multiomics pipeline (Swarup Lab, UCI). It is intended for
power users on HPC clusters who want to extend or rebuild the containers from
scratch.

> **Canonical build recipes** live in [`docs/defs/`](./defs/) — five Singularity
> definition files, one per container. Build with `singularity build --fakeroot
> <name>.sif docs/defs/<name>.def`, or use the `hpc_defs/BUILD_ON_HPC.sh` wrapper
> to build them all with logging and a SHA256 manifest.
>
> Every version pin, GitHub commit, and pitfall in this doc is reproduced from
> those `.def` files and the v3.4 build logs (`scgpu_build.log`,
> `seurat_build.log`, `next_build.log`). The earlier sandbox-based builder
> (`mac_build_containers.sh`) is retained as a Mac convenience wrapper but no
> longer the authoritative recipe.

## Container overview

| Container                  | Size  | Base image                                   | Purpose                                                        |
| -------------------------- | ----- | -------------------------------------------- | -------------------------------------------------------------- |
| `scgpu_extended.sif`       | 3.7 G | `ghcr.io/scverse/scvi-tools:py3.11-cu12-base`| scVI / scANVI, CellBender, CellTypist, scrublet, MOFA+, muon   |
| `snapatac_extended.sif`    | 4.7 G | `python:3.10-slim`                           | SnapATAC2, scPrinter, scATAnno, MACS3, cupy/rmm, deeptools     |
| `seurat_extended.sif`      | 1.8 G | `rocker/r-ver:4.4.3`                         | Seurat 5, hdWGCNA, CellChat, MAST, WGCNA, zellkonverter (+ baked basilisk env) |
| `cicero.sif`               | 1.9 G | `condaforge/mambaforge:latest`               | R + Cicero (via Monocle3), Bioconductor, rtracklayer, Gviz     |
| `scenicplus.sif`           | 1.8 G | `python:3.11.8-slim`                         | SCENIC+, pycisTopic, pySCENIC, Mallet, graph-tool              |

Sizes above are the actual `.sif` sizes from `sif_output/` for the v3.4 set,
**except `seurat_extended.sif`** which lists the projected v3.5 size after
the basilisk-env bake (was 1.1 G in v3.4; +500–800 MB for the baked conda
env). Build dates of the canonical artifacts: cicero 2026-03-12, seurat
2026-03-13 (v3.4) / pending rebuild (v3.5), scgpu 2026-03-17 (CellBender
fix), scenicplus 2026-03-30, snapatac 2026-04-08 (scATAnno addition).

## Build system

Each container has a standalone Singularity definition file under
[`docs/defs/`](./defs/). The recommended path is:

```bash
# On an HPC compute node (NOT a login node)
srun --partition=free --time=04:00:00 --mem=24G --cpus-per-task=8 --pty bash
module load singularity     # or: module load apptainer
cd /path/to/forge

# Build all five with --fakeroot, log per container, summarize at the end:
bash hpc_defs/BUILD_ON_HPC.sh all

# Or a single container:
bash hpc_defs/BUILD_ON_HPC.sh seurat_extended
```

For each container, the wrapper runs `singularity build --fakeroot
<name>.sif docs/defs/<name>.def` (falling back to an unprivileged build if
`--fakeroot` is unavailable), captures stdout+stderr to
`singularity_cache/build_logs/`, and at the end prints a results table plus
a SHA256 manifest you can diff against the expected hashes (see
"Verifying a build" below).

A Mac builder is also available — `mac_build_containers.sh` provisions a
Lima VM and runs the same `.def` builds inside it. This is convenient for
preparing `.sif` files on a laptop, then `scp`-ing to HPC. The recipes are
identical; only the host orchestration differs.

### Build order

Cicero is built first in the `all` target because its mamba solve is the
longest single step (and is most likely to expose disk-space or networking
issues early):

```
cicero → scgpu_extended → snapatac_extended → seurat_extended → scenicplus
```

There are no inter-container dependencies — any one container can be rebuilt
in isolation. The order is purely a matter of failing-fast on the slowest
solve.

### Build prerequisites

- macOS (any architecture; Apple Silicon uses Rosetta x86_64 emulation
  automatically) **or** a Linux x86_64 build host with `apptainer >= 1.4`
  and `--fakeroot` support.
- Homebrew and Lima on macOS (the script installs both if missing).
- ~80 GB free disk for the Lima VM image plus another ~15 GB for the
  combined `.sif` output.
- Network access to: `ghcr.io`, `quay.io`, `docker.io`, `cran.r-project.org`,
  `bioconductor.org`, `pypi.org`, `download.pytorch.org`, `github.com`,
  `conda-forge`.

### Build times (observed)

These are wall-clock times on an Apple Silicon Mac (M-series, 16 GB host
RAM, Rosetta x86_64 emulation):

| Container                | First build | Rebuild (cached base) |
| ------------------------ | ----------- | --------------------- |
| `cicero.sif`             | ~30–45 min  | ~25 min               |
| `scgpu_extended.sif`     | ~20–30 min  | ~15 min               |
| `snapatac_extended.sif`  | ~25–35 min  | ~20 min               |
| `seurat_extended.sif`    | ~50–70 min  | ~45 min               |
| `scenicplus.sif`         | ~20–30 min  | ~15 min               |

Native x86_64 Linux builds are ~2–3× faster than Apple Silicon under
Rosetta.

---

## 1. `scgpu_extended.sif`

GPU Python container: scVI/scANVI, CellBender, CellTypist, scrublet, MOFA+
(Python side), muon, scanpy.

### Resolved versions (from `scgpu_build.log`, 2026-03-13)

| Package         | Version            |
| --------------- | ------------------ |
| scvi-tools      | 1.4.2              |
| cellbender      | 0.3.2 (commit `4334e89`) |
| muon            | 0.1.7              |
| mudata          | 0.3.3              |
| mofapy2         | 0.7.3              |
| mofax           | (PyPI latest)      |
| celltypist      | 1.7.1              |
| scrublet        | 0.2.3              |
| scanpy          | 1.11.5             |
| anndata         | 0.12.10            |
| torch           | 2.4.0+cu121        |
| numpy           | 2.1.0 (from base)  |

### Definition file

The full recipe is in [`docs/defs/scgpu_extended.def`](./defs/scgpu_extended.def). Build with:

```bash
singularity build --fakeroot scgpu_extended.sif docs/defs/scgpu_extended.def
```

### Pitfalls

- **CellBender 0.3.2 from PyPI is broken on PyTorch 2.x.** Pinning to
  GitHub commit `4334e89` is mandatory. Do not relax this to "latest" — the
  upstream PyPI release has not been re-rolled.
- **The base image already ships PyTorch 2.4.0+cu121 and numpy 2.1.0.** Do
  not downgrade either; scvi-tools 1.4 and the rest of the stack work with
  numpy 2.x in this container. (This is the opposite of `snapatac_extended`,
  where numpy must be held at < 2.0.)
- **CellTypist models are baked in at build time** so the container works
  offline on compute nodes with no outbound network. They live at
  `~/.celltypist/data/models` inside the container; the pipeline sets
  `HOME=/tmp` at runtime so the models are re-discovered via
  `CELLTYPIST_FOLDER` if you override the path.

### Runtime requirements

- `--nv` on the `singularity exec`/`run` invocation to expose host
  `libcuda.so.1`. The pipeline sets this for `scgpu_*` processes via
  `containerOptions = '--nv'` in `nextflow.config`.

---

## 2. `snapatac_extended.sif`

ATAC-side Python container: SnapATAC2, scPrinter, scATAnno, MACS3, plus
GPU-accelerated chromVAR via cupy/rmm.

### Resolved versions

From the v3.4 build (inferred from the pinned constraints + the actual SIF):

| Package         | Version constraint   |
| --------------- | -------------------- |
| snapatac2       | latest (currently 2.7.x) |
| macs3           | latest               |
| numpy           | `< 2.0` (forced after deeptools install) |
| cupy-cuda12x    | `>= 13.0, < 14`      |
| rmm-cu12        | `>= 24.0, < 25`      |
| cuda-bindings   | `12.9.4`             |
| torch           | CUDA 12.8 build      |
| kaleido         | `0.2.1`              |
| scATAnno        | latest               |
| scPrinter       | git HEAD (no PyPI release) |

### Definition file

The full recipe is in [`docs/defs/snapatac_extended.def`](./defs/snapatac_extended.def). Build with:

```bash
singularity build --fakeroot snapatac_extended.sif docs/defs/snapatac_extended.def
```

### Pitfalls

- **The numpy 1.x / 2.x boundary is the single biggest source of build
  failures.** snapatac2 requires numpy < 2.0. cupy 14.x and rmm 25.x are
  numpy-2.x-only. deeptools transitively pulls numpy 2.x. The cure is to
  re-pin `"numpy<2.0"` *after* the offending install and to pin
  `cupy-cuda12x>=13,<14` and `rmm-cu12>=24,<25`. The `%test` block
  asserts numpy is still 1.x to catch silent regressions.
- **kaleido must be pinned to `0.2.1`.** The 1.x rewrite ships a different
  API and does not work with the plotly version in this stack. The v3.2
  patch (`container_rebuild_fix.patch`) was added because kaleido failed
  silently in earlier builds; the `%test` block now exports a PNG to
  confirm it actually works end-to-end.
- **scPrinter is installed with `--no-deps`** because its `setup.py`
  over-specifies versions that conflict with the chosen torch / numpy
  pins. All scPrinter runtime deps are installed manually in the
  preceding step.
- **MACS3 has no pip wheels — it compiles from source.** Hence the
  `gcc/g++/make` system deps. This is fine on x86_64 Debian but is the
  reason the build cannot run on a base image without compilers.

---

## 3. `seurat_extended.sif`

R container: Seurat 5, hdWGCNA, CellChat, MAST, WGCNA.

### Resolved versions (from `seurat_build.log`)

| Package | Version |
| ------- | ------- |
| R       | 4.4.3 (from `rocker/r-ver:4.4.3`) |
| Seurat  | 5.4.0   |
| ggrepel | 0.9.6 (pinned — 0.9.7+ requires R ≥ 4.5) |
| magick (R) | 2.9.1 (rebuilt from source against HDRI ImageMagick) |
| BiocManager | 3.20 |
| WGCNA, MAST, hdWGCNA, CellChat | git HEAD at build time |

### Definition file

The full recipe is in [`docs/defs/seurat_extended.def`](./defs/seurat_extended.def). Build with:

```bash
singularity build --fakeroot seurat_extended.sif docs/defs/seurat_extended.def
```

> **Note on the request.** The original task list mentioned `DESeq2` and
> `clusterProfiler`. These are **not** installed in the canonical
> `seurat_extended.sif`. Instead the container ships `MAST` for DE testing
> and `enrichR` + `GSVA` + `UCell` + `GeneOverlap` for pathway / signature
> analysis. If you need DESeq2 + clusterProfiler, add them to the
> Bioconductor block above; they have no special build requirements
> beyond what is already installed.

### Pitfalls

- **`ggrepel` must be pinned to 0.9.6.** The CRAN current version requires
  R ≥ 4.5 and will fail to install against R 4.4.3.
- **`magick` from CRAN binary links against the non-HDRI ImageMagick**,
  which Ubuntu/Debian no longer ships. The container rebuilds magick from
  source after creating a `Magick++.pc → Magick++-6.Q16HDRI.pc` symlink so
  pkg-config can find it.
- **hdWGCNA via `install_local` silently fails.** Use `install_github`. The
  earlier `hpc_defs/seurat_extended.def` used `install_local` with a
  pre-downloaded tarball — that path is buggy and was abandoned in v3.1.
- **`R_LIBS_USER=/dev/null` at runtime.** When the pipeline mounts the user's
  HPC home directory, R picks up a personal library that often contains
  ABI-incompatible binaries from a different conda env. The container
  bakes an `Renviron.site` lockdown in `cicero.sif` (see below); for
  `seurat_extended.sif` the lockdown happens via the runtime env var,
  configured in `nextflow.config` via `containerOptions = '--env R_LIBS_USER=/dev/null'`.
- **basilisk env must be baked at `/opt/basilisk_cache` (v3.5).** zellkonverter
  uses basilisk to manage its anndata conda env. basilisk's default cache
  path is derived from `tools::R_user_dir("basilisk", "cache")`, which under
  the pipeline runtime (`--contain --bind /tmp --env HOME=/tmp
  --env XDG_CACHE_HOME=/tmp/cache`) resolves to somewhere under `/tmp`. The
  `--bind /tmp` mount then shadows whatever was at that path in the image,
  so basilisk effectively sees an empty cache and tries to bootstrap a fresh
  conda env on every job — which (a) requires outbound network access, often
  unavailable on compute nodes, and (b) creates a node-local env in the
  host's `/tmp` that disappears when `/tmp` is cleaned. Symptom: CellChat
  works on some nodes and not others, depending on whether a prior job
  happened to populate that host's `/tmp`.

  The v3.5 fix invokes `basilisk::obtainEnvironmentPath()` during `%post`
  with `BASILISK_EXTERNAL_DIR=/opt/basilisk_cache` set, so the conda env
  materializes inside the image at `/opt/basilisk_cache/` (not shadowed by
  `--bind /tmp`). `%environment` then exports `BASILISK_EXTERNAL_DIR` so
  basilisk uses the baked env at runtime. No more node lottery.

  For belt-and-suspenders, also set this in `nextflow.config` for the
  CellChat / r_seurat / r_cellchat processes:

  ```groovy
  containerOptions = '--env R_LIBS_USER=/dev/null --env BASILISK_EXTERNAL_DIR=/opt/basilisk_cache'
  ```

  (The container's `%environment` already sets this, but an explicit
  `--env` at runtime is harmless and protects against `--cleanenv` /
  `--no-init` flags that strip `%environment` from being applied.)

---

## 4. `cicero.sif`

R container for chromatin co-accessibility (Cicero on Monocle3 backend) plus
Bioconductor genomic-range utilities.

### Definition file

The cicero build is the most fragile of the five because mamba's solver
crashes mid-run on a libxml2 self-upgrade. The crash is *expected* and the
R install is fine afterwards. Translating the build to a clean `.def` is
possible but the imperative form below is what we actually ship:

The full recipe is in [`docs/defs/cicero.def`](./defs/cicero.def). Build with:

```bash
singularity build --fakeroot cicero.sif docs/defs/cicero.def
```

### Pitfalls

- **Mamba's exit code is unreliable** during the big solve because of the
  libxml2 self-upgrade. The build script (`mac_build_containers.sh`)
  swallows the non-zero exit and verifies R is functional with an explicit
  `Rscript -e 'library(rtracklayer); library(Gviz)'` check. Reproduce
  that pattern if you re-tool the build.
- **HPC personal R libraries leak in via `$HOME`.** The fix is twofold:
  (1) bake `Renviron.site` with empty `R_LIBS*` (done at build time —
  see above), and (2) set `R_LIBS_USER=/dev/null` at runtime
  (`nextflow.config: containerOptions = '--env R_LIBS_USER=/dev/null'`).
  Either one alone is insufficient because `Renviron.site` can be
  overridden by `~/.Renviron` if the user happens to have one.
- **`r-ggrastr` is included** even though it's not strictly required by
  Cicero, because several downstream Cicero-derived plots use it.
- **BPCells must be installed before Monocle3** because newer Monocle3
  optionally uses BPCells for sparse matrix backing.

---

## 5. `scenicplus.sif`

SCENIC+ container: pycisTopic + pySCENIC + Mallet LDA + graph-tool.

### Resolved versions (from `next_build.log`)

| Package        | Version              |
| -------------- | -------------------- |
| Python         | 3.11.8 (pinned — scenicplus requires ≤ 3.11.8) |
| scenicplus     | 1.0a2 (commit `840dab8`) |
| pycisTopic     | 2.0a0 (git HEAD)     |
| pyscenic       | 0.12.1+11.g06bafba   |
| pycistarget    | 1.1 (commit `03e886e`) |
| LoomXpy        | 0.4.2 (commit `61995ff`) |
| numpy          | 1.26.4               |
| pandas         | 1.5.0 (pinned by scenicplus) |
| scipy          | 1.12.0               |
| polars         | 0.20.13              |
| Mallet         | 202108 (baked into `/opt/Mallet-202108`) |
| Java JRE       | OpenJDK 17 (default-jre-headless) |
| graph-tool     | via miniforge env at `/opt/miniforge3/envs/gt` |

### Definition file

The full recipe is in [`docs/defs/scenicplus.def`](./defs/scenicplus.def). Build with:

```bash
singularity build --fakeroot scenicplus.sif docs/defs/scenicplus.def
```

### Pitfalls

- **scenicplus pins Python ≤ 3.11.8.** Do not bump the base image to
  `python:3.12-slim`.
- **Two patches must be applied after `pip install`**: pycisTopic's
  gene_annotation.py (BioMart str cast) and gensim's matutils.py
  (scipy.linalg.triu → numpy.triu). Both are baked into the build script;
  the `%test` block verifies the patches survived.
- **graph-tool cannot be pip-installed.** It is a C++ library with Boost/CGAL
  bindings, only distributed via conda-forge. The build installs miniforge
  *without* adding it to PATH, then creates a dedicated `gt` env. The
  system Python (3.11) and the gt env Python (3.11 with numpy 2.4) are
  intentionally separate — pipeline scripts that need graph-tool must use
  `/opt/miniforge3/envs/gt/bin/python` explicitly.
- **Earlier v3.2 included `pygraphistry` + RAPIDS cu12.** This was dropped
  in v3.3 because it conflicted with scenicplus's pinned pandas 1.5 /
  scipy 1.12 / dask 2024.2.1. If you need network visualization, use
  graph-tool's drawing API rather than pygraphistry.
- **Mallet is baked into the image** at `/opt/Mallet-202108` with a
  `/usr/local/bin/mallet` symlink so pycisTopic does not need to download
  it at runtime on compute nodes without outbound network.

---

## Build commands

### From a `.def` file on HPC (canonical path)

The repository's `hpc_defs/BUILD_ON_HPC.sh` is a thin wrapper around
`singularity build --fakeroot docs/defs/<name>.def`. It builds in the
recommended order, captures per-build logs to
`singularity_cache/build_logs/`, falls back to an unprivileged build if
`--fakeroot` is unavailable, and emits a results table plus a SHA256
manifest at the end.

```bash
# Allocate a compute node (do NOT build on a login node)
srun --partition=free --time=04:00:00 --mem=24G --cpus-per-task=8 --pty bash

module load singularity      # or: module load apptainer
cd /path/to/forge

bash hpc_defs/BUILD_ON_HPC.sh all                    # build all five
bash hpc_defs/BUILD_ON_HPC.sh seurat_extended        # one container
bash hpc_defs/BUILD_ON_HPC.sh seurat_extended --rebuild   # force rebuild
bash hpc_defs/BUILD_ON_HPC.sh all --no-test          # skip %test blocks
```

By default the wrapper writes `.sif` files to
`$REPO_DIR/singularity_cache/`. Override with `FORGE_SIF_DIR=/some/path
bash hpc_defs/BUILD_ON_HPC.sh all`. The script skips any container whose
`.sif` already exists; pass `--rebuild` to force.

If `--fakeroot` is unavailable on your cluster, the unprivileged fallback
will still work for most steps but may fail on `apt-get` calls; in that
case build on a workstation with `--fakeroot` and `scp` the `.sif` over.

You can also build one container by invoking `singularity build` directly:

```bash
singularity build --fakeroot singularity_cache/seurat_extended.sif \
    docs/defs/seurat_extended.def
```

### From the Mac builder (convenience, not authoritative)

`mac_build_containers.sh` provisions a Lima VM (Apptainer 1.4.5, Rosetta
x86_64 emulation on Apple Silicon) and runs the same `.def` builds inside
it. Use this when preparing `.sif` files on a laptop:

```bash
chmod +x mac_build_containers.sh

./mac_build_containers.sh                  # build all five (~2-3 hours total)
./mac_build_containers.sh scgpu            # build just scgpu_extended
./mac_build_containers.sh snapatac         # build just snapatac_extended
./mac_build_containers.sh seurat           # build just seurat_extended
./mac_build_containers.sh cicero           # build just cicero
./mac_build_containers.sh scenicplus       # build just scenicplus
```

Output lands in `./sif_output/`; `scp` over to HPC when done.

If a build fails mid-step, the sandbox at
`/tmp/container_builds/<name>_sandbox` (inside the Lima VM) is preserved so
you can re-enter and continue manually:

```bash
limactl shell apptainer-builder
apptainer exec --fakeroot --writable /tmp/container_builds/<name>_sandbox <fix command>
apptainer build sif_output/<name>.sif /tmp/container_builds/<name>_sandbox
```

---

## Post-build steps

### Cache population

| Container             | Bake-time cache                                      | Runtime cache writes                                     |
| --------------------- | ---------------------------------------------------- | -------------------------------------------------------- |
| `scgpu_extended.sif`  | CellTypist models (~500 MB) at `~/.celltypist/data/` | numba JIT (`/tmp/numba_cache`), matplotlib font cache    |
| `snapatac_extended.sif` | None                                                | numba JIT, cupy compile cache (`/tmp/cupy_cache`)        |
| `seurat_extended.sif` | basilisk conda env at `/opt/basilisk_cache/` (v3.5; ~500–800 MB) | None                                                     |
| `cicero.sif`          | None                                                | None                                                     |
| `scenicplus.sif`      | Mallet 202108 at `/opt/Mallet-202108`                | cistarget databases (NOT baked — user supplies)          |

**scPrinter dispersion models** are not baked into `snapatac_extended.sif` —
they are downloaded on first use into `$XDG_CACHE_HOME` (which the pipeline
points at `/tmp/cache`). On a network-isolated compute node you must
pre-populate this cache by running scPrinter once on a node with outbound
network, then copying the cache to a shared path bound into `/tmp/cache`
at runtime.

**cisTarget databases** (`*.feather`, ~5–25 GB per genome) are not in the
scenicplus container; they live outside on `/dfs7` and are bound in at
runtime.

### Runtime bind paths

The pipeline (`nextflow.config`) uses these `singularity.runOptions`:

```
--contain
--bind /dfs7
--bind /tmp
--env NUMBA_CACHE_DIR=/tmp/numba_cache
--env MPLCONFIGDIR=/tmp/matplotlib
--env XDG_CACHE_HOME=/tmp/cache
--env CUPY_CACHE_DIR=/tmp/cupy_cache
--env HOME=/tmp
```

`--contain` suppresses the implicit `$HOME` and `/tmp` mounts plus any
`bind path` directives from `/etc/singularity/singularity.conf`. Everything
the pipeline needs is then explicitly re-added. Critically:

- **`/dfs7`** — the lab's Swarup Lab share on UCI HPC3. Replace with
  your cluster's data path. **You must add it** or the pipeline cannot
  see input or write output.
- **`/tmp`** — required so the cache env vars above resolve to writable
  paths. The pipeline does heavy numba JIT compilation and cupy kernel
  caching here; budget ≥ 20 GB scratch per parallel task on large runs.
- **`HOME=/tmp`** — prevents R/Python from picking up the user's HPC
  home (R personal library, Python user-site, plotly/celltypist configs).

Per-process overrides (also in `nextflow.config`):

- `scgpu_*` processes get `containerOptions = '--nv'` to expose host
  NVIDIA drivers (`libcuda.so.1`).
- All R processes (`r_cicero`, `r_seurat`, `r_cellchat`) get
  `containerOptions = '--env R_LIBS_USER=/dev/null'` to harden the
  HPC-personal-library lockdown.

### Verifying a build

There are three layers of validation, from cheapest to most expensive.

**1. SHA256 manifest.** Every `.def` build inside `BUILD_ON_HPC.sh` ends
with a SHA256 table. Two builds from the same `.def` on the same builder
should produce bit-identical hashes (assuming no upstream package updates
in the meantime). The reference hashes for the published v3.4 builds are:

| File                      |       Size | SHA256                                                             |
| ------------------------- | ---------: | ------------------------------------------------------------------ |
| `cicero.sif`              | 2038841344 | `dc11280027ce23d2fe696cb23e05736330938b0cf79a3618ecaa638fc567f165` |
| `scgpu_extended.sif`      | 3906072576 | `ec11ed947e1d2c4645e5b65a5725627bf37f2527fe89fea81020c0b7cb03cf4d` |
| `snapatac_extended.sif`   | 5101576192 | `2679b4cfc2843d253a34cd6d1822fb612d40588001f03b5efb476f7016bb3e13` |
| `seurat_extended.sif`     | 1131126784 | `c840bbf5f292aca8d1690d32c0aeabab07945468535aca4eeb8c59caa45af06a` |
| `scenicplus.sif`          | 1954689024 | `392e8c77766652a6995741224849053ddf5924aeb2d4c481785399f094991bb7` |

To check on HPC at any time:

```bash
cd /path/to/forge/singularity_cache
for f in cicero.sif scgpu_extended.sif snapatac_extended.sif seurat_extended.sif scenicplus.sif; do
  [ -f "$f" ] || { echo "MISSING: $f"; continue; }
  printf '%-26s  %12s  %s\n' "$f" "$(stat -c %s "$f")" "$(sha256sum "$f" | awk '{print $1}')"
done
```

Sizes should match exactly; hashes will match only if the upstream packages
have not been republished since the reference build. If sizes differ, the
package set has drifted — use step 2 to find out what.

**2. Package manifest.** If hashes don't match, dump installed-package
versions and diff against the resolved-version tables above:

```bash
cd /path/to/forge/singularity_cache
{ for f in scgpu_extended.sif snapatac_extended.sif scenicplus.sif; do
    echo "===== $f ====="
    singularity exec --contain --bind /tmp "$f" python -c \
      "import importlib.metadata as m; [print(f'{d.name}=={d.version}') for d in sorted(m.distributions(), key=lambda x: x.name.lower())]"
  done
  for f in seurat_extended.sif cicero.sif; do
    echo "===== $f ====="
    singularity exec --contain --bind /tmp "$f" Rscript -e \
      'ip <- installed.packages()[,c("Package","Version")]; ip <- ip[order(tolower(ip[,1])),]; for (i in seq_len(nrow(ip))) cat(sprintf("%s==%s\n", ip[i,1], ip[i,2]))'
  done; } > forge_container_audit_$(date +%Y%m%d).txt
```

**3. Functional smoke test.** Confirms the major libraries import and (for
the GPU containers) CUDA is visible:

```bash
module load singularity

singularity exec singularity_cache/scgpu_extended.sif \
    python -c 'import scvi, cellbender, celltypist; print("scgpu OK")'
singularity exec singularity_cache/snapatac_extended.sif \
    python -c 'import snapatac2, scprinter, scATAnno; print("snapatac OK")'
singularity exec singularity_cache/seurat_extended.sif \
    Rscript -e 'library(Seurat); library(hdWGCNA); library(CellChat); cat("seurat OK\n")'
singularity exec singularity_cache/cicero.sif \
    Rscript -e 'library(cicero); library(monocle3); cat("cicero OK\n")'
singularity exec singularity_cache/scenicplus.sif \
    python -c 'import scenicplus, pycisTopic, pyscenic; print("scenicplus OK")'

# GPU containers — confirm CUDA on a GPU node:
srun --gres=gpu:1 --pty singularity exec --nv \
    singularity_cache/scgpu_extended.sif \
    python -c 'import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))'
```

The `%test` block in each `.def` already runs at build time and asserts the
same imports plus the pinning invariants (numpy 1.x in `snapatac_extended`,
`ggrepel == 0.9.6` in `seurat_extended`, both `sed` patches survived in
`scenicplus`, etc.). Re-run with `singularity test <name>.sif` if needed.

---

## Changelog summary (from `mac_build_containers.sh` header)

| Version | Date       | Key change                                                                     |
| ------- | ---------- | ------------------------------------------------------------------------------ |
| v3.0    | 2026-03-11 | Consolidated rebuild; +procps in all containers; scenicplus added to `all`     |
| v3.1    | 2026-03-12 | snapatac: pin cupy 13.x / rmm 24.x; seurat: rebuild magick from source (HDRI)  |
| v3.2    | 2026-03-19 | scgpu: +mofax/scrublet/CellBender@4334e89; snapatac: +pygenometracks/deeptools; scenicplus: +Mallet, +pygraphistry |
| v3.3    | 2026-03-29 | scenicplus: dropped RAPIDS/pygraphistry; +graph-tool via miniforge; pycisTopic + gensim patches |
| v3.4    | 2026-04-08 | snapatac: +scATAnno + harmonypy + leidenalg + python-igraph                    |
| v3.5    | 2026-05-30 | seurat: bake basilisk env at `/opt/basilisk_cache` (fixes `/tmp` shadow lottery for CellChat / zellkonverter) |
