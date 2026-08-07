#!/usr/bin/env python3
"""
Phase 3 Step 0 — does FORGE's GPU stack fall back to CPU?

Mirrors the exact call patterns FORGE uses (train() with NO accelerator
argument) on a tiny synthetic dataset. Run inside scgpu_extended.sif WITHOUT
--nv so no GPU is visible; if these succeed, the same code paths work CPU-only.
"""
import sys
import numpy as np

print("=" * 62)
print("torch / CUDA visibility")
print("=" * 62)
import torch
print(f"  torch                 : {torch.__version__}")
print(f"  torch.cuda.is_available: {torch.cuda.is_available()}   <-- expect False")

import anndata as ad
import scanpy as sc
import scipy.sparse as sp

rng = np.random.default_rng(0)
N_CELLS, N_GENES, N_PEAKS = 300, 200, 400


def make_rna():
    X = sp.csr_matrix(rng.poisson(0.5, size=(N_CELLS, N_GENES)).astype("float32"))
    a = ad.AnnData(X)
    a.obs_names = [f"cell_{i}" for i in range(N_CELLS)]
    a.var_names = [f"gene_{j}" for j in range(N_GENES)]
    a.obs["batch"] = ["b0" if i % 2 else "b1" for i in range(N_CELLS)]
    a.obs["cell_type"] = ["A" if i % 3 else "B" for i in range(N_CELLS)]
    a.layers["counts"] = a.X.copy()
    return a


results = {}

# ---- 1. scVI --------------------------------------------------------------
print("\n" + "=" * 62)
print("1. scVI  (FORGE: bin/train_scvi.py -> scvi_ref.train(max_epochs=...))")
print("=" * 62)
try:
    import scvi
    print(f"  scvi-tools: {scvi.__version__}")
    a = make_rna()
    scvi.model.SCVI.setup_anndata(a, layer="counts", batch_key="batch")
    m = scvi.model.SCVI(a, n_latent=5)
    m.train(max_epochs=2)                      # no accelerator= , exactly as FORGE
    _ = m.get_latent_representation()
    print("  RESULT: scVI trained on CPU  [PASS]")
    results["scVI"] = "PASS"
except Exception as e:
    print(f"  RESULT: FAILED -> {type(e).__name__}: {e}")
    results["scVI"] = "FAIL"

# ---- 2. scANVI ------------------------------------------------------------
print("\n" + "=" * 62)
print("2. scANVI  (FORGE: bin/train_scanvi.py)")
print("=" * 62)
try:
    import scvi
    a = make_rna()
    scvi.model.SCVI.setup_anndata(a, layer="counts", batch_key="batch")
    base = scvi.model.SCVI(a, n_latent=5)
    base.train(max_epochs=2)
    sca = scvi.model.SCANVI.from_scvi_model(base, unlabeled_category="Unknown",
                                            labels_key="cell_type")
    sca.train(max_epochs=2)
    print("  RESULT: scANVI trained on CPU  [PASS]")
    results["scANVI"] = "PASS"
except Exception as e:
    print(f"  RESULT: FAILED -> {type(e).__name__}: {e}")
    results["scANVI"] = "FAIL"

# ---- 3. MultiVI -----------------------------------------------------------
print("\n" + "=" * 62)
print("3. MultiVI  (FORGE: bin/run_multivi_integration.py -> mvi.train(...))")
print("=" * 62)
try:
    import scvi
    n_tot = N_GENES + N_PEAKS
    X = sp.csr_matrix(rng.poisson(0.4, size=(N_CELLS, n_tot)).astype("float32"))
    a = ad.AnnData(X)
    a.obs_names = [f"cell_{i}" for i in range(N_CELLS)]
    a.var_names = ([f"gene_{j}" for j in range(N_GENES)]
                   + [f"peak_{k}" for k in range(N_PEAKS)])
    a.var["modality"] = (["Gene Expression"] * N_GENES
                         + ["Peaks"] * N_PEAKS)
    a.obs["sample_id"] = ["s0" if i % 2 else "s1" for i in range(N_CELLS)]
    scvi.model.MULTIVI.setup_anndata(a, batch_key="sample_id")
    mvi = scvi.model.MULTIVI(a, n_genes=N_GENES, n_regions=N_PEAKS, n_latent=5)
    mvi.train(max_epochs=2)
    _ = mvi.get_latent_representation()
    print("  RESULT: MultiVI trained on CPU  [PASS]")
    results["MultiVI"] = "PASS"
except Exception as e:
    print(f"  RESULT: FAILED -> {type(e).__name__}: {e}")
    results["MultiVI"] = "FAIL"

# ---- 4. MOFA+ -------------------------------------------------------------
print("\n" + "=" * 62)
print("4. MOFA+  (FORGE: bin/mofa_int.py)")
print("=" * 62)
try:
    from mofapy2.run.entry_point import entry_point
    ep = entry_point()
    d = [[rng.normal(size=(N_CELLS, N_GENES)).astype("float32")],
         [rng.normal(size=(N_CELLS, N_PEAKS)).astype("float32")]]
    ep.set_data_options(scale_views=False)
    ep.set_data_matrix(d, likelihoods=["gaussian", "gaussian"])
    ep.set_model_options(factors=3)
    ep.set_train_options(iter=3, convergence_mode="fast", gpu_mode=False, verbose=False)
    ep.build()
    ep.run()
    print("  RESULT: MOFA+ ran with gpu_mode=False  [PASS]")
    results["MOFA+"] = "PASS"
except Exception as e:
    print(f"  RESULT: FAILED -> {type(e).__name__}: {e}")
    results["MOFA+"] = "FAIL"

# ---- 5. CellBender CLI reachable + CPU default ---------------------------
print("\n" + "=" * 62)
print("5. CellBender  (FORGE: bin/cellbender_wrapper.py -> base_cli.main)")
print("=" * 62)
try:
    import cellbender
    from cellbender.base_cli import main  # noqa: F401
    print(f"  cellbender import OK ({getattr(cellbender, '__version__', 'n/a')})")
    print("  NOTE: CellBender uses CPU unless --cuda is passed; FORGE's wrapper")
    print("        does not pass it. Treated as CPU-capable.")
    results["CellBender"] = "PASS (import)"
except Exception as e:
    print(f"  RESULT: FAILED -> {type(e).__name__}: {e}")
    results["CellBender"] = "FAIL"

print("\n" + "=" * 62)
print("SUMMARY — CPU-only viability")
print("=" * 62)
for k, v in results.items():
    print(f"  {k:<12} {v}")
bad = [k for k, v in results.items() if v.startswith("FAIL")]
print("\n" + ("ALL CLEAR: T2 can run CPU-only." if not bad
              else f"BLOCKERS: {bad}"))
sys.exit(1 if bad else 0)
