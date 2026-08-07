#!/usr/bin/env python3
"""
Phase 3 Step 0c — confirm explicit accelerator='cpu' / 'auto' work.

FORGE's bin/train_scvi.py, train_scanvi.py and run_multivi_integration.py now
take --accelerator and pass it to .train(accelerator=...). This asserts that both
'cpu' and 'auto' are accepted by scvi-tools 1.4.2 for all three model classes, so
the CPU-only tutorial run is explicit rather than relying on auto-detection.

Run inside scgpu_extended.sif WITHOUT --nv.
"""
import sys
import numpy as np
import anndata as ad
import mudata as mu
import scipy.sparse as sp
import scvi
import torch

print(f"scvi-tools {scvi.__version__} | torch {torch.__version__} | "
      f"cuda={torch.cuda.is_available()}")

rng = np.random.default_rng(0)
N, G, P = 200, 120, 200
fails = []


def rna_adata():
    a = ad.AnnData(sp.csr_matrix(rng.poisson(0.5, size=(N, G)).astype("float32")))
    a.obs_names = [f"c{i}" for i in range(N)]
    a.var_names = [f"g{j}" for j in range(G)]
    a.obs["batch"] = ["b0" if i % 2 else "b1" for i in range(N)]
    a.obs["cell_type"] = ["A" if i % 3 else "B" for i in range(N)]
    a.layers["counts"] = a.X.copy()
    return a


for acc in ("cpu", "auto"):
    print(f"\n{'=' * 58}\naccelerator = {acc!r}\n{'=' * 58}")

    # scVI
    try:
        a = rna_adata()
        scvi.model.SCVI.setup_anndata(a, layer="counts", batch_key="batch")
        m = scvi.model.SCVI(a, n_latent=5)
        m.train(max_epochs=2, accelerator=acc)
        print(f"  scVI    [PASS]")
    except Exception as e:
        print(f"  scVI    [FAIL] {type(e).__name__}: {e}")
        fails.append(f"scVI/{acc}")

    # scANVI
    try:
        a = rna_adata()
        scvi.model.SCVI.setup_anndata(a, layer="counts", batch_key="batch")
        base = scvi.model.SCVI(a, n_latent=5)
        base.train(max_epochs=2, accelerator=acc)
        sca = scvi.model.SCANVI.from_scvi_model(
            base, unlabeled_category="Unknown", labels_key="cell_type")
        sca.train(max_epochs=2, n_samples_per_label=10, accelerator=acc)
        print(f"  scANVI  [PASS]")
    except Exception as e:
        print(f"  scANVI  [FAIL] {type(e).__name__}: {e}")
        fails.append(f"scANVI/{acc}")

    # MultiVI — FORGE's setup_mudata path
    try:
        r = rna_adata()
        at = ad.AnnData(sp.csr_matrix(rng.poisson(0.3, size=(N, P)).astype("float32")))
        at.obs_names = r.obs_names.copy()
        at.var_names = [f"p{k}" for k in range(P)]
        at.obs["sample_id"] = ["s0" if i % 2 else "s1" for i in range(N)]
        r.obs["sample_id"] = at.obs["sample_id"].values
        md = mu.MuData({"rna": r, "atac": at})
        scvi.model.MULTIVI.setup_mudata(
            md, rna_layer=None, atac_layer=None, protein_layer=None,
            batch_key="sample_id",
            modalities={"rna_layer": "rna", "atac_layer": "atac",
                        "protein_layer": None, "batch_key": "rna"},
        )
        mvi = scvi.model.MULTIVI(md, n_latent=5)
        mvi.train(max_epochs=2, adversarial_mixing=True, accelerator=acc)
        print(f"  MultiVI [PASS]")
    except Exception as e:
        print(f"  MultiVI [FAIL] {type(e).__name__}: {e}")
        fails.append(f"MultiVI/{acc}")

print(f"\n{'=' * 58}")
print("ALL PASS — explicit accelerator works for cpu and auto."
      if not fails else f"FAILURES: {fails}")
sys.exit(1 if fails else 0)
