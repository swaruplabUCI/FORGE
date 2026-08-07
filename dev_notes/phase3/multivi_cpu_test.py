#!/usr/bin/env python3
"""
Phase 3 Step 0b — MultiVI CPU fallback, mirroring FORGE exactly.

The first attempt (cpu_fallback_test.py) used the DEPRECATED
MULTIVI.setup_anndata and failed with "Please set up your AnnData with
MULTIVI.setup_anndata first" — a bug in the test, not a CPU problem.

FORGE's bin/run_multivi_integration.py uses MULTIVI.setup_mudata on a MuData with
'rna' and 'atac' modalities, then MULTIVI(...) and .train(adversarial_mixing=True)
with NO accelerator= argument. This reproduces that path verbatim on tiny
synthetic data.

Run inside scgpu_extended.sif WITHOUT --nv so no GPU is visible.
"""
import sys
import numpy as np
import anndata as ad
import mudata as mu
import scipy.sparse as sp
import torch
from scvi.model import MULTIVI

print(f"torch {torch.__version__} | cuda.is_available={torch.cuda.is_available()}  <-- expect False")

rng = np.random.default_rng(0)
N, G, P = 300, 200, 400

# Counts in .X for both modalities (FORGE passes rna_layer=None/atac_layer=None,
# meaning "use .X as count data").
rna = ad.AnnData(sp.csr_matrix(rng.poisson(0.5, size=(N, G)).astype("float32")))
rna.obs_names = [f"cell_{i}" for i in range(N)]
rna.var_names = [f"gene_{j}" for j in range(G)]
rna.obs["sample_id"] = ["s0" if i % 2 else "s1" for i in range(N)]

atac = ad.AnnData(sp.csr_matrix(rng.poisson(0.3, size=(N, P)).astype("float32")))
atac.obs_names = rna.obs_names.copy()
atac.var_names = [f"peak_{k}" for k in range(P)]
atac.obs["sample_id"] = rna.obs["sample_id"].values

mdata = mu.MuData({"rna": rna, "atac": atac})

try:
    MULTIVI.setup_mudata(
        mdata,
        rna_layer=None,
        atac_layer=None,
        protein_layer=None,
        batch_key="sample_id",
        size_factor_key=None,
        categorical_covariate_keys=None,
        continuous_covariate_keys=None,
        idx_layer=None,
        modalities={
            "rna_layer": "rna",
            "atac_layer": "atac",
            "protein_layer": None,
            "batch_key": "rna",
        },
    )
    mvi = MULTIVI(
        mdata,
        modality_weights="equal",
        modality_penalty="Jeffreys",
        n_latent=5,
        region_factors=True,
        gene_likelihood="zinb",
        dispersion="gene",
        use_batch_norm="none",
        use_layer_norm="both",
        latent_distribution="normal",
    )
    mvi.train(max_epochs=2, adversarial_mixing=True)   # no accelerator=, as FORGE
    latent = mvi.get_latent_representation()
    print(f"\n  latent shape: {latent.shape}")
    print("  RESULT: MultiVI trained on CPU  [PASS]")
    sys.exit(0)
except Exception as e:
    print(f"\n  RESULT: FAILED -> {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
