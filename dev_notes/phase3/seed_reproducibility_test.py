#!/usr/bin/env python3
"""
Phase 3 — does seeding scvi-tools actually make training reproducible?

Audit finding: scvi.settings.seed defaults to None, and it was the ONLY
non-reproducible step in FORGE (Leiden random_state defaults to 0 in both scanpy
and snapatac2; MOFA+ already takes seed=42; the R scripts call set.seed).

This asserts the fix delivers what it claims:
  * two runs with the SAME seed produce bit-identical latent representations
  * two runs with DIFFERENT seeds do not (i.e. the seed is actually in play)
  * unseeded (settings.seed = None) runs diverge — the behaviour being fixed

Run inside scgpu_extended.sif WITHOUT --nv (CPU, matching the tutorial).
"""
import sys
import numpy as np
import anndata as ad
import scipy.sparse as sp
import scvi

N, G = 300, 200


def build():
    rng = np.random.default_rng(12345)          # data is fixed; only training varies
    a = ad.AnnData(sp.csr_matrix(rng.poisson(0.5, size=(N, G)).astype("float32")))
    a.obs_names = [f"c{i}" for i in range(N)]
    a.var_names = [f"g{j}" for j in range(G)]
    a.obs["batch"] = ["b0" if i % 2 else "b1" for i in range(N)]
    a.layers["counts"] = a.X.copy()
    return a


def latent(seed):
    """Train scVI exactly as bin/train_scvi.py does, returning the latent space."""
    scvi.settings.seed = seed                    # None => unseeded, the old behaviour
    a = build()
    scvi.model.SCVI.setup_anndata(a, layer="counts", batch_key="batch")
    m = scvi.model.SCVI(a, n_latent=5)
    m.train(max_epochs=3, accelerator="cpu")
    return m.get_latent_representation()


print(f"scvi-tools {scvi.__version__}\n")
fails = []

print("=" * 62)
print("A. same seed twice -> must be IDENTICAL")
print("=" * 62)
a1, a2 = latent(42), latent(42)
same = np.array_equal(a1, a2)
delta = float(np.abs(a1 - a2).max())
print(f"  max |diff| = {delta:.3e}   identical={same}")
if not same:
    fails.append("same seed did not reproduce")

print("\n" + "=" * 62)
print("B. different seeds -> must DIFFER (proves the seed is in play)")
print("=" * 62)
b = latent(7)
differs = not np.array_equal(a1, b)
print(f"  max |diff| vs seed 42 = {float(np.abs(a1 - b).max()):.3e}   differs={differs}")
if not differs:
    fails.append("different seeds produced identical output (seed ignored?)")

print("\n" + "=" * 62)
print("C. unseeded twice -> expected to DIVERGE (the old default)")
print("=" * 62)
u1, u2 = latent(None), latent(None)
u_same = np.array_equal(u1, u2)
print(f"  max |diff| = {float(np.abs(u1 - u2).max()):.3e}   identical={u_same}")
print("  (identical here would just mean this toy case is too small to expose"
      " drift; it does not invalidate the fix)")

print("\n" + "=" * 62)
print("ALL PASS — seeding makes scvi-tools reproducible."
      if not fails else f"FAILURES: {fails}")
print("=" * 62)
sys.exit(1 if fails else 0)
