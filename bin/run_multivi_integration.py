#!/usr/bin/env python

import argparse
import muon as mu
import numpy as np
import scvi
from scvi.model import MULTIVI


def parse_args():
    p = argparse.ArgumentParser(
        description="Run scvi-tools MultiVI integration on a MuData object",
    )
    p.add_argument("--mudata", required=True, help="Input MuData .h5mu file")
    p.add_argument("--out_model", required=True, help="Output directory for saved MultiVI model")
    p.add_argument("--out_mudata", required=True, help="Output MuData .h5mu with MultiVI annotations")
    p.add_argument("--n_epochs", type=int, default=200, help="Number of training epochs")
    p.add_argument("--batch_key", type=str, default="sample_id", help="obs batch key in MuData")
    p.add_argument(
        "--modality_weights",
        type=str,
        default="equal",
        choices=["equal", "cell", "universal"],
        help="MultiVI modality weighting scheme",
    )
    p.add_argument(
        "--modality_penalty",
        type=str,
        default="Jeffreys",
        choices=["Jeffreys", "MMD", "None"],
        help="MultiVI modality alignment penalty",
    )
    p.add_argument(
        "--latent_dim",
        type=int,
        default=20,
        help="Dimensionality of the latent space",
    )
    return p.parse_args()


def main():
    args = parse_args()

    # Load MuData produced by BUILD_MUDATA
    mdata = mu.read_h5mu(args.mudata)

    # Sanity checks on modalities
    if "rna" not in mdata.mod or "atac" not in mdata.mod:
        raise ValueError(
            "Expected MuData modalities 'rna' and 'atac' but found %s" % list(mdata.mod.keys())
        )

    # Register MuData with MultiVI using setup_mudata
    # rna_layer / atac_layer None => use .X for those modalities as count data
    MULTIVI.setup_mudata(
        mdata,
        rna_layer=None,
        atac_layer=None,
        protein_layer=None,
        batch_key=args.batch_key,
        size_factor_key=None,
        categorical_covariate_keys=None,
        continuous_covariate_keys=None,
        idx_layer=None,
        modalities={
            "rna_layer": "rna",
            "atac_layer": "atac",
            "protein_layer": None,
            "batch_key": "rna",  # batch lives in obs of the 'rna' modality
        },
    )

    # Instantiate MultiVI model with MuData
    mvi = MULTIVI(
        mdata,
        modality_weights=args.modality_weights,
        modality_penalty=args.modality_penalty,
        n_latent=args.latent_dim,
        region_factors=True,
        gene_likelihood="zinb",
        dispersion="gene",
        use_batch_norm="none",
        use_layer_norm="both",
        latent_distribution="normal",
    )

    # Train the model with adversarial mixing to align modalities
    mvi.train(
        max_epochs=args.n_epochs,
        adversarial_mixing=True,
    )

    # Store joint latent representation back into MuData
    mdata.obsm["X_MultiVI"] = mvi.get_latent_representation(modality="joint")

    # Normalized expression for RNA modality
    expr = mvi.get_normalized_expression(return_numpy=True)
    if expr.shape[0] != mdata.mod["rna"].n_obs:
        raise RuntimeError("Normalized expression cells do not match RNA cells")
    mdata.mod["rna"].layers["multivi_norm"] = expr

    # Normalized accessibility probabilities for ATAC modality
    acc = mvi.get_normalized_accessibility(return_numpy=True)
    if acc.shape[0] != mdata.mod["atac"].n_obs:
        raise RuntimeError("Normalized accessibility cells do not match ATAC cells")
    mdata.mod["atac"].layers["multivi_accessibility"] = acc

    # Optional: store library size factors in .obs
    lib_sizes = mvi.get_library_size_factors()
    if lib_sizes is not None:
        mdata.obs["multivi_libsize_expr"] = lib_sizes["expression"]
        mdata.obs["multivi_libsize_acc"] = lib_sizes["accessibility"]

    # Save model and updated MuData
    mvi.save(args.out_model, overwrite=True)
    mdata.write_h5mu(args.out_mudata)


if __name__ == "__main__":
    main()
