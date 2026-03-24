#!/usr/bin/env python3
"""Train scVI on prepared reference"""

import argparse
import scanpy as sc
import scvi
import os
import torch
torch.set_float32_matmul_precision('medium')

def train_refmap_scvi_model(adata, number_of_epochs, batch_key, results_dir):
    """
    Trains a scVI model with proper batch key handling
    """
    SCVI_LATENT_KEY = f"X_scVI_{number_of_epochs}e"
    
    # Verify batch key exists
    if batch_key not in adata.obs.columns:
        raise ValueError(
            f"Batch key '{batch_key}' not found in adata.obs. "
            f"Available columns: {adata.obs.columns.tolist()}"
        )
    
    print(f"Using batch_key='{batch_key}' with {adata.obs[batch_key].nunique()} unique batches")
    print(f"Batch distribution:\n{adata.obs[batch_key].value_counts()}")
    
    # Ensure counts layer exists
    if 'counts' not in adata.layers:
        if 'UMIs' in adata.layers:
            adata.layers['counts'] = adata.layers['UMIs'].copy()
            print("Created 'counts' layer from 'UMIs'")
        else:
            adata.layers['counts'] = adata.X.copy()
            print("Created 'counts' layer from X matrix")
    
    # Setup SCVI
    print(f"Setting up scVI with batch_key='{batch_key}' and layer='counts'")
    scvi.model.SCVI.setup_anndata(
        adata, 
        batch_key=batch_key,
        layer="counts"
    )
    
    # Initialize SCVI model
    print("Initializing scVI model...")
    scvi_ref = scvi.model.SCVI(
        adata,
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
        n_latent=10,
        gene_likelihood='zinb'
    )
    
    # Train SCVI model
    print(f"Training scVI for {number_of_epochs} epochs...")
    scvi_ref.train(max_epochs=number_of_epochs)
    
    # Save the SCVI model
    scvi_model_path = os.path.join(results_dir, f"scvi_{number_of_epochs}e")
    print(f"Saving scVI model to {scvi_model_path}", flush=True)
    scvi_ref.save(scvi_model_path, overwrite=True, save_anndata=False)
    print(f"  [OK] Model weights saved (anndata excluded from model dir)", flush=True)

    # Store latent representation
    print("Computing latent representation for all cells...", flush=True)
    adata.obsm[SCVI_LATENT_KEY] = scvi_ref.get_latent_representation()
    print(f"  [OK] Latent representation stored in adata.obsm['{SCVI_LATENT_KEY}']", flush=True)

    # Compute neighbors and UMAP
    print("Computing neighbors using scVI latent representation...", flush=True)
    sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
    print("  [OK] Neighbors computed", flush=True)

    print("Computing UMAP...", flush=True)
    sc.tl.umap(adata)
    print("  [OK] UMAP computed", flush=True)

    # Optional: Leiden clustering at multiple resolutions
    resolutions = [0.2, 0.5, 1.0]
    for res in resolutions:
        key = f"leiden_scvi_{res}".replace(".", "_")
        print(f"Computing leiden clustering at resolution {res}...", flush=True)
        sc.tl.leiden(adata, key_added=key, resolution=res)
        n_clusters = len(adata.obs[key].unique())
        print(f"  [OK] {n_clusters} clusters at resolution {res}", flush=True)

    return adata

def main():
    parser = argparse.ArgumentParser(description="Train scVI on prepared reference")
    parser.add_argument('--prepared_ref', required=True, help='Prepared reference h5ad')
    parser.add_argument('--results_dir', required=True)
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--batch_key', default='sample', help='Batch key (default: sample)')
    parser.add_argument('--save_adata', action='store_true', help='Save updated adata with embeddings')
    
    args = parser.parse_args()
    
    os.makedirs(args.results_dir, exist_ok=True)
    
    # Load prepared reference
    print(f"Loading prepared reference from {args.prepared_ref}")
    ref = sc.read_h5ad(args.prepared_ref)
    print(f"Loaded {ref.n_obs} cells and {ref.n_vars} genes")
    
    # Train scVI model
    ref_trained = train_refmap_scvi_model(
        adata=ref,
        number_of_epochs=args.epochs,
        batch_key=args.batch_key,
        results_dir=args.results_dir
    )
    
    # Optionally save reference with SCVI embedding
    if args.save_adata:
        ref_out = os.path.join(args.results_dir, f'reference_scvi_{args.epochs}e.h5ad')
        print(f"Writing reference h5ad with embeddings to {ref_out}...", flush=True)
        ref_trained.write(ref_out)
        print(f"  [OK] Saved reference with SCVI embedding: {ref_out}", flush=True)
    
    print("\n" + "="*60)
    print("SCVI TRAINING COMPLETE!")
    print("="*60)
    print(f"Model saved: {args.results_dir}/scvi_{args.epochs}e")
    print(f"Latent representation: X_scVI_{args.epochs}e")

if __name__ == '__main__':
    main()
