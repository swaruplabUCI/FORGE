#!/usr/bin/env python3
"""
Train scANVI using pre-trained SCVI model  v2.2.0

FIX-33c: Changed early_stopping_monitor from 'elbo_train' to 'elbo_validation'.
  Monitoring training loss defeats early stopping (train loss always decreases).
  Validation loss correctly detects overfitting and stops at the right epoch.
"""

import argparse
import scanpy as sc
import scvi
import os
import torch
from h5ad_compat import sanitize_adata
torch.set_float32_matmul_precision('medium')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prepared_ref', required=True, help='Prepared reference h5ad')
    parser.add_argument('--prepared_query', required=True, help='Prepared query h5ad')
    parser.add_argument('--label_key_file', required=True, help='File containing label key')
    parser.add_argument('--scvi_model_dir', required=True, help='Directory containing trained SCVI model')
    parser.add_argument('--results_dir', required=True)
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--seed', type=int, default=42,
                        help="Random seed for scvi-tools training; scvi.settings.seed "
                             "is None (unseeded) unless set. Seeds torch/numpy/lightning.")
    parser.add_argument('--accelerator', default='auto',
                        choices=['auto', 'gpu', 'cpu'],
                        help="Device for scvi-tools training. 'auto' (default) uses a GPU "
                             "when one is visible and falls back to CPU. 'cpu' forces CPU, "
                             "which is what CPU-only tutorial runs use. Passed straight "
                             "through to SCANVI.train(accelerator=...).")
    
    args = parser.parse_args()

    # Seed before any model construction or training (see --seed help).
    scvi.settings.seed = args.seed
    print(f"scvi.settings.seed = {args.seed}")

    os.makedirs(args.results_dir, exist_ok=True)

    # Read label key
    with open(args.label_key_file) as f:
        SCANVI_LABELS_KEY = f.read().strip()
    
    print(f"Using label key: {SCANVI_LABELS_KEY}")
    
    # Load prepared datasets
    print("Loading prepared reference...")
    ref = sc.read_h5ad(args.prepared_ref)
    
    print("Loading prepared query...")
    adata_query = sc.read_h5ad(args.prepared_query)
    
    # ============================================
    # CRITICAL: LOAD PRE-TRAINED SCVI MODEL
    # ============================================
    print("\n" + "="*60)
    print("STEP 1: LOADING PRE-TRAINED SCVI MODEL")
    print("="*60)
    
    print(f"Loading SCVI model from {args.scvi_model_dir}")
    scvi_ref = scvi.model.SCVI.load(args.scvi_model_dir, adata=ref)
    
    print("Successfully loaded pre-trained SCVI model")
    
    # Add latent representation to reference
    ref.obsm[f'X_scVI_{args.epochs}e'] = scvi_ref.get_latent_representation()
    
    # ============================================
    # STEP 2: Train SCANVI from loaded SCVI model
    # ============================================
    print("\n" + "="*60)
    print("STEP 2: TRAINING SCANVI FROM PRE-TRAINED SCVI")
    print("="*60)
    
    scanvi_ref = scvi.model.SCANVI.from_scvi_model(
        scvi_ref,
        unlabeled_category="Unknown",
        labels_key=SCANVI_LABELS_KEY
    )
    
    print(f"Training SCANVI for {min(20, args.epochs)} epochs "
          f"(accelerator={args.accelerator})...")
    scanvi_ref.train(max_epochs=min(20, args.epochs), n_samples_per_label=100,
                     accelerator=args.accelerator)
    
    scanvi_model_path = os.path.join(args.results_dir, f"scanvi_{SCANVI_LABELS_KEY}_{args.epochs}e")
    scanvi_ref.save(scanvi_model_path, overwrite=True)
    print(f"Saved SCANVI model: {scanvi_model_path}")
    
    # ============================================
    # STEP 3: Query SCANVI with samples
    # ============================================
    print("\n" + "="*60)
    print("STEP 3: QUERYING SCANVI FOR PREDICTIONS")
    print("="*60)
    
    adata_query.obs[SCANVI_LABELS_KEY] = "Unknown"
    
    scanvi_query = scvi.model.SCANVI.load_query_data(
        adata_query,
        scanvi_model_path
    )
    
    print(f"Fine-tuning on query for {args.epochs} epochs...")
    # FIX-33c: Monitor validation loss for early stopping (not training loss).
    train_kwargs = {
        "early_stopping": True,
        "early_stopping_monitor": "elbo_validation",
        "early_stopping_patience": 10,
        "early_stopping_min_delta": 0.001,
        "plan_kwargs": {"weight_decay": 0.0},
    }
    
    scanvi_query.train(
        max_epochs=args.epochs,
        check_val_every_n_epoch=10,
        accelerator=args.accelerator,
        **train_kwargs
    )
    
    # Add predictions
    adata_query.obs['scanvi_prediction'] = scanvi_query.predict()
    adata_query.obsm[f'X_scanvi_{args.epochs}e'] = scanvi_query.get_latent_representation()
    
    output_path = os.path.join(args.results_dir, f'annotated_{args.epochs}e.h5ad')
    sanitize_adata(adata_query, output_path)
    adata_query.write(output_path)
    
    print(f"\nSaved annotated query: {output_path}")
    print(f"Predicted cell types: {adata_query.obs['scanvi_prediction'].nunique()}")
    
    # Print cell type distribution
    print("\nCell type distribution:")
    print(adata_query.obs['scanvi_prediction'].value_counts())

if __name__ == '__main__':
    main()
