#!/usr/bin/env python3
"""
Barcode Matching Utilities for Multiome Integration
Handles cross-modality barcode alignment for RNA+ATAC pairing
"""

import pandas as pd
import numpy as np
import re
from pathlib import Path
import sys

def standardize_barcode(bc: str, sample_prefix: str = None) -> str:
    """
    Standardize barcode format across modalities
    
    Examples:
        RNA:  'AAACGAAAGAACCATA-1' or '90plus_L1_SMK_june_filtered_S16:AAACGAAAGAACCATA-1'
        ATAC: 'AAACGAAAGAACCATA-1' or 'S16_june:AAACGAAAGAACCATA-1'
        
    Returns: 'SAMPLE:AAACGAAAGAACCATA-1' format
    """
    # Remove any existing sample prefix
    if ':' in bc:
        parts = bc.split(':')
        bc_only = parts[-1]
    else:
        bc_only = bc
    
    # Remove trailing -1 if present
    bc_clean = re.sub(r'-\d+$', '', bc_only)
    
    # Add sample prefix if provided
    if sample_prefix:
        return f"{sample_prefix}:{bc_clean}"
    else:
        return bc_clean

def extract_sample_id_from_path(file_path: str, modality: str = 'rna') -> str:
    """
    Extract sample ID from file path
    
    RNA example: '90plus_L1_SMK_june_filtered_S16.h5ad' → 'S16_june'
    ATAC example: 'S16_june.h5ad' → 'S16_june'
    """
    stem = Path(file_path).stem
    
    if modality == 'rna':
        # Pattern: *_filtered_S16 or *_filtered_L1_Donor_0
        match = re.search(r'_filtered_(S\d+|L\d+_Donor_\d+)$', stem)
        if match:
            sample_id = match.group(1)
            # Detect batch
            if 'june' in stem.lower():
                return f"{sample_id}_june"
            elif 'july' in stem.lower():
                return f"{sample_id}_july"
            elif 'nov' in stem.lower():
                return f"{sample_id}_nov"
            return sample_id
    
    elif modality == 'atac':
        # Pattern: S16_june or L1_Donor_0_july
        if re.match(r'(S\d+|L\d+_Donor_\d+)_(june|july|nov)', stem):
            return stem
    
    # Fallback
    return stem

def match_rna_atac_barcodes(rna_barcodes: list, atac_barcodes: list, 
                            sample_id: str, method: str = 'direct') -> tuple:
    """
    Match RNA and ATAC barcodes for the same sample
    
    Args:
        rna_barcodes: List of RNA cell barcodes
        atac_barcodes: List of ATAC cell barcodes
        sample_id: Sample identifier
        method: 'direct' (10x barcodes) or 'positional' (BD Rhapsody)
        
    Returns:
        (common_barcodes, rna_indices, atac_indices)
    """
    if method == 'direct':
        # Standard 10x: barcodes should match directly
        rna_clean = [standardize_barcode(bc) for bc in rna_barcodes]
        atac_clean = [standardize_barcode(bc) for bc in atac_barcodes]
        
        common = sorted(list(set(rna_clean) & set(atac_clean)))
        
        if len(common) == 0:
            print(f"  WARNING: No direct barcode overlap for {sample_id}", file=sys.stderr)
            return [], [], []
        
        # Map back to original indices
        rna_map = {standardize_barcode(bc): i for i, bc in enumerate(rna_barcodes)}
        atac_map = {standardize_barcode(bc): i for i, bc in enumerate(atac_barcodes)}
        
        rna_idx = [rna_map[bc] for bc in common]
        atac_idx = [atac_map[bc] for bc in common]
        
        # Create globally unique barcodes
        global_barcodes = [f"{sample_id}:{bc}" for bc in common]
        
        return global_barcodes, rna_idx, atac_idx
    
    elif method == 'positional':
        # BD Rhapsody: cells are matched by position
        n_cells = min(len(rna_barcodes), len(atac_barcodes))
        
        global_barcodes = [f"{sample_id}:cell_{i:06d}" for i in range(n_cells)]
        rna_idx = list(range(n_cells))
        atac_idx = list(range(n_cells))
        
        return global_barcodes, rna_idx, atac_idx
    
    else:
        raise ValueError(f"Unknown matching method: {method}")

def validate_multiome_overlap(rna_n: int, atac_n: int, overlap_n: int, 
                              min_retention: float = 0.5) -> bool:
    """
    Validate that RNA/ATAC overlap meets minimum threshold
    """
    if overlap_n == 0:
        return False
    
    retention = overlap_n / min(rna_n, atac_n)
    return retention >= min_retention

def main():
    """Test barcode matching utilities"""
    # Test standardization
    test_barcodes = [
        "AAACGAAAGAACCATA-1",
        "S16_june:AAACGAAAGAACCATA-1",
        "90plus_L1_SMK_june_filtered_S16:AAACGAAAGAACCATA-1"
    ]
    
    print("Barcode Standardization Test:")
    for bc in test_barcodes:
        std = standardize_barcode(bc, sample_prefix="S16_june")
        print(f"  {bc} → {std}")
    
    # Test sample ID extraction
    test_paths = [
        ("90plus_L1_SMK_june_filtered_S16.h5ad", "rna"),
        ("S16_june.h5ad", "atac"),
        ("90plus_WTA_10plex_L1_nov_filtered_L1_Donor_0.h5ad", "rna")
    ]
    
    print("\nSample ID Extraction Test:")
    for path, modality in test_paths:
        sid = extract_sample_id_from_path(path, modality)
        print(f"  {path} ({modality}) → {sid}")

if __name__ == "__main__":
    main()
