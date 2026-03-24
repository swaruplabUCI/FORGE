#!/usr/bin/env python3
"""
gene_coordinate_resolver.py

Resolves gene symbols to genomic coordinates using GTF annotations.
Supports both human (hg38) and mouse (mm10) references.
"""

import argparse
import json
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import re
from typing import Dict, List, Tuple, Optional

def parse_gtf_attributes(attr_string: str) -> dict:
    """Parse GTF attribute string into dictionary."""
    attrs = {}
    for match in re.finditer(r'(\w+) "([^"]+)"', attr_string):
        attrs[match.group(1)] = match.group(2)
    return attrs

def load_gtf_genes(gtf_path: str) -> pd.DataFrame:
    """Load gene information from GTF file."""
    print(f"Loading GTF from {gtf_path}...")
    
    genes = []
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            if parts[2] == 'gene':
                attrs = parse_gtf_attributes(parts[8])
                
                # Extract gene symbol (try multiple attribute names)
                gene_symbol = None
                for attr_name in ['gene_name', 'gene_symbol', 'Name']:
                    if attr_name in attrs:
                        gene_symbol = attrs[attr_name]
                        break
                
                if not gene_symbol:
                    continue
                
                genes.append({
                    'gene_symbol': gene_symbol,
                    'gene_id': attrs.get('gene_id', ''),
                    'chr': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'gene_type': attrs.get('gene_type', attrs.get('gene_biotype', 'unknown'))
                })
    
    df = pd.DataFrame(genes)
    print(f"Loaded {len(df)} genes from GTF")
    
    # Prioritize protein-coding genes
    df['priority'] = df['gene_type'].apply(
        lambda x: 0 if x == 'protein_coding' else 1
    )
    
    return df

def get_tss(row: pd.Series) -> int:
    """Get TSS position based on strand."""
    if row['strand'] == '+':
        return row['start']
    else:
        return row['end']

def resolve_gene_coordinates(
    gene_symbols: List[str],
    gtf_df: pd.DataFrame,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 2000,
    enhancer_max_distance: int = 500000
) -> Dict:
    """Resolve gene symbols to genomic coordinates."""
    
    results = {}
    
    for gene in gene_symbols:
        # Find gene in GTF (case-insensitive)
        matches = gtf_df[gtf_df['gene_symbol'].str.upper() == gene.upper()]
        
        if matches.empty:
            print(f"WARNING: Gene '{gene}' not found in GTF")
            results[gene] = {
                'found': False,
                'error': f"Gene not found in reference"
            }
            continue
        
        # If multiple matches, prefer protein-coding
        matches = matches.sort_values('priority')
        best_match = matches.iloc[0]
        
        tss = get_tss(best_match)
        
        # Define promoter region
        if best_match['strand'] == '+':
            promoter_start = max(0, tss - promoter_upstream)
            promoter_end = tss + promoter_downstream
        else:
            promoter_start = max(0, tss - promoter_downstream)
            promoter_end = tss + promoter_upstream
        
        # Define enhancer search region
        enhancer_start = max(0, tss - enhancer_max_distance)
        enhancer_end = tss + enhancer_max_distance
        
        results[gene] = {
            'found': True,
            'gene_id': best_match['gene_id'],
            'chr': best_match['chr'],
            'strand': best_match['strand'],
            'gene_type': best_match['gene_type'],
            'tss': tss,
            'gene_start': best_match['start'],
            'gene_end': best_match['end'],
            'promoter': {
                'chr': best_match['chr'],
                'start': promoter_start,
                'end': promoter_end
            },
            'enhancer_search': {
                'chr': best_match['chr'],
                'start': enhancer_start,
                'end': enhancer_end
            }
        }
        
        if len(matches) > 1:
            results[gene]['warning'] = f"Multiple matches found, using {best_match['gene_type']} gene"
    
    return results

def merge_with_manual_coordinates(auto_coords: Dict, manual_coords: Dict) -> Dict:
    """Merge automatic coordinates with manual overrides."""
    
    merged = auto_coords.copy()
    
    for gene, coords in manual_coords.items():
        if gene in merged and merged[gene]['found']:
            # Override with manual coordinates
            if 'promoter' in coords:
                merged[gene]['promoter'].update(coords['promoter'])
                merged[gene]['manual_override'] = True
            if 'enhancer_search' in coords:
                merged[gene]['enhancer_search'].update(coords['enhancer_search'])
        else:
            # Add manual entry for gene not in GTF
            merged[gene] = {
                'found': True,
                'manual_entry': True,
                'promoter': coords.get('promoter', {}),
                'enhancer_search': coords.get('enhancer_search', {})
            }
    
    return merged

def convert_numpy_types(obj):
    """Convert numpy types to native Python types for JSON serialization"""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    else:
        return obj


def main():
    parser = argparse.ArgumentParser(
        description="Resolve gene symbols to genomic coordinates"
    )
    parser.add_argument(
        '--genes',
        nargs='+',
        required=True,
        help='Gene symbols to resolve'
    )
    parser.add_argument(
        '--species',
        required=True,
        choices=['human', 'mouse'],
        help='Species'
    )
    parser.add_argument(
        '--gtf-human',
        default=None,
        help='Path to human GTF'
    )
    parser.add_argument(
        '--gtf-mouse', 
        default=None,
        help='Path to mouse GTF'
    )
    parser.add_argument(
        '--manual-coords',
        help='JSON file with manual coordinate overrides'
    )
    parser.add_argument(
        '--promoter-upstream',
        type=int,
        default=2000,
        help='Upstream distance from TSS for promoter (default: 2000)'
    )
    parser.add_argument(
        '--promoter-downstream',
        type=int,
        default=2000,
        help='Downstream distance from TSS for promoter (default: 2000)'
    )
    parser.add_argument(
        '--output',
        default='gene_coordinates.json',
        help='Output JSON file'
    )
    parser.add_argument(
        '--format',
        choices=['json', 'tsv', 'both'],
        default='both',
        help='Output format'
    )
    
    args = parser.parse_args()
    
    # Load appropriate GTF
    gtf_path = args.gtf_human if args.species == 'human' else args.gtf_mouse
    if not Path(gtf_path).exists():
        print(f"ERROR: GTF file not found: {gtf_path}")
        sys.exit(1)
    
    gtf_df = load_gtf_genes(gtf_path)
    
    # Resolve coordinates
    auto_coords = resolve_gene_coordinates(
        args.genes,
        gtf_df,
        args.promoter_upstream,
        args.promoter_downstream
    )
    
    # Load manual overrides if provided
    if args.manual_coords:
        with open(args.manual_coords, 'r') as f:
            manual_coords = json.load(f)
        coords = merge_with_manual_coordinates(auto_coords, manual_coords)
    else:
        coords = auto_coords
    
    # Report results
    print("\n" + "="*60)
    print("GENE COORDINATE RESOLUTION RESULTS")
    print("="*60)
    
    for gene, info in coords.items():
        if info['found']:
            print(f"\n{gene}:")
            if 'manual_override' in info:
                print("  [MANUAL OVERRIDE APPLIED]")
            if 'manual_entry' in info:
                print("  [MANUAL ENTRY]")
            else:
                print(f"  Gene ID: {info.get('gene_id', 'N/A')}")
                print(f"  Type: {info.get('gene_type', 'N/A')}")
                print(f"  Strand: {info.get('strand', 'N/A')}")
                print(f"  TSS: {info['chr']}:{info.get('tss', 'N/A')}")
            
            if 'promoter' in info:
                p = info['promoter']
                print(f"  Promoter: {p['chr']}:{p['start']}-{p['end']}")
                print(f"    Length: {p['end'] - p['start']} bp")
            
            if 'warning' in info:
                print(f"  ⚠ {info['warning']}")
        else:
            print(f"\n{gene}: NOT FOUND - {info.get('error', '')}")
    
    # Save results
    if args.format in ['json', 'both']:
        with open(args.output, 'w') as f:
            coords_converted = convert_numpy_types(coords)
            json.dump(coords_converted, f, indent=2)
        print(f"\nCoordinates saved to {args.output}")
    
    if args.format in ['tsv', 'both']:
        tsv_output = args.output.replace('.json', '.tsv')
        rows = []
        for gene, info in coords.items():
            if info['found'] and 'promoter' in info:
                rows.append({
                    'gene': gene,
                    'chr': info['promoter']['chr'],
                    'start': info['promoter']['start'],
                    'end': info['promoter']['end'],
                    'strand': info.get('strand', ''),
                    'tss': info.get('tss', ''),
                    'gene_type': info.get('gene_type', '')
                })
        
        if rows:
            df = pd.DataFrame(rows)
            df.to_csv(tsv_output, sep='\t', index=False)
            print(f"TSV saved to {tsv_output}")

if __name__ == '__main__':
    main()
