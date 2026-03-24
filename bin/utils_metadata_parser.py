#!/usr/bin/env python3
"""
Smart metadata parser - detects modalities, conditions, and study design
"""

import pandas as pd
import sys
from pathlib import Path

class MetadataParser:
    def __init__(self, metadata_path):
        self.df = pd.read_csv(metadata_path)
        self.study_design = self._detect_study_design()
    
    def _detect_study_design(self):
        """Detect available modalities and experimental design"""
        design = {
            'has_rna': False,
            'has_atac': False,
            'has_conditions': False,
            'is_comparative': False,
            'is_multiome': False,
            'modality': None,
            'condition_key': None,
            'batch_key': None,
            'cell_type_key': None
        }
        
        # Detect modalities
        if 'rna_file' in self.df.columns or 'wta_file' in self.df.columns:
            design['has_rna'] = True
        
        if 'fragment_file' in self.df.columns or 'atac_file' in self.df.columns:
            design['has_atac'] = True
        
        # Determine modality mode
        if design['has_rna'] and design['has_atac']:
            design['is_multiome'] = True
            design['modality'] = 'multiome'
        elif design['has_rna']:
            design['modality'] = 'rna_only'
        elif design['has_atac']:
            design['modality'] = 'atac_only'
        else:
            raise ValueError("No RNA or ATAC files detected in metadata!")
        
        # Detect experimental conditions (for comparative analysis)
        condition_candidates = ['condition', 'group', 'treatment', 'genotype', 'disease_state']
        for col in condition_candidates:
            if col in self.df.columns:
                unique_vals = self.df[col].nunique()
                if unique_vals >= 2:  # At least 2 groups for comparison
                    design['has_conditions'] = True
                    design['is_comparative'] = True
                    design['condition_key'] = col
                    break
        
        # Detect batch key
        batch_candidates = ['batch', 'sample_batch', 'sequencing_batch', 'lane']
        for col in batch_candidates:
            if col in self.df.columns:
                design['batch_key'] = col
                break
        
        # Detect cell type annotations (if pre-existing)
        celltype_candidates = ['cell_type', 'celltype', 'cluster', 'annotation']
        for col in celltype_candidates:
            if col in self.df.columns:
                design['cell_type_key'] = col
                break
        
        return design
    
    def get_workflow_mode(self):
        """Return recommended workflow configuration"""
        d = self.study_design
        
        mode = {
            'run_rna': d['has_rna'],
            'run_atac': d['has_atac'],
            'run_multiome_integration': d['is_multiome'],
            'run_differential_atac': d['is_comparative'] and d['has_atac'],
            'run_differential_rna': d['is_comparative'] and d['has_rna'],
            'run_cicero': d['has_atac'],
            'run_chromvar': d['has_atac'],
            'run_scprint': d['has_atac'],
            'run_cellchat': d['has_rna'],
            'run_hdwgcna': d['has_rna'],
            'run_mofa': d['is_multiome'],
            'condition_key': d['condition_key'],
            'batch_key': d['batch_key'] or 'sample',
            'modality': d['modality']
        }
        
        return mode
    
    def validate_files_exist(self, base_dirs):
        """Validate that files referenced in metadata actually exist"""
        missing = []
        
        for idx, row in self.df.iterrows():
            for col in ['rna_file', 'fragment_file', 'atac_file', 'wta_file']:
                if col in row and pd.notna(row[col]):
                    file_found = False
                    for base_dir in base_dirs:
                        full_path = Path(base_dir) / row[col]
                        if full_path.exists():
                            file_found = True
                            break
                    
                    if not file_found:
                        missing.append(f"Row {idx}: {col}={row[col]}")
        
        return missing

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python metadata_parser.py <metadata.csv>")
        sys.exit(1)
    
    parser = MetadataParser(sys.argv[1])
    mode = parser.get_workflow_mode()
    
    import json
    print(json.dumps(mode, indent=2))
