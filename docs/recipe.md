# Enhancer-Resolved Multiscale Footprinting: Agent Recipe

## Purpose

Three protocols for characterizing transcription factor binding at putative enhancer elements using multiscale footprinting (scprinter/PRINT/seq2PRINT). Recipe A uses ATAC data only. Recipe B extends into multiome (RNA+ATAC) where paired expression provides additional evidence for enhancer-gene linkage, TF activity, and biological interpretation. Recipe C integrates cell-cell communication inference (CellChat) with the regulatory footprinting framework, tracing intercellular signals from ligand-receptor expression through to chromatin-level evidence of signal reception in target cells.

All recipes are biology-agnostic. They yield structured results — per-cell-state, per-TF, per-enhancer binding evidence — regardless of tissue or disease context.

---

## Shared Prerequisites

```
# Python environment
scprinter (scp), scanpy (sc), anndata, pandas, numpy, matplotlib, pyBigWig

# R environment (Recipe B extends into R for hdWGCNA)
Seurat, hdWGCNA, WGCNA, Cicero/Signac (if using R-based CCAN calling)

# Data
- fragments.tsv.gz (10x cellranger-atac or cellranger-arc output)
- barcodes with cell-type/cluster annotations
- genome reference: scp.genome.hg38 or scp.genome.mm10
```

---

## RECIPE A: ATAC-Only Enhancer Footprinting

### Conceptual Flow

```
Fragments
  -> Peak Calling
  -> chromVAR (TF motif enrichment per cell)
  -> Cicero CCANs (co-accessibility networks)
  -> Extract distal enhancer peaks from CCANs
  -> Motif scan enhancer peaks for chromVAR-nominated TFs
  -> scprinter footprinting at motif-containing enhancers
  -> Interpret: per-cell-state binding evidence at enhancers
```

### Step 1: Import Fragments and Call Peaks

```python
import scprinter as scp
import pandas as pd
import numpy as np

genome = scp.genome.hg38  # or mm10

printer = scp.pp.import_fragments(
    path_to_frags='fragments.tsv.gz',
    barcodes=barcodes,           # array of valid barcode strings
    savename='project.h5ad',
    genome=genome,
    sorted_by_barcode=False
)

# Two peak sets: one for seq2PRINT training, one for chromVAR
scp.pp.call_peaks(printer=printer, preset='seq2PRINT', group_names=['all'])
scp.pp.call_peaks(printer=printer, preset='chromvar', group_names=['chromvar_all'])
```

**Output:** printer object with peak catalogs stored on disk.

### Step 2: chromVAR — Identify Candidate TFs per Cell State

```python
import scanpy as sc

# Build peak-by-cell matrix
adata_peaks = scp.pp.make_peak_matrix(printer, regions='chromvar_peaks.bed')

# Compute GC bias and sample background peaks
scp.chromvar.get_bias(adata_peaks, genome)
scp.chromvar.sample_bg_peaks(adata_peaks, method='chromvar', niterations=250)

# Scan peaks for motif occurrences
motifs = scp.motifs.FigR_Human_Motifs(genome)  # or custom motif set
motifs.chromvar_scan(adata_peaks)

# Compute deviation z-scores (GPU accelerated)
chromvar_results = scp.chromvar.compute_deviations(adata_peaks, device='cuda')

# Cluster on chromVAR deviations
sc.tl.pca(chromvar_results)
sc.pp.neighbors(chromvar_results)
sc.tl.leiden(chromvar_results, resolution=0.5)
sc.tl.umap(chromvar_results)
```

**Output:** `chromvar_results` AnnData where `.X` contains per-cell TF deviation z-scores, `.obs['leiden']` contains cluster labels.

**Decision point:** Identify the top N TFs enriched per cluster. These are your candidate TFs for footprinting.

```python
# Example: extract top 5 TFs per cluster by mean deviation
import pandas as pd

cluster_tf_means = pd.DataFrame(
    chromvar_results.X,
    index=chromvar_results.obs['leiden'],
    columns=chromvar_results.var_names
).groupby(level=0).mean()

top_tfs_per_cluster = {}
for cluster in cluster_tf_means.index:
    top_tfs_per_cluster[cluster] = (
        cluster_tf_means.loc[cluster]
        .sort_values(ascending=False)
        .head(5).index.tolist()
    )
```

### Step 3: Cicero CCANs — Identify Putative Enhancer-Promoter Links

This step runs in R (Cicero/Signac) or can use ArchR. The output needed is a table of co-accessible peak pairs with correlation scores and CCAN group assignments.

```r
# --- R code (Cicero via Signac or standalone) ---
library(Signac)
library(Cicero)  # or use cicero through monocle3

# Input: peak-by-cell matrix, cell UMAP coordinates
# Run Cicero to get co-accessibility scores
conns <- run_cicero(cds, genomic_coords)  # returns peak-peak connections with coaccess scores

# Filter for significant connections
sig_conns <- conns[conns$coaccess > 0.25, ]

# Generate CCANs (groups of co-accessible peaks)
ccans <- generate_ccans(sig_conns, coaccess_cutoff = 0.25)

# Export: save connections and CCAN assignments
write.csv(sig_conns, 'cicero_connections.csv')
write.csv(ccans, 'cicero_ccans.csv')
```

**Output:** `cicero_connections.csv` (peak1, peak2, coaccess score) and `cicero_ccans.csv` (peak, CCAN_id).

### Step 4: Extract Enhancer Regions from CCANs

Back in Python. The logic: for each CCAN, identify which peak overlaps a TSS (promoter anchor), then collect the remaining peaks as candidate enhancers linked to that gene.

```python
import pybedtools

# Load CCAN assignments and a TSS reference
ccans = pd.read_csv('cicero_ccans.csv')
tss = pd.read_csv('tss_reference.bed', sep='\t', header=None,
                   names=['chr', 'start', 'end', 'gene'])

# Parse peak coordinates from CCAN table
ccans[['chr', 'start', 'end']] = ccans['peak'].str.split('[:-]', expand=True)
ccans['start'] = ccans['start'].astype(int)
ccans['end'] = ccans['end'].astype(int)

# Intersect with TSS to label promoter peaks
peaks_bed = pybedtools.BedTool.from_dataframe(ccans[['chr', 'start', 'end', 'peak', 'CCAN']])
tss_bed = pybedtools.BedTool.from_dataframe(tss)
promoter_peaks = peaks_bed.intersect(tss_bed, wa=True, wb=True).to_dataframe()

# For each CCAN with a promoter anchor, extract the non-promoter (enhancer) peaks
promoter_peak_ids = set(promoter_peaks['name'])  # 'name' = peak ID from col 4
enhancer_peaks = ccans[~ccans['peak'].isin(promoter_peak_ids)].copy()

# Add the linked gene(s) via CCAN membership
ccan_to_gene = promoter_peaks.groupby('score')['thickEnd'].apply(list).to_dict()
# (column names depend on bedtools output; adjust accordingly)
enhancer_peaks['linked_gene'] = enhancer_peaks['CCAN'].map(ccan_to_gene)
```

**Output:** DataFrame of enhancer peaks with columns `[chr, start, end, CCAN, linked_gene]`.

### Step 5: Motif Scan Enhancer Peaks for Candidate TFs

Scan the enhancer regions for occurrences of your chromVAR-nominated TF motifs. This tells you *which enhancers harbor which TF motifs*.

```python
# Option A: Use scprinter's built-in motif scanning
# Create an AnnData of just the enhancer regions
enhancer_regions_bed = enhancer_peaks[['chr', 'start', 'end']]
enhancer_regions_bed.to_csv('enhancer_regions.bed', sep='\t', header=False, index=False)

adata_enh = scp.pp.make_peak_matrix(printer, regions='enhancer_regions.bed')
motifs.chromvar_scan(adata_enh)
# adata_enh.varm['motif_scan'] now contains a (peaks x motifs) binary match matrix

# Option B: External FIMO scan (more control over p-value thresholds)
# fimo --thresh 1e-4 motifs.meme enhancer_regions.fa > fimo_results.tsv
```

**Output:** A mapping of (enhancer_region, TF_motif, motif_position_within_region).

### Step 6: Construct Targeted Region Sets for Footprinting

For each TF of interest, collect the enhancer regions containing its motif. Optionally re-center on the motif occurrence for cleaner footprint visualization.

```python
# Example: build region sets for top TFs
# Assume motif_scan is a binary DataFrame (regions x TFs) from Step 5

def get_tf_enhancer_regions(tf_name, enhancer_peaks_df, motif_scan_df):
    """Return enhancer regions containing a given TF motif."""
    matching_idx = motif_scan_df.index[motif_scan_df[tf_name] > 0]
    return enhancer_peaks_df.loc[matching_idx, ['chr', 'start', 'end']].copy()

# Build a dict of {tf_name: regions_dataframe}
tf_region_sets = {}
for tf in candidate_tfs:
    regions_df = get_tf_enhancer_regions(tf, enhancer_peaks, motif_scan)
    if len(regions_df) >= 10:  # minimum region count for meaningful footprinting
        regions_df.columns = ['Chromosome', 'Start', 'End']
        regions_df = regions_df.reset_index(drop=True)
        tf_region_sets[tf] = regions_df
```

**Output:** Dictionary mapping each TF to a scprinter-compatible regions DataFrame.

### Step 7: Define Pseudobulk Cell Groupings

Use the chromVAR-derived clusters (or any biologically meaningful grouping) to define barcode groups for scprinter.

```python
# Map cluster labels to barcode groups
barcode_groups = pd.DataFrame({
    'barcode': chromvar_results.obs.index,
    'group': chromvar_results.obs['leiden'].values
})

grouping, uniq_groups = scp.utils.df2cell_grouping(printer, barcode_groups)
```

### Step 8: seq2PRINT Model Training (Optional but Recommended)

If you want sequence-informed binding predictions (not just statistical footprints), train the base model and LoRA-adapt per pseudobulk.

```python
# Train base model on all cells (requires GPU, ~1-4 hours)
model_config = scp.tl.seq_model_config(
    printer, model_name='Base_Bulk', preset='seq2PRINT'
)
scp.tl.launch_seq2print(model_config)

# LoRA fine-tune per pseudobulk group
lora_config = scp.tl.seq_lora_model_config(
    printer,
    pretrain_model='Base_Bulk_model.pt',
    model_name='LoRA_per_group'
)
scp.tl.launch_seq2print(lora_config)
```

**Output:** Trained seq2PRINT models (base + per-group LoRA).

### Step 9: Compute Binding Scores and Multiscale Footprints

Run footprinting on each TF's enhancer region set.

```python
# Load models
printer.load_disp_model()
printer.load_bindingscore_model("TF", scp.datasets.pretrained_TFBS_model)
printer.load_bindingscore_model("Nuc", scp.datasets.pretrained_NucBS_model)

# For each TF's enhancer region set:
for tf_name, regions in tf_region_sets.items():

    # TFBS binding scores
    scp.tl.get_binding_score(
        printer, grouping, uniq_groups, regions,
        model_key='TF', n_jobs=32, contextRadius=100,
        save_key=f"enhancer_{tf_name}_TFBS", backed=True
    )

    # Multiscale footprints (scales 2-100 bp)
    scp.tl.get_footprint_score(
        printer, grouping, uniq_groups, regions,
        modes=np.arange(2, 101), n_jobs=32,
        save_key=f"enhancer_{tf_name}_footprint", backed=True
    )
```

**Output:** Per-TF, per-cell-state binding score matrices and multiscale footprint tensors stored in the printer h5ad.

### Step 10: Visualization and Interpretation

```python
# Heatmap: binding score across cell types at a specific enhancer region
scp.pl.plot_binding_score(
    printer,
    save_key=f'enhancer_{tf_name}_TFBS',
    group_names=uniq_groups[order],
    kind='heatmap',
    region=regions.iloc[0],
    row_label=cell_type_labels[order]
)

# Multiscale footprint at a specific enhancer
scp.pl.plot_footprints(
    printer,
    save_key=f'enhancer_{tf_name}_footprint',
    group_names=uniq_groups[:5],
    region=regions.iloc[0],
    cmap='Blues', vmin=0.5, vmax=2.0
)

# Stacked footprint heatmap across all cell types at one region
scp.pl.plot_footprints(
    printer,
    save_key=f'enhancer_{tf_name}_footprint',
    group_names=uniq_groups[order],
    row_label=cell_type_labels[order],
    region=regions.iloc[0],
    stack=True, scales=[20, 50, 100],
    cmap='Blues', vmin=0.5, vmax=2.0
)
```

### Recipe A Summary Table

| Step | Tool | Input | Output | Modality |
|------|------|-------|--------|----------|
| 1 | scprinter | fragments.tsv.gz | printer object + peaks | ATAC |
| 2 | scprinter/chromVAR | peaks + barcodes | TF deviation z-scores per cell | ATAC |
| 3 | Cicero | peak matrix + UMAP | CCANs + co-accessibility scores | ATAC |
| 4 | pybedtools | CCANs + TSS ref | enhancer peaks linked to genes | ATAC |
| 5 | scprinter motifs / FIMO | enhancer peaks + motif DB | (enhancer, TF) match matrix | ATAC |
| 6 | pandas | match matrix + candidate TFs | per-TF region sets | ATAC |
| 7 | scprinter | cluster labels + barcodes | pseudobulk groupings | ATAC |
| 8 | seq2PRINT | all fragments + GPU | base + LoRA models | ATAC |
| 9 | scprinter | region sets + groupings | binding scores + footprints | ATAC |
| 10 | scprinter.pl | scores + footprints | heatmaps, footprint plots | ATAC |

---

## RECIPE B: Multiome-Informed Enhancer Footprinting

Recipe B extends Recipe A with three upgrades that require paired RNA+ATAC in the same cells:

1. **Replace Cicero CCANs with DORCs** — peak-gene links based on expression correlation (stronger evidence than co-accessibility)
2. **Replace chromVAR clusters with RNA-informed clusters** — transcriptionally defined cell states for pseudobulk grouping
3. **Add SCENIC+ eRegulons** — pre-computed TF -> enhancer -> gene triplets as an alternative or complementary entry point
4. **Cross-modal validation** — correlate footprint depth with TF expression and target gene expression

### Conceptual Flow

```
Paired RNA + ATAC (muData / matched AnnData objects)
  |
  +-> MultiVI joint embedding -> RNA-informed clusters (Step B1)
  |
  +-> DORC analysis: peak-gene expression correlation (Step B2)
  |     Output: DORC genes + their linked enhancer peaks
  |
  +-> SCENIC+ eGRN inference (Step B3)
  |     Output: eRegulons = (TF, enhancer set, target gene set)
  |
  +-> Motif scan DORC/eRegulon enhancer regions (Step B4)
  |
  +-> scprinter footprinting with RNA-informed groupings (Step B5)
  |
  +-> Cross-modal validation: footprint vs expression (Step B6)
  |
  +-> hdWGCNA integration (Step B7)
```

### Step B1: RNA-Informed Cell State Definitions

Use the RNA modality (or a joint embedding like MultiVI) to define cell states at higher resolution than ATAC clustering allows.

```python
# Assuming mudata object with mudata['rna'] and mudata['atac']
# Or two matched AnnData objects with aligned barcodes

import scanpy as sc

# Standard RNA processing
adata_rna = mudata['rna'].copy()
sc.pp.normalize_total(adata_rna, target_sum=1e4)
sc.pp.log1p(adata_rna)
sc.pp.highly_variable_genes(adata_rna, n_top_genes=3000)
sc.tl.pca(adata_rna)
sc.pp.neighbors(adata_rna)
sc.tl.leiden(adata_rna, resolution=0.8)
sc.tl.umap(adata_rna)

# These RNA-defined clusters become the pseudobulk groupings for scprinter
rna_barcode_groups = pd.DataFrame({
    'barcode': adata_rna.obs.index,
    'group': adata_rna.obs['leiden'].values
})

# Convert to scprinter grouping format
grouping_rna, uniq_groups_rna = scp.utils.df2cell_grouping(printer, rna_barcode_groups)
```

**Key advantage over Recipe A:** RNA clusters distinguish cell states that share accessibility landscapes but differ transcriptionally (e.g., naive vs memory T cells, subtypes within a tumor).

### Step B2: DORC Analysis — Expression-Correlated Enhancer Peaks

DORCs replace Cicero's co-accessibility links with a stronger evidence standard: peaks whose accessibility correlates with gene expression in the same cell.

```python
adata_atac = mudata['atac'].copy()
adata_rna_matched = mudata['rna'].copy()

# Align barcodes (should already be aligned in muData)
shared = list(set(adata_rna_matched.obs.index) & set(adata_atac.obs.index))
adata_rna_matched = adata_rna_matched[shared].copy()
adata_atac = adata_atac[shared].copy()

# Normalize RNA
sc.pp.normalize_total(adata_rna_matched, target_sum=1e4)
sc.pp.log1p(adata_rna_matched)

# Compute peak-gene correlations (50kb window around TSS)
dorc_all = scp.dorc.fast_gene_peak_corr(
    adata_atac, adata_rna_matched,
    genome=scp.genome.hg38,
    tss_df=scp.datasets.FigR_hg38TSSRanges,
    window_pad_size=50000,
    n_jobs=32, n_bg=100, pos_only=True
)

# Filter significant associations
dorc_sig = dorc_all[dorc_all['pvalZ'] <= 0.05]

# Identify DORC genes (>= 7 significantly correlated peaks)
dorc_gene_list = scp.dorc.dorc_j_plot(
    dorc_sig, cutoff=7, label_top=25, return_gene_list=True
)

# Extract the enhancer peaks linked to DORC genes
# dorc_sig contains columns: ['peak', 'gene', 'corr', 'pvalZ']
dorc_enhancer_peaks = dorc_sig[dorc_sig['gene'].isin(dorc_gene_list)].copy()
# Parse peak coordinates
dorc_enhancer_peaks[['chr', 'start', 'end']] = (
    dorc_enhancer_peaks['peak'].str.split('[:-]', expand=True)
)
dorc_enhancer_peaks['start'] = dorc_enhancer_peaks['start'].astype(int)
dorc_enhancer_peaks['end'] = dorc_enhancer_peaks['end'].astype(int)
```

**Output:** `dorc_enhancer_peaks` DataFrame with columns `[chr, start, end, gene, corr, pvalZ]` — enhancer peaks with known target genes and quantified regulatory strength.

**Key advantage over Cicero CCANs:** Each enhancer-gene link is supported by expression correlation in the same cell, not just co-accessibility. Genes with many such peaks (DORCs) are enriched for lineage-defining regulators.

### Step B3: SCENIC+ eRegulon Inference (Alternative/Complementary Entry Point)

SCENIC+ provides a complete TF -> enhancer -> gene mapping using motif enrichment + GRNBoost2 co-expression. This can serve as the primary source of region sets, or as independent validation of DORC-derived regions.

```python
# SCENIC+ runs as its own pipeline (scenicplus package)
# Key outputs needed for footprinting:
#   1. eRegulon table: TF, list of target enhancers, list of target genes, direction (+/-)
#   2. Region-to-gene links with correlation scores
#
# After running SCENIC+, extract enhancer region sets per eRegulon:

# Pseudocode (actual API depends on scenicplus version):
for eregulon in scenicplus_results.eRegulons:
    tf_name = eregulon.TF
    enhancer_regions = eregulon.target_regions  # list of (chr, start, end)
    target_genes = eregulon.target_genes
    direction = eregulon.direction  # activating (+) or repressing (-)

    # Format for scprinter
    regions_df = pd.DataFrame(enhancer_regions, columns=['Chromosome', 'Start', 'End'])
    eregulon_region_sets[tf_name] = regions_df
```

**Key advantage:** SCENIC+ provides the TF identity, the enhancer set, AND the target gene set as a pre-computed triplet. No need to separately run chromVAR + Cicero/DORC + motif scanning. However, footprinting adds the physical binding evidence layer that SCENIC+ infers only statistically.

### Step B4: Motif Scan DORC Enhancer Peaks

Same as Recipe A Step 5, but applied to DORC-derived or eRegulon-derived enhancer regions.

```python
# For DORC-derived enhancers, you still need to determine which TFs might bind there.
# Two approaches:
#
# Approach 1: chromVAR nomination (same as Recipe A)
#   Use chromVAR top TFs + motif scan of DORC enhancer peaks
#
# Approach 2: SCENIC+ eRegulon nomination
#   The TF is already specified by the eRegulon; skip motif scanning
#
# Approach 3: Both — use DORC regions but restrict to TFs nominated by SCENIC+ eRegulons
#   This gives you the intersection: enhancers correlated with expression (DORC)
#   that are also part of an inferred regulatory network (SCENIC+)

# Construct the combined region sets:
for tf_name in eregulon_region_sets:
    # Option: intersect eRegulon enhancers with DORC-linked enhancers
    ereg_regions = eregulon_region_sets[tf_name]
    dorc_regions = dorc_enhancer_peaks[['chr', 'start', 'end']].drop_duplicates()

    # BedTools intersection for high-confidence enhancer set
    ereg_bed = pybedtools.BedTool.from_dataframe(ereg_regions)
    dorc_bed = pybedtools.BedTool.from_dataframe(dorc_regions)
    overlap = ereg_bed.intersect(dorc_bed, wa=True).to_dataframe()

    if len(overlap) >= 10:
        overlap.columns = ['Chromosome', 'Start', 'End']
        high_confidence_regions[tf_name] = overlap
```

### Step B5: scprinter Footprinting with RNA-Informed Groupings

Identical to Recipe A Steps 8-9, but using `grouping_rna` / `uniq_groups_rna` from Step B1.

```python
printer.load_disp_model()
printer.load_bindingscore_model("TF", scp.datasets.pretrained_TFBS_model)
printer.load_bindingscore_model("Nuc", scp.datasets.pretrained_NucBS_model)

for tf_name, regions in high_confidence_regions.items():

    scp.tl.get_binding_score(
        printer, grouping_rna, uniq_groups_rna, regions,
        model_key='TF', n_jobs=32, contextRadius=100,
        save_key=f"multiome_enh_{tf_name}_TFBS", backed=True
    )

    scp.tl.get_footprint_score(
        printer, grouping_rna, uniq_groups_rna, regions,
        modes=np.arange(2, 101), n_jobs=32,
        save_key=f"multiome_enh_{tf_name}_FP", backed=True
    )
```

### Step B6: Cross-Modal Validation (Multiome-Exclusive)

This is impossible without paired data. For each (TF, enhancer, target gene) triplet, verify that footprint depth, TF expression, and target gene expression are concordant across cell states.

```python
# For each TF:
#   1. Extract per-group mean footprint score at its enhancer regions
#   2. Extract per-group mean TF expression (RNA)
#   3. Extract per-group mean target gene expression (RNA)
#   4. Correlate all three

results = []
for tf_name in high_confidence_regions:

    # Footprint scores (from scprinter, aggregated across enhancer regions)
    # ... extract from printer h5ad ...

    # TF expression per group
    tf_expr_per_group = []
    for group_label in uniq_groups_rna:
        cells = adata_rna.obs[adata_rna.obs['leiden'] == group_label].index
        tf_expr_per_group.append(adata_rna[cells, tf_name].X.mean())

    # Target gene expression per group (from DORC or eRegulon)
    target_genes = dorc_enhancer_peaks[
        dorc_enhancer_peaks['gene'].isin(dorc_gene_list)
    ]['gene'].unique()  # or from eRegulon

    for gene in target_genes:
        gene_expr_per_group = []
        for group_label in uniq_groups_rna:
            cells = adata_rna.obs[adata_rna.obs['leiden'] == group_label].index
            gene_expr_per_group.append(adata_rna[cells, gene].X.mean())

        # Correlate: footprint ~ TF expression, footprint ~ gene expression
        results.append({
            'TF': tf_name,
            'target_gene': gene,
            'corr_footprint_TFexpr': np.corrcoef(footprint_scores, tf_expr_per_group)[0,1],
            'corr_footprint_geneexpr': np.corrcoef(footprint_scores, gene_expr_per_group)[0,1],
            'corr_TFexpr_geneexpr': np.corrcoef(tf_expr_per_group, gene_expr_per_group)[0,1],
        })

validation_df = pd.DataFrame(results)
```

**Interpretation guide:**

| footprint ~ TF expr | footprint ~ gene expr | Meaning |
|---|---|---|
| High positive | High positive | Classic activator: TF is expressed, bound, and target is on |
| High positive | High negative | Candidate repressor: TF is expressed and bound, target is off |
| Low / none | High positive | Indirect regulation: target is on but not through this TF's direct binding |
| Discordant | — | Post-transcriptional regulation, chromatin memory, or false positive |

### Step B7: hdWGCNA Integration (Optional Extension)

Connect footprinting results to co-expression modules to ask: do genes regulated by the same TF (per footprinting) fall into the same co-expression module?

```r
# In R:
# 1. Run hdWGCNA on the RNA modality (see hdWGCNA core workflow)
# 2. Get module assignments for DORC genes
# 3. Test: are DORC genes linked to the same TF (by footprinting)
#    enriched in the same module?

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# For each TF's target gene set (from DORC or eRegulon):
#   Compute overlap with each hdWGCNA module using Fisher's exact test
#   Significant overlap = the TF's regulatory program maps onto a
#   co-expression module, providing independent RNA-level validation
```

---

## RECIPE C: CellChat-Guided Regulatory Footprinting

### Rationale

CellChat infers intercellular signaling from ligand-receptor co-expression in scRNA-seq data. But a predicted signal is only biologically meaningful if the receiver cell actually transduces it — and transduction ultimately means a downstream TF binds chromatin and activates target genes. This recipe traces the full signal transduction chain from sender to receiver chromatin by combining CellChat with the footprinting framework established in Recipes A and B.

**Prerequisite:** Recipe B completed (SCENIC+ eRegulons, DORCs, chromVAR deviations, and scprinter printer object all available).

### Conceptual Flow

```
CellChat (RNA) -> Significant signaling pathways with sender/receiver cell types
    |
    v
Curated pathway-to-TF dictionary -> Candidate downstream TFs per pathway
    |
    v
chromVAR validation -> Confirm downstream TF motif enrichment in receiver cells
    |
    v
eRegulon / DORC lookup -> Extract enhancer targets of validated downstream TFs
    |
    v
scprinter footprinting -> Physical binding evidence in receiver-cell pseudobulks
    |
    v
Cross-modal chain correlation -> Ligand (sender) ~ Footprint (receiver) ~ Target gene (receiver)
```

### Step C1: CellChat — Infer Signaling Networks

CellChat runs in R on the RNA modality. The key outputs are: significant ligand-receptor interactions, pathway classifications, and sender/receiver centrality scores per cell type.

```r
library(CellChat)
library(patchwork)

# Create CellChat object from normalized RNA counts + cell labels
cellchat <- createCellChat(object = seurat_obj, group.by = "cell_type")

# Set ligand-receptor database (human or mouse)
CellChatDB <- CellChatDB.human  # or CellChatDB.mouse
cellchat@DB <- CellChatDB

# Standard CellChat pipeline
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Network analysis: centrality scores identify senders, receivers, mediators
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Extract significant interactions as a table
sig_interactions <- subsetCommunication(cellchat)
# Columns: source, target, ligand, receptor, pathway_name, prob, pval

# Export for Python
write.csv(sig_interactions, 'cellchat_significant_interactions.csv', row.names = FALSE)

# Also export the pathway-level summary
pathway_summary <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(pathway_summary, 'cellchat_pathway_summary.csv', row.names = FALSE)
```

**Output:** `cellchat_significant_interactions.csv` — each row is a (sender, receiver, ligand, receptor, pathway, probability, p-value) tuple.

### Step C2: Curated Pathway-to-Downstream-TF Dictionary

Map each CellChat signaling pathway to its canonical downstream transcription factors. This dictionary is biology-agnostic in that these pathways are conserved across tissues; the tissue-specific part is *which pathways are active*, which CellChat already determined.

```python
import pandas as pd

# Curated mapping: CellChat pathway name -> list of canonical downstream TFs
# These are the TFs activated in RECEIVER cells upon signal transduction.
# Sources: KEGG, Reactome, primary literature.
# NOTE: Expand this dictionary as needed for your specific dataset.

PATHWAY_TO_TFS = {
    # TGF-beta superfamily
    'TGFb':     ['SMAD2', 'SMAD3', 'SMAD4'],
    'BMP':      ['SMAD1', 'SMAD5', 'SMAD9', 'SMAD4'],
    'ACTIVIN':  ['SMAD2', 'SMAD3'],
    'GDF':      ['SMAD1', 'SMAD5'],
    'NODAL':    ['SMAD2', 'SMAD3', 'FOXH1'],

    # WNT family
    'WNT':      ['TCF7', 'LEF1', 'TCF7L1', 'TCF7L2'],
    'ncWNT':    ['JUN', 'FOS', 'ATF2', 'NFAT5'],  # non-canonical: PCP/Ca2+ arms

    # Notch
    'NOTCH':    ['RBPJ', 'HES1', 'HEY1', 'HEY2'],

    # JAK-STAT cytokines
    'IFN-I':    ['STAT1', 'STAT2', 'IRF9'],
    'IFN-II':   ['STAT1'],
    'IL6':      ['STAT3'],
    'IL2':      ['STAT5A', 'STAT5B'],
    'IL4':      ['STAT6'],
    'IL10':     ['STAT3'],
    'IL1':      ['NFKB1', 'RELA', 'JUN'],
    'OSM':      ['STAT3'],
    'LIF':      ['STAT3'],
    'CSF':      ['STAT5A', 'SPI1'],
    'CSF3':     ['STAT3', 'CEBPE'],

    # TNF / NF-kB axis
    'TNF':      ['NFKB1', 'RELA', 'NFKB2', 'RELB'],
    'LIGHT':    ['NFKB1', 'RELA'],
    'RANKL':    ['NFKB1', 'FOS', 'NFATC1'],
    'TWEAK':    ['NFKB1'],
    'CD40':     ['NFKB1', 'RELA', 'NFKB2'],
    'BAFF':     ['NFKB2', 'RELB'],

    # Receptor tyrosine kinase / MAPK
    'EGF':      ['FOS', 'JUN', 'ELK1', 'MYC', 'ETS1'],
    'FGF':      ['ETS1', 'ETS2', 'FOS', 'JUN'],
    'PDGF':     ['FOS', 'JUN', 'ETS1', 'SRF'],
    'HGF':      ['FOS', 'JUN', 'ETS1'],
    'VEGF':     ['ETS1', 'FOS', 'JUN'],
    'IGF':      ['FOXO1', 'FOXO3', 'MYC'],
    'KIT':      ['MITF', 'FOS', 'JUN'],
    'NGF':      ['CREB1', 'FOS', 'JUN'],

    # Hedgehog
    'SHH':      ['GLI1', 'GLI2', 'GLI3'],
    'IHH':      ['GLI1', 'GLI2'],

    # Hippo / YAP
    'HIPPO':    ['TEAD1', 'TEAD2', 'TEAD3', 'TEAD4'],

    # Retinoic acid
    'RA':       ['RARA', 'RARB', 'RARG', 'RXRA'],

    # Chemokines
    'CXCL':     ['NFKB1', 'STAT3', 'FOS', 'JUN'],
    'CCL':      ['NFKB1', 'STAT3', 'FOS'],

    # Cell adhesion / contact-based
    'CDH':      ['CTNNB1', 'TCF7', 'LEF1'],  # cadherin -> beta-catenin
    'JAM':      ['CEBPB'],
    'NECTIN':   ['FOS', 'JUN'],
    'EPHA':     ['FOS', 'JUN'],
    'EPHB':     ['FOS', 'JUN'],

    # ECM-receptor
    'COLLAGEN': ['SRF', 'FOS', 'JUN'],
    'LAMININ':  ['FOS', 'JUN'],
    'FN1':      ['SRF', 'FOS', 'JUN'],
    'SPP1':     ['NFKB1', 'FOS', 'JUN'],
    'THBS':     ['NFKB1'],
    'TENASCIN': ['FOS', 'JUN'],

    # Others
    'SEMA3':    ['NFAT5', 'FOS'],
    'GAS':      ['STAT3'],
    'PROS':     ['AXL'],  # TAM receptor -> limited TF data
    'ANNEXIN':  ['NFKB1'],
    'MK':       ['STAT3', 'FOS'],
    'PTN':      ['FOS', 'JUN'],
    'GALECTIN': ['NFKB1', 'NFAT5'],
    'ANGPT':    ['STAT3', 'FOXO1'],
    'APELIN':   ['KLF2', 'NRF2'],
    'EDN':      ['FOS', 'JUN', 'NFAT5'],
    'MIF':      ['NFKB1', 'AP1'],
}

# Load CellChat results
sig_interactions = pd.read_csv('cellchat_significant_interactions.csv')

# Map each interaction to downstream TFs
sig_interactions['downstream_TFs'] = (
    sig_interactions['pathway_name'].map(PATHWAY_TO_TFS)
)

# Drop pathways without a mapping (unmapped or novel)
mapped = sig_interactions.dropna(subset=['downstream_TFs']).copy()

# Build the working table: (receiver_celltype, pathway, downstream_TF)
receiver_tf_table = []
for _, row in mapped.iterrows():
    for tf in row['downstream_TFs']:
        receiver_tf_table.append({
            'sender': row['source'],
            'receiver': row['target'],
            'pathway': row['pathway_name'],
            'ligand': row['ligand'],
            'receptor': row['receptor'],
            'downstream_TF': tf,
            'interaction_prob': row['prob'],
        })

receiver_tf_df = pd.DataFrame(receiver_tf_table)
```

**Output:** `receiver_tf_df` — each row is a (sender, receiver, pathway, ligand, receptor, downstream_TF, probability) record. This is the hypothesis table: for each predicted signaling event, we now have a specific TF to look for in the receiver cell's chromatin.

### Step C3: chromVAR Validation in Receiver Cells

Confirm that the predicted downstream TFs actually show elevated motif accessibility in receiver cells. This filters out pathway-TF mappings that are canonically correct but not active in this dataset.

```python
# chromvar_results from Recipe A/B Step 2
# receiver_tf_df from Step C2

validated_pairs = []

for _, row in receiver_tf_df.drop_duplicates(
    subset=['receiver', 'downstream_TF']
).iterrows():

    tf = row['downstream_TF']
    receiver_type = row['receiver']

    # Check if this TF motif exists in the chromVAR results
    if tf not in chromvar_results.var_names:
        continue

    # Get chromVAR deviation scores for this TF in receiver vs all other cells
    receiver_mask = chromvar_results.obs['cell_type'] == receiver_type
    receiver_scores = chromvar_results[receiver_mask, tf].X.flatten()
    other_scores = chromvar_results[~receiver_mask, tf].X.flatten()

    # Statistical test: is the motif enriched in receiver cells?
    from scipy.stats import mannwhitneyu
    stat, pval = mannwhitneyu(receiver_scores, other_scores, alternative='greater')
    mean_diff = receiver_scores.mean() - other_scores.mean()

    if pval < 0.05 and mean_diff > 0:
        validated_pairs.append({
            **row.to_dict(),
            'chromvar_pval': pval,
            'chromvar_mean_diff': mean_diff,
        })

validated_df = pd.DataFrame(validated_pairs)
```

**Output:** `validated_df` — subset of `receiver_tf_df` where the downstream TF's motif is significantly enriched in chromatin accessibility in the predicted receiver cell type. These are the signal transduction events with chromatin-level support.

### Step C4: Retrieve Enhancer Targets from eRegulons / DORCs

For each validated (receiver, downstream_TF) pair, look up the TF's target enhancers from SCENIC+ eRegulons computed in Recipe B. Cross-reference with DORCs for additional confidence.

```python
# eregulon_region_sets from Recipe B Step B3
# dorc_enhancer_peaks from Recipe B Step B2
# high_confidence_regions from Recipe B Step B4 (eRegulon ∩ DORC overlap)

import pybedtools

signaling_footprint_targets = {}

for _, row in validated_df.drop_duplicates(
    subset=['receiver', 'downstream_TF']
).iterrows():

    tf = row['downstream_TF']
    receiver = row['receiver']
    pathway = row['pathway']
    key = f"{pathway}__{receiver}__{tf}"

    # Primary: use eRegulon enhancer regions for this TF
    if tf in eregulon_region_sets:
        regions = eregulon_region_sets[tf].copy()
    # Fallback: use DORC peaks that contain this TF's motif (from Recipe B Step B4)
    elif tf in high_confidence_regions:
        regions = high_confidence_regions[tf].copy()
    else:
        # Last resort: scan all DORC enhancer peaks for this TF's motif
        # (requires motif_scan matrix from Recipe A/B Step 5)
        if tf in motif_scan.columns:
            matching = motif_scan.index[motif_scan[tf] > 0]
            if len(matching) >= 10:
                regions = dorc_enhancer_peaks.loc[matching, ['chr', 'start', 'end']].copy()
                regions.columns = ['Chromosome', 'Start', 'End']
                regions = regions.reset_index(drop=True)
            else:
                continue
        else:
            continue

    if len(regions) >= 10:
        signaling_footprint_targets[key] = {
            'regions': regions,
            'tf': tf,
            'receiver': receiver,
            'pathway': pathway,
            'sender': row['sender'],
            'ligand': row['ligand'],
            'receptor': row['receptor'],
        }
```

**Output:** `signaling_footprint_targets` — dictionary keyed by `pathway__receiver__TF`, each containing a scprinter-compatible regions DataFrame and metadata about the signaling chain.

### Step C5: Footprint Downstream TFs in Receiver Cell Pseudobulks

Construct pseudobulks from receiver cell type barcodes and footprint.

```python
# Build receiver-cell-type pseudobulks for scprinter
# Use RNA-informed cell type labels (from Recipe B Step B1)

printer.load_disp_model()
printer.load_bindingscore_model("TF", scp.datasets.pretrained_TFBS_model)
printer.load_bindingscore_model("Nuc", scp.datasets.pretrained_NucBS_model)

for key, target_info in signaling_footprint_targets.items():

    regions = target_info['regions']
    tf = target_info['tf']

    # Footprint using the full cell-type groupings (all types)
    # so we can compare receiver vs non-receiver
    scp.tl.get_binding_score(
        printer, grouping_rna, uniq_groups_rna, regions,
        model_key='TF', n_jobs=32, contextRadius=100,
        save_key=f"signal_{key}_TFBS", backed=True
    )

    scp.tl.get_footprint_score(
        printer, grouping_rna, uniq_groups_rna, regions,
        modes=np.arange(2, 101), n_jobs=32,
        save_key=f"signal_{key}_FP", backed=True
    )
```

**Future expansion:** Instead of cell-type-level pseudobulks, split receiver cells by receptor expression quantile (high vs. low) using KNN-based pseudobulk sampling. This would test whether receptor expression level within a cell type correlates with downstream TF binding depth, providing within-cell-type dose-response evidence. This requires sufficient cell numbers per quantile bin (recommended: >= 500 cells per bin for adequate fragment depth).

### Step C6: Signal Transduction Chain Correlation

The final validation: across cell states, correlate the full chain — ligand expression in sender, footprint depth of downstream TF in receiver, and target gene expression in receiver.

```python
chain_results = []

for key, target_info in signaling_footprint_targets.items():

    tf = target_info['tf']
    sender_type = target_info['sender']
    receiver_type = target_info['receiver']
    pathway = target_info['pathway']
    ligand = target_info['ligand']
    receptor = target_info['receptor']

    # --- Ligand expression in sender (per sample/condition if available) ---
    # If dataset has multiple samples or conditions, compute per-condition means.
    # If single-sample, compute per-cluster means as a proxy.
    ligand_expr = []
    for group_label in uniq_groups_rna:
        cells = adata_rna.obs[adata_rna.obs['leiden'] == group_label].index
        if ligand in adata_rna.var_names:
            ligand_expr.append(float(adata_rna[cells, ligand].X.mean()))
        else:
            ligand_expr.append(np.nan)

    # --- Receptor expression in receiver cells per group ---
    receptor_expr = []
    for group_label in uniq_groups_rna:
        cells = adata_rna.obs[adata_rna.obs['leiden'] == group_label].index
        if receptor in adata_rna.var_names:
            receptor_expr.append(float(adata_rna[cells, receptor].X.mean()))
        else:
            receptor_expr.append(np.nan)

    # --- Downstream TF footprint depth in receiver per group ---
    # Extract from scprinter results (aggregate binding score across enhancer regions)
    # footprint_scores = ... (extract from printer h5ad, save_key=f"signal_{key}_TFBS")

    # --- Target gene expression per group ---
    # Get target genes from eRegulon or DORC for this TF
    if tf in eregulon_region_sets:
        # Look up target genes from SCENIC+ eRegulon
        target_genes = [g for g in scenicplus_target_genes.get(tf, [])
                        if g in adata_rna.var_names]
    else:
        # Fall back to DORC genes linked to peaks containing this TF's motif
        target_genes = dorc_enhancer_peaks[
            dorc_enhancer_peaks['peak'].isin(
                motif_scan.index[motif_scan.get(tf, pd.Series(dtype=float)) > 0]
            )
        ]['gene'].unique().tolist()

    for gene in target_genes[:20]:  # cap at 20 to limit compute
        gene_expr = []
        for group_label in uniq_groups_rna:
            cells = adata_rna.obs[adata_rna.obs['leiden'] == group_label].index
            gene_expr.append(float(adata_rna[cells, gene].X.mean()))

        chain_results.append({
            'pathway': pathway,
            'sender': sender_type,
            'receiver': receiver_type,
            'ligand': ligand,
            'receptor': receptor,
            'downstream_TF': tf,
            'target_gene': gene,
            # Pairwise correlations across cell groups:
            # 'corr_ligand_footprint': np.corrcoef(ligand_expr, footprint_scores)[0,1],
            # 'corr_footprint_target': np.corrcoef(footprint_scores, gene_expr)[0,1],
            # 'corr_receptor_footprint': np.corrcoef(receptor_expr, footprint_scores)[0,1],
            # 'corr_ligand_target': np.corrcoef(ligand_expr, gene_expr)[0,1],
        })

chain_df = pd.DataFrame(chain_results)
```

**Output:** `chain_df` — for each signaling chain (pathway → sender → ligand → receptor → receiver → TF → target gene), quantified correlations at every link.

### Interpretation Guide: Signal Transduction Evidence Tiers

| Evidence Tier | Criteria | Interpretation |
|---|---|---|
| **Tier 1: Predicted** | CellChat significant (prob > threshold, p < 0.05) | Ligand-receptor co-expression suggests signaling |
| **Tier 2: Chromatin-validated** | Tier 1 + downstream TF motif enriched in receiver chromVAR | Signal is associated with open chromatin at TF motifs |
| **Tier 3: Footprint-confirmed** | Tier 2 + scprinter binding score elevated at TF's enhancer targets in receiver | Physical evidence of TF binding at regulatory elements |
| **Tier 4: Functionally complete** | Tier 3 + target gene expression correlated with footprint depth | Full chain from signal to transcriptional output confirmed |

Interactions reaching Tier 4 represent the strongest possible computational evidence for functional intercellular signaling without perturbation experiments.

### Recipe C Summary Table

| Step | Tool | Input | Output |
|------|------|-------|--------|
| C1 | CellChat (R) | RNA counts + cell labels | Significant L-R interactions with pathway labels |
| C2 | Curated dictionary | CellChat pathways | (receiver, pathway, downstream_TF) hypothesis table |
| C3 | chromVAR | Hypothesis table + deviation scores | Validated (receiver, TF) pairs with chromatin support |
| C4 | eRegulon / DORC lookup | Validated TFs | Enhancer region sets per downstream TF |
| C5 | scprinter | Region sets + receiver pseudobulks | Binding scores + multiscale footprints |
| C6 | Cross-modal correlation | Footprints + RNA expression | Full signal transduction chain evidence table |

---

## Comparison: What Each Layer Adds

| Evidence Layer | Recipe A (ATAC-only) | Recipe B (Multiome) | Recipe C (CellChat + Multiome) |
|---|---|---|---|
| **Enhancer identification** | Cicero CCANs (co-accessibility) | DORCs (expression correlation) | eRegulon/DORC enhancers of downstream TFs |
| **Enhancer-gene linkage** | Statistical co-accessibility | Direct expression correlation in same cell | Inherited from Recipe B |
| **TF nomination** | chromVAR deviation scores | chromVAR + SCENIC+ eRegulons | Pathway-to-TF dictionary + chromVAR validation |
| **TF -> enhancer mapping** | Motif scan (sequence match) | Motif enrichment + GRNBoost2 correlation | eRegulon lookup (pre-computed) |
| **Cell state definition** | ATAC clusters (chromVAR) | RNA-informed clusters (transcriptional) | RNA-informed + sender/receiver identity |
| **Physical binding evidence** | scprinter footprinting | scprinter footprinting (identical) | scprinter footprinting at signal-response loci |
| **Target gene validation** | Not possible | Footprint ~ expression correlation | Full chain: ligand ~ footprint ~ target gene |
| **Regulatory direction** | Not determinable | Activator vs repressor from correlation sign | Activator/repressor + signal origin (sender) |
| **Module-level context** | Not available | hdWGCNA module overlap | hdWGCNA + CellChat communication patterns |
| **Intercellular context** | Not available | Not available | Sender-receiver-pathway attribution |

---

## Critical API Notes for Agent Execution

1. **scprinter regions format:** Always `pd.DataFrame` with columns `['Chromosome', 'Start', 'End']`. Any BED-like source works.
2. **scprinter is disk-backed:** Always call `printer.close()` when done. Use `printer.load_printer()` to reopen.
3. **chromVAR GPU:** Requires CUDA. Initialize RMM before compute: `rmm.reinitialize(managed_memory=True, pool_allocator=True)`.
4. **seq2PRINT training:** GPU-intensive. Base model ~1-4 hours. LoRA ~80x faster than training separate models.
5. **DORC analysis requires paired data:** `fast_gene_peak_corr` needs matched RNA + ATAC AnnData with identical barcode indices.
6. **SCENIC+ is a separate pipeline:** Run independently, then import eRegulon tables for region set construction.
7. **Footprint modes:** `np.arange(2, 101)` covers TF-scale (2-20bp) through nucleosome-scale (50-100bp+). The multiscale decomposition is what distinguishes direct TF binding from nucleosome positioning.
8. **Minimum region count:** Aim for >= 10 regions per TF for meaningful aggregate footprints. Fewer regions = noisy signal.
