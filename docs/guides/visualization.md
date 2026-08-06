# Visualization

FORGE renders figures as pipeline outputs rather than as a separate manual step.
Most visualization blocks are opt-in, because they depend on upstream results that
are themselves expensive.

---

## What renders by default

With a standard run you get QC plots, UMAPs, annotation summaries, MOFA/MultiVI
latent-space plots, and Cicero target plots — no extra configuration.

Cicero target plots render at three window/threshold tiers so you can judge how
much a link depends on the threshold:

```groovy
cicero {
    target_genes    = ['APOE', 'TREM2']
    plot_windows    = "100000,250000,500000"
    plot_thresholds = "0.10,0.05,0.025"
    plot_labels     = "strict,permissive,exploratory"
}
```

---

## BigWig export

Coverage tracks, and the prerequisite for genome-browser panels:

```groovy
enhancer_viz {
    run           = true
    bin_size      = 10
    normalization = 'RPKM'      // 'RPKM' | 'CPM' | 'BPM' | 'None'
    min_cells     = 100         // skip groups below this
    blacklist     = null
    cell_type_col = null        // null → resolved chain; 'cell_type_broad' for fewer tracks
    condition_col = null        // set → also export per-(CT × condition) tracks
}
```

`condition_col` is required for any differential browser mode: without it there is
only one track per cell type, so there is nothing to contrast. Set it to your
`condition_group` key when running differential analyses.

## Genome browser

```groovy
browser_viz {
    enabled    = true
    mode       = 'absolute'    // 'absolute' | 'differential'
    window     = 15000         // bp padding around the gene body
    n_bins     = 700
    max_value  = null          // null → auto 99.5th percentile
    cell_types = null          // null → all cell types
}
```

Matplotlib-native accessibility tracks per gene. `differential` mode requires
bigWigs exported with `condition_col` set.

Leaving `max_value = null` auto-scales to the 99.5th percentile, which is usually
what you want — a single extreme bin otherwise flattens every track.

## Composite panels

```groovy
enhancer_viz {
    target_genes = ['APOE', 'TREM2']
    window_kb    = 100
    dpi          = 200
    composite_filter {
        top_n_per_gene_per_ct = 3
        fdr_cutoff            = 0.05
        lfc_cutoff            = 0.5
        min_motif_sites       = 1
    }
}
```

`composite_filter` selects (gene, TF, cell type) triples worth plotting. Without
it, the full cartesian product runs to tens of thousands of panels — the filter is
what makes this stage useful rather than merely large.

## Footprint strips

Requires footprinting. Two gates:

```groovy
enhancer_footprinting { msfp_enabled = true }    // outer — computes footprints
msfp_strip {
    enabled       = true                          // inner — renders strips
    mode          = 'all_three'                   // 'absolute' | 'differential' | 'all_three'
    context_bp    = 500000
    top_n_regions = 100
    top_k_genes   = 5
}
promoter_overlay { enabled = true }
```

Enhancer-strip genes are discovered per (cell type, TF); promoter-strip genes come
from `scprinter.target_genes`. See
[Regulatory analysis](regulatory.md#enhancer-strips).

---

## Re-rendering without recomputing

Figures usually need several iterations. Do not re-run the pipeline for them —
`VIZ_ONLY` rebuilds plots from persisted artifacts:

```bash
nextflow run main.nf -entry VIZ_ONLY -c my_study.config
```

```groovy
params.viz_only {
    peak_matrix_h5ad   = 'results/atac/final/peak_matrix.h5ad'
    cicero_connections = 'results/cicero/connections.csv'
    cicero_ccan        = 'results/cicero/CCAN_assignments.tsv.gz'
    cicero_cds         = 'results/cicero/input_cds_ordered.rds'
    target_genes       = 'APOE,TREM2,CST7'
}
```

The three Cicero keys are required together.

---

## Editable figures for publication

FORGE writes PNGs. For figures you intend to edit in Illustrator or Inkscape,
re-export as SVG with text kept as text rather than outlines:

```python
import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'   # keeps text selectable
```

Two practical notes from preparing the published figures:

- Use **SVG**, not PDF. A rasterized-plus-alpha layer in a PDF becomes an
  SMask, which many viewers render as a red X.
- Rasterize only the dense data layer (`rasterized=True` on the heavy artist) and
  leave axes, labels, and legends as vector. You get an editable figure without a
  200 MB file.

---

## Outputs

| Path | Contents |
|---|---|
| `results/enhancer_viz/tracks/` | BigWigs and track configs |
| `results/enhancer_viz/composites/` | Composite panels |
| `results/atac/final/` | QC plots, UMAPs |
| `results/multiome/*/visualizations/` | MOFA and MultiVI plots |
| `logs/nextflow/report.html` | Per-process runtime and peak memory |

## See also

- [Regulatory analysis](regulatory.md) — what these figures are built from
- [On-ramps & resuming](../onramps.md)
