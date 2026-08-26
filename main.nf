#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================================
// REFACTOR5: Unified Multiomics Pipeline v3.0.0
// Merges BD multi-sample and PBMC single-sample pipelines into one.
// Handles both single-sample (10x-style) and multi-sample (BD-style)
// datasets through config parameterization.
//
// Key changes from v2.x:
//   - Generic batch_dirs map replaces hardcoded june/july/nov lookups
//   - cell_type_key computed ONCE at pipeline start
//   - On-ramp entry points for pre-computed intermediates
//   - Fragment files passed as channel input to REGULATORY_ANALYSIS
//   - Both PBMC and BD enhancer footprinting recipes included
//   - BD multiome join diagnostics preserved (fail-fast with orphan warnings)
//   - MAP_TF_TO_TARGET_GENES included
//   - SCENICPLUS_VISUALIZE included
// ============================================================================

// ============================================================================
// RNA MODULES
// ============================================================================
include { CELLBENDER } from './modules/rna/cellbender'
include { RNA_QC } from './modules/rna/qc'
include { CONCAT_BATCHES } from './modules/rna/concat'

// ============================================================================
// RNA INTEGRATION MODULES
// ============================================================================
include { PREPARE_REFERENCE } from './modules/integration/prepare_reference'
include { TRAIN_SCVI } from './modules/integration/scvi'
include { TRAIN_SCANVI } from './modules/integration/scanvi'
include { RUN_CELLTYPIST } from './modules/cellannotator/celltypist'
include { RUN_MARKER_ANNOTATION } from './modules/cellannotator/marker_annotation'

// ============================================================================
// VISUALIZATION MODULES
// ============================================================================
include { PLOT_POST_SCANVI } from './modules/visualization/plot_post_scanvi'
include { RUN_CELLCHAT } from './modules/cellchat/cellchat'
include { CELLCHAT_PER_CONDITION } from './modules/cellchat/cellchat_compare'
include { CELLCHAT_COMPARE } from './modules/cellchat/cellchat_compare'
include { CONVERT_H5AD_TO_SEURAT } from './modules/conversion/h5ad_to_seurat'
include { HDWGCNA_PER_CELLTYPE } from './modules/hdwgcna/hdwgcna.nf'
include { HDWGCNA_DIFFERENTIAL } from './modules/hdwgcna/hdwgcna_differential'
include { HDWGCNA_ENRICHMENT } from './modules/hdwgcna/hdwgcna_differential'

// ============================================================================
// RNA DIFFERENTIAL EXPRESSION MODULES
// ============================================================================
include { ASSIGN_TEST_GROUPS } from './modules/rna/differential_expression'
include { CONVERT_H5AD_FOR_MAST } from './modules/rna/differential_expression'
include { EXTRACT_CELL_TYPES_FOR_MAST } from './modules/rna/differential_expression'
include { RUN_MAST_DE } from './modules/rna/differential_expression'
include { CREATE_VOLCANO_PLOTS } from './modules/rna/differential_expression'
include { RUN_GO_ENRICHMENT } from './modules/rna/differential_expression'

// ============================================================================
// ATAC MODULES (UNIFIED)
// ============================================================================
include { ATAC_INITIAL_QC; ATAC_MAKE_THRESHOLDS } from './modules/atac/atac_initial_qc'
include { ATAC_FINAL_PIPELINE } from './modules/atac/consolidated_pipeline'

// Cell-type annotation for ATAC
include { ATAC_CELLTYPE_ANNOTATION } from './modules/atac/celltype_annotation'
include { ATAC_SCATANNO } from './modules/atac/atac_scatanno'
include { MERGE_ANNOTATIONS } from './modules/atac/merge_annotations'

// Differential ATAC analysis
include { EXTRACT_ATAC_CELL_TYPES } from './modules/atac/snapatac_diff'
include { SNAPATAC_DIFFERENTIAL } from './modules/atac/snapatac_diff'

// ============================================================================
// REGULATORY ANALYSIS MODULES
// ============================================================================

// Cicero modules
include { CICERO_TRIPLETS } from './modules/cicero/cicero_triplets'
include { CICERO_ESTIMATE_DP } from './modules/cicero/cicero_estimate_dp'
include { CICERO_ESTIMATE_DP as CICERO_ESTIMATE_DP_CTRL } from './modules/cicero/cicero_estimate_dp'
include { CICERO_ESTIMATE_DP as CICERO_ESTIMATE_DP_TRT  } from './modules/cicero/cicero_estimate_dp'
include { CICERO_FULL_CHROM  } from './modules/cicero/cicero_full_chrom'
include { CICERO_FULL_CHROM as CICERO_FULL_CHROM_CTRL } from './modules/cicero/cicero_full_chrom'
include { CICERO_FULL_CHROM as CICERO_FULL_CHROM_TRT  } from './modules/cicero/cicero_full_chrom'
include { CICERO_JOIN        } from './modules/cicero/cicero_join'
include { CICERO_JOIN as CICERO_JOIN_CTRL } from './modules/cicero/cicero_join'
include { CICERO_JOIN as CICERO_JOIN_TRT  } from './modules/cicero/cicero_join'
include { CICERO_TARGET_PLOTS } from './modules/cicero/cicero_target_plots'

// ChromVAR modules
include { GPU_CHROMVAR } from './modules/chromvar/gpu_chromvar'
include { VIS_CHROMVAR } from './modules/chromvar/vischromvar'
include { EXTRACT_CHROMVAR_MOTIFS } from './modules/chromvar/extract_motifs'
include { MAP_TF_TO_TARGET_GENES } from './modules/chromvar/map_tf_targets'
// Shi et al. TF accessibility tester (differential or descriptive mode)
include { DIFFERENTIAL_TF_ACCESSIBILITY } from './modules/chromvar/differential_tf_accessibility'

// Shi et al. stratified Cicero + co-accessibility comparison (absorbed from SDas)
// Stratified legs reuse the per-chromosome path via CICERO_ESTIMATE_DP_CTRL/TRT
// + CICERO_FULL_CHROM_CTRL/TRT + CICERO_JOIN_CTRL/TRT (aliased above).
// The monolithic CICERO_FULL_STRATIFIED process in cicero_stratified.nf is
// intentionally NOT included; kept as dead code for backward reference.
include { CICERO_TRIPLETS_STRATIFIED } from './modules/cicero/cicero_stratified'
include { CICERO_TRIPLETS_STRATIFIED as CICERO_TRIPLETS_STRATIFIED_TRT } from './modules/cicero/cicero_stratified'
include { COMPARE_COACCESSIBILITY } from './modules/cicero/cicero_stratified'
include { CICERO_TRIPLETS_PER_CT    } from './modules/cicero/cicero_per_ct'
include { CICERO_ESTIMATE_DP_PER_CT } from './modules/cicero/cicero_per_ct'
include { CICERO_FULL_CHROM_PER_CT  } from './modules/cicero/cicero_per_ct'
include { CICERO_JOIN_PER_CT        } from './modules/cicero/cicero_per_ct'
include { BUILD_CT_ANNOTATION       } from './modules/cicero/build_ct_annotation'

// scPRINT modules
include { SCPRINTER_BARCODES } from './modules/scprint/barcodes'
include { RESOLVE_GENE_COORDINATES                              } from './modules/scprint/resolve_coordinates'
include { RESOLVE_GENE_COORDINATES as RESOLVE_OVERLAY_COORDINATES } from './modules/scprint/resolve_coordinates'
include { SCPRINTER_BUILD_PRINTER } from './modules/scprint/build_printer'
include { SCPRINTER_FOOTPRINTING } from './modules/scprint/footprinter'
include { SCPRINTER_FOOTPRINTING as SCPRINTER_FOOTPRINTING_DIFF } from './modules/scprint/footprinter'
include { SCPRINTER_MOTIF_SCAN } from './modules/scprint/motif_scan'

// ============================================================================
// MULTIOME INTEGRATION MODULES
// ============================================================================
include { BUILD_MUDATA } from './modules/multiome/build_mudata'

// MOFA integration - supports both modes
include { MOFA_INTEGRATE } from './modules/multiome/mofa_integrate'
include { MOFA_VISUALIZE } from './modules/multiome/mofa_visualize'

// Bootstrap MOFA components (preserved for low-memory systems)
include { BOOTSTRAP_MOFA_INTEGRATION } from './modules/multiome/bootstrap_mofa'
include { CONSENSUS_ANALYSIS } from './modules/multiome/bootstrap_mofa'
include { ANALYZE_MEMORY_LOG } from './modules/multiome/bootstrap_mofa'

// Optional cellismo conversion
include { CONVERT_CELLISMO; CONCAT_CELLISMO } from './modules/multiome/convert_cellismo'

// MultiVI integration
include { MULTIVI_INTEGRATE } from './modules/multiome/multivi_integrate'
include { MULTIVI_VISUALIZE } from './modules/multiome/multivi_visualize'

// MultiVI downstream analysis
include { MULTIVI_MASKING_SWEEP_ONE; MULTIVI_MASKING_SWEEP_AGGREGATE } from './modules/multiome/multivi_masking_sweep'
include { MULTIVI_DRIVER_FACTORS } from './modules/multiome/multivi_driver_factors'
include { MULTIVI_GAP_FILL } from './modules/multiome/multivi_gap_fill'
include { MULTIVI_VALIDATE } from './modules/multiome/multivi_validate'

// ============================================================================
// PYcistopic / SCENIC+ / DORC MODULES
// ============================================================================
include { PYCISTOPIC_PHASE1         } from './modules/multiome/pycistopic_phase1'
include { PYCISTOPIC_ATAC_PREPARE   } from './modules/multiome/pycistopic_prepare_atac'
include { PYCISTOPIC_PER_GROUP      } from './modules/multiome/pycistopic_per_group'
include { PYCISTOPIC_MERGE_OBJECTS  } from './modules/multiome/pycistopic_merge_objects'
include { PYCISTOPIC_RUN_LDA        } from './modules/multiome/pycistopic_run_lda'
include { PYCISTOPIC_FINALIZE_LDA   } from './modules/multiome/pycistopic_finalize_lda'
include { SCENICPLUS_RUN       } from './modules/multiome/scenicplus_run'
include { SCENICPLUS_VISUALIZE } from './modules/multiome/scenicplus_visualize'
include { SCENICPLUS_GRN_VIZ  } from './modules/multiome/scenicplus_grn_viz'
include { EXPORT_MUDATA_RNA } from './modules/multiome/export_mudata_rna'
include { SCPRINTER_DORC       } from './modules/multiome/scprinter_dorc'

// ============================================================================
// ENHANCER FOOTPRINTING RECIPE MODULES (Recipes A/B/C/D)
// ============================================================================

// Phase 1: ATAC-Only Enhancer Footprinting (Recipe A)
include { EXTRACT_CCAN_ENHANCERS } from './modules/scprint/extract_ccan_enhancers'
include { MOTIF_SCAN_ENHANCERS   } from './modules/scprint/motif_scan_enhancers'
include { ENHANCER_FOOTPRINTING  } from './modules/scprint/enhancer_footprinting'
// Phase 3: ATAC-only TF-gene regulatory network (Shi et al. TF_Net equivalent,
// continuous scPrinter binding × Cicero co-accessibility). Absorbed from SDas.
include { BUILD_TF_GENE_NETWORK } from './modules/scprint/build_tf_gene_network'
include { PLOT_TF_GENE_NETWORK  } from './modules/scprint/plot_tf_gene_network'

// Phase 2: Multiome Integration (Recipe B)
include { EXTRACT_EREGULON_REGIONS } from './modules/multiome/extract_eregulon_regions'
include { CROSS_MODAL_VALIDATION   } from './modules/multiome/cross_modal_validation'

// Phase 3: CellChat-Guided Footprinting (Recipe C)
include { CELLCHAT_TO_TF_HYPOTHESES } from './modules/cellchat/cellchat_to_tf_hypotheses'
include { EXTRACT_SIGNALING_TARGETS } from './modules/cellchat/extract_signaling_targets'
include { SIGNAL_CHAIN_CORRELATION  } from './modules/cellchat/signal_chain_correlation'

// Phase 4: Composite Enhancer Visualization (Recipe D)
include { PREPARE_ENHANCER_VIZ_TRACKS } from './modules/visualization/enhancer_viz'
include { COMPOSITE_ENHANCER_VIZ      } from './modules/visualization/enhancer_viz'

// D1b: B-tier modules absorbed from SDas_nf
include { POST_QC_REPORT               } from './modules/atac/post_qc_report'
include { ATAC_DESCRIPTIVE_REPORT      } from './modules/atac/descriptive_report'
include { ENHANCER_FOOTPRINTING_PER_CT        } from './modules/scprint/enhancer_footprinting_per_ct'
include { ENHANCER_FOOTPRINTING_PER_CT_STRIP  } from './modules/scprint/enhancer_footprinting_per_ct_strip'
include { RANK_ENHANCER_STRIP_GENES           } from './modules/scprint/rank_enhancer_strip_genes'
include { BUILD_VIZ_CANDIDATES         } from './modules/visualization/build_viz_candidates'
include { AGGREGATE_FP_STATS           } from './modules/visualization/aggregate_fp_stats'
include { EXPORT_ATAC_BIGWIGS          } from './modules/visualization/export_atac_bigwigs'
include { PROMOTER_MSFP_PER_CT         } from './modules/scprint/promoter_msfp_per_ct'
include { RENDER_PROMOTER_MSFP_OVERLAY } from './modules/visualization/render_promoter_msfp_overlay'
// Cis-rewiring (2026-05-06) — per-condition stratified CCAN extraction +
// union-peakset motif scan + per-TF gained-enhancer motif-presence panels.
include { EXTRACT_CCAN_ENHANCERS as EXTRACT_CCAN_ENHANCERS_CTRL } from './modules/scprint/extract_ccan_enhancers'
include { EXTRACT_CCAN_ENHANCERS as EXTRACT_CCAN_ENHANCERS_TRT  } from './modules/scprint/extract_ccan_enhancers'
include { BUILD_UNION_ENHANCER_PEAKS    } from './modules/scprint/build_union_enhancer_peaks'
include { MOTIF_SCAN_ENHANCERS_UNION    } from './modules/scprint/motif_scan_enhancers_union'
include { MOTIF_IN_GAINED_CCANS           } from './modules/visualization/motif_in_gained_ccans'
include { RENDER_CIS_REWIRING_MOTIF_STACK } from './modules/visualization/render_cis_rewiring_motif_stack'
include { RENDER_MSFP_PROMOTER_STRIP      } from './modules/visualization/render_msfp_promoter_strip'
include { RENDER_MSFP_ENHANCER_STRIP      } from './modules/visualization/render_msfp_enhancer_strip'
include { RENDER_GENOME_BROWSER           } from './modules/visualization/render_genome_browser'
include { RENDER_CICERO_LOLLIPOP          } from './modules/visualization/render_cicero_lollipop'

// SHI_FIGURES — Shi et al. 2025 figure equivalents (1E, 2B-E, 4B-E, 5A-F)
// Tier A (single-condition compatible) + Tier B (require >=2 conditions).
include { ANNOTATE_PEAK_TYPES       } from './modules/visualization/annotate_peak_types'
include { NMF_ENHANCER_PROGRAMS     } from './modules/visualization/nmf_enhancer_programs'
include { MARKER_COVERAGE_TRACKS    } from './modules/visualization/marker_coverage'
include { SELECT_SHI_CANDIDATES     } from './modules/visualization/select_shi_candidates'
include { COACC_CORRELATION_MATRIX  } from './modules/visualization/coacc_correlation_matrix'
include { DA_PEAK_BREAKDOWN         } from './modules/visualization/da_peak_breakdown'
include { DA_LOG2FC_HEATMAPS        } from './modules/visualization/da_log2fc_heatmaps'
include { TF_DIFFERENTIAL_VOLCANO   } from './modules/visualization/tf_differential_volcano'
include { CURATED_TF_NETWORKS       } from './modules/visualization/curated_tf_networks'
include { LOCUS_TF_BINDING          } from './modules/visualization/locus_tf_binding'


// ============================================================================
// GLOBAL: Compute cell_type_key ONCE
// ============================================================================
def has_reference = (params.species == 'human' && params.ref_dir_human_integrated) ||
                    (params.species == 'mouse' && params.ref_dir_mouse_integrated)
// Unified canonical cell-type column. BUILD_MUDATA and PLOT_POST_SCANVI write
// 'cell_type' by precedence (scanvi > celltypist > marker) and stamp provenance
// in 'cell_type_source'. Downstream modules use this single key regardless of
// which annotation tool ran.
//
// PLOT_POST_SCANVI runs BEFORE BUILD_MUDATA unifies, so the key it receives
// reflects the column the current annotation tool wrote on the h5ad.
def cell_type_key = (params.rna?.annotation_method == 'markers')
    ? 'cell_type_marker' : 'cell_type'

// ATAC cell type column: depends on annotation mode
// marker_file → 'cell_type', scATAnno → 'cell_type_prediction', CellTypist → 'celltypist_prediction'
def atac_cell_type_key = params.atac.marker_file ? 'cell_type' :
    (params.atac.annotation_method == 'scatanno' ? 'cell_type_prediction' : 'celltypist_prediction')

// Broad ATAC cell type column: condensed re-mapping written by MERGE_ANNOTATIONS
// (bin/merge_annotations.py applies CELLTYPIST_BROAD_MAP / scatanno_broad_map and
// stores the result under this column on peak_matrix_annotated.h5ad). Consumed by
// SHI Tier B (NMF_ENHANCER_PROGRAMS) which pseudobulks across broad classes.
//
// TODO (QOL refactor): replace with channel-driven contract — have MERGE_ANNOTATIONS
// emit a sidecar (e.g. cell_type_broad_col.txt) carrying the column name it just
// wrote, then thread that as a value channel through SHI_FIGURES.take. Producer
// becomes the single source of truth and atlas/column-name swaps need no Groovy edit.
def atac_broad_cell_type_key = params.atac?.broad_cell_type_key ?: 'cell_type_broad'


// ============================================================================
// HELPER: Resolve directory for a manifest row using generic batch_dirs map
// ============================================================================
// For RNA data_dir lookups
def resolveRnaDir(row) {
    if (row.data_dir) return row.data_dir
    def dir = params.batch_dirs?.get(row.batch, null)
    if (!dir) error "No directory configured for batch '${row.batch}'. Set batch_dirs.${row.batch} in config."
    // Support sub-directory pattern (e.g., nov batch with lane subdirs)
    if (row.original_lane_id && params.batch_dirs_use_lane_subdir?.contains(row.batch)) {
        return "${dir}/${row.original_lane_id}"
    }
    return dir
}

// For ATAC fragment directory lookups (barcode-sorted demux fragments)
def resolveAtacDir(row) {
    if (row.data_dir) return row.data_dir
    def dir = params.atac_batch_dirs?.get(row.batch, null)
    if (!dir) error "No ATAC directory configured for batch '${row.batch}'. Set atac_batch_dirs.${row.batch} in config."
    return dir
}

// For ATAC coord-sorted fragment directory lookups
def resolveAtacCoordDir(row) {
    if (row.coord_data_dir) return row.coord_data_dir
    def dir = params.atac_coord_batch_dirs?.get(row.batch, null)
    if (!dir) error "No ATAC coord-sorted directory configured for batch '${row.batch}'. Set atac_coord_batch_dirs.${row.batch} in config."
    return dir
}

// HELPER: Resolve per-sample demux files from generic config maps
// Returns file('NO_FILE') when no demux is configured for this sample's batch.
def resolveDemuxMetadata(String sample) {
    def match = params.demux_metadata?.find { batch, path -> sample.contains(batch) }
    return match ? file(match.value) : file('NO_FILE_METADATA')
}
def resolveDemuxSouporcell(String sample) {
    def match = params.demux_souporcell_dirs?.find { batch, path -> sample.contains(batch) }
    return match ? file(match.value) : file('NO_FILE_SOUPORCELL')
}


// ============================================================================
// FIX-P0-7 (A2): Normalize sample_type for case-insensitive matching
// ============================================================================
def normalizeSampleType(String raw) {
    if (!raw) return raw
    def v = raw.trim().toLowerCase()
    // Accept common aliases
    if (v in ['rna', 'lane'])  return 'lane'
    if (v in ['atac', 'demux']) return 'demux'
    return v
}

def isLane(row)  { normalizeSampleType(row.sample_type) == 'lane' }
def isDemux(row) { normalizeSampleType(row.sample_type) == 'demux' }

// FIX-A5: Trim whitespace from all CSV fields to prevent path resolution failures
// FIX-A6: Filter predicate to skip empty/blank rows from trailing newlines
def trimRow(row) { row.collectEntries { k, v -> [k, v?.trim()] } }
def isNonEmptyRow(row) { row.sample_id?.trim() }

// ============================================================================
// HELPER: Levenshtein distance for fuzzy column name matching (FIX-A1)
// ============================================================================
def levenshteinClose(String a, String b, int maxDist = 3) {
    a = a.toLowerCase(); b = b.toLowerCase()
    if (a == b) return false  // exact match is not a "close" suggestion
    int n = a.length(), m = b.length()
    if (Math.abs(n - m) > maxDist) return false
    int[][] d = new int[n + 1][m + 1]
    (0..n).each { d[it][0] = it }
    (0..m).each { d[0][it] = it }
    (1..n).each { i -> (1..m).each { j ->
        d[i][j] = [d[i-1][j] + 1, d[i][j-1] + 1,
                    d[i-1][j-1] + (a[i-1] == b[j-1] ? 0 : 1)].min()
    }}
    return d[n][m] <= maxDist && d[n][m] > 0
}

// ============================================================================
// FIX-P0: STARTUP VALIDATION (PRE-FLIGHT CHECKLIST)
// Catches genome build mismatches, species inconsistencies, manifest issues,
// file existence, parameter divergence, and cross-parameter consistency
// BEFORE any compute starts.
// ============================================================================
def validateStartupParams() {
    def errors = []
    def warnings = []
    def checks_passed = []
    def rna_rows_found = false   // MIS-14: set true when manifest contains a lane row with rna_file

    // --- P0-7 (A2): sample_type normalization is handled in-line at each
    //     splitCsv filter (see below). Here we warn if non-standard values exist.
    if (params.metadata_file) {
        def csv = file(params.metadata_file)
        if (csv.exists()) {
            def lines = csv.readLines()
            if (lines.size() > 1) {
                def header = lines[0].split(',').collect { it.trim() }
                def stIdx = header.findIndexOf { it == 'sample_type' }
                def sidIdx = header.findIndexOf { it == 'sample_id' }

                // FIX-A1: Full manifest CSV schema validation with fuzzy suggestions
                def required_base = ['sample_id', 'sample_type']
                def required_lane = ['rna_file']
                def required_demux = ['fragment_file']
                def all_known = ['sample_id', 'sample_type', 'data_dir', 'rna_file', 'fragment_file',
                                 'batch', 'condition_group', 'original_lane_id', 'coord_data_dir']
                def missing = required_base.findAll { col -> !header.contains(col) }
                if (missing) {
                    // Fuzzy match: suggest close column names
                    def suggestions = missing.collect { req ->
                        def close = header.findAll { h -> levenshteinClose(h, req) }
                        close ? "${req} (did you mean: ${close.join(', ')}?)" : req
                    }
                    errors << "Manifest CSV missing required columns: ${suggestions}. Found: ${header}"
                }

                // FIX-A1: Check lane-specific and demux-specific required columns
                def dataRows = lines[1..-1].findAll { it.trim() }
                if (stIdx >= 0) {
                    def types = dataRows.collect { it.split(',')[stIdx]?.trim() }.findAll { it }.unique()
                    def hasLane = types.any { normalizeSampleType(it) == 'lane' }
                    def hasDemux = types.any { normalizeSampleType(it) == 'demux' }

                    if (hasLane) {
                        def missingLane = required_lane.findAll { col -> !header.contains(col) }
                        if (missingLane) {
                            def sugg = missingLane.collect { req ->
                                def close = header.findAll { h -> levenshteinClose(h, req) }
                                close ? "${req} (did you mean: ${close.join(', ')}?)" : req
                            }
                            errors << "Manifest has 'lane' rows but missing required columns: ${sugg}. Found: ${header}"
                        }
                    }
                    if (hasDemux) {
                        def missingDemux = required_demux.findAll { col -> !header.contains(col) }
                        if (missingDemux) {
                            def sugg = missingDemux.collect { req ->
                                def close = header.findAll { h -> levenshteinClose(h, req) }
                                close ? "${req} (did you mean: ${close.join(', ')}?)" : req
                            }
                            errors << "Manifest has 'demux' rows but missing required columns: ${sugg}. Found: ${header}"
                        }
                    }

                    // P0-7 (A2): Warn on non-standard sample_type values
                    def nonStandard = types.findAll { !(it in ['lane', 'demux']) }
                    if (nonStandard) {
                        def couldBe = nonStandard.findAll { it.toLowerCase() in ['lane', 'demux', 'rna', 'atac'] }
                        if (couldBe) {
                            errors << "Manifest sample_type contains non-standard values: ${nonStandard}. " +
                                      "Expected 'lane' or 'demux' (case-sensitive). Did you mean: ${couldBe.collect { it.toLowerCase() == 'rna' ? 'lane' : it.toLowerCase() }}?"
                        }
                    }
                }

                // P0-3 (A3-partial): Check sample_id uniqueness per sample_type
                if (sidIdx >= 0 && stIdx >= 0) {
                    def sampleIds = dataRows.collect { line ->
                        def cols = line.split(',')
                        [cols[sidIdx]?.trim(), cols[stIdx]?.trim()]
                    }
                    def laneIds = sampleIds.findAll { it[1] == 'lane' }.collect { it[0] }
                    def demuxIds = sampleIds.findAll { it[1] == 'demux' }.collect { it[0] }
                    def laneDupes = laneIds.countBy { it }.findAll { k, v -> v > 1 }.keySet()
                    def demuxDupes = demuxIds.countBy { it }.findAll { k, v -> v > 1 }.keySet()
                    if (laneDupes) errors << "Duplicate sample_id in lane rows: ${laneDupes}"
                    if (demuxDupes) errors << "Duplicate sample_id in demux rows: ${demuxDupes}"
                }

                // P0-11 (A4) + MIS-15 (2026-05-04): Require condition_group column for ANY
                // workflow that consumes a condition axis (ATAC-diff, RNA-diff, TF-diff in
                // differential mode, stratified Cicero, SHI Tier B, disease-stratified enhancer).
                def cgIdx = header.findIndexOf { it == 'condition_group' }
                def needsConditionGroup = (params.differential?.run ?: false) ||
                    (params.differential_rna?.run ?: false) ||
                    (params.differential_tf?.run && params.differential_tf?.mode == 'differential') ||
                    (params.cicero?.stratified ?: false) ||
                    (params.shi_figures?.enabled && params.shi_figures?.treatment && params.shi_figures?.control) ||
                    (params.enhancer_footprinting?.disease_stratified ?: false)
                if (cgIdx < 0 && needsConditionGroup) {
                    errors << "Manifest CSV missing 'condition_group' column but a condition-aware " +
                              "workflow is enabled (differential / differential_rna / differential_tf / " +
                              "cicero.stratified / shi_figures Tier B / disease_stratified). " +
                              "Add a 'condition_group' column to the manifest or disable these workflows."
                } else if (cgIdx >= 0) {
                    // Check for empty/missing condition_group values
                    def emptyRows = dataRows.findAll { line ->
                        def cols = line.split(',')
                        cgIdx >= cols.size() || !cols[cgIdx]?.trim()
                    }
                    if (emptyRows) {
                        warnings << "Found ${emptyRows.size()} manifest row(s) with empty condition_group. " +
                                    "These will default to 'Control' in differential analysis."
                    }
                }

                // FIX-A7 + FIX-D1: Validate file paths exist (works in both preview and production mode)
                def ddIdx = header.findIndexOf { it == 'data_dir' }
                def rfIdx = header.findIndexOf { it == 'rna_file' }
                def ffIdx = header.findIndexOf { it == 'fragment_file' }
                def btIdx = header.findIndexOf { it == 'batch' }

                dataRows.each { line ->
                    def cols = line.split(',').collect { it.trim() }
                    def stype = stIdx >= 0 && stIdx < cols.size() ? normalizeSampleType(cols[stIdx]) : null
                    def sid = sidIdx >= 0 && sidIdx < cols.size() ? cols[sidIdx] : 'unknown'

                    if (stype == 'lane' && rfIdx >= 0) {
                        // Resolve RNA file path
                        def dataDir = (ddIdx >= 0 && ddIdx < cols.size() && cols[ddIdx]) ? cols[ddIdx] : null
                        def batch = (btIdx >= 0 && btIdx < cols.size() && cols[btIdx]) ? cols[btIdx] : null
                        if (!dataDir && batch) dataDir = params.batch_dirs?.get(batch, null)
                        def rnaFname = rfIdx < cols.size() ? cols[rfIdx] : null
                        if (rnaFname) { rna_rows_found = true }   // MIS-14
                        if (dataDir && rnaFname) {
                            def rnaPath = file("${dataDir}/${rnaFname}")
                            if (!rnaPath.exists()) {
                                errors << "RNA file not found for sample '${sid}': ${rnaPath}"
                            } else {
                                // FIX-D1: If it's a directory (MEX), check required files
                                if (rnaPath.isDirectory()) {
                                    def mexRequired = ['matrix.mtx.gz', 'barcodes.tsv.gz', 'features.tsv.gz']
                                    def mexMissing = mexRequired.findAll { !file("${rnaPath}/${it}").exists() }
                                    if (mexMissing) {
                                        errors << "MEX directory for sample '${sid}' missing files: ${mexMissing} in ${rnaPath}"
                                    }
                                }
                                // FIX-D2: Probe h5 file format (race-safe: open, check, close immediately)
                                else if (rnaFname.endsWith('.h5')) {
                                    try {
                                        // Quick header probe — just check the file is readable and non-empty
                                        if (rnaPath.size() < 100) {
                                            warnings << "RNA h5 file for sample '${sid}' is suspiciously small (${rnaPath.size()} bytes): ${rnaPath}"
                                        }
                                    } catch (Exception e) {
                                        warnings << "Could not probe RNA h5 file for sample '${sid}': ${e.message}"
                                    }
                                }
                            }
                        }
                    }

                    if (stype == 'demux' && ffIdx >= 0) {
                        // Resolve fragment file path
                        def dataDir = (ddIdx >= 0 && ddIdx < cols.size() && cols[ddIdx]) ? cols[ddIdx] : null
                        def batch = (btIdx >= 0 && btIdx < cols.size() && cols[btIdx]) ? cols[btIdx] : null
                        if (!dataDir && batch) dataDir = params.atac_fragment_dirs_bc?.get(batch, null)
                        def fragFname = ffIdx < cols.size() ? cols[ffIdx] : null
                        if (dataDir && fragFname) {
                            def fragFull = fragFname.contains('.') ? fragFname : "${fragFname}.bed.gz"
                            def fragPath = file("${dataDir}/${fragFull}")
                            if (!fragPath.exists()) {
                                errors << "ATAC fragment file not found for sample '${sid}': ${fragPath}"
                            } else {
                                // FIX-C5: Validate fragment file column format (first non-comment line)
                                try {
                                    def firstLine = null
                                    def is = fragPath.newInputStream()
                                    def stream = fragFull.endsWith('.gz') ?
                                        new java.util.zip.GZIPInputStream(is) : is
                                    def br = new java.io.BufferedReader(new java.io.InputStreamReader(stream))
                                    try {
                                        def l
                                        while ((l = br.readLine()) != null) {
                                            if (!l.startsWith('#')) { firstLine = l; break }
                                        }
                                    } finally {
                                        br.close()
                                    }
                                    if (firstLine) {
                                        def fields = firstLine.split('\t')
                                        if (fields.size() < 4) {
                                            errors << "Fragment file for sample '${sid}' has ${fields.size()} tab-separated columns " +
                                                      "(expected >= 4: chr, start, end, barcode[, count]): ${fragPath}"
                                        }
                                        // FIX-C1: Check position-sorted (start should be numeric)
                                        try {
                                            Long.parseLong(fields[1])
                                        } catch (NumberFormatException e) {
                                            errors << "Fragment file for sample '${sid}' column 2 is not numeric " +
                                                      "(expected position-sorted BED format): '${fields[1]}' in ${fragPath}"
                                        }
                                    }
                                } catch (Exception e) {
                                    warnings << "Could not probe fragment file for sample '${sid}': ${e.message}"
                                }
                            }
                        }
                    }
                }

                checks_passed << "Manifest schema (${dataRows.size()} rows)"
            }
        } else {
            errors << "Manifest CSV not found: ${params.metadata_file}"
        }
    }

    // FIX-C3: ATAC sample_metadata consistency with main manifest
    if (params.atac?.sample_metadata && params.metadata_file) {
        def atacMeta = file(params.atac.sample_metadata)
        def mainCsv = file(params.metadata_file)
        if (atacMeta.exists() && mainCsv.exists()) {
            try {
                def mainLines = mainCsv.readLines()
                def mainHeader = mainLines[0].split(',').collect { it.trim() }
                def mainStIdx = mainHeader.findIndexOf { it == 'sample_type' }
                def mainSidIdx = mainHeader.findIndexOf { it == 'sample_id' }
                if (mainStIdx >= 0 && mainSidIdx >= 0) {
                    def demuxSids = mainLines[1..-1].findAll { it.trim() }
                        .collect { it.split(',') }
                        .findAll { cols -> mainStIdx < cols.size() && normalizeSampleType(cols[mainStIdx]?.trim()) == 'demux' }
                        .collect { cols -> mainSidIdx < cols.size() ? cols[mainSidIdx]?.trim() : null }
                        .findAll { it }
                        .toSet()

                    def atacLines = atacMeta.readLines()
                    def atacHeader = atacLines[0].split(',').collect { it.trim() }
                    def atacSidIdx = atacHeader.findIndexOf { it == 'sample_id' } ?:
                                     atacHeader.findIndexOf { it == 'sample' }
                    if (atacSidIdx >= 0 && atacLines.size() > 1) {
                        def atacSids = atacLines[1..-1].findAll { it.trim() }
                            .collect { it.split(',') }
                            .collect { cols -> atacSidIdx < cols.size() ? cols[atacSidIdx]?.trim() : null }
                            .findAll { it }
                            .toSet()
                        def inAtacNotMain = atacSids - demuxSids
                        def inMainNotAtac = demuxSids - atacSids
                        if (inAtacNotMain) {
                            warnings << "ATAC sample_metadata has ${inAtacNotMain.size()} sample(s) not in main manifest demux rows: " +
                                        "${inAtacNotMain.take(5)}${inAtacNotMain.size() > 5 ? '...' : ''}"
                        }
                        if (inMainNotAtac) {
                            warnings << "Main manifest has ${inMainNotAtac.size()} demux sample(s) not in ATAC sample_metadata: " +
                                        "${inMainNotAtac.take(5)}${inMainNotAtac.size() > 5 ? '...' : ''}"
                        }
                        if (!inAtacNotMain && !inMainNotAtac) {
                            checks_passed << "ATAC sample_metadata consistency (${atacSids.size()} samples match)"
                        }
                    }
                }
            } catch (Exception e) {
                warnings << "Could not validate ATAC sample_metadata consistency: ${e.message}"
            }
        }
    }

    // --- MIS-01 (2026-05-04): species must be explicitly set
    def allowedSpecies = ['human', 'mouse']
    if (!params.species) {
        errors << "params.species is unset. Set to 'human' or 'mouse' in your dataset " +
                  "config (configs/datasets/<dataset>.config). The silent 'human' default " +
                  "was removed because it produced misannotated mouse runs."
    } else if (!(params.species in allowedSpecies)) {
        errors << "params.species='${params.species}' is invalid. Allowed: ${allowedSpecies}."
    }

    // --- P0-6 (B2) + P0-8 (E1/E2): Species consistency
    def speciesMap = [human: 'hsapiens', mouse: 'mmusculus']
    def expectedPycisSpecies = speciesMap[params.species]
    if (expectedPycisSpecies && params.pycistopic?.species) {
        if (params.pycistopic.species != expectedPycisSpecies) {
            errors << "Species mismatch: params.species='${params.species}' implies " +
                      "pycistopic.species='${expectedPycisSpecies}' but got '${params.pycistopic.species}'. " +
                      "This will cause pycisTopic to use wrong genome annotations."
        }
    }

    // P0-4 (B1): GTF genome build consistency check
    // Check that GTF path contains species-consistent build string
    def gtfPath = params.species == 'human' ? params.gtf_human_full : params.gtf_mouse_full
    if (gtfPath) {
        def gtfStr = gtfPath.toString().toLowerCase()
        if (params.species == 'human' && (gtfStr.contains('mm10') || gtfStr.contains('mm39') || gtfStr.contains('grcm'))) {
            errors << "GTF genome build mismatch: params.species='human' but GTF path contains mouse build: ${gtfPath}"
        }
        if (params.species == 'mouse' && (gtfStr.contains('hg38') || gtfStr.contains('hg19') || gtfStr.contains('grch'))) {
            errors << "GTF genome build mismatch: params.species='mouse' but GTF path contains human build: ${gtfPath}"
        }
    }

    // P0-4 (B1): pycisTopic GTF check
    if (params.pycistopic?.gtf) {
        def ptGtf = params.pycistopic.gtf.toString().toLowerCase()
        if (params.species == 'human' && (ptGtf.contains('mm10') || ptGtf.contains('mm39') || ptGtf.contains('grcm'))) {
            errors << "pycisTopic GTF contains mouse build string but params.species='human': ${params.pycistopic.gtf}"
        }
        if (params.species == 'mouse' && (ptGtf.contains('hg38') || ptGtf.contains('hg19') || ptGtf.contains('grch'))) {
            errors << "pycisTopic GTF contains human build string but params.species='mouse': ${params.pycistopic.gtf}"
        }
    }

    // P1-20 (R1-23): pycistopic.gtf must be set when pycistopic.run is enabled
    if (params.pycistopic?.run && !params.pycistopic?.gtf) {
        errors << "pycistopic.run=true but pycistopic.gtf is null. " +
                  "Provide a Gencode GTF path in your dataset config (e.g., pycistopic { gtf = '/path/to/gencode.gtf' }). " +
                  "Without it, pycisTopic falls back to BioMart which fails due to pyarrow incompatibility."
    }

    // P0-9 (B4): SCENIC+ cisTarget species check
    if (params.scenicplus?.cistarget_rankings) {
        def ctStr = params.scenicplus.cistarget_rankings.toString().toLowerCase()
        if (params.species == 'human' && (ctStr.contains('mm10') || ctStr.contains('mm9'))) {
            errors << "SCENIC+ cisTarget rankings appear to be for mouse but params.species='human': ${params.scenicplus.cistarget_rankings}"
        }
        if (params.species == 'mouse' && (ctStr.contains('hg38') || ctStr.contains('hg19'))) {
            errors << "SCENIC+ cisTarget rankings appear to be for human but params.species='mouse': ${params.scenicplus.cistarget_rankings}"
        }
    }

    // P0-13 (B5): Blacklist BED genome build check
    def blacklist = params.blacklist_bed ?: params.pycistopic?.blacklist_bed
    if (blacklist) {
        def blStr = blacklist.toString().toLowerCase()
        if (params.species == 'human' && (blStr.contains('mm10') || blStr.contains('mm9'))) {
            errors << "Blacklist BED appears to be for mouse but params.species='human': ${blacklist}"
        }
        if (params.species == 'mouse' && (blStr.contains('hg38') || blStr.contains('hg19'))) {
            errors << "Blacklist BED appears to be for human but params.species='mouse': ${blacklist}"
        }
    }

    // P0-14 (E3): CellBender expected_cells sanity check
    if (params.cellbender?.expected_cells) {
        def ec = params.cellbender.expected_cells as int
        if (ec > 50000) {
            warnings << "CellBender expected_cells=${ec} is very high (>50000). " +
                        "Verify this matches your data or ambient RNA removal may be degraded."
        }
        if (ec < 100) {
            warnings << "CellBender expected_cells=${ec} is very low (<100). " +
                        "Verify this matches your data or ambient RNA removal may be degraded."
        }
    }
    checks_passed << "Species/genome consistency"

    // R3-4: GTF file validation for all GTF-consuming workflows
    // Checks: (1) not the literal string "null", (2) file exists on disk
    def gtfChecks = [:]  // label → path
    def speciesGtf = params.species == 'human' ? params.gtf_human_full : params.gtf_mouse_full

    // scPRINTER / enhancer / chromvar workflows all use scprinter.gtf_human/mouse
    def scprinterGtf = params.species == 'human' ? params.scprinter?.gtf_human : params.scprinter?.gtf_mouse
    if (params.scprinter?.run || params.enhancer_footprinting?.run || params.enhancer_viz?.run || params.chromvar?.run) {
        gtfChecks["scprinter.gtf_${params.species} (scPRINTER/enhancer/chromVAR)"] = scprinterGtf
    }

    // Cicero
    if (params.cicero?.run) {
        gtfChecks["cicero.gtf_full"] = params.cicero?.gtf_full
        if (params.cicero?.target_genes) {
            gtfChecks["cicero.gtf_plot"] = params.cicero?.gtf_plot
        }
    }

    // pycisTopic (already has null check at line 520, add file-existence here)
    if (params.pycistopic?.run && params.pycistopic?.gtf) {
        gtfChecks["pycistopic.gtf"] = params.pycistopic.gtf
    }

    // SCENIC+
    if (params.scenicplus?.run && params.scenicplus?.gtf) {
        gtfChecks["scenicplus.gtf"] = params.scenicplus.gtf
    }

    // Validate each GTF
    gtfChecks.each { label, gtfVal ->
        def gtfString = gtfVal?.toString()
        if (!gtfString || gtfString == 'null') {
            errors << "${label} resolves to '${gtfString}'. " +
                      "Set the species-appropriate GTF path in your dataset config " +
                      "(e.g., gtf_human_full = '/path/to/gencode.v44.annotation.gtf')."
        } else if (!file(gtfString).exists()) {
            errors << "${label} file does not exist: ${gtfString}. " +
                      "Verify the path and ensure the GTF is accessible from the execution host."
        }
    }
    def gtfErrors = gtfChecks.count { label, gtfVal ->
        def s = gtfVal?.toString(); !s || s == 'null' || !file(s).exists()
    }
    if (gtfChecks && gtfErrors == 0) {
        checks_passed << "GTF files (${gtfChecks.size()} paths validated)"
    }

    // FIX-B3: Reference atlas directory must contain h5ad files if reference path is set
    def refDir = params.species == 'human' ? params.ref_dir_human_integrated : params.ref_dir_mouse_integrated
    if (refDir) {
        def refPath = file(refDir)
        if (!refPath.exists()) {
            errors << "Reference atlas directory does not exist: ${refDir}"
        } else if (refPath.isDirectory()) {
            def h5adFiles = []
            refPath.eachFileMatch(~/.*\.h5ad/) { h5adFiles << it.name }
            if (!h5adFiles) {
                errors << "Reference atlas directory contains no .h5ad files: ${refDir}. " +
                          "scANVI annotation requires at least one reference h5ad file."
            } else {
                checks_passed << "Reference atlas (${h5adFiles.size()} h5ad files)"
            }
        }
    }

    // MIS-17 (2026-05-04): scATAnno reference atlas must be set + exist when used.
    // Previously instance configs defaulted reference_atlas to mouse_kidney_*.h5ad,
    // which silently misannotated other tissues. Validator now requires explicit
    // path and warns on tissue_type ↔ filename mismatch.
    if (params.atac?.run && params.atac?.auto_annotate &&
        params.atac?.annotation_method == 'scatanno' && !params.atac?.marker_file) {
        def atlas = params.scatanno?.reference_atlas
        if (!atlas) {
            errors << "atac.annotation_method='scatanno' requires params.scatanno.reference_atlas " +
                      "(path to a .h5ad reference). Set this explicitly in your dataset config."
        } else if (!file(atlas).exists()) {
            errors << "scATAnno reference atlas not found: ${atlas}"
        } else {
            checks_passed << "scATAnno reference atlas (${file(atlas).name})"
            // Tissue sanity check: warn (don't error) on obvious tissue keyword mismatch.
            def atlasName = atlas.toString().toLowerCase()
            def tissue = params.atac?.tissue_type?.toString()?.toLowerCase()
            def tissueKeywords = ['kidney', 'brain', 'pbmc', 'liver', 'lung', 'heart', 'gut', 'skin']
            def atlasTissues = tissueKeywords.findAll { atlasName.contains(it) }
            if (tissue && atlasTissues && !atlasTissues.contains(tissue)) {
                warnings << "scATAnno atlas filename suggests ${atlasTissues} but " +
                            "params.atac.tissue_type='${tissue}'. Verify this is the right atlas."
            }
        }
    }

    // FIX-E8: Validate mofa.mode against allowed values before compute
    if (params.mofa?.run) {
        def allowedModes = ['high_memory', 'bootstrap']
        if (!(params.mofa.mode in allowedModes)) {
            errors << "Invalid mofa.mode='${params.mofa.mode}'. Allowed values: ${allowedModes}. " +
                      "Check for typos in your dataset config."
        } else {
            checks_passed << "MOFA mode (${params.mofa.mode})"
        }
    }

    // FIX-F12: Validate scprinter.genome matches pipeline-wide species/genome
    if (params.scprinter?.run) {
        def expectedGenomes = [human: ['hg38', 'hg19', 'grch38', 'grch37'],
                               mouse: ['mm10', 'mm39', 'grcm38', 'grcm39']]
        def validGenomes = expectedGenomes[params.species] ?: []
        if (params.scprinter.genome && !(params.scprinter.genome.toLowerCase() in validGenomes)) {
            errors << "scprinter.genome='${params.scprinter.genome}' does not match params.species='${params.species}'. " +
                      "Expected one of: ${validGenomes}"
        } else {
            checks_passed << "scPRINTER genome (${params.scprinter.genome})"
        }
    }

    // FIX-B6: Log CellTypist model being used
    if (params.celltypist?.enabled) {
        def model = params.celltypist.model ?: 'default'
        checks_passed << "CellTypist model: ${model}"
        // Optional tissue-type consistency hint
        def brainModels = ['Developing_Human_Brain', 'Pan_Fetal_Human', 'Adult_Human_PrefrontalCortex']
        def immuneModels = ['Immune_All_Low', 'Immune_All_High', 'Pan_Immune_CellTypist']
        def tissue = params.atac?.tissue_type ?: 'unknown'
        if (tissue == 'brain' && immuneModels.any { model.contains(it) }) {
            warnings << "CellTypist model '${model}' appears immune-focused but tissue_type='brain'. " +
                        "Consider a brain-specific model (e.g., Developing_Human_Brain.pkl)."
        }
        if (tissue == 'pbmc' && brainModels.any { model.contains(it) }) {
            warnings << "CellTypist model '${model}' appears brain-focused but tissue_type='pbmc'. " +
                        "Consider an immune-specific model (e.g., Immune_All_Low.pkl)."
        }
    }

    // FIX-R1-7: Validate container files exist (derived from params.containers)
    // Skipped under -preview, which constructs the process graph but never
    // launches a task, so no container is ever entered. Demoting this to a
    // warning is what lets someone validate a fresh clone, or their own config,
    // before spending anything on image pulls.
    if (params.containers) {
        def uniqueSifs = params.containers.values().collect { it.toString() }.unique()
        def containerMissing = []
        uniqueSifs.each { sifPath ->
            def sif = file(sifPath)
            if (!sif.exists()) {
                containerMissing << sif.name
            }
        }
        if (containerMissing) {
            if (workflow.preview) {
                warnings << "Missing container files: ${containerMissing}. " +
                            "Not required for -preview, but a real run needs them in singularity_cache/."
            } else {
                errors << "Missing container files: ${containerMissing}. " +
                          "Expected in singularity_cache/. Run container build/pull first."
            }
        } else {
            checks_passed << "Containers (${uniqueSifs.size()} SIF files)"
        }
    }

    // --- MIS-11 / MIS-25 (2026-05-04): resource_tier typo guard.
    // 'auto' is a documented alias for 'small'; case-sensitive match catches
    // typos like 'Medium' that previously fell through silently to small.
    // 'test'     — minimal tier used by -profile test for -preview (1 CPU / 1 GB).
    // 'tutorial' — CPU-only tier for the tier-2 tutorial dataset; real work, no GPU.
    def allowedTiers = ['small', 'medium', 'large', 'auto', 'test', 'tutorial']
    if (params.resource_tier && !(params.resource_tier in allowedTiers)) {
        errors << "params.resource_tier='${params.resource_tier}' is invalid (case-sensitive). " +
                  "Allowed: ${allowedTiers}. Check for typos like capital letters."
    } else if (params.resource_tier) {
        checks_passed << "Resource tier (${params.resource_tier})"
    }

    // --- MIS-07 (2026-05-04): differential.run=true requires non-empty comparisons
    if (params.differential?.run) {
        def cmps = params.differential?.comparisons
        if (!cmps || (cmps instanceof List && cmps.isEmpty())) {
            errors << "differential.run=true but differential.comparisons=[]. " +
                      "Set comparisons to a list of [treatment, control] pairs " +
                      "(e.g. [['TG','WT']]) or set differential.run=false."
        }
    }

    // --- MIS-14 (2026-05-04): rna.run=true requires manifest with at least one RNA row.
    // rna_rows_found is set during manifest parsing above (lane rows with non-null rna_file).
    // Skip when no metadata_file is set — that's a different mode (covered by other checks).
    if (params.rna?.run && params.metadata_file && !rna_rows_found) {
        errors << "rna.run=true but the manifest contains no rows with a non-null rna_file. " +
                  "Either populate rna_file in the manifest, or set rna.run=false for ATAC-only runs."
    }

    // --- MIS-04 (2026-05-04): onramp wiring guard for the 6 forward-declared
    // keys with no consumer anywhere. After the D2.1 back-port from PBMC,
    // 14 keys are wired (3 top-level + 10 inside REGULATORY_ANALYSIS + 1 mudata).
    // The keys below remain placeholders and would silently no-op if set.
    def neverWiredOnrampKeys = ['cistopic_obj_pkl', 'seurat_rds', 'da_peaks_dir',
                                'cicero_connections_ctrl', 'cicero_connections_trt',
                                'cicero_ccan_ctrl', 'cicero_ccan_trt']
    if (params.onramp) {
        params.onramp.each { k, v ->
            if (v && k in neverWiredOnrampKeys) {
                errors << "params.onramp.${k} is set but is forward-declared only — no " +
                          "consumer exists in any FORGE main.nf. Unset this key; the " +
                          "producer step will run regardless."
            }
        }
    }

    // --- MIS-04 partial-set validators on the multi-key bundles.
    // Cicero global onramp: connections + ccan + cds are an all-or-none triple.
    def cicero_om_keys = [params.onramp?.cicero_connections,
                          params.onramp?.cicero_ccan,
                          params.onramp?.cicero_cds]
    def cicero_om_set  = cicero_om_keys.count { it }
    if (cicero_om_set in [1, 2]) {
        errors << "Cicero onramp is an all-or-none triple: cicero_connections + " +
                  "cicero_ccan + cicero_cds. Got ${cicero_om_set}/3 set."
    }
    // ChromVAR onramp: deviations + raw are a pair.
    def chromvar_om_keys = [params.onramp?.chromvar_deviations,
                            params.onramp?.chromvar_raw]
    def chromvar_om_set = chromvar_om_keys.count { it }
    if (chromvar_om_set == 1) {
        errors << "ChromVAR onramp is a pair: chromvar_deviations + chromvar_raw. " +
                  "Got 1/2 set."
    }
    // ATAC onramp side keys must be set when atac_peak_matrix_h5ad is used and
    // the consuming downstream is enabled.
    if (params.onramp?.atac_peak_matrix_h5ad) {
        if (params.scprinter?.run && !params.onramp?.atac_individual_samples_dir) {
            warnings << "atac_peak_matrix_h5ad onramp set without atac_individual_samples_dir " +
                       "— scPRINTer cannot build per-sample fragment binding."
        }
        if (params.enhancer_footprinting?.run && !params.onramp?.atac_anndataset) {
            warnings << "atac_peak_matrix_h5ad onramp set without atac_anndataset — " +
                       "ENHANCER_FOOTPRINTING_RECIPES will be unable to bind snapatac2 anndataset."
        }
    }
    // RNA per-sample dir is required for multiome when rna_integrated is onramped.
    if (params.onramp?.rna_integrated_h5ad && params.run_multiome_integration &&
        !params.onramp?.rna_per_sample_h5ads_dir && !params.onramp?.mudata_h5mu) {
        errors << "rna_integrated_h5ad onramp + run_multiome_integration=true requires " +
                  "either rna_per_sample_h5ads_dir (for fresh multiome) or mudata_h5mu " +
                  "(skip multiome). Set one, or set run_multiome_integration=false."
    }

    // --- MIS-22 (2026-05-04): disease_stratified=true requires resolved condition labels
    if (params.enhancer_footprinting?.disease_stratified) {
        def ctrl = params.enhancer_footprinting?.control_condition ?:
                   params.differential?.control_condition
        def trt  = params.enhancer_footprinting?.treatment_condition ?:
                   params.differential?.treatment_condition
        if (!ctrl || !trt) {
            errors << "enhancer_footprinting.disease_stratified=true but condition labels " +
                      "are unset (control='${ctrl}', treatment='${trt}'). Set " +
                      "enhancer_footprinting.{control,treatment}_condition or fall back to " +
                      "differential.{control,treatment}_condition. Otherwise, output filenames " +
                      "would resolve to literal 'null' strings."
        }
    }

    // --- Emit warnings ---
    warnings.each { log.warn it }

    // Fail hard on validation errors
    if (errors) {
        def msg = "\n" + "="*80 + "\n" +
                  "PRE-FLIGHT CHECKLIST FAILED (${errors.size()} error(s)):\n" +
                  "="*80 + "\n" +
                  errors.withIndex().collect { err, i -> "  ${i+1}. ${err}" }.join("\n") + "\n" +
                  "="*80
        error msg
    }

    // FIX-16: Enhanced startup banner with pre-flight summary
    log.info """
    PRE-FLIGHT CHECKLIST PASSED (${checks_passed.size()} checks):
      ${checks_passed.collect { "  [OK] ${it}" }.join('\n      ')}
    ${warnings ? "  Warnings: ${warnings.size()} (see above)" : "  No warnings."}
    """
}


// ============================================================================
// RNA WORKFLOW
// ============================================================================
workflow RNA {
    main:
    log.info """
    RNA PROCESSING WORKFLOW
    Cell type key: ${cell_type_key}
    """

    // ========================================================================
    // SOURCE-OF-TRUTH MANIFEST PARSING (RNA)
    // ========================================================================
    ch_all_files = Channel.fromPath(params.metadata_file)
        .splitCsv(header: true)
        .filter { isNonEmptyRow(it) }  // FIX-A6: skip phantom rows
        .map { trimRow(it) }           // FIX-A5: trim whitespace from fields
        .filter { isLane(it) }         // FIX-P0-7: case-insensitive sample_type
        .map { row ->
            def rna_dir = resolveRnaDir(row)
            def rna_fname = row.rna_file
            def rna_file = file("${rna_dir}/${rna_fname}")
            // FIX-A7: Validate paths in both production and preview mode
            if (!rna_file.exists()) {
                error "Manifest validation failed: RNA file not found -> ${rna_file}"
            }
            tuple(row.sample_id, rna_file)
        }

    // Step 1 - Run CellBender
    log.info "Step 1: Running CellBender ambient RNA removal..."
    CELLBENDER(ch_all_files)

    // Step 2 - Run QC on CellBender-filtered data
    // Demux files resolved per-sample from generic maps (demux_metadata, demux_souporcell_dirs).
    // For datasets without demultiplexing (e.g., 10x PBMC), both resolve to NO_FILE.
    log.info "Step 2: Running RNA QC and demultiplexing..."
    ch_qc_input = CELLBENDER.out.filtered_h5.map { sample, h5 ->
        tuple(sample, h5, resolveDemuxMetadata(sample), resolveDemuxSouporcell(sample))
    }
    RNA_QC(ch_qc_input)

    // Step 3: Flatten and collect all demultiplexed samples
    log.info "Step 3: Concatenating RNA batches..."
    ch_all_qc = RNA_QC.out.filtered_h5ad
        .flatMap { sample, files -> files }
        .collect()

    CONCAT_BATCHES(ch_all_qc)

    if (has_reference) {
        // ============================================================
        // PATH A: Reference atlas provided -> scANVI annotation
        // ============================================================
        log.info "Reference atlas detected -- using scANVI annotation path"

        log.info "Step 4: Preparing reference for integration..."
        PREPARE_REFERENCE(
            CONCAT_BATCHES.out.concatenated,
            params.species
        )

        log.info "Step 5: Training scVI foundation model..."
        TRAIN_SCVI(
            PREPARE_REFERENCE.out.prepared_ref
        )

        log.info "Step 6: Training scANVI model for cell-type annotation..."
        TRAIN_SCANVI(
            PREPARE_REFERENCE.out.prepared_ref,
            PREPARE_REFERENCE.out.prepared_query,
            PREPARE_REFERENCE.out.label_key,
            TRAIN_SCVI.out.scvi_model_dir
        )

        if (params.rna.annotation_method == 'markers') {
            if (!params.rna.marker_file) {
                error "params.rna.annotation_method='markers' requires params.rna.marker_file"
            }
            log.info "Step 6b: Running marker-gene annotation..."
            RUN_MARKER_ANNOTATION(
                TRAIN_SCANVI.out.annotated,
                Channel.value(file(params.rna.marker_file))
            )
            log.info "Step 7: Generating post-integration visualizations..."
            PLOT_POST_SCANVI(RUN_MARKER_ANNOTATION.out.annotated_h5ad, cell_type_key)
        } else {
            log.info "Step 6b: Running CellTypist annotation..."
            RUN_CELLTYPIST(TRAIN_SCANVI.out.annotated)

            log.info "Step 7: Generating post-integration visualizations..."
            PLOT_POST_SCANVI(RUN_CELLTYPIST.out.annotated_h5ad, cell_type_key)
        }

    } else {
        // ============================================================
        // PATH B: No reference atlas -> CellTypist (direct)
        // ============================================================
        log.info "No reference atlas provided -- using CellTypist (direct) path"

        if (params.rna.annotation_method == 'markers') {
            if (!params.rna.marker_file) {
                error "params.rna.annotation_method='markers' requires params.rna.marker_file"
            }
            log.info "Step 5-alt: Running marker-gene annotation..."
            RUN_MARKER_ANNOTATION(
                CONCAT_BATCHES.out.concatenated,
                Channel.value(file(params.rna.marker_file))
            )
            log.info "Step 7: Generating post-integration visualizations..."
            PLOT_POST_SCANVI(RUN_MARKER_ANNOTATION.out.annotated_h5ad, cell_type_key)
        } else {
            log.info "Step 5-alt: Running CellTypist annotation..."
            RUN_CELLTYPIST(CONCAT_BATCHES.out.concatenated)

            log.info "Step 7: Generating post-integration visualizations..."
            PLOT_POST_SCANVI(RUN_CELLTYPIST.out.annotated_h5ad, cell_type_key)
        }
    }

    // Debug: View what was emitted
    PLOT_POST_SCANVI.out.annotated_updated.view { "H5AD file emitted: $it" }

    // ====================================================================
    // Step 8: CellChat L-R Analysis (global-first + optional differential)
    // ====================================================================
    if (params.cellchat.run) {
        log.info "Step 8: Running CellChat ligand-receptor analysis..."

        // --- Tier 1: Global CellChat (condition-agnostic) ---
        log.info "  Step 8a: Global CellChat (all cells, no condition split)..."
        ch_cellchat_input = PLOT_POST_SCANVI.out.annotated_updated
            .map { file -> tuple("integrated", file) }

        RUN_CELLCHAT(
            ch_cellchat_input,
            cell_type_key,
            "none",
            params.species
        )

        // --- Tier 2: Per-condition CellChat + comparison ---
        def cellchat_has_conditions = (params.cellchat.condition_key &&
                                       params.cellchat.condition_key != "none" &&
                                       params.cellchat.conditions &&
                                       !params.cellchat.conditions.isEmpty())

        if (cellchat_has_conditions) {
            log.info "  Step 8b: Per-condition CellChat comparative analysis..."
            log.info "  Conditions: ${params.cellchat.conditions}"

            ch_conditions = Channel.from(params.cellchat.conditions)

            CELLCHAT_PER_CONDITION(
                PLOT_POST_SCANVI.out.annotated_updated,
                ch_conditions,
                cell_type_key,
                params.cellchat.condition_key,
                params.species,
                file(params.cellchat.group_mapping ?: 'NO_FILE')
            )

            CELLCHAT_COMPARE(
                CELLCHAT_PER_CONDITION.out.cellchat_rds.collect(),
                CELLCHAT_PER_CONDITION.out.condition_label.collect(),
                params.species
            )
        } else {
            log.info "  Skipping CellChat comparative analysis (no conditions defined)"
        }
    }

    // ====================================================================
    // Step 9: hdWGCNA Co-expression Network Analysis
    // ====================================================================
    if (params.hdwgcna.run) {
        log.info "Step 9: Running hdWGCNA co-expression network analysis..."

        // --- Step 9a: Convert h5ad to Seurat RDS ---
        CONVERT_H5AD_TO_SEURAT(
            PLOT_POST_SCANVI.out.annotated_updated,
            cell_type_key
        )

        // --- Step 9b: Create per-cell-type channel ---
        ch_hdwgcna_input = CONVERT_H5AD_TO_SEURAT.out.seurat_rds
            .map { seurat_file ->
                def cell_types_file = file(seurat_file.toString().replace('.rds', '_celltypes.txt'))
                def cell_types = cell_types_file.readLines()
                cell_types.collect { ct -> tuple(seurat_file, ct.trim()) }
            }
            .flatMap()

        // --- Tier 1: Global per-cell-type network construction ---
        log.info "  Step 9c: Tier 1 -- per-cell-type network construction..."
        HDWGCNA_PER_CELLTYPE(
            ch_hdwgcna_input,
            cell_type_key,
            params.hdwgcna.metadata ?: 'NO_FILE'
        )

        // --- Tier 1b: Enrichment + Network Visualization ---
        // out.seurat_obj emits *_seurat_with_wgcna.rds (a Seurat object).
        // out.results emits *_complete_results.rds (a list wrapper) — wrong type for downstream R scripts.
        log.info "  Step 9d: Tier 1 -- enrichment & network visualization..."
        HDWGCNA_ENRICHMENT(
            HDWGCNA_PER_CELLTYPE.out.seurat_obj.map { ct, rds -> rds },
            HDWGCNA_PER_CELLTYPE.out.seurat_obj.map { ct, rds -> ct },
            params.species
        )

        // --- Tier 2: Differential DME analysis ---
        def hdwgcna_has_conditions = (params.hdwgcna.condition_key &&
                                      params.hdwgcna.condition_key != "none" &&
                                      params.hdwgcna.control_condition &&
                                      params.hdwgcna.treatment_condition)

        if (hdwgcna_has_conditions) {
            log.info "  Step 9e: Tier 2 -- differential DME analysis..."
            log.info "  ${params.hdwgcna.control_condition} vs ${params.hdwgcna.treatment_condition}"

            def traits_str = params.hdwgcna.traits ? params.hdwgcna.traits.join(',') : ""

            HDWGCNA_DIFFERENTIAL(
                HDWGCNA_PER_CELLTYPE.out.seurat_obj.map { ct, rds -> rds },
                HDWGCNA_PER_CELLTYPE.out.seurat_obj.map { ct, rds -> ct },
                cell_type_key,
                params.hdwgcna.condition_key,
                params.hdwgcna.control_condition,
                params.hdwgcna.treatment_condition,
                traits_str,
                file(params.hdwgcna?.group_mapping ?: params.differential_rna?.group_mapping ?: params.cellchat?.group_mapping ?: 'NO_FILE')
            )
        } else {
            log.info "  Skipping hdWGCNA differential analysis (no conditions defined)"
        }
    }

    emit:
    integrated_rna = PLOT_POST_SCANVI.out.annotated_updated
    qc_h5ads = RNA_QC.out.filtered_h5ad
    pre_qc_plots = CONCAT_BATCHES.out.plots
    post_qc_plots = PLOT_POST_SCANVI.out.plots
    cellbender_reports = CELLBENDER.out.report

    // CellChat outputs
    cellchat_results    = params.cellchat.run ? RUN_CELLCHAT.out.cellchat_rds : Channel.empty()
    cellchat_plots      = params.cellchat.run ? RUN_CELLCHAT.out.plots : Channel.empty()
    cellchat_csv        = params.cellchat.run ? RUN_CELLCHAT.out.csv : Channel.empty()
    cellchat_comparison = (params.cellchat.run && params.cellchat.conditions && !params.cellchat.conditions.isEmpty()) ?
        CELLCHAT_COMPARE.out.comparison_rds : Channel.empty()
    cellchat_comp_plots = (params.cellchat.run && params.cellchat.conditions && !params.cellchat.conditions.isEmpty()) ?
        CELLCHAT_COMPARE.out.plots : Channel.empty()

    // hdWGCNA outputs
    hdwgcna_results     = params.hdwgcna.run ? HDWGCNA_PER_CELLTYPE.out.results : Channel.empty()
    hdwgcna_figures     = params.hdwgcna.run ? HDWGCNA_PER_CELLTYPE.out.figures : Channel.empty()
    hdwgcna_logs        = params.hdwgcna.run ? HDWGCNA_PER_CELLTYPE.out.log : Channel.empty()
    hdwgcna_enrichment  = params.hdwgcna.run ? HDWGCNA_ENRICHMENT.out.enrichr_plots : Channel.empty()
    hdwgcna_networks    = params.hdwgcna.run ? HDWGCNA_ENRICHMENT.out.network_plots : Channel.empty()
    hdwgcna_dme         = (params.hdwgcna.run && params.hdwgcna.condition_key && params.hdwgcna.condition_key != "none") ?
        HDWGCNA_DIFFERENTIAL.out.dme_results : Channel.empty()
    hdwgcna_dme_plots   = (params.hdwgcna.run && params.hdwgcna.condition_key && params.hdwgcna.condition_key != "none") ?
        HDWGCNA_DIFFERENTIAL.out.dme_plots : Channel.empty()
    hdwgcna_trait_plots = (params.hdwgcna.run && params.hdwgcna.condition_key && params.hdwgcna.condition_key != "none") ?
        HDWGCNA_DIFFERENTIAL.out.module_trait_plots : Channel.empty()
}

// ============================================================================
// RNA DIFFERENTIAL EXPRESSION WORKFLOW
// ============================================================================
workflow RNA_DIFFERENTIAL {

    take:
        annotated_h5ad

    main:
        log.info """
        RNA DIFFERENTIAL EXPRESSION (MAST)
        Comparisons: ${params.differential_rna.comparisons.size()}
        """

        // Step 1: Assign test groups
        log.info "Step 1: Assigning test groups to samples..."
        ASSIGN_TEST_GROUPS(
            annotated_h5ad,
            file(params.differential_rna.group_mapping)
        )

        // Step 2: Convert to Seurat RDS for MAST
        log.info "Step 2: Converting h5ad to Seurat format for MAST..."
        CONVERT_H5AD_FOR_MAST(
            ASSIGN_TEST_GROUPS.out.h5ad
        )

        // Step 2b: Extract valid cell types for MAST
        log.info "Step 2b: Extracting valid cell types for MAST..."

        EXTRACT_CELL_TYPES_FOR_MAST(
            ASSIGN_TEST_GROUPS.out.h5ad
        )

        def cell_types_file_ch = EXTRACT_CELL_TYPES_FOR_MAST.out.cell_types

        def cell_types_list_ch = cell_types_file_ch.map { f ->
            if( !f.exists() ) {
                throw new IllegalStateException("Cell types file not found: $f")
            }
            def lst = f.readLines()
                   .collect { it.trim() }
                   .findAll { it }
            log.info "Valid cell types for MAST (${lst.size()}): ${lst.join(', ')}"
            return lst
        }

        def mast_cell_types = cell_types_list_ch.flatten()

        // Step 3: Set up comparisons from config
        log.info "Step 3: Setting up comparisons from params.differential_rna.comparisons..."

        if (!params.differential_rna.comparisons || params.differential_rna.comparisons.size() == 0) {
            log.warn "No comparisons defined in params.differential_rna.comparisons — RNA DE will be skipped."
        }

        def comparisons_ch = Channel.fromList(params.differential_rna.comparisons)
            .map { it -> tuple(it[0], it[1]) }

        def mast_de_params = comparisons_ch
            .combine(mast_cell_types)
            .map { group1, group2, cell_type ->
                tuple(cell_type, group1, group2)
            }

        // Step 4: Run MAST differential expression
        log.info "Step 4: Running MAST differential expression..."

        RUN_MAST_DE(
            CONVERT_H5AD_FOR_MAST.out.seurat_rds,
            mast_de_params
        )

        // Step 5: Create volcano plots
        log.info "Step 5: Generating volcano plots..."
        CREATE_VOLCANO_PLOTS(
            RUN_MAST_DE.out.de_results
        )

        // Step 6: Run GO enrichment
        if( params.differential_rna.run_go_enrichment ) {
            log.info "Step 6: Running GO enrichment analysis..."
            RUN_GO_ENRICHMENT(
                RUN_MAST_DE.out.de_results
            )
        }

    emit:
        de_results    = RUN_MAST_DE.out.de_results
        volcano_plots = CREATE_VOLCANO_PLOTS.out.plots
        go_results    = params.differential_rna.run_go_enrichment ?
                        RUN_GO_ENRICHMENT.out.enrichment_results :
                        Channel.empty()
}

// ============================================================================
// ATAC WORKFLOW (UNIFIED TWO-STAGE QC)
// ============================================================================
workflow ATAC_INITIAL {
    main:
    log.info """
    ATAC INITIAL QC (Stage 1: Uniform Thresholds)
    """

    // ========================================================================
    // SOURCE-OF-TRUTH MANIFEST PARSING (ATAC)
    // ========================================================================
    ch_all_fragments = Channel.fromPath(params.metadata_file)
        .splitCsv(header: true)
        .filter { isNonEmptyRow(it) }  // FIX-A6: skip phantom rows
        .map { trimRow(it) }           // FIX-A5: trim whitespace from fields
        .map { row ->
            def atac_dir = resolveAtacDir(row)
            def frag_fname = row.fragment_file.contains('.') ? row.fragment_file : "${row.fragment_file}.bed.gz"
            def fragment_file = file("${atac_dir}/${frag_fname}")
            // FIX-A7: Validate paths in both production and preview mode
            if (!fragment_file.exists()) {
                error "Manifest validation failed: ATAC fragment file not found -> ${fragment_file}"
            }
            fragment_file
        }
        .collect()

    // Debug: count files
    ch_all_fragments.view { files -> "Total ATAC fragment files: ${files.size()}" }

    // Run initial QC with uniform thresholds
    ATAC_INITIAL_QC(
        ch_all_fragments,
        file(params.atac.sample_metadata)
    )

    // Build sample-specific thresholds from QC statistics
    ATAC_MAKE_THRESHOLDS(
        ATAC_INITIAL_QC.out.sample_stats
    )

    emit:
    anndataset         = ATAC_INITIAL_QC.out.anndataset
    qc_plots           = ATAC_INITIAL_QC.out.qc_plots
    summary            = ATAC_INITIAL_QC.out.summary
    sample_stats       = ATAC_INITIAL_QC.out.sample_stats
    thresholds         = ATAC_MAKE_THRESHOLDS.out.thresholds_file
    individual_samples = ATAC_INITIAL_QC.out.individual_samples
}

// ============================================================================
// ATAC FINAL WORKFLOW (Stage 2: Sample-Specific Thresholds + Annotation)
// ============================================================================
workflow ATAC_FINAL {
    take:
    thresholds_file

    main:
    log.info """
    ATAC FINAL QC (Stage 2: Sample-Specific Filtering)
    """

    // ========================================================================
    // SOURCE-OF-TRUTH MANIFEST PARSING (ATAC FINAL)
    // ========================================================================
    ch_demux_fragments = Channel.fromPath(params.metadata_file)
        .splitCsv(header: true)
        .filter { isNonEmptyRow(it) }  // FIX-A6: skip phantom rows
        .map { trimRow(it) }           // FIX-A5: trim whitespace from fields
        .map { row ->
            def atac_dir = resolveAtacDir(row)
            def frag_fname = row.fragment_file.contains('.') ? row.fragment_file : "${row.fragment_file}.bed.gz"
            def fragment_file = file("${atac_dir}/${frag_fname}")
            fragment_file
        }

    // Run final ATAC pipeline with sample-specific thresholds
    ATAC_FINAL_PIPELINE(
        ch_demux_fragments.collect(),
        file(params.atac.sample_metadata),
        thresholds_file
    )

    // Automated cell-type annotation
    if (params.atac.auto_annotate) {
        if (params.atac.marker_file) {
            // Super-user mode: marker-based annotation → 'cell_type' column
            log.info "ATAC annotation: marker-based (super-user mode, marker_file provided)"
            ATAC_CELLTYPE_ANNOTATION(
                ATAC_FINAL_PIPELINE.out.cluster_scores
            )
            MERGE_ANNOTATIONS(
                ATAC_FINAL_PIPELINE.out.peak_matrix,
                ATAC_CELLTYPE_ANNOTATION.out.annotations,
                file(params.atac.sample_metadata),
                'marker'
            )
            peak_matrix_annotated = MERGE_ANNOTATIONS.out.peak_matrix
        } else if (params.atac.annotation_method == 'scatanno') {
            // scATAnno mode: reference-based peak annotation → 'cell_type_prediction' column
            log.info "ATAC annotation: scATAnno reference-based (annotation_method = scatanno)"
            ATAC_SCATANNO(
                ATAC_FINAL_PIPELINE.out.anndataset,
                file(params.scatanno.reference_atlas)
            )
            MERGE_ANNOTATIONS(
                ATAC_FINAL_PIPELINE.out.peak_matrix,
                ATAC_SCATANNO.out.annotations,
                file(params.atac.sample_metadata),
                'scatanno'
            )
            peak_matrix_annotated = MERGE_ANNOTATIONS.out.peak_matrix
        } else {
            error "ATAC auto_annotate is enabled but no annotation method was selected. Set params.atac.marker_file for super-user marker mode, or params.atac.annotation_method = 'scatanno' for reference-based annotation. ATAC_CELLTYPIST on gene activity has been removed — use scATAnno instead."
        }
    } else {
        peak_matrix_annotated = ATAC_FINAL_PIPELINE.out.peak_matrix
    }

    emit:
    anndataset         = ATAC_FINAL_PIPELINE.out.anndataset
    peak_matrix        = peak_matrix_annotated
    peak_matrix_raw    = ATAC_FINAL_PIPELINE.out.peak_matrix
    gene_matrix        = ATAC_FINAL_PIPELINE.out.gene_matrix
    qc_plots           = ATAC_FINAL_PIPELINE.out.qc_plots
    summary            = ATAC_FINAL_PIPELINE.out.summary
    individual_samples = ATAC_FINAL_PIPELINE.out.individual_samples
    cluster_scores     = ATAC_FINAL_PIPELINE.out.cluster_scores
    atac_annotations   = (params.atac.auto_annotate && params.atac.marker_file) ? ATAC_CELLTYPE_ANNOTATION.out.annotations : Channel.empty()
}

// ============================================================================
// DIFFERENTIAL ATAC ANALYSIS
// ============================================================================
workflow ATAC_DIFFERENTIAL {
    take:
    peak_matrix
    metadata

    main:

    // FIX-E7: If cell_types list is empty, auto-discover from annotated peak matrix.
    // This reinforces wiring from celltypist/scanvi predictions — users don't need
    // to know cell types before running the pipeline.
    if (params.differential.cell_types && params.differential.cell_types.size() > 0) {
        ch_cell_types = Channel.from(params.differential.cell_types)
        log.info """
        DIFFERENTIAL ACCESSIBILITY ANALYSIS
        Comparisons: ${params.differential.comparisons.size()}
        Cell types: ${params.differential.cell_types.size()} (from config)
        Total tests: ${params.differential.comparisons.size() * params.differential.cell_types.size()}
        """
    } else {
        log.info """
        DIFFERENTIAL ACCESSIBILITY ANALYSIS
        Comparisons: ${params.differential.comparisons.size()}
        Cell types: auto-discovering from annotated peak matrix (params.differential.cell_types was empty)
        """
        EXTRACT_ATAC_CELL_TYPES(peak_matrix)
        ch_cell_types = EXTRACT_ATAC_CELL_TYPES.out.cell_types
            .splitText()
            .map { it.trim() }
            .filter { it }
    }

    ch_comparisons = Channel.fromList(params.differential.comparisons)
        .map { it -> tuple(it[0], it[1]) }

    ch_tasks = ch_comparisons
        .combine(ch_cell_types)

    // C4: snapatac_diff.nf now takes a single 3-tuple (cell_type, treatment, control).
    // ch_tasks emits (treatment, control, cell_type) from ch_comparisons.combine(ch_cell_types).
    SNAPATAC_DIFFERENTIAL(
        peak_matrix,
        metadata,
        ch_tasks.map { tuple(it[2], it[0], it[1]) }
    )

    emit:
    da_peaks = SNAPATAC_DIFFERENTIAL.out.da_peaks
    plots = SNAPATAC_DIFFERENTIAL.out.plots
}

// ============================================================================
// REGULATORY ANALYSIS WORKFLOW
// Uses PBMC's fragment_files input pattern for SCPRINTER_BUILD_PRINTER.
// ============================================================================
workflow REGULATORY_ANALYSIS {
    take:
    peak_matrix
    individual_samples
    da_peaks_optional
    metadata
    fragment_files    // Fragment files for SCPRINTER_BUILD_PRINTER (channel input)
    cell_type_col     // FIX-46: obs column name for cell types in peak matrix

    main:

    def has_da_peaks = (params.differential.run ?: false)

    def is_discovery_mode = !params.scprinter.target_genes || params.scprinter.target_genes.isEmpty()

    def use_chromvar_for_cicero = (params.cicero.use_chromvar_targets ?: false) && is_discovery_mode && params.chromvar.run

    log.info """
    REGULATORY ANALYSIS WORKFLOW (parallel architecture)
    - Cicero:   Co-accessibility networks ${params.cicero.run ? '(ENABLED)' : '(disabled)'}
    - ChromVAR: TF motif enrichment (global -- no conditions required)
    - scPRINT:  TF footprinting (${is_discovery_mode ? 'DISCOVERY' : 'TARGETED'} mode)
    - Mode:     ${is_discovery_mode ? 'Per-cell-type ChromVAR discovery -> fan-out' : 'User target genes -> single call'}
    - Cicero targets: ${use_chromvar_for_cicero ? 'ChromVAR-driven (data-driven)' : 'config list'}
    """

    // ================================================================
    // PARALLEL LEG A: Cicero Co-Accessibility (genome-wide)
    // 2026-04-30 onramp (D2.1 back-port 2026-05-04): bind ch_cicero_{connections,
    // ccan,cds} from either params.onramp.cicero_* triple OR run the per-chrom
    // fan-out. Downstream consumers (CICERO_TARGET_PLOTS, scPRINTer FP,
    // MAP_TF_TO_TARGET_GENES, emit) all use the channel variables.
    // ================================================================
    def cicero_om = (params.onramp?.cicero_connections &&
                     params.onramp?.cicero_ccan &&
                     params.onramp?.cicero_cds)
    def ch_cicero_connections, ch_cicero_ccan, ch_cicero_cds

    if (params.cicero.run) {
        if (cicero_om) {
            log.info "ON-RAMP: Cicero (connections + CCAN + CDS)"
            ch_cicero_connections = Channel.value(file(params.onramp.cicero_connections))
            ch_cicero_ccan        = Channel.value(file(params.onramp.cicero_ccan))
            ch_cicero_cds         = Channel.value(file(params.onramp.cicero_cds))
        } else {
            log.info "Running Cicero co-accessibility network analysis (per-chrom fan-out)..."

            CICERO_TRIPLETS(peak_matrix)

            // Phase 3 rewrite: global distance_parameter -> per-chrom run_cicero -> rbind.
            // Validated rho=1.0 vs single-process baseline in tests/cicero_parallel_test/.
            CICERO_ESTIMATE_DP(CICERO_TRIPLETS.out.triplets)

            // Mouse: chr1..19 + X/Y/M; human: chr1..22 + X/Y/M. chrUn_* excluded
            // upstream in the R scripts.
            def cicero_chroms = (params.species == 'mouse' ? (1..19) : (1..22))
                                    .collect { "chr${it}" } + ['chrX', 'chrY', 'chrM']
            chroms_ch = Channel.fromList(cicero_chroms)

            per_chrom_in = chroms_ch
                .combine(CICERO_ESTIMATE_DP.out.cicero_cds)
                .combine(CICERO_ESTIMATE_DP.out.gene_ann)
                .combine(CICERO_ESTIMATE_DP.out.dp)
            CICERO_FULL_CHROM(per_chrom_in)

            CICERO_JOIN(
                CICERO_FULL_CHROM.out.chrom_conns.map { it[1] }.collect(),
                CICERO_ESTIMATE_DP.out.ordered_cds,
                params.cicero.gtf_full,
                ""
            )

            ch_cicero_connections = CICERO_JOIN.out.connections
            ch_cicero_ccan        = CICERO_JOIN.out.ccan
            ch_cicero_cds         = CICERO_JOIN.out.cds
        }

        if (!use_chromvar_for_cicero && params.cicero.target_genes && !params.cicero.target_genes.isEmpty()) {
            log.info "Rendering Cicero target plots with static gene list: ${params.cicero.target_genes}"
            CICERO_TARGET_PLOTS(
                ch_cicero_connections,
                ch_cicero_ccan,
                ch_cicero_cds,
                params.cicero.gtf_plot,
                params.cicero.target_genes
            )
        }
    } else {
        ch_cicero_connections = Channel.empty()
        ch_cicero_ccan        = Channel.empty()
        ch_cicero_cds         = Channel.empty()
    }

    // ================================================================
    // PARALLEL LEG B: GPU ChromVAR TF Motif Enrichment
    // 2026-04-30 onramp (D2.1 back-port 2026-05-04): bind ch_chromvar_dev /
    // ch_chromvar_raw from either params.onramp.{chromvar_deviations,
    // chromvar_raw} OR GPU_CHROMVAR output.
    // ================================================================
    def chromvar_om = (params.onramp?.chromvar_deviations &&
                       params.onramp?.chromvar_raw)
    def ch_chromvar_dev, ch_chromvar_raw

    if (params.chromvar.run) {
        if (chromvar_om) {
            log.info "ON-RAMP: ChromVAR (deviations + raw)"
            ch_chromvar_dev = Channel.value(file(params.onramp.chromvar_deviations))
            ch_chromvar_raw = Channel.value(file(params.onramp.chromvar_raw))
        } else {
            log.info "Running GPU-accelerated ChromVAR analysis..."

            def da_peaks_collected = has_da_peaks ?
                da_peaks_optional.collect() :
                Channel.value(file('NO_FILE_DA_PEAKS'))

            GPU_CHROMVAR(
                peak_matrix,
                Channel.value(file('NO_FILE_METADATA')),
                params.scprinter.cache_dir,
                params.scprinter.pfms ?: '',
                params.scprinter.genome,
                da_peaks_collected
            )
            ch_chromvar_dev = GPU_CHROMVAR.out.chromvar_dev
            ch_chromvar_raw = GPU_CHROMVAR.out.chromvar_raw
        }

        if (has_da_peaks) {
            VIS_CHROMVAR(ch_chromvar_dev)
        } else {
            log.info "Skipping VIS_CHROMVAR (requires differential conditions for permutation tests)"
        }

        if (is_discovery_mode && params.scprinter.run) {
            EXTRACT_CHROMVAR_MOTIFS(
                ch_chromvar_dev,
                params.chromvar.top_n_per_celltype,
                params.chromvar.min_motif_zscore,
                atac_cell_type_key,
                params.chromvar?.global_top_n ?: 0
            )

            EXTRACT_CHROMVAR_MOTIFS.out.report.view {
                "\n==== PER-CELL-TYPE CHROMVAR MOTIFS ====\n${it.text}\n======================================="
            }
        }
    } else {
        ch_chromvar_dev = Channel.empty()
        ch_chromvar_raw = Channel.empty()
    }

    // ================================================================
    // ChromVAR-driven Cicero Target Plots
    // ================================================================
    if (use_chromvar_for_cicero && params.cicero.run && params.scprinter.run) {
        log.info "Rendering Cicero target plots with ChromVAR-discovered TFs..."

        ch_chromvar_target_genes = EXTRACT_CHROMVAR_MOTIFS.out.motif_list
            .map { json_file ->
                new groovy.json.JsonSlurper().parseText(json_file.text).all_unique_tfs
            }

        CICERO_TARGET_PLOTS(
            ch_cicero_connections,
            ch_cicero_ccan,
            ch_cicero_cds,
            params.cicero.gtf_plot,
            ch_chromvar_target_genes
        )
    }

    // ================================================================
    // STEP 3: scPRINTer TF Footprinting
    // 2026-04-30 onramp (D2.1 back-port 2026-05-04): bind ch_printer from
    // params.onramp.printer_h5ad OR build it via SCPRINTER_BARCODES + BUILD.
    // ================================================================
    def printer_om = params.onramp?.printer_h5ad
    def ch_printer

    if (params.scprinter.run) {
        if (printer_om) {
            log.info "ON-RAMP: scPRINTer printer h5ad: ${printer_om}"
            ch_printer = Channel.value(file(printer_om))
        } else {
            log.info "Running scPRINTer TF footprinting workflow (build printer)..."

            // ---- Shared setup: barcodes + printer (both modes) ----
            ch_samples_with_names = individual_samples.flatten()
                .map { h5ad_file ->
                    def sample_name = h5ad_file.baseName.replaceAll('.h5ad$', '')
                    tuple(h5ad_file, sample_name)
            }

            SCPRINTER_BARCODES(
                ch_samples_with_names.map { it[0] }.collect(),
                ch_samples_with_names.map { it[1] }.collect()
            )

            // PBMC-style: fragments come via channel input
            SCPRINTER_BUILD_PRINTER(
                SCPRINTER_BARCODES.out.barcodes,
                fragment_files
            )

            ch_printer = SCPRINTER_BUILD_PRINTER.out.printer
        }

        // Manual coordinate overrides (both modes)
        def manual_coords = params.scprinter.gene_coordinates ?
            file(params.scprinter.gene_coordinates) :
            file('NO_FILE')

        // ============================================================
        // DISCOVERY MODE: per-cell-type ChromVAR TFs -> fan-out
        // ============================================================
        if (is_discovery_mode && params.chromvar.run) {
            log.info "DISCOVERY MODE -- mapping ChromVAR TFs to target genes via motif-peak-CCAN linkage"

            // FIX-R3-3: Use species-appropriate GTF
            def tf_map_gtf = params.species == 'human' ? params.scprinter.gtf_human : params.scprinter.gtf_mouse
            MAP_TF_TO_TARGET_GENES(
                ch_chromvar_raw,
                EXTRACT_CHROMVAR_MOTIFS.out.motif_list,
                ch_cicero_ccan,
                tf_map_gtf
            )

            MAP_TF_TO_TARGET_GENES.out.report.view {
                "\n==== TF -> TARGET GENE MAPPING ====\n${it.text}\n===================================="
            }

            ch_per_celltype = MAP_TF_TO_TARGET_GENES.out.tf_targets
                .flatMap { json_file ->
                    def data = new groovy.json.JsonSlurper().parseText(json_file.text)
                    data.per_celltype_targets.collect { ct_name, target_list ->
                        tuple(ct_name, target_list)
                    }
                }

            ch_all_targets = MAP_TF_TO_TARGET_GENES.out.tf_targets
                .map { json_file ->
                    new groovy.json.JsonSlurper().parseText(json_file.text).all_target_genes
                }

            RESOLVE_GENE_COORDINATES(
                params.species,
                ch_all_targets,
                manual_coords
            )

            RESOLVE_GENE_COORDINATES.out.report.view {
                "\n==== GENE COORDINATE RESOLUTION ====\n${it.text}\n===================================="
            }

            ch_fp = ch_per_celltype
                .combine(RESOLVE_GENE_COORDINATES.out.coordinates)

            log.info "Fan-out footprinting at TF target gene promoters..."

            def cicero_conns_for_fp = params.cicero.run ?
                ch_cicero_connections.ifEmpty(file('NO_FILE')) :
                Channel.value(file('NO_FILE'))
            def pfm_for_fp = params.scprinter.pfms ?
                Channel.value(file(params.scprinter.pfms)) :
                Channel.value(file('NO_FILE'))

            SCPRINTER_FOOTPRINTING(
                peak_matrix,
                metadata,
                ch_printer,
                ch_fp.map { it[0] },
                ch_fp.map { it[1] },
                ch_fp.map { it[2] },
                '',
                '',
                MAP_TF_TO_TARGET_GENES.out.tf_targets,
                cicero_conns_for_fp,
                pfm_for_fp,
                cell_type_col
            )

            if (has_da_peaks) {
                log.info "DISCOVERY differential footprinting per cell type..."

                // Skip DIFF for cell types with no target genes — the script
                // exits 0 without producing footprints_*.h5ad (required output),
                // and the 780 MB printer copy in .command.sh is wasted. Non-DIFF
                // SCPRINTER_FOOTPRINTING above is intentionally left unfiltered.
                ch_fp_diff = ch_fp.filter { t -> t[1] && t[1].size() > 0 }

                SCPRINTER_FOOTPRINTING_DIFF(
                    peak_matrix,
                    metadata,
                    ch_printer,
                    ch_fp_diff.map { it[0] },
                    ch_fp_diff.map { it[1] },
                    ch_fp_diff.map { it[2] },
                    params.differential.control_condition,
                    params.differential.treatment_condition,
                    MAP_TF_TO_TARGET_GENES.out.tf_targets,
                    cicero_conns_for_fp,
                    pfm_for_fp,
                    cell_type_col
                )

                SCPRINTER_MOTIF_SCAN(
                    peak_matrix,
                    da_peaks_optional.collect(),
                    ch_printer,
                    ch_fp_diff.map { it[0] },
                    SCPRINTER_FOOTPRINTING_DIFF.out.footprints.collect()
                )
            }

        // ============================================================
        // TARGETED MODE: user genes -> single call + optional differential
        // ============================================================
        } else if (!is_discovery_mode) {
            log.info "TARGETED MODE -- user-specified genes: ${params.scprinter.target_genes}"

            RESOLVE_GENE_COORDINATES(
                params.species,
                Channel.value(params.scprinter.target_genes),
                manual_coords
            )

            RESOLVE_GENE_COORDINATES.out.report.view {
                "\n==== GENE COORDINATE RESOLUTION ====\n${it.text}\n===================================="
            }

            // Per-ct fan-out mirrors SCPRINTER_FOOTPRINTING_DIFF below: a literal
            // 'targeted' as cell_type fails downstream because the python filters
            // adata.obs[cell_type_col] == 'targeted' and gets 0 cells.
            ch_ct_global = Channel.from(params.differential.cell_types)

            SCPRINTER_FOOTPRINTING(
                peak_matrix,
                metadata,
                ch_printer,
                ch_ct_global,
                Channel.value(params.scprinter.target_genes),
                RESOLVE_GENE_COORDINATES.out.coordinates,
                '',
                '',
                file('NO_FILE'),
                file('NO_FILE'),
                params.scprinter.pfms ? file(params.scprinter.pfms) : file('NO_FILE'),
                cell_type_col
            )

            if (has_da_peaks && params.differential.cell_types && !params.differential.cell_types.isEmpty()) {
                log.info "TARGETED differential footprinting per cell type..."

                ch_ct_diff = Channel.from(params.differential.cell_types)

                SCPRINTER_FOOTPRINTING_DIFF(
                    peak_matrix,
                    metadata,
                    ch_printer,
                    ch_ct_diff,
                    Channel.value(params.scprinter.target_genes),
                    RESOLVE_GENE_COORDINATES.out.coordinates,
                    params.differential.control_condition,
                    params.differential.treatment_condition,
                    file('NO_FILE'),
                    file('NO_FILE'),
                    params.scprinter.pfms ? file(params.scprinter.pfms) : file('NO_FILE'),
                    cell_type_col
                )

                SCPRINTER_MOTIF_SCAN(
                    peak_matrix,
                    da_peaks_optional.collect(),
                    ch_printer,
                    ch_ct_diff,
                    SCPRINTER_FOOTPRINTING_DIFF.out.footprints.collect()
                )
            } else if (has_da_peaks) {
                log.warn """
                WARNING: differential.run = true but differential.cell_types is empty.
                Cannot fan-out differential footprinting without cell types.
                Set differential.cell_types = ['Astrocytes', 'Microglia', ...] to enable.
                """
            }

        // ============================================================
        // EDGE CASE: discovery mode but chromvar.run = false
        // ============================================================
        } else {
            error """
            ERROR: Discovery mode (empty target_genes) requires chromvar.run = true.
            Either:
              1. Set chromvar.run = true for discovery mode, OR
              2. Provide scprinter.target_genes = ['GENE1', 'GENE2', ...] for targeted mode.
            """
        }
    }

    // ================================================================
    // SHI-STYLE TF ACCESSIBILITY TESTING (differential or descriptive)
    //   Absorbed from SDas. Gated by params.differential_tf.run.
    //   mode='differential' — per (cell_type, trt, ctrl); needs 2+ conditions
    //   mode='descriptive'  — per cell_type vs rest-of-cells (single-cond OK)
    // ================================================================
    if ((params.differential_tf?.run ?: false) && params.chromvar.run) {
        def tf_mode = params.differential_tf.mode ?: 'descriptive'
        def tf_cts  = params.differential_tf?.cell_types ?: []
        if (!tf_cts) {
            error "params.differential_tf.cell_types must be non-empty when differential_tf.run=true"
        }
        def tf_tasks = []
        if (tf_mode == 'differential') {
            def cmps = params.differential_tf?.comparisons ?: []
            if (!cmps) {
                error "differential_tf.mode='differential' requires differential_tf.comparisons=[[trt,ctrl], ...]"
            }
            tf_cts.each { ct -> cmps.each { c -> tf_tasks << tuple(ct as String, c[0] as String, c[1] as String) } }
        } else {
            // descriptive — 3rd tuple value is an unused placeholder (module ignores it)
            def filters = params.differential_tf?.conditions ?: [null]
            tf_cts.each { ct -> filters.each { f -> tf_tasks << tuple(ct as String, (f ?: 'all') as String, 'rest' as String) } }
        }
        log.info "TF accessibility (mode=${tf_mode}): ${tf_tasks.size()} tasks"
        def ch_tf_tasks = Channel.from(tf_tasks)
        DIFFERENTIAL_TF_ACCESSIBILITY(
            ch_chromvar_dev,
            peak_matrix,
            ch_tf_tasks
        )
    }

    // ================================================================
    // SHI-STYLE STRATIFIED CICERO + CO-ACCESSIBILITY COMPARISON
    //   Absorbed from SDas. Runs Cicero SEPARATELY per condition and
    //   compares the two co-accessibility maps. Runs IN ADDITION TO
    //   the default global Cicero (both useful).
    //
    //   D3: gate now auto-activates whenever differential.run=true AND a
    //   condition_key is declared (the disease axis is intrinsically
    //   tied to differential mode). params.cicero.stratified=true remains
    //   a manual override for non-differential workflows that still want
    //   per-condition Cicero.
    //
    //   TODO (separate effort): per-cell-type Cicero fan-out as the
    //   default architecture (always-on; orthogonal to disease axis).
    //   That's a new architecture not present in SDas — needs design
    //   doc + emit-shape rewrite + downstream consumer updates in
    //   ENHANCER_FOOTPRINTING_RECIPES. Not part of the SDas absorption
    //   series.
    // ================================================================
    def cicero_strat_auto   = (params.differential?.run ?: false) &&
                              params.differential?.condition_key
    def cicero_strat_manual = params.cicero?.stratified ?: false
    // 2026-05-06: workflow-scope channels for cis-rewiring; populated inside
    // the stratified-cicero block, default empty so emit gates are well-formed
    // when stratified Cicero doesn't run.
    def ch_strat_ctrl_links       = Channel.empty()
    def ch_strat_trt_links        = Channel.empty()
    def ch_strat_ctrl_peaks       = Channel.empty()
    def ch_strat_trt_peaks        = Channel.empty()
    def ch_strat_ctrl_connections = Channel.empty()
    def ch_strat_trt_connections  = Channel.empty()
    if ((cicero_strat_auto || cicero_strat_manual) && params.cicero.run) {
        def strat_ctrl = params.cicero?.control_condition ?:
                         params.differential?.control_condition
        def strat_trt  = params.cicero?.treatment_condition ?:
                         params.differential?.treatment_condition
        def strat_key  = params.cicero?.condition_key ?:
                         (params.differential?.condition_key ?: 'condition')
        if (!strat_ctrl || !strat_trt) {
            error "Stratified Cicero requires control_condition and treatment_condition (set params.cicero.{control,treatment}_condition or params.differential.{control,treatment}_condition)"
        }
        log.info "Stratified Cicero: ${strat_ctrl} vs ${strat_trt} (condition_key=${strat_key}) — per-chromosome fan-out per leg"

        def cicero_chroms_strat = (params.species == 'mouse' ? (1..19) : (1..22))
                                      .collect { "chr${it}" } + ['chrX', 'chrY', 'chrM']

        // --- Control leg: triplets → estimate_dp → per-chrom → join ---
        CICERO_TRIPLETS_STRATIFIED(peak_matrix, strat_key, strat_ctrl)
        CICERO_ESTIMATE_DP_CTRL(CICERO_TRIPLETS_STRATIFIED.out.triplets)
        CICERO_FULL_CHROM_CTRL(
            Channel.fromList(cicero_chroms_strat)
                .combine(CICERO_ESTIMATE_DP_CTRL.out.cicero_cds)
                .combine(CICERO_ESTIMATE_DP_CTRL.out.gene_ann)
                .combine(CICERO_ESTIMATE_DP_CTRL.out.dp)
        )
        CICERO_JOIN_CTRL(
            CICERO_FULL_CHROM_CTRL.out.chrom_conns.map { it[1] }.collect(),
            CICERO_ESTIMATE_DP_CTRL.out.ordered_cds,
            params.cicero.gtf_full,
            "stratified/${strat_ctrl}"
        )

        // --- Treatment leg: symmetric ---
        CICERO_TRIPLETS_STRATIFIED_TRT(peak_matrix, strat_key, strat_trt)
        CICERO_ESTIMATE_DP_TRT(CICERO_TRIPLETS_STRATIFIED_TRT.out.triplets)
        CICERO_FULL_CHROM_TRT(
            Channel.fromList(cicero_chroms_strat)
                .combine(CICERO_ESTIMATE_DP_TRT.out.cicero_cds)
                .combine(CICERO_ESTIMATE_DP_TRT.out.gene_ann)
                .combine(CICERO_ESTIMATE_DP_TRT.out.dp)
        )
        CICERO_JOIN_TRT(
            CICERO_FULL_CHROM_TRT.out.chrom_conns.map { it[1] }.collect(),
            CICERO_ESTIMATE_DP_TRT.out.ordered_cds,
            params.cicero.gtf_full,
            "stratified/${strat_trt}"
        )

        COMPARE_COACCESSIBILITY(
            CICERO_JOIN_CTRL.out.connections,
            CICERO_JOIN_TRT.out.connections,
            strat_ctrl,
            strat_trt
        )

        // 2026-05-06: per-condition CCAN→gene-link extraction. Powers the
        // cis-rewiring directionality proxy (gained enhancers per gene) and
        // fills the SHI-expected gene_links_{ctrl,trt} paths. enhancer_gtf
        // resolved later in REGULATORY_ANALYSIS — recompute here too.
        def strat_enhancer_gtf = params.species == 'human' ? params.scprinter.gtf_human : params.scprinter.gtf_mouse
        EXTRACT_CCAN_ENHANCERS_CTRL(
            CICERO_JOIN_CTRL.out.connections,
            CICERO_JOIN_CTRL.out.ccan,
            strat_enhancer_gtf,
            strat_ctrl
        )
        EXTRACT_CCAN_ENHANCERS_TRT(
            CICERO_JOIN_TRT.out.connections,
            CICERO_JOIN_TRT.out.ccan,
            strat_enhancer_gtf,
            strat_trt
        )
        ch_strat_ctrl_links       = EXTRACT_CCAN_ENHANCERS_CTRL.out.gene_links
        ch_strat_trt_links        = EXTRACT_CCAN_ENHANCERS_TRT.out.gene_links
        ch_strat_ctrl_peaks       = EXTRACT_CCAN_ENHANCERS_CTRL.out.enhancer_peaks
        ch_strat_trt_peaks        = EXTRACT_CCAN_ENHANCERS_TRT.out.enhancer_peaks
        ch_strat_ctrl_connections = CICERO_JOIN_CTRL.out.connections
        ch_strat_trt_connections  = CICERO_JOIN_TRT.out.connections
    }

    // ================================================================
    // EMIT
    // ================================================================
    emit:
    cicero_connections   = params.cicero.run ? ch_cicero_connections : Channel.empty()
    cicero_ccan          = params.cicero.run ? ch_cicero_ccan : Channel.empty()
    chromvar_deviations  = params.chromvar.run ? ch_chromvar_dev : Channel.empty()
    chromvar_per_ct      = (is_discovery_mode && params.chromvar.run && params.scprinter.run) ?
        EXTRACT_CHROMVAR_MOTIFS.out.motif_list : Channel.empty()
    scprinter_printer    = params.scprinter.run ?
        ch_printer : Channel.empty()
    scprinter_footprints = params.scprinter.run ?
        SCPRINTER_FOOTPRINTING.out.footprints : Channel.empty()
    // 2026-05-06: emit gate must match invocation gate. SCPRINTER_FOOTPRINTING_DIFF
    // fires only when (DISCOVERY mode + has_da_peaks) OR (TARGETED mode + has_da_peaks
    // + differential.cell_types non-empty). Prior gate accessed .out unconditionally
    // when (has_da_peaks && scprinter.run), throwing access-without-invocation in
    // TARGETED mode with empty cell_types.
    scprinter_diff       = (has_da_peaks && params.scprinter.run && (
            is_discovery_mode ||
            (params.differential?.cell_types && !params.differential.cell_types.isEmpty())
        )) ? SCPRINTER_FOOTPRINTING_DIFF.out.footprints : Channel.empty()
    // tf_target_genes.json — input for BUILD_TF_GENE_NETWORK (Phase 3, ATAC-only GRN)
    tf_targets           = (is_discovery_mode && params.chromvar.run) ?
        MAP_TF_TO_TARGET_GENES.out.tf_targets : Channel.empty()
    // 2026-05-04: tf_differential CSVs from DIFFERENTIAL_TF_ACCESSIBILITY,
    // piped through to ENHANCER_FOOTPRINTING_RECIPES (DSL2 cross-workflow
    // scope fix). Empty when differential_tf.run=false; consumers must guard
    // their compute on .ifEmpty([]) or equivalent.
    tf_diff              = ((params.differential_tf?.run ?: false) && params.chromvar.run) ?
        DIFFERENTIAL_TF_ACCESSIBILITY.out.tf_diff : Channel.empty()
    // 2026-05-06: gene_coordinates.json from RESOLVE_GENE_COORDINATES — needed
    // by ENHANCER_FOOTPRINTING_RECIPES.PROMOTER_MSFP_PER_CT (DSL2 cross-workflow
    // scope fix per feedback_dsl2_xworkflow_scope). Always populated when
    // scprinter.run=true since RESOLVE_GENE_COORDINATES fires in both discovery
    // and targeted modes.
    gene_coordinates     = params.scprinter.run ?
        RESOLVE_GENE_COORDINATES.out.coordinates : Channel.empty()
    // 2026-05-06: per-condition CCAN→gene-link + enhancer-peak channels for
    // cis-rewiring. Empty when stratified Cicero didn't run (single-condition
    // datasets, params.differential.run=false, etc.).
    cicero_strat_ctrl_links       = ch_strat_ctrl_links
    cicero_strat_trt_links        = ch_strat_trt_links
    cicero_strat_ctrl_peaks       = ch_strat_ctrl_peaks
    cicero_strat_trt_peaks        = ch_strat_trt_peaks
    cicero_strat_ctrl_connections = ch_strat_ctrl_connections
    cicero_strat_trt_connections  = ch_strat_trt_connections
}

// ============================================================================
// MULTIOME INTEGRATION WORKFLOW
// ============================================================================
workflow MULTIOME_INTEGRATION {
    take:
    rna_h5ad_files
    atac_h5ad_files
    scanvi_predictions
    metadata_csv
    sample_map
    atac_peak_matrix_annotated

    main:
    log.info """
    MULTIOME INTEGRATION WORKFLOW
    """

    // Step 1: Build unified MuData object
    log.info "Building MuData with sample-matched RNA+ATAC pairing..."

    BUILD_MUDATA(
        rna_h5ad_files.flatten().collect(),
        atac_h5ad_files.flatten().collect(),
        scanvi_predictions,
        metadata_csv,
        sample_map,
        atac_peak_matrix_annotated
    )

    // Export RNA modality from MuData for DORC / downstream use
    EXPORT_MUDATA_RNA(
        BUILD_MUDATA.out.mudata
    )

    // Step 2: MOFA integration (mode-dependent)
    if (params.mofa.run) {

        if (params.mofa.mode == 'high_memory') {
            log.info """
            MOFA+ HIGH-MEMORY MODE (Single-Shot Integration)
            - Using ALL cells (no subsampling)
            - Factors: ${params.mofa.n_factors}
            """.stripIndent()

            MOFA_INTEGRATE(
                BUILD_MUDATA.out.mudata
            )

            MOFA_VISUALIZE(
                MOFA_INTEGRATE.out.model,
                BUILD_MUDATA.out.mudata,
                MOFA_INTEGRATE.out.metadata
            )

            integrated_output = MOFA_VISUALIZE.out.integrated_mudata
            mofa_model = MOFA_INTEGRATE.out.model
            mofa_plots = MOFA_VISUALIZE.out.plots
            bootstrap_outputs = Channel.empty()

        } else if (params.mofa.mode == 'bootstrap') {
            log.info """
            BOOTSTRAP MOFA INTEGRATION (Low-Memory Mode)
            - Iterations: ${params.mofa.bootstrap.n_iterations}
            - Sample fraction: ${params.mofa.bootstrap.sample_fraction}
            - Memory will be logged to memory_log.jsonl
            """.stripIndent()

            BOOTSTRAP_MOFA_INTEGRATION(
                BUILD_MUDATA.out.mudata,
                params.mofa.bootstrap.n_iterations,
                params.mofa.bootstrap.sample_fraction,
                params.mofa.n_factors
            )

            ANALYZE_MEMORY_LOG(
                BOOTSTRAP_MOFA_INTEGRATION.out.memory_log
            )

            CONSENSUS_ANALYSIS(
                BOOTSTRAP_MOFA_INTEGRATION.out.results_dir
            )

            integrated_output = BUILD_MUDATA.out.mudata
            mofa_model = BOOTSTRAP_MOFA_INTEGRATION.out.models
            mofa_plots = CONSENSUS_ANALYSIS.out.stability_plot
            bootstrap_outputs = BOOTSTRAP_MOFA_INTEGRATION.out.results_dir

        } else {
            error "Invalid MOFA mode: ${params.mofa.mode}. Use 'high_memory' or 'bootstrap'"
        }

    } else {
        integrated_output = BUILD_MUDATA.out.mudata
        mofa_model = Channel.empty()
        mofa_plots = Channel.empty()
        bootstrap_outputs = Channel.empty()
    }

    // Step 3: MultiVI integration (always run when enabled)
    if (params.multivi.run) {
        log.info "Starting MultiVI integration on MuData (rna+atac)"

        MULTIVI_INTEGRATE(
            BUILD_MUDATA.out.mudata
        )

        log.info "Generating MultiVI visualizations and metrics..."
        MULTIVI_VISUALIZE(
            MULTIVI_INTEGRATE.out.integrated,
            MULTIVI_INTEGRATE.out.model
        )

        // Step 3b: MultiVI masking sweep (parallel with driver factors)
        // Fan out (fraction, seed) pairs into independent tasks so each MultiVI
        // fit runs in a fresh Python process — the scvi-tools host-RAM leak
        // cannot accumulate across fits this way. AGGREGATE joins the per-fit
        // outputs into the canonical results CSV + summary JSON + plots.
        if (params.multivi.masking_sweep.run) {
            log.info "Running MultiVI masking sweep benchmark (fan-out)..."
            masking_fractions = params.multivi.masking_sweep.fractions.toString().tokenize(',').collect { it.trim() }
            masking_seeds     = params.multivi.masking_sweep.seeds.toString().tokenize(',').collect { it.trim() }
            masking_sweep_ch  = Channel.fromList(masking_fractions)
                                       .combine(Channel.fromList(masking_seeds))
                                       .combine(BUILD_MUDATA.out.mudata)
            MULTIVI_MASKING_SWEEP_ONE(masking_sweep_ch)
            masking_fit_outputs = MULTIVI_MASKING_SWEEP_ONE.out.rna
                .mix(MULTIVI_MASKING_SWEEP_ONE.out.atac)
                .mix(MULTIVI_MASKING_SWEEP_ONE.out.integration)
                .collect()
            MULTIVI_MASKING_SWEEP_AGGREGATE(masking_fit_outputs)
        }

        // Step 3c: MultiVI driver factor analysis
        if (params.multivi.driver_factors.run) {
            log.info "Running MultiVI driver factor analysis..."
            // Pass MOFA-integrated MuData for comparison if available
            def mofa_input = (params.mofa.run && params.multivi.driver_factors.compare_mofa) ?
                MOFA_VISUALIZE.out.integrated_mudata : Channel.fromPath("NO_MOFA")
            MULTIVI_DRIVER_FACTORS(
                MULTIVI_INTEGRATE.out.integrated,
                MULTIVI_INTEGRATE.out.model,
                mofa_input
            )
        }

        // Step 3d: MultiVI gap-filling
        // Uses the pre-intersection RNA (concatenated) and ATAC (peak matrix)
        // to identify cells lost during BUILD_MUDATA inner join
        if (params.multivi.gap_fill.run) {
            log.info "Running MultiVI gap-filling for QC-filtered cells..."
            MULTIVI_GAP_FILL(
                BUILD_MUDATA.out.mudata,
                Channel.fromPath(params.multivi.gap_fill.rna_h5ad),
                Channel.fromPath(params.multivi.gap_fill.atac_h5ad)
            )

            // Step 3e: Biological validation of gap-filled cells
            log.info "Validating gap-filled cells..."
            MULTIVI_VALIDATE(
                MULTIVI_GAP_FILL.out.gap_filled,
                Channel.fromPath(params.multivi.gap_fill.rna_h5ad),
                Channel.fromPath(params.multivi.gap_fill.atac_h5ad)
            )
        }

    } else {
        log.info "Skipping MultiVI integration (disabled in config)"
    }

    emit:
    mudata = BUILD_MUDATA.out.mudata
    integrated = integrated_output
    stats = BUILD_MUDATA.out.stats
    mofa_model = mofa_model
    mofa_plots = mofa_plots
    rna_for_dorc = EXPORT_MUDATA_RNA.out

    // MultiVI outputs
    multivi_model = params.multivi.run ? MULTIVI_INTEGRATE.out.model : Channel.empty()
    multivi_output = params.multivi.run ? MULTIVI_INTEGRATE.out.integrated : Channel.empty()
    multivi_plots = params.multivi.run ? MULTIVI_VISUALIZE.out.plots : Channel.empty()
    multivi_metrics = params.multivi.run ? MULTIVI_VISUALIZE.out.metrics : Channel.empty()
    multivi_report = params.multivi.run ? MULTIVI_VISUALIZE.out.report : Channel.empty()

    // Bootstrap-specific outputs for debugging
    bootstrap_summary = (params.mofa.mode == 'bootstrap') ?
        BOOTSTRAP_MOFA_INTEGRATION.out.summary : Channel.empty()
    memory_log = (params.mofa.mode == 'bootstrap') ?
        ANALYZE_MEMORY_LOG.out.memory_stats : Channel.empty()
    consensus_results = (params.mofa.mode == 'bootstrap') ?
        CONSENSUS_ANALYSIS.out.results_dir : Channel.empty()
}

// ============================================================================
// MULTIOME GRN WORKFLOW (pycisTopic + SCENIC+ + DORC)
// ============================================================================
workflow MULTIOME_GRN {
    take:
    rna_h5ad
    atac_peak_matrix
    metadata_csv
    rna_for_dorc
    mudata_stats
    blacklist_bed
    cell_type_key

    main:
    log.info """
    MULTIOME GRN WORKFLOW
    - pycisTopic: topics, DARs, gene activity
    - SCENIC+: eRegulons & AUCell
    - DORC: peak-gene associations and DORC scores
    """

    // pycisTopic preparation — 2 modes (same Phase 2-3, differ only in Phase 1):
    //   atac_only=false (default): build cell metadata from RNA h5ad (multiome)
    //   atac_only=true           : build cell metadata from ATAC peak matrix directly
    if (params.pycistopic.run) {
        log.info "Running pycisTopic preparation (atac_only=${params.pycistopic.atac_only ?: false})..."

        def _pyc_species = params.pycistopic.species
            ?: [human: 'hsapiens', mouse: 'mmusculus'].get(params.species, params.species)

        // ── Phase 1 (mode-dependent) ─────────────────────────────────────────
        def _phase1_dir_ch        = null
        def _phase1_group_list_ch = null
        def _phase1_cell_meta_ch  = null
        def _phase1_qc_dir_ch     = null
        def _phase1_blacklist_ch  = null
        def _phase1_bigwigs_ch    = null

        def _pyc_cond_col  = params.pycistopic.condition_col ?: 'condition'
        def _pyc_min_cells = params.pycistopic.min_cells     ?: 100

        if (params.pycistopic.atac_only) {
            // ATAC-only: cell types come from ATAC peak matrix obs (no RNA needed)
            PYCISTOPIC_ATAC_PREPARE(
                atac_peak_matrix,
                metadata_csv,
                _pyc_species,
                blacklist_bed,
                file(params.pycistopic.gtf),
                cell_type_key,
                _pyc_cond_col,
                _pyc_min_cells
            )
            _phase1_dir_ch        = PYCISTOPIC_ATAC_PREPARE.out.phase1_dir
            _phase1_group_list_ch = PYCISTOPIC_ATAC_PREPARE.out.group_list
            _phase1_cell_meta_ch  = PYCISTOPIC_ATAC_PREPARE.out.cell_metadata
            _phase1_qc_dir_ch     = PYCISTOPIC_ATAC_PREPARE.out.qc_dir
            _phase1_blacklist_ch  = PYCISTOPIC_ATAC_PREPARE.out.blacklist
            _phase1_bigwigs_ch    = PYCISTOPIC_ATAC_PREPARE.out.pseudobulk_bigwigs
        } else {
            // Multiome (default): cell types come from RNA h5ad annotations
            PYCISTOPIC_PHASE1(
                metadata_csv,
                rna_h5ad,
                _pyc_species,
                blacklist_bed,
                file(params.pycistopic.gtf),
                cell_type_key,
                _pyc_cond_col,
                _pyc_min_cells
            )
            _phase1_dir_ch        = PYCISTOPIC_PHASE1.out.phase1_dir
            _phase1_group_list_ch = PYCISTOPIC_PHASE1.out.group_list
            _phase1_cell_meta_ch  = PYCISTOPIC_PHASE1.out.cell_metadata
            _phase1_qc_dir_ch     = PYCISTOPIC_PHASE1.out.qc_dir
            _phase1_blacklist_ch  = PYCISTOPIC_PHASE1.out.blacklist
            _phase1_bigwigs_ch    = PYCISTOPIC_PHASE1.out.pseudobulk_bigwigs
        }

        // ── Phase 2: one CistopicObject per CT×condition (parallel fan-out) ──
        def _groups_ch = _phase1_group_list_ch.splitCsv(header: true, sep: '\t')
        PYCISTOPIC_PER_GROUP(_groups_ch.combine(_phase1_dir_ch))

        // ── Phase 3a: merge all per-group PKLs into merged_cistopic.pkl ──────
        PYCISTOPIC_MERGE_OBJECTS(
            PYCISTOPIC_PER_GROUP.out.pkl.collect(),
            _phase1_cell_meta_ch
        )

        // ── Phase 3b: one LDA job per topic count (parallel fan-out) ─────────
        def _topics_list = (params.pycistopic.topics ?: '10,20,30')
            .split(',').collect { it.trim().toInteger() }
        PYCISTOPIC_RUN_LDA(
            Channel.fromList(_topics_list),
            PYCISTOPIC_MERGE_OBJECTS.out.merged_pkl.first()
        )

        // ── Phase 3c: evaluate models, binarize, DARs, region_sets ───────────
        PYCISTOPIC_FINALIZE_LDA(
            PYCISTOPIC_RUN_LDA.out.topic_pkl.collect(),
            PYCISTOPIC_MERGE_OBJECTS.out.merged_pkl,
            _phase1_cell_meta_ch,
            _phase1_qc_dir_ch,
            _phase1_blacklist_ch,
            _pyc_species
        )
    }

    // SCENIC+ via Snakemake
    if (params.scenicplus.run && params.pycistopic.run) {
        log.info "Running SCENIC+ Snakemake pipeline..."

        SCENICPLUS_RUN(
            PYCISTOPIC_FINALIZE_LDA.out.cistopic_obj,
            rna_for_dorc,
            PYCISTOPIC_FINALIZE_LDA.out.region_sets,
            params.scenicplus.ctx_rankings,
            params.scenicplus.ctx_scores,
            params.scenicplus.motif_annotations,
            params.scenicplus.bc_transform_func
        )

        // SCENIC+ visualization
        SCENICPLUS_VISUALIZE(
            SCENICPLUS_RUN.out.scplus_mudata,
            SCENICPLUS_RUN.out.aucell_direct,
            SCENICPLUS_RUN.out.ereg_direct,
            cell_type_key
        )

        // SCENIC+ GRN network visualization (graph-tool)
        SCENICPLUS_GRN_VIZ(
            SCENICPLUS_RUN.out.ereg_direct,
            SCENICPLUS_RUN.out.aucell_direct,
            SCENICPLUS_VISUALIZE.out.rss_csv,
            SCENICPLUS_RUN.out.tf2g,
            cell_type_key
        )
    }

    // DORC analysis with scPrinter
    if (params.dorc.run) {
        log.info "Running scPrinter DORC analysis..."

        SCPRINTER_DORC(
            atac_peak_matrix,
            rna_for_dorc,
            "none"
        )
    }

    emit:
    // pycisTopic outputs
    cistopic_obj       = params.pycistopic.run ? PYCISTOPIC_FINALIZE_LDA.out.cistopic_obj   : Channel.empty()
    region_sets        = params.pycistopic.run ? PYCISTOPIC_FINALIZE_LDA.out.region_sets     : Channel.empty()
    gene_activity      = params.pycistopic.run ? PYCISTOPIC_FINALIZE_LDA.out.gene_activity   : Channel.empty()
    pseudobulk_bigwigs = params.pycistopic.run
        ? (params.pycistopic.atac_only
            ? PYCISTOPIC_ATAC_PREPARE.out.pseudobulk_bigwigs
            : PYCISTOPIC_PHASE1.out.pseudobulk_bigwigs)
        : Channel.empty()

    // SCENIC+ outputs
    scplus_mudata      = (params.scenicplus.run && params.pycistopic.run) ? SCENICPLUS_RUN.out.scplus_mudata      : Channel.empty()
    scplus_acc_gex     = (params.scenicplus.run && params.pycistopic.run) ? SCENICPLUS_RUN.out.acc_gex_mudata     : Channel.empty()
    scplus_ereg_direct = (params.scenicplus.run && params.pycistopic.run) ? SCENICPLUS_RUN.out.ereg_direct        : Channel.empty()
    scplus_r2g         = (params.scenicplus.run && params.pycistopic.run) ? SCENICPLUS_RUN.out.r2g                : Channel.empty()
    scplus_ereg_ext    = (params.scenicplus.run && params.pycistopic.run) ? SCENICPLUS_RUN.out.ereg_extended      : Channel.empty()
    scplus_aucell_dir  = (params.scenicplus.run && params.pycistopic.run) ? SCENICPLUS_RUN.out.aucell_direct      : Channel.empty()
    scplus_aucell_ext  = (params.scenicplus.run && params.pycistopic.run) ? SCENICPLUS_RUN.out.aucell_extended    : Channel.empty()

    // DORC outputs
    dorc_all      = params.dorc.run ? SCPRINTER_DORC.out.dorc_all    : Channel.empty()
    dorc_sig      = params.dorc.run ? SCPRINTER_DORC.out.dorc_sig    : Channel.empty()
    dorc_scores   = params.dorc.run ? SCPRINTER_DORC.out.dorc_scores : Channel.empty()
}

// ============================================================================
// ENHANCER FOOTPRINTING RECIPES (A/B/C/D)
// ============================================================================
workflow ENHANCER_FOOTPRINTING_RECIPES {
    take:
    peak_matrix
    printer
    chromvar_dev_ch
    chromvar_motifs_ch
    cicero_conns_ch
    cicero_ccan_ch
    ereg_direct_ch
    r2g_ch
    dorc_sig_ch
    rna_h5ad_ch
    cellchat_csv_ch
    pseudobulk_bigwigs_ch
    cell_type_col             // FIX-43: obs column name for cell types
    tf_targets_ch             // tf_target_genes.json from MAP_TF_TO_TARGET_GENES (Phase 3)
    anndataset_ch             // D1b: ATAC AnnDataSet for EXPORT_ATAC_BIGWIGS
    scprinter_footprints_ch   // 2026-04-30: pass-through sync channel; replaces
                              //   direct SCPRINTER_FOOTPRINTING.out access (DSL2
                              //   cross-workflow scope violation fix)
    tf_diff_ch                // 2026-05-04: pass-through for DIFFERENTIAL_TF_ACCESSIBILITY
                              //   .out.tf_diff (DSL2 cross-workflow scope fix +
                              //   correct gate alignment with differential_tf.run).
                              //   Empty when differential_tf.run=false.
    gene_coordinates_ch       // 2026-05-06: pass-through for RESOLVE_GENE_COORDINATES
                              //   .out.coordinates — feeds PROMOTER_MSFP_PER_CT
                              //   (DSL2 cross-workflow scope fix). Empty when
                              //   scprinter.run=false.
    cicero_strat_ctrl_links_ch        // 2026-05-06: per-condition CCAN→gene-link TSVs
    cicero_strat_trt_links_ch         //   from EXTRACT_CCAN_ENHANCERS_{CTRL,TRT}.
    cicero_strat_ctrl_peaks_ch        //   Powers cis-rewiring; empty when stratified
    cicero_strat_trt_peaks_ch         //   Cicero didn't run.
    cicero_strat_ctrl_connections_ch  // 2026-06-01: raw stratified Cicero connections
    cicero_strat_trt_connections_ch   //   (Peak1/Peak2/coaccess .tsv.gz) for lollipop.

    main:

    // ================================================================
    // PHASE 1: ATAC-Only Enhancer Footprinting (Recipe A)
    // ================================================================
    log.info "ENHANCER FOOTPRINTING RECIPES: Phase 1 (ATAC-only)"

    // SHI Tier A wiring (2026-05-04): declare bigwig channels at workflow
    // scope so the emit: block can expose them regardless of which Phase 4
    // branch fires. Populated below inside the EXPORT_ATAC_BIGWIGS branch.
    def shi_bw_manifest_ch = Channel.empty()
    def shi_bw_dir_ch      = Channel.empty()

    // FIX-R3-3: Use species-appropriate GTF for enhancer extraction
    def enhancer_gtf = params.species == 'human' ? params.scprinter.gtf_human : params.scprinter.gtf_mouse

    EXTRACT_CCAN_ENHANCERS(
        cicero_conns_ch,
        cicero_ccan_ch,
        enhancer_gtf,
        ''  // condition_label: empty for default run
    )

    MOTIF_SCAN_ENHANCERS(
        EXTRACT_CCAN_ENHANCERS.out.enhancer_peaks,
        chromvar_motifs_ch
    )

    // D1b: condition-aware control/treatment from differential params.
    // Replaces D1a empty-string placeholders.
    def enh_ctrl = (params.differential?.run ?: false) ?
        (params.differential.control_condition   ?: '') : ''
    def enh_trt  = (params.differential?.run ?: false) ?
        (params.differential.treatment_condition ?: '') : ''

    // D1b: per-CT vs sharded toggle. Both paths emit the same
    // {summary, footprints, binding_scores, global_plots} structure;
    // bind enh_fp_* references so downstream wiring is arch-agnostic.
    def enh_use_per_ct = params.enhancer_footprinting.use_per_ct ?: false
    // 2026-05-26: MSFP compute gate. false → skip the expensive
    // ENHANCER_FOOTPRINTING_PER_CT / ENHANCER_FOOTPRINTING wall-time and all
    // downstream consumers. MOTIF_SCAN_ENHANCERS + Cicero + browser/lollipop
    // still run. Instance configs flip true for full MSFP runs.
    def msfp_enabled = params.enhancer_footprinting.msfp_enabled ?: false
    if (!msfp_enabled) {
        log.info "ENHANCER FOOTPRINTING RECIPES: MSFP compute DISABLED (msfp_enabled=false)"
    }

    def enh_fp_summary        = Channel.empty()
    def enh_fp_footprints     = Channel.empty()
    def enh_fp_binding_scores = Channel.empty()
    def enh_fp_global_plots   = Channel.empty()
    // Hoisted so ENHANCER_FOOTPRINTING_PER_CT_STRIP (the strip re-run, further
    // down in its own msfp_strip gate) can join the same per-CT manifest+beds
    // channel the first pass used. Declared with `def` inside the gate below, it
    // would not be visible there.
    def ch_per_ct_input       = Channel.empty()

    if (msfp_enabled && enh_use_per_ct) {
        // Per-CT fan-out: one task per ct loads printer/peak matrix once
        // and iterates all TFs. Cuts task count ~657 → ~10/33.
        def per_ct_manifest_dir = file("${workflow.workDir}/per_ct_manifests")
        per_ct_manifest_dir.mkdirs()

        ch_per_ct_input = MOTIF_SCAN_ENHANCERS.out.manifest
            .flatMap { manifest_file ->
                def data = new groovy.json.JsonSlurper().parseText(manifest_file.text)
                data.region_sets.collect { entry ->
                    tuple(entry.cell_type, [tf: entry.tf, bed_file: entry.bed_file])
                }
            }
            .groupTuple()
            .combine(MOTIF_SCAN_ENHANCERS.out.region_sets)
            .map { ct, entries, region_sets_dir ->
                def safe_ct = ct.replaceAll(/[\/\s\(\)]+/, '_')
                def manifest_path = file("${per_ct_manifest_dir}/manifest_${safe_ct}.json")
                def newContent = groovy.json.JsonOutput.toJson(entries)
                if (!manifest_path.exists() || manifest_path.text != newContent) {
                    manifest_path.text = newContent
                }
                def beds = entries.collect { e -> file("${region_sets_dir}/${e.bed_file}") }
                tuple(ct, manifest_path, beds)
            }

        ENHANCER_FOOTPRINTING_PER_CT(
            ch_per_ct_input,
            printer,
            peak_matrix,
            cell_type_col,
            enh_ctrl,
            enh_trt
        )

        enh_fp_summary        = ENHANCER_FOOTPRINTING_PER_CT.out.summary
        enh_fp_footprints     = ENHANCER_FOOTPRINTING_PER_CT.out.footprints
        enh_fp_binding_scores = ENHANCER_FOOTPRINTING_PER_CT.out.binding_scores
        enh_fp_global_plots   = ENHANCER_FOOTPRINTING_PER_CT.out.global_plots

    } else if (msfp_enabled) {
        // Sharded fan-out: one task per (cell_type, TF) pair (legacy).
        ch_enhancer_tasks = MOTIF_SCAN_ENHANCERS.out.manifest
            .flatMap { manifest_file ->
                def data = new groovy.json.JsonSlurper().parseText(manifest_file.text)
                data.region_sets.collect { entry ->
                    tuple(entry.cell_type, entry.tf, entry.bed_file)
                }
            }

        ch_enhancer_fp = ch_enhancer_tasks
            .combine(MOTIF_SCAN_ENHANCERS.out.region_sets)
            .map { ct, tf, bed_fname, region_sets_dir ->
                tuple(file("${region_sets_dir}/${bed_fname}"), ct, tf)
            }

        ENHANCER_FOOTPRINTING(
            ch_enhancer_fp.map { it[0] },
            printer,
            peak_matrix,
            ch_enhancer_fp.map { it[1] },
            ch_enhancer_fp.map { it[2] },
            cicero_conns_ch.ifEmpty(file('NO_CICERO_CONNS')),
            cell_type_col,
            enh_ctrl,
            enh_trt
        )

        enh_fp_summary        = ENHANCER_FOOTPRINTING.out.summary
        enh_fp_footprints     = ENHANCER_FOOTPRINTING.out.footprints
        enh_fp_binding_scores = ENHANCER_FOOTPRINTING.out.binding_scores
        enh_fp_global_plots   = ENHANCER_FOOTPRINTING.out.global_plots
    }

    // ================================================================
    // PHASE 3: TF-Gene Regulatory Network (ATAC-only, Shi et al. TF_Net equivalent)
    //   Uses continuous scPrinter binding scores × Cicero co-accessibility.
    //   Absorbed from SDas_nf. No RNA dependency.
    // ================================================================
    if (msfp_enabled && (params.enhancer_footprinting.build_network ?: false)) {
        log.info "ENHANCER FOOTPRINTING RECIPES: Phase 3 (TF-gene regulatory network)"
        def tf_gtf = params.species == 'human' ?
            params.scprinter.gtf_human : params.scprinter.gtf_mouse
        BUILD_TF_GENE_NETWORK(
            tf_targets_ch,
            EXTRACT_CCAN_ENHANCERS.out.gene_links,
            cicero_conns_ch.ifEmpty(file('NO_CICERO_CONNS')).first(),
            enh_fp_summary.collect(),
            tf_gtf
        )
        PLOT_TF_GENE_NETWORK(BUILD_TF_GENE_NETWORK.out.adjacency)
    }

    // ================================================================
    // PHASE 2: Multiome Integration (Recipe B)
    // ================================================================
    def has_scenic = (params.scenicplus.run ?: false) && (params.pycistopic.run ?: false)
    def has_dorc = (params.dorc.run ?: false)

    if (has_scenic && msfp_enabled) {
        log.info "ENHANCER FOOTPRINTING RECIPES: Phase 2 (Multiome integration)"

        def dorc_sig_file = has_dorc ?
            dorc_sig_ch : Channel.value(file('NO_FILE_dorc'))

        EXTRACT_EREGULON_REGIONS(
            ereg_direct_ch,
            r2g_ch,
            dorc_sig_file
        )

        CROSS_MODAL_VALIDATION(
            enh_fp_footprints.collect(),
            rna_h5ad_ch,
            EXTRACT_EREGULON_REGIONS.out.target_genes
        )
    }

    // ================================================================
    // PHASE 3: CellChat-Guided Footprinting (Recipe C)
    // ================================================================
    if (params.enhancer_recipe_c.run) {
        log.info "ENHANCER FOOTPRINTING RECIPES: Phase 3 (CellChat-guided)"

        // FIX-R1-17/F9: Validate pathway_to_tfs.json at point of use (pipeline-bundled file).
        // Fail early with a clear message instead of a cryptic JSON parse error.
        def ptfJson = file("${projectDir}/data/pathway_to_tfs.json")
        if (!ptfJson.exists()) {
            error "Recipe C requires data/pathway_to_tfs.json but the file is missing from the project directory: ${ptfJson}"
        }
        try {
            def ptfContent = new groovy.json.JsonSlurper().parseText(ptfJson.text)
            if (!(ptfContent instanceof Map) || ptfContent.size() == 0) {
                error "data/pathway_to_tfs.json is empty or not a JSON object. Expected {pathway: [TF1, TF2, ...], ...}"
            }
            // Spot-check: at least one entry should map to a non-empty list
            def hasValidEntry = ptfContent.any { k, v -> v instanceof List && v.size() > 0 }
            if (!hasValidEntry) {
                error "data/pathway_to_tfs.json has no valid pathway->TF mappings. Each key should map to a non-empty list of TF gene names."
            }
            log.info "  pathway_to_tfs.json validated: ${ptfContent.size()} pathways"
        } catch (groovy.json.JsonException e) {
            error "data/pathway_to_tfs.json is malformed JSON: ${e.message}"
        }

        CELLCHAT_TO_TF_HYPOTHESES(
            cellchat_csv_ch,
            chromvar_dev_ch,
            ptfJson
        )

        def ereg_regions_ch = (has_scenic && msfp_enabled) ?
            EXTRACT_EREGULON_REGIONS.out.region_sets :
            Channel.value(file('NO_FILE_ereg_regions'))
        def ereg_manifest_ch = (has_scenic && msfp_enabled) ?
            EXTRACT_EREGULON_REGIONS.out.manifest :
            Channel.value(file('NO_FILE_ereg_manifest'))
        def ccan_regions_ch = MOTIF_SCAN_ENHANCERS.out.region_sets
        def ccan_manifest_ch = MOTIF_SCAN_ENHANCERS.out.manifest

        EXTRACT_SIGNALING_TARGETS(
            CELLCHAT_TO_TF_HYPOTHESES.out.validated_pairs,
            ereg_regions_ch,
            ereg_manifest_ch,
            ccan_regions_ch,
            ccan_manifest_ch
        )

        if (msfp_enabled) {
            SIGNAL_CHAIN_CORRELATION(
                enh_fp_footprints.collect(),
                rna_h5ad_ch,
                EXTRACT_SIGNALING_TARGETS.out.metadata
            )
        }
    }

    // ================================================================
    // PHASE 4: Composite Enhancer Visualization (Recipe D)
    // ================================================================
    if (params.enhancer_viz.run) {
        log.info "ENHANCER FOOTPRINTING RECIPES: Phase 4 (Composite Visualization)"

        def viz_cell_type_col = params.enhancer_viz.cell_type_col ?: cell_type_col
        def viz_min_cells     = params.enhancer_viz.min_cells     ?: 100

        // D1b: ATAC-only pseudobulk bigWigs via snap.ex.export_coverage.
        //
        // EXPORT_ATAC_BIGWIGS always runs: it produces BROAD-CLASS bigwigs
        // (cell_type_broad from scATAnno) for SHI MARKER_COVERAGE_TRACKS.
        //
        // For SCENIC+ datasets pycistopic already produces fine-grained bigwigs
        // for PREPARE_ENHANCER_VIZ_TRACKS (fine resolution, condition-aware).
        // Those two purposes are distinct: SHI needs broad classes; enhancer viz
        // tracks need fine labels.  Previously gating EXPORT_ATAC_BIGWIGS with
        // !has_scenic caused MARKER_COVERAGE_TRACKS to silently skip on all
        // SCENIC+ instances (AD, Brain_Mm_BD, etc).  2026-05-24 fix: always run.
        // 2026-05-26: pass condition_col so export_atac_bigwigs.py also writes
        // per-CT × condition BigWigs into manifest.by_condition. Required for
        // RENDER_GENOME_BROWSER in differential mode. When enhancer_viz.condition_col
        // is null, the script skips the second pass — backwards compatible.
        def bw_condition_col = params.enhancer_viz?.condition_col ?: 'none'

        EXPORT_ATAC_BIGWIGS(
            anndataset_ch,
            peak_matrix,
            viz_cell_type_col,
            viz_min_cells,
            bw_condition_col
        )
        // SHI Tier A wiring: broad-class manifest always available now.
        shi_bw_manifest_ch = EXPORT_ATAC_BIGWIGS.out.manifest
        shi_bw_dir_ch      = EXPORT_ATAC_BIGWIGS.out.bigwigs

        // 2026-05-26: RENDER_GENOME_BROWSER — per-gene matplotlib browser tracks.
        // Gated on browser_viz.enabled; target_genes shared with enhancer_viz.
        if (params.browser_viz?.enabled ?: false) {
            def browser_mode      = params.browser_viz?.mode ?: 'absolute'
            def browser_genes     = params.enhancer_viz?.target_genes ?: []
            def browser_cell_types = params.browser_viz?.cell_types ?:
                                     params.enhancer_viz?.target_genes ? null : null
            def browser_gtf       = file(params.species == 'human' ?
                                    params.scprinter.gtf_human : params.scprinter.gtf_mouse)
            def browser_ctrl      = enh_ctrl ?: 'none'
            def browser_trt       = enh_trt  ?: 'none'
            def no_da             = file('NO_FILE')

            if (browser_genes && !browser_genes.isEmpty()) {
                log.info "GENOME BROWSER: enabled for ${browser_genes.size()} genes, mode=${browser_mode}"

                // FIX (attempt 3): combine() always splatts the top-level List from the
                // right-hand channel into the combined tuple element-by-element.
                // Removing the no-op .map didn't help — combine(ch_browser_bws) where
                // ch_browser_bws emits [bw1,...,bwN] still produces a flat N+2-element tuple.
                //
                // Correct approach: use flatMap to embed bws INSIDE a [gene, bws] parent tuple
                // BEFORE any combine. combine() only splatts the top-level emission; nested
                // list elements inside a tuple slot are preserved.
                //   flatMap → [gene, [bw1,...,bwN]]       (bws nested as element[1])
                //   .combine(manifest) → [gene, [bw1,...,bwN], manifest]  ← bws intact ✓
                def ch_browser_bws      = EXPORT_ATAC_BIGWIGS.out.bigwigs.collect()
                def ch_browser_manifest = EXPORT_ATAC_BIGWIGS.out.manifest

                def ch_browser_in = ch_browser_bws
                    .flatMap { bws_raw ->
                        // Normalize: glob output emits List<Path>; .collect() on a 1-item
                        // channel may double-wrap to [[bw1,...,bwN]]. flatten() collapses
                        // one nesting level safely whether single- or double-wrapped.
                        def bws = bws_raw instanceof List ? bws_raw.flatten() : [bws_raw]
                        browser_genes.collect { gene -> tuple(gene, bws) }
                    }
                    .combine(ch_browser_manifest)
                    .map { gene, bws, manifest ->
                        def ct_str = (browser_cell_types instanceof List)
                            ? browser_cell_types.join(',')
                            : (browser_cell_types ?: '')
                        tuple(gene, manifest, bws, ct_str, browser_mode)
                    }

                RENDER_GENOME_BROWSER(
                    ch_browser_in,
                    browser_gtf,
                    browser_ctrl,
                    browser_trt,
                    no_da
                )
            } else {
                log.warn "GENOME BROWSER: enabled but no target genes configured (set enhancer_viz.target_genes)"
            }
        }

        // PREPARE_ENHANCER_VIZ_TRACKS gets fine-grained bigwigs from pycistopic
        // (SCENIC+ path) or from EXPORT_ATAC_BIGWIGS (non-SCENIC+ path).
        def bigwig_dir_ch = has_scenic ?
            pseudobulk_bigwigs_ch.ifEmpty(file('NO_BIGWIGS')).collect() :
            EXPORT_ATAC_BIGWIGS.out.bigwigs.collect()

        PREPARE_ENHANCER_VIZ_TRACKS(
            cicero_conns_ch,
            EXTRACT_CCAN_ENHANCERS.out.enhancer_peaks,
            bigwig_dir_ch,
            file(params.species == 'human' ? params.scprinter.gtf_human : params.scprinter.gtf_mouse),
            params.enhancer_viz.target_genes,
            MOTIF_SCAN_ENHANCERS.out.manifest
        )

        // D1b: Swarup-chain (gene, TF, ct) candidate filter. Replaces the old
        // genes × TFs cartesian fan-out for COMPOSITE_ENHANCER_VIZ.
        // 2026-05-04: gate aligned to differential_tf.run (the param that
        // actually invokes DIFFERENTIAL_TF_ACCESSIBILITY at main.nf:1884) and
        // input piped via REGULATORY_ANALYSIS.out.tf_diff to satisfy DSL2
        // cross-workflow scope rules. AD's first run hit the prior bug because
        // differential.run=true triggered the ternary's true branch even
        // though differential_tf.run=false meant DTA was never invoked.
        def diff_csvs_ch = tf_diff_ch.collect().ifEmpty([])

        BUILD_VIZ_CANDIDATES(
            PREPARE_ENHANCER_VIZ_TRACKS.out.track_manifest,
            MOTIF_SCAN_ENHANCERS.out.manifest,
            MOTIF_SCAN_ENHANCERS.out.region_sets,
            diff_csvs_ch
        )

        // D1b: per-(ct, TF) scalar metric rollup; gated on msfp_enabled + build_network
        // so BUILD_TF_GENE_NETWORK.out.adjacency is available.
        if (msfp_enabled && (params.enhancer_footprinting.build_network ?: false)) {
            def promoter_fp_dir_ch = scprinter_footprints_ch
                .collect()
                .map { _files -> file("${params.outdir}/scprinter/footprints") }
                .first()

            AGGREGATE_FP_STATS(
                enh_fp_footprints.collect(),
                enh_fp_binding_scores.collect(),
                diff_csvs_ch,
                BUILD_TF_GENE_NETWORK.out.adjacency,
                MOTIF_SCAN_ENHANCERS.out.region_sets,
                MOTIF_SCAN_ENHANCERS.out.manifest,
                PREPARE_ENHANCER_VIZ_TRACKS.out.track_manifest,
                promoter_fp_dir_ch
            )
        }

        // D1b: per-(gene, TF, ct) fan-out replaces D1a placeholders that
        // passed '' for cell_type and NO_FILE for promoter dir.
        def ch_viz_tasks = BUILD_VIZ_CANDIDATES.out.candidates
            .splitCsv(header: true)
            .map { row -> tuple(row.gene as String, row.tf as String, row.cell_type as String) }

        // D1b: collapse fp PNG inputs to a single dir symlink. .collect() acts
        // as the synchronization barrier; per-ct architecture publishes to a
        // different root, so pick by use_per_ct.
        // 2026-05-26: .ifEmpty([]) ensures the channel fires immediately when
        // msfp_enabled=false so COMPOSITE_ENHANCER_VIZ is not blocked.
        def fp_root = enh_use_per_ct ?
            file("${params.outdir}/enhancer_footprinting_per_ct") :
            file("${params.outdir}/enhancer_footprinting/footprints")
        def fp_dir_ch = enh_fp_global_plots
            .collect()
            .ifEmpty([])
            .map { _files -> fp_root }
            .first()

        // D1b: Tier-1 promoter MSFP staging — same dir-symlink trick.
        def promoter_dir_ch = scprinter_footprints_ch
            .collect()
            .map { _files -> file("${params.outdir}/scprinter/footprints") }
            .first()

        COMPOSITE_ENHANCER_VIZ(
            PREPARE_ENHANCER_VIZ_TRACKS.out.track_manifest.first(),
            PREPARE_ENHANCER_VIZ_TRACKS.out.track_inis.collect(),
            MOTIF_SCAN_ENHANCERS.out.manifest.first(),
            ch_viz_tasks.map { it[0] },
            ch_viz_tasks.map { it[1] },
            ch_viz_tasks.map { it[2] },
            fp_dir_ch,
            promoter_dir_ch
        )
    }

    // ================================================================
    // PHASE 5: Promoter MSFP + Motif Overlay (per-condition, TF-binding scale)
    //   Per-(cell_type, gene) compute -> render PNG with stacked per-condition
    //   heatmaps at TF-binding scales (<=30 bp by default) above a JASPAR
    //   motif track for the top-N TFs ranked by AGGREGATE_FP_STATS.out
    //   .triple_csv ((n_sites_in_gene_window desc, bind_dip_depth_mean desc)).
    //   Differential-only by gate; single-condition runs auto-skip.
    //   CT enumeration reuses MOTIF_SCAN_ENHANCERS.out.manifest, which is
    //   already filtered by the chromvar resolution-floor gate
    //   (feedback_celltype_resolution_floor).
    // ================================================================
    def overlay_enabled = (params.promoter_overlay?.enabled ?: false)
    def overlay_trt = params.promoter_overlay?.treatment ?:
        params.differential?.treatment_condition
    def overlay_ctrl = params.promoter_overlay?.control ?:
        params.differential?.control_condition
    def overlay_cond_col = params.promoter_overlay?.condition_col ?:
        params.differential?.condition_key
    // 2026-05-06: prefer promoter_overlay.target_genes; fall back to
    // scprinter.target_genes for back-compat. Decoupled because setting
    // scprinter.target_genes non-empty flips the pipeline into TARGETED mode,
    // which short-circuits MOTIF_SCAN_ENHANCERS (DISCOVERY-only) and starves
    // the overlay's CT-enumeration channel.
    def overlay_target_genes = (params.promoter_overlay?.target_genes ?: params.scprinter?.target_genes) ?: []
    def overlay_top_n = (params.promoter_overlay?.top_n_tfs ?: 4) as Integer
    def overlay_palette = ['#d62728','#1f77b4','#2ca02c','#9467bd',
                           '#ff7f0e','#17becf','#8c564b','#e377c2']

    if (msfp_enabled &&
        overlay_enabled &&
        overlay_trt && overlay_ctrl &&
        overlay_target_genes && !overlay_target_genes.isEmpty() &&
        (params.enhancer_footprinting.build_network ?: false)) {
        log.info "PROMOTER MSFP OVERLAY: enabled (cts via MOTIF_SCAN_ENHANCERS, " +
                 "TFs via AGGREGATE_FP_STATS triple_csv, top_n=${overlay_top_n})"

        def overlay_genes_csv = overlay_target_genes.join(',')
        def overlay_palette_csv = overlay_palette.take(overlay_top_n).join(',')

        // Resolve coords for the overlay target genes using a 4000bp window
        // (upstream=2000, downstream=2000 via withName:RESOLVE_OVERLAY_COORDINATES).
        // Separate from gene_coordinates_ch which uses the scPrinter footprinting
        // window and covers a different gene set.
        RESOLVE_OVERLAY_COORDINATES(
            params.species,
            Channel.value(overlay_target_genes),
            file('NO_FILE')
        )

        def ch_overlay_compute_in = MOTIF_SCAN_ENHANCERS.out.manifest
            .flatMap { manifest_file ->
                def data = new groovy.json.JsonSlurper().parseText(manifest_file.text)
                data.region_sets.collect { it.cell_type }.unique()
            }
            .map { ct -> tuple(ct, overlay_genes_csv) }

        PROMOTER_MSFP_PER_CT(
            ch_overlay_compute_in,
            printer,
            peak_matrix,
            RESOLVE_OVERLAY_COORDINATES.out.coordinates,
            cell_type_col,
            overlay_cond_col,
            overlay_ctrl,
            overlay_trt
        )

        def ch_fp_per_ct_gene = PROMOTER_MSFP_PER_CT.out.fp
            .transpose()
            .map { ct, f ->
                def stem = f.baseName
                def gene = stem.replaceFirst(/^promoter_fp_/, '').split('__')[0]
                tuple(ct, gene, f)
            }

        def ch_tf_per_ct_gene = AGGREGATE_FP_STATS.out.triple_csv
            .splitCsv(header: true)
            .filter { row -> row.gene && (row.gene in overlay_target_genes) }
            .map { row ->
                def ns = (row.n_sites_in_gene_window ?: '0') as Double
                def bd = (row.bind_dip_depth_mean    ?: '0') as Double
                tuple(row.cell_type as String, row.gene as String,
                      row.tf as String, ns, bd)
            }
            .groupTuple(by: [0,1])
            .map { ct, gene, tfs, ns_list, bd_list ->
                def ranked = [tfs, ns_list, bd_list].transpose()
                    .sort { -((it[1] as Double) * 1000.0 + (it[2] as Double)) }
                    .take(overlay_top_n)
                tuple(ct, gene, ranked.collect { it[0] }.join(','))
            }

        def ch_render_in = ch_fp_per_ct_gene
            .join(ch_tf_per_ct_gene, by: [0,1])
            .map { ct, gene, f, tfs ->
                tuple(ct, gene, f, tfs, overlay_palette_csv)
            }

        RENDER_PROMOTER_MSFP_OVERLAY(ch_render_in, overlay_ctrl, overlay_trt)

        // 2026-05-26: new-style strip render (zoom + logo + ref seq).
        // Gated on msfp_strip.enabled; reuses the same PROMOTER_MSFP_PER_CT h5ads.
        if (params.msfp_strip?.enabled ?: false) {
            def strip_mode = params.msfp_strip?.mode ?: 'all_three'
            log.info "RENDER_MSFP_PROMOTER_STRIP: mode=${strip_mode}"

            // Collapse ch_tf_per_ct_gene (already ranked per (ct, gene)) to per-ct.
            // overlay_palette contains hex colors, not TF names — the script
            // needs real TF names to scan JASPAR motifs.
            def ch_tfs_per_ct = ch_tf_per_ct_gene
                .map { ct, gene, tfs_csv -> tuple(ct, tfs_csv) }
                .groupTuple(by: [0])
                .map { ct, tfs_csvs ->
                    def seen = new LinkedHashSet()
                    tfs_csvs.each { csv -> csv.split(',').each { seen << it.trim() } }
                    tuple(ct, seen.take(overlay_top_n).join(','))
                }

            // Group all h5ads per (ct, TF set) and stage into scan_dir
            def ch_strip_in = PROMOTER_MSFP_PER_CT.out.fp
                .join(ch_tfs_per_ct, by: [0])
                .map { ct, h5ads, tfs ->
                    tuple(ct, tfs, overlay_genes_csv,
                          h5ads instanceof List ? h5ads : [h5ads],
                          strip_mode)
                }

            RENDER_MSFP_PROMOTER_STRIP(ch_strip_in, overlay_ctrl, overlay_trt)
        }
    }

    // 2026-05-30: MSFP enhancer strip — discovery architecture.
    // Phase 2.5: RANK_ENHANCER_STRIP_GENES discovers per-(CT, TF) target genes
    //   from tfbs binding scores × Cicero co-accessibility (ATAC-only, no RNA).
    // Phase 2.6: ENHANCER_FOOTPRINTING_PER_CT_STRIP re-runs with per-CT gene lists
    //   to populate obsm in the h5ads (required for render_msfp_enhancer_strip.py).
    // Phase 5:   RENDER_MSFP_ENHANCER_STRIP renders strips per (CT, TF, gene) using
    //   h5ads from the strip re-run, joined to per-(CT,TF) gene lists.
    //
    // Supersedes the params-driven target_genes list: genes are now discovered
    // per (CT, TF) rather than applied uniformly, so every rendered strip has
    // co-accessibility evidence behind it.
    //
    // Gating: only fires when msfp_enabled + msfp_strip.enabled.
    // Cell-count gate is implicit: RANK only receives h5ads from CTs that already
    // passed the upstream min-cell / ChromVAR z-score gates via ENHANCER_FOOTPRINTING_PER_CT.
    // CTs with no Cicero-ranked genes produce empty gene lists → filter excludes them
    // from the STRIP re-run and RENDER, avoiding null-input tasks.
    //
    // TODO: pycisTopic ATAC-only split from SCENIC+ — when available, pipe
    // pycistopic_topics_ch into RANK_ENHANCER_STRIP_GENES for stronger evidence tier.
    if (msfp_enabled && (params.msfp_strip?.enabled ?: false)) {
        def _strip_mode_raw = params.msfp_strip?.mode ?: 'absolute'
        // enhancer strip only supports absolute|differential; map all_three → differential
        def enh_strip_mode  = (_strip_mode_raw == 'all_three') ? 'differential' : _strip_mode_raw
        def enh_strip_gtf   = file(params.species == 'human' ?
            params.scprinter.gtf_human : params.scprinter.gtf_mouse)

        // ── Phase 2.5: rank per-(CT, TF) target genes ──────────────────────
        RANK_ENHANCER_STRIP_GENES(
            enh_fp_binding_scores.collect(),
            cicero_conns_ch.ifEmpty(file('NO_CICERO')).first(),
            enh_strip_gtf,
            params.msfp_strip?.top_n_regions ?: 100,
            params.msfp_strip?.top_k_genes   ?: 5
        )

        // ── Phase 2.6: re-run footprinting with per-CT strip gene lists ─────
        // Build (ct, strip_genes_csv) channel from per_ct_genes.csv.
        // Filter out CTs with empty gene lists before joining to avoid null-input tasks.
        def ch_ct_strip_genes = RANK_ENHANCER_STRIP_GENES.out.per_ct_genes
            .splitCsv(header: true)
            .filter { row -> row.strip_target_genes && row.strip_target_genes.trim() }
            .map    { row -> tuple(row.cell_type as String,
                                   row.strip_target_genes as String) }

        // Join per-CT footprinting inputs with their strip gene list.
        // ch_per_ct_input = (ct, manifest_json, beds) — same channel used for first pass.
        def ch_strip_per_ct_input = ch_per_ct_input
            .join(ch_ct_strip_genes, by: 0)  // (ct, manifest, beds, strip_genes_csv)

        // Pass strip_genes inside the per-CT tuple so it stays synchronized.
        ENHANCER_FOOTPRINTING_PER_CT_STRIP(
            ch_strip_per_ct_input,   // (ct, manifest, beds, strip_genes_csv) — 4-element tuple
            printer,
            peak_matrix,
            cell_type_col,
            enh_ctrl,
            enh_trt
        )

        // ── Phase 5: render strips per (CT, TF, gene) ──────────────────────
        // Build (ct, tf, gene, h5ad) from strip h5ads joined to per-(CT,TF) gene lists.
        // Both channels are process outputs within this workflow → multicast, no exhaustion.
        def ch_strip_fp_by_name = ENHANCER_FOOTPRINTING_PER_CT_STRIP.out.footprints
            .flatMap { files ->
                (files instanceof List ? files : [files]).collect { f -> tuple(f.name, f) }
            }

        // Filename → (ct, tf) lookup from the motif scan manifest. Reproduces the
        // Python sanitize rule from run_enhancer_footprinting.py, which names the
        // h5ads with a lax replace-only rule: safe = name.replace("/", "_")
        // (slashes ONLY — spaces are preserved).
        def ch_strip_name_to_ct_tf = MOTIF_SCAN_ENHANCERS.out.manifest
            .flatMap { manifest_file ->
                def data = new groovy.json.JsonSlurper().parseText(manifest_file.text)
                data.region_sets.collect { entry ->
                    def safe_ct = (entry.cell_type as String).replace('/', '_')
                    def safe_tf = (entry.tf as String).replace('/', '_')
                    tuple("enhancer_footprints_${safe_ct}_${safe_tf}.h5ad",
                          entry.cell_type as String,
                          entry.tf as String)
                }
            }

        // per_ct_tf_genes.json: {ct: {tf: [genes]}} → (ct, tf, [genes]) channel
        def ch_ct_tf_gene_lists = RANK_ENHANCER_STRIP_GENES.out.per_ct_tf_genes
            .flatMap { json_file ->
                def data = new groovy.json.JsonSlurper().parseText(json_file.text)
                def results = []
                data.each { ct, tf_map ->
                    tf_map.each { tf, genes ->
                        if (genes) results << tuple(ct as String, tf as String, genes as List)
                    }
                }
                results
            }

        def ch_enh_strip_in = ch_strip_fp_by_name
            .join(ch_strip_name_to_ct_tf, by: 0)         // (fname, h5ad, ct, tf)
            .map { fname, h5ad, ct, tf -> tuple(ct, tf, h5ad) }
            .join(ch_ct_tf_gene_lists, by: [0, 1])        // (ct, tf, h5ad, [genes])
            .flatMap { ct, tf, h5ad, genes ->
                genes.collect { gene ->
                    tuple(ct, tf, gene, h5ad, file('NO_FILE'), enh_strip_mode)
                }
            }

        log.info "RENDER_MSFP_ENHANCER_STRIP: discovery mode, mode=${enh_strip_mode}, " +
                 "top_n_regions=${params.msfp_strip?.top_n_regions ?: 100}, " +
                 "top_k_genes=${params.msfp_strip?.top_k_genes ?: 5}"

        RENDER_MSFP_ENHANCER_STRIP(
            ch_enh_strip_in,
            enh_strip_gtf,
            enh_ctrl,
            enh_trt
        )
    }

    // ================================================================
    // PHASE 6: Cis-Rewiring (per-TF gained-CCAN motif-presence panels)
    //   Build union enhancer peakset (ctrl ∪ trt) → motif scan on union →
    //   for top-N TFs (ranked from AGGREGATE_FP_STATS triple_csv): split each
    //   gene's gained enhancer→gene links by whether the enhancer hosts the
    //   TF motif → render directionality bar plot. The union scan is the
    //   methodological fix for the SDas SREBF1/Cers2 caveat — scanning only
    //   the control peakset misses motifs at treatment-only enhancers.
    //   Gates: differential.run + cis_rewiring.enabled + per-condition links
    //   from EXTRACT_CCAN_ENHANCERS_{CTRL,TRT} (i.e. stratified Cicero ran).
    // ================================================================
    def cis_enabled = (params.cis_rewiring?.enabled ?: false)
    def cis_top_n_tfs   = (params.cis_rewiring?.top_n_tfs ?: 10) as Integer
    def cis_top_n_genes = (params.cis_rewiring?.top_n_genes ?: 25) as Integer
    def cis_min_delta   = (params.cis_rewiring?.min_delta ?: 1) as Integer
    // 2026-05-26: cis_rewiring gate decoupled from build_network.
    // When msfp_enabled=false or build_network=false, TFs are sourced from
    // cis_rewiring.target_tfs (explicit list in config) instead of
    // AGGREGATE_FP_STATS ranking. Empty target_tfs with build_network=false
    // means the motif_in_gained / motif_stack blocks fire but lollipop skips.
    if (cis_enabled && overlay_trt && overlay_ctrl) {
        log.info "CIS-REWIRING: enabled (top_n_tfs=${cis_top_n_tfs}, top_n_genes=${cis_top_n_genes}, min_delta=${cis_min_delta})"

        // Build union enhancer peakset (ctrl ∪ trt) — methodological fix
        BUILD_UNION_ENHANCER_PEAKS(
            cicero_strat_ctrl_peaks_ch,
            cicero_strat_trt_peaks_ch
        )

        // Reuse the existing chromvar motif list (DISCOVERY mode); same TFs
        // as the global scan but evaluated against the union peakset.
        MOTIF_SCAN_ENHANCERS_UNION(
            BUILD_UNION_ENHANCER_PEAKS.out.peaks_union,
            chromvar_motifs_ch
        )

        // Top-N TF selection: prefer AGGREGATE_FP_STATS ranking (requires
        // msfp_enabled + build_network); fall back to cis_rewiring.target_tfs.
        def ch_top_tfs
        if (msfp_enabled && (params.enhancer_footprinting.build_network ?: false)) {
            ch_top_tfs = AGGREGATE_FP_STATS.out.triple_csv
                .splitCsv(header: true)
                .map { row ->
                    def ns = (row.n_sites_in_gene_window ?: '0') as Double
                    def bd = (row.bind_dip_depth_mean ?: '0') as Double
                    tuple(row.tf as String, ns, bd)
                }
                .groupTuple(by: 0)
                .map { tf, ns_list, bd_list ->
                    def total_ns = ns_list.collect { it as Double }.sum() ?: 0.0
                    def mean_bd  = bd_list ?
                        (bd_list.collect { it as Double }.sum() / bd_list.size()) : 0.0
                    tuple(tf, total_ns, mean_bd)
                }
                .toSortedList { a, b -> (b[1] <=> a[1]) ?: (b[2] <=> a[2]) }
                .flatMap { sorted -> sorted.take(cis_top_n_tfs) }
                .map { tf, ns, bd -> tf }
        } else {
            def explicit_tfs = params.cis_rewiring?.target_tfs ?: []
            if (explicit_tfs) {
                log.info "CIS-REWIRING: using cis_rewiring.target_tfs (${explicit_tfs.size()} TFs) — build_network/msfp_enabled not active"
                ch_top_tfs = Channel.fromList(explicit_tfs.take(cis_top_n_tfs))
            } else {
                log.warn "CIS-REWIRING: msfp_enabled=false and cis_rewiring.target_tfs is empty — motif_in_gained/motif_stack will fire but lollipop will not"
                ch_top_tfs = Channel.empty()
            }
        }

        // Per-TF motif BED extraction from the union scan's region_sets dir.
        // region_sets manifest names them as <safe_ct>_<TF>.bed; for cis-rewiring
        // we want a single TF-wide BED across all CTs (any enhancer hosting the
        // TF motif anywhere). Manifest carries per-(ct,tf) BED paths; collect
        // and union per TF below.
        def ch_motif_per_tf = MOTIF_SCAN_ENHANCERS_UNION.out.manifest
            .flatMap { manifest_file ->
                def data = new groovy.json.JsonSlurper().parseText(manifest_file.text)
                data.region_sets.collect { rs ->
                    tuple(rs.tf as String, rs.bed_file as String)
                }
            }
            .combine(MOTIF_SCAN_ENHANCERS_UNION.out.region_sets)
            .map { tf, bed_fname, region_sets_dir ->
                tuple(tf, file("${region_sets_dir}/${bed_fname}"))
            }
            .groupTuple(by: 0)
            .map { tf, beds ->
                tuple(tf, beds)
            }

        // Join top-N TFs against per-TF BEDs, attach per-condition links.
        // Module concatenates per-CT BEDs into a TF-wide motif BED internally.
        def ch_motif_in_gained = ch_top_tfs
            .map { tf -> tuple(tf, true) }
            .join(ch_motif_per_tf)
            .map { tf, _flag, beds ->
                tuple(tf, beds.findAll { it.size() > 0 })
            }
            .filter { tf, beds -> beds && !beds.isEmpty() }
            .combine(cicero_strat_ctrl_links_ch)
            .combine(cicero_strat_trt_links_ch)
            .map { tf, beds, links_ctrl, links_trt ->
                tuple(tf, beds, links_ctrl, links_trt)
            }

        if (params.cis_rewiring?.motif_in_gained_enabled ?: true) {
            MOTIF_IN_GAINED_CCANS(ch_motif_in_gained, overlay_ctrl, overlay_trt)
        }
        if (params.cis_rewiring?.motif_stack_enabled ?: true) {
            RENDER_CIS_REWIRING_MOTIF_STACK(ch_motif_in_gained, overlay_ctrl, overlay_trt)
        }

        // 2026-05-26: RENDER_CICERO_LOLLIPOP — per-(CT, TF) lollipop chart
        // of Δ CCAN arc counts comparing ctrl vs trt.
        // Gated on cis_rewiring.lollipop_enabled (default true).
        // Requires: ch_top_tfs (non-empty) + per-CT motif BEDs from
        // MOTIF_SCAN_ENHANCERS_UNION + per-condition cicero gz files.
        if ((params.cis_rewiring?.lollipop_enabled ?: true) &&
            (cicero_strat_ctrl_connections_ch || cicero_strat_trt_connections_ch)) {
            def lollipop_gtf = file(params.species == 'human' ?
                params.scprinter.gtf_human : params.scprinter.gtf_mouse)

            // Per-CT motif BED for each top-N TF (across all CTs in union scan)
            def ch_lollipop_motif = MOTIF_SCAN_ENHANCERS_UNION.out.manifest
                .flatMap { manifest_file ->
                    def data = new groovy.json.JsonSlurper().parseText(manifest_file.text)
                    data.region_sets.collect { rs ->
                        tuple(rs.tf as String, rs.cell_type as String, rs.bed_file as String)
                    }
                }
                .combine(MOTIF_SCAN_ENHANCERS_UNION.out.region_sets)
                .map { tf, ct, bed_fname, region_sets_dir ->
                    tuple(tf, ct, file("${region_sets_dir}/${bed_fname}"))
                }
                // For lollipop we want one BED per (CT, TF); use the first CT
                // with a motif BED that passes. The lollipop script receives
                // the per-CT CCAN base so CT is the join key.
                .groupTuple(by: [0, 1])
                .map { tf, ct, beds -> tuple(tf, ct, beds[0]) }

            // Restrict to top-N TFs
            def ch_lollipop_in = ch_top_tfs
                .map { tf -> tuple(tf, true) }
                .join(
                    ch_lollipop_motif.map { tf, ct, bed -> tuple(tf, ct, bed) }
                        .groupTuple(by: 0)
                        .map { tf, cts, beds -> tuple(tf, [cts, beds].transpose()) },
                    by: 0
                )
                .flatMap { tf, _flag, ct_bed_pairs ->
                    ct_bed_pairs.collect { ct, bed -> tuple(ct, tf, bed) }
                }
                // Add per-condition cicero connections (raw Peak1/Peak2/coaccess gz)
                .combine(cicero_strat_ctrl_connections_ch)
                .combine(cicero_strat_trt_connections_ch)
                .map { ct, tf, bed, ctrl_gz, trt_gz ->
                    tuple(ct, tf, bed, ctrl_gz, trt_gz)
                }

            RENDER_CICERO_LOLLIPOP(
                ch_lollipop_in,
                lollipop_gtf,
                overlay_ctrl,
                overlay_trt
            )
        }
    }

    emit:
    enhancer_peaks     = EXTRACT_CCAN_ENHANCERS.out.enhancer_peaks
    motif_scan         = MOTIF_SCAN_ENHANCERS.out.motif_scan
    region_sets        = MOTIF_SCAN_ENHANCERS.out.region_sets
    enhancer_fps       = enh_fp_footprints
    cross_modal        = (has_scenic && msfp_enabled) ? CROSS_MODAL_VALIDATION.out.validation_table : Channel.empty()
    evidence_tiers     = (params.enhancer_recipe_c.run && msfp_enabled) ? SIGNAL_CHAIN_CORRELATION.out.evidence_tiers : Channel.empty()
    enhancer_viz       = (params.enhancer_viz.run ?: false) ? COMPOSITE_ENHANCER_VIZ.out.composite_png : Channel.empty()
    shi_bw_manifest    = shi_bw_manifest_ch
    shi_bw_dir         = shi_bw_dir_ch
}

// ============================================================================
// MAIN WORKFLOW (UNIFIED ENTRY POINT)
// ============================================================================
workflow {

    // FIX-P0: Run startup validation before any compute
    validateStartupParams()

    log.info """

    REFACTOR5: Unified Multiomics Pipeline v3.0.0
    RNA + ATAC Integration & Regulatory Analysis

    Species:       ${params.species}
    Output:        ${params.outdir}
    Metadata:      ${params.metadata_file}
    Cell type key: ${cell_type_key}

    On-ramps (wired in src; D2.1 back-port from PBMC instance is deferred):
      rna_integrated_h5ad:    ${params.onramp?.rna_integrated_h5ad ?: 'none'}
      atac_peak_matrix_h5ad:  ${params.onramp?.atac_peak_matrix_h5ad ?: 'none'}
      mudata_h5mu:            ${params.onramp?.mudata_h5mu ?: 'none'}

    """

    // ========================================================================
    // RNA PROCESSING (with granular on-ramps)
    // 2026-04-30 refactor (D2.1 back-port 2026-05-04): expose ch_rna_qc_h5ads +
    // ch_rna_cellchat_csv as channel variables so downstream consumers
    // (MULTIOME_INTEGRATION, ENHANCER_FOOTPRINTING_RECIPES) can pull from
    // either RNA.out.* or onramp.
    // ========================================================================
    def ch_integrated_rna
    def ch_rna_qc_h5ads
    def ch_rna_cellchat_csv

    if (params.onramp?.rna_integrated_h5ad) {
        log.info "ON-RAMP: Pre-integrated RNA: ${params.onramp.rna_integrated_h5ad}"
        ch_integrated_rna = Channel.value(file(params.onramp.rna_integrated_h5ad))

        // Per-sample QC'd h5ads — required by MULTIOME_INTEGRATION (sample-level join with ATAC)
        if (params.onramp?.rna_per_sample_h5ads_dir) {
            log.info "ON-RAMP: RNA per-sample h5ads from: ${params.onramp.rna_per_sample_h5ads_dir}"
            // Emit (sample_id, file) tuples matching RNA.out.qc_h5ads shape.
            // Sample id is inferred from basename, stripping common QC suffixes.
            ch_rna_qc_h5ads = Channel.fromPath("${params.onramp.rna_per_sample_h5ads_dir}/*.h5ad")
                .map { f ->
                    def sid = f.baseName.replaceAll(/_filtered_All$/, '')
                                        .replaceAll(/_filtered$/, '')
                    tuple(sid, f)
                }
        } else {
            log.warn "RNA integrated onramp set without rna_per_sample_h5ads_dir — MULTIOME_INTEGRATION will be unable to build MuData."
            ch_rna_qc_h5ads = Channel.empty()
        }

        // CellChat CSV — optional, only consumed by ENHANCER_FOOTPRINTING_RECIPES Phase 3
        if (params.onramp?.rna_cellchat_csv) {
            ch_rna_cellchat_csv = Channel.value(file(params.onramp.rna_cellchat_csv))
        } else {
            ch_rna_cellchat_csv = Channel.empty()
        }

        rna_completed = true
        rna_from_onramp = true
    } else if (params.rna.run) {
        log.info "Starting RNA workflow..."
        RNA()
        ch_integrated_rna   = RNA.out.integrated_rna
        ch_rna_qc_h5ads     = RNA.out.qc_h5ads
        ch_rna_cellchat_csv = RNA.out.cellchat_csv
        rna_completed = true
        rna_from_onramp = false
    } else {
        log.info "Skipping RNA workflow (disabled in config)"
        ch_integrated_rna   = Channel.empty()
        ch_rna_qc_h5ads     = Channel.empty()
        ch_rna_cellchat_csv = Channel.empty()
        rna_completed = false
        rna_from_onramp = false
    }

    // ========================================================================
    // RNA DIFFERENTIAL EXPRESSION
    // ========================================================================
    if (params.differential_rna.run && rna_completed) {
        log.info """
        RNA DIFFERENTIAL EXPRESSION ENABLED
        Comparisons: ${params.differential_rna.comparisons.size()}
        """

        RNA_DIFFERENTIAL(
            ch_integrated_rna
        )
    }

    // ========================================================================
    // ATAC PROCESSING (with granular on-ramps)
    // 2026-04-30 refactor (D2.1 back-port 2026-05-04): onramp source produces
    // the same channels as ATAC_INITIAL/ATAC_FINAL so REGULATORY_ANALYSIS,
    // ENHANCER_FOOTPRINTING, and downstream consumers run identically
    // regardless of source. Each major artifact is independently onrampable.
    // ========================================================================
    def ch_atac_peak_matrix
    def ch_atac_individual_samples
    def ch_atac_anndataset
    def da_peaks_ch

    if (params.onramp?.atac_peak_matrix_h5ad) {
        log.info "ON-RAMP: ATAC peak matrix: ${params.onramp.atac_peak_matrix_h5ad}"
        ch_atac_peak_matrix = Channel.value(file(params.onramp.atac_peak_matrix_h5ad))

        // Per-sample h5ads — required by SCPRINTER (per-sample fragment binding)
        if (params.onramp?.atac_individual_samples_dir) {
            log.info "ON-RAMP: ATAC per-sample h5ads from: ${params.onramp.atac_individual_samples_dir}"
            ch_atac_individual_samples = Channel.fromPath("${params.onramp.atac_individual_samples_dir}/*.h5ad").collect()
        } else {
            log.warn "ATAC peak_matrix onramp set without atac_individual_samples_dir — scPRINTer will be unable to build per-sample fragment binding."
            ch_atac_individual_samples = Channel.empty()
        }

        // AnnDataSet (snapatac2 binding object) — required by ENHANCER_FOOTPRINTING_RECIPES
        if (params.onramp?.atac_anndataset) {
            log.info "ON-RAMP: ATAC anndataset: ${params.onramp.atac_anndataset}"
            ch_atac_anndataset = Channel.value(file(params.onramp.atac_anndataset))
        } else {
            ch_atac_anndataset = Channel.empty()
        }

        atac_completed = true
        atac_from_onramp = true
        da_peaks_ch = Channel.empty()  // descriptive + DA only valid in non-onramp ATAC mode
    } else if (params.atac.run && params.atac.run_initial_qc) {
        log.info "Starting ATAC workflow..."

        // Stage 1: Initial QC
        ATAC_INITIAL()

        // Stage 2: Final QC with sample-specific thresholds
        ATAC_FINAL(ATAC_INITIAL.out.thresholds)

        ch_atac_peak_matrix         = ATAC_FINAL.out.peak_matrix
        ch_atac_individual_samples  = ATAC_FINAL.out.individual_samples
        ch_atac_anndataset          = ATAC_FINAL.out.anndataset
        atac_completed = true
        atac_from_onramp = false

        // ====================================================================
        // ATAC Descriptive Report (always-on after ATAC_FINAL)
        // D1b: SDas_nf line 214 — unconditional cell-count summaries
        // ====================================================================
        ATAC_DESCRIPTIVE_REPORT(
            ATAC_FINAL.out.peak_matrix,
            file(params.atac.sample_metadata)
        )

        // ====================================================================
        // ATAC Post-hoc QC Report (opt-in via params.qc.post_hoc_report)
        // D1b: re-renders FORGE-style TSSE marginal scatter from per-sample
        // ATAC_INITIAL_QC outputs. Only valid in non-onramp ATAC mode.
        // ====================================================================
        if (params.qc?.post_hoc_report ?: false) {
            POST_QC_REPORT(
                ATAC_INITIAL.out.individual_samples.collect(),
                params.atac.min_tsse        ?: params.atac.initial_min_tsse   ?: 6,
                params.atac.min_counts      ?: params.atac.initial_min_counts ?: 5000,
                params.atac.max_counts      ?: params.atac.initial_max_counts ?: 100000
            )
        }

        // Differential ATAC Analysis (requires condition_key)
        if (params.differential.run && params.differential.condition_key) {
            log.info """
            DIFFERENTIAL ATAC ANALYSIS ENABLED
            Condition key: ${params.differential.condition_key}
            Control:   ${params.differential.control_condition}
            Treatment: ${params.differential.treatment_condition}
            """

            ATAC_DIFFERENTIAL(
                ATAC_FINAL.out.peak_matrix,
                file(params.atac.sample_metadata)
            )

            da_peaks_ch = ATAC_DIFFERENTIAL.out.da_peaks
        } else {
            log.info "Skipping differential ATAC (no condition_key in metadata)"
            da_peaks_ch = Channel.empty()
        }
    } else {
        log.info "Skipping ATAC workflow (disabled in config)"
        atac_completed = false
        atac_from_onramp = false
        ch_atac_peak_matrix         = Channel.empty()
        ch_atac_individual_samples  = Channel.empty()
        ch_atac_anndataset          = Channel.empty()
        da_peaks_ch = Channel.empty()
    }

    // ========================================================================
    // REGULATORY ANALYSIS (Cicero, ChromVAR, scPRINTer)
    // 2026-04-30: Hoisted out of ATAC elif so it runs from onramp branch too.
    // Per-process onramps (cicero/chromvar/printer) are gated INSIDE the
    // workflow on params.onramp.{cicero_*, chromvar_*, printer_h5ad}.
    // ========================================================================
    if (atac_completed && params.atac.run) {
        log.info "Starting regulatory analysis workflow..."

        // Parse fragment files from manifest for SCPRINTER_BUILD_PRINTER
        ch_reg_fragments = Channel.fromPath(params.metadata_file)
            .splitCsv(header: true)
            .filter { isNonEmptyRow(it) }  // FIX-A6
            .map { trimRow(it) }           // FIX-A5
            .map { row ->
                def atac_dir = resolveAtacDir(row)
                def frag_fname = row.fragment_file.contains('.') ? row.fragment_file : "${row.fragment_file}.bed.gz"
                file("${atac_dir}/${frag_fname}")
            }
            .collect()

        REGULATORY_ANALYSIS(
            ch_atac_peak_matrix,
            ch_atac_individual_samples,
            da_peaks_ch,
            file(params.atac.sample_metadata),
            ch_reg_fragments,
            atac_cell_type_key
        )
    }

    // ========================================================================
    // MULTIOME INTEGRATION (requires both RNA and ATAC, OR mudata onramp)
    // 2026-04-30 refactor (D2.1 back-port 2026-05-04): gate uses can_build_mudata,
    // allowing onramped RNA per-sample h5ads + onramped ATAC per-sample h5ads
    // to be joined just like fresh outputs.
    // ========================================================================
    def can_build_mudata = rna_completed && atac_completed &&
        (!rna_from_onramp || params.onramp?.rna_per_sample_h5ads_dir) &&
        (!atac_from_onramp || params.onramp?.atac_individual_samples_dir)
    def mudata_from_onramp = false

    if (params.onramp?.mudata_h5mu) {
        log.info "ON-RAMP: Using pre-computed MuData: ${params.onramp.mudata_h5mu}"
        mudata_from_onramp = true
        // Skip MULTIOME_INTEGRATION entirely -- downstream GRN workflows
        // should read from the on-ramp mudata directly.
    } else if (params.run_multiome_integration && can_build_mudata) {
        log.info """
        MULTIOME INTEGRATION
        Waiting for RNA and ATAC workflows to complete...
        """

        // Explicit dependency: wait for ATAC samples to materialize (fresh or onramp)
        ch_atac_individual_samples
            .collect()
            .subscribe { atac_files ->
                log.info "ATAC samples ready: ${atac_files.size()} files"
            }

        // ========================================================================
        // FAIL-FAST MODALITY JOINING (RNA + ATAC)
        // BD multiome join diagnostics: logs all channel keys, warns on orphans,
        // reports ALL mismatches before erroring.
        // ========================================================================
        ch_expected_samples = Channel.fromPath(params.metadata_file)
            .splitCsv(header: true)
            .filter { isNonEmptyRow(it) }  // FIX-A6
            .map { trimRow(it) }           // FIX-A5
            .map { row ->
                def condition = row.condition_group?.trim()
                if (!condition) {
                    log.warn "Sample '${row.sample_id}' has no condition_group — defaulting to 'Control'"
                    condition = 'Control'
                }
                def meta = [id: row.sample_id, batch: row.batch, condition: condition]
                tuple(row.sample_id, meta)
            }

        // Map flat RNA outputs back to tuples with keys
        // For single-file (non-demux) lanes, use lane_sample_id directly.
        // For demux (BD brain): reconstruct donor + batch from filename.
        ch_rna_mapped = ch_rna_qc_h5ads
            .flatMap { sample_id, files ->
                def file_list = files instanceof List ? files : [files]
                if (file_list.size() == 1) {
                    return [tuple(sample_id, file_list[0])]
                }
                file_list.collect { f ->
                    def donor = f.baseName.replaceAll(/^.*_filtered_/, '')
                    def batch = sample_id.tokenize('_').last()
                    tuple("${donor}_${batch}", f)
                }
            }

        // Map flat ATAC outputs back to tuples with keys
        // FIX-45b: Filter out auxiliary h5ads (peak_matrix, gene_matrix, atac_complete,
        // atac_celltypist_annotations, peak_matrix_annotated, cluster_avg_gene_scores)
        // that are not per-sample files.
        def atac_auxiliary = ['peak_matrix', 'gene_matrix', 'atac_complete',
                              'atac_celltypist_annotations', 'scatanno_annotations',
                              'peak_matrix_annotated',
                              'cluster_avg_gene_scores'] as Set
        ch_atac_mapped = ch_atac_individual_samples
            .flatMap { it instanceof List ? it : [it] }
            .filter { file -> !atac_auxiliary.contains(file.baseName) }
            .map { file ->
                def sample_id = file.baseName
                tuple(sample_id, file)
            }

        // Log available keys from each channel for debugging
        ch_rna_mapped.map { it[0] }.collect().subscribe { keys ->
            log.info "MULTIOME JOIN: RNA keys (${keys.size()}): ${keys.sort().take(10)}${keys.size() > 10 ? '... (+' + (keys.size()-10) + ' more)' : ''}"
        }
        ch_atac_mapped.map { it[0] }.collect().subscribe { keys ->
            log.info "MULTIOME JOIN: ATAC keys (${keys.size()}): ${keys.sort().take(10)}${keys.size() > 10 ? '... (+' + (keys.size()-10) + ' more)' : ''}"
        }
        ch_expected_samples.map { it[0] }.collect().subscribe { keys ->
            log.info "MULTIOME JOIN: Expected keys (${keys.size()}): ${keys.sort().take(10)}${keys.size() > 10 ? '... (+' + (keys.size()-10) + ' more)' : ''}"
        }

        // Join with remainder to detect mismatches, then validate
        ch_joined_raw = ch_expected_samples
            .join(ch_rna_mapped, by: 0, remainder: true)
            .join(ch_atac_mapped, by: 0, remainder: true)

        // Separate orphans (keys not in manifest) from expected samples
        ch_joined_raw
            .filter { it[1] == null }
            .map { it[0] }
            .collect()
            .subscribe { orphans ->
                if (orphans) {
                    log.warn "MULTIOME JOIN: ${orphans.size()} orphan key(s) found in RNA/ATAC but not in manifest: ${orphans.sort()}"
                    log.warn "These samples will be EXCLUDED from multiome integration. If unexpected, check manifest sample_id values."
                }
            }

        // Collect all validation errors for expected samples, then report
        ch_joined_raw
            .filter { it[1] != null }
            .collect(flat: false)
            .subscribe { items ->
                if (items == null || items.isEmpty()) {
                    log.warn "MULTIOME JOIN: No expected samples found in join -- check manifest and channel emissions."
                    return
                }
                def missing_rna = items.findAll { it.size() < 3 || it[2] == null }.collect { it[0] }
                def missing_atac = items.findAll { it.size() < 4 || it[3] == null }.collect { it[0] }
                if (missing_rna) {
                    log.error "MULTIOME JOIN: ${missing_rna.size()} sample(s) missing RNA: ${missing_rna.sort()}"
                }
                if (missing_atac) {
                    log.error "MULTIOME JOIN: ${missing_atac.size()} sample(s) missing ATAC: ${missing_atac.sort()}"
                }
            }

        ch_paired_multiome = ch_joined_raw
            .filter { it[1] != null }
            .map { it ->
                def sample_id = it[0]
                def meta = it[1]
                def rna = it[2]
                def atac = it[3]
                if (rna == null) error "Fail-fast: Missing RNA for sample ${sample_id}. Check RNA_QC output filenames vs manifest sample_id."
                if (atac == null) error "Fail-fast: Missing ATAC for sample ${sample_id}. Check ATAC h5ad basenames vs manifest sample_id."

                return [meta, rna, atac]
            }

        ch_synced_rna = ch_paired_multiome.map { meta, rna, atac -> rna }.collect()
        ch_synced_atac = ch_paired_multiome.map { meta, rna, atac -> atac }.collect()

        // Generate sample pairing map so BUILD_MUDATA knows which files go together
        ch_sample_map_file = ch_paired_multiome
            .map { meta, rna, atac -> "${meta.id},${rna.name},${atac.name}" }
            .collectFile(name: 'sample_map.csv', seed: 'sample_id,rna_file,atac_file', newLine: true)

        MULTIOME_INTEGRATION(
            ch_synced_rna,
            ch_synced_atac,
            ch_integrated_rna,
            file(params.metadata_file),
            ch_sample_map_file,
            ch_atac_peak_matrix
        )
    } else if (params.run_multiome_integration && (!rna_completed || !atac_completed)) {
        log.warn """
        MULTIOME INTEGRATION SKIPPED
        Reason: Both RNA and ATAC workflows must complete successfully
        RNA completed: ${rna_completed}
        ATAC completed: ${atac_completed}
        """
    }

    // ========================================================================
    // MULTIOME GRN (pycistopic + SCENIC+ + DORC)
    // 2026-04-30 (D2.1 back-port 2026-05-04): can_run_multiome_grn replaces
    // ad-hoc onramp gates; MULTIOME_INTEGRATION.out is only available when
    // multiome was actually built (not when mudata_from_onramp).
    // ========================================================================
    def can_run_multiome_grn = can_build_mudata && !mudata_from_onramp &&
                               params.run_multiome_integration &&
                               (params.pycistopic.run || params.scenicplus.run || params.dorc.run)

    if (can_run_multiome_grn) {
        log.info """
        MULTIOME GRN ANALYSIS
        Running pycisTopic, SCENIC+ and/or DORC as requested in config
        """
        MULTIOME_GRN(
            ch_integrated_rna,
            ch_atac_peak_matrix,
            file(params.metadata_file),
            MULTIOME_INTEGRATION.out.rna_for_dorc,
            MULTIOME_INTEGRATION.out.stats,
            file(params.pycistopic.blacklist_bed),
            cell_type_key
        )
    }

    // ========================================================================
    // ENHANCER FOOTPRINTING RECIPES (A/B/C/D)
    // Requires: ATAC completed + REGULATORY_ANALYSIS (Cicero, ChromVAR, scPRINTer).
    // 2026-04-30 (D2.1 back-port 2026-05-04): atac_from_onramp gate removed —
    // EFR now runs from onramp branches when atac side keys are populated.
    // ========================================================================
    if (params.enhancer_footprinting.run && atac_completed &&
        params.cicero.run && params.chromvar.run && params.scprinter.run) {
        log.info """
        ENHANCER FOOTPRINTING RECIPES
        Phase 1 (ATAC-only):     always
        Phase 2 (Multiome):      ${(params.scenicplus.run ?: false) ? 'enabled' : 'disabled'}
        Phase 3 (CellChat):      ${(params.enhancer_recipe_c.run ?: false) ? 'enabled' : 'disabled'}
        Phase 4 (Visualization): ${(params.enhancer_viz.run ?: false) ? 'enabled' : 'disabled'}
        """

        def has_grn = can_run_multiome_grn

        def ereg_direct_ch = (has_grn && params.scenicplus.run && params.pycistopic.run) ?
            MULTIOME_GRN.out.scplus_ereg_direct : Channel.empty()
        def r2g_ch = (has_grn && params.scenicplus.run && params.pycistopic.run) ?
            MULTIOME_GRN.out.scplus_r2g : Channel.empty()
        def dorc_sig_ch = (has_grn && params.dorc.run) ?
            MULTIOME_GRN.out.dorc_sig : Channel.empty()
        def rna_ch = rna_completed ? ch_integrated_rna : Channel.empty()
        // CellChat CSV: ch_rna_cellchat_csv is populated when RNA ran fresh
        // OR when params.onramp.rna_cellchat_csv was supplied; otherwise empty.
        def cellchat_csv_ch = (rna_completed && params.cellchat.run) ?
            ch_rna_cellchat_csv : Channel.empty()
        def bigwigs_ch = (has_grn && params.pycistopic.run) ?
            MULTIOME_GRN.out.pseudobulk_bigwigs : Channel.empty()

        ENHANCER_FOOTPRINTING_RECIPES(
            ch_atac_peak_matrix,
            REGULATORY_ANALYSIS.out.scprinter_printer,
            REGULATORY_ANALYSIS.out.chromvar_deviations,
            REGULATORY_ANALYSIS.out.chromvar_per_ct,
            REGULATORY_ANALYSIS.out.cicero_connections,
            REGULATORY_ANALYSIS.out.cicero_ccan,
            ereg_direct_ch,
            r2g_ch,
            dorc_sig_ch,
            rna_ch,
            cellchat_csv_ch,
            bigwigs_ch,
            atac_cell_type_key,
            REGULATORY_ANALYSIS.out.tf_targets,
            ch_atac_anndataset,                            // D1b
            REGULATORY_ANALYSIS.out.scprinter_footprints,  // 2026-04-30: pass-through (no cross-workflow access)
            REGULATORY_ANALYSIS.out.tf_diff,               // 2026-05-04: pass-through (DSL2 scope fix; gate-aligned w/ differential_tf.run)
            REGULATORY_ANALYSIS.out.gene_coordinates,      // 2026-05-06: pass-through for PROMOTER_MSFP_PER_CT
            REGULATORY_ANALYSIS.out.cicero_strat_ctrl_links,        // 2026-05-06: cis-rewiring
            REGULATORY_ANALYSIS.out.cicero_strat_trt_links,
            REGULATORY_ANALYSIS.out.cicero_strat_ctrl_peaks,
            REGULATORY_ANALYSIS.out.cicero_strat_trt_peaks,
            REGULATORY_ANALYSIS.out.cicero_strat_ctrl_connections,   // 2026-06-01: lollipop fix
            REGULATORY_ANALYSIS.out.cicero_strat_trt_connections
        )
    }

    // ========================================================================
    // PER-CT × CONDITION CICERO
    // Single-pass: BUILD_CT_ANNOTATION runs after peak_matrix_annotated.h5ad
    // is ready, then fans out to (CT, condition) strata in the same run.
    // Onramp override: set params.cicero_per_ct.annotation_csv to skip
    // BUILD_CT_ANNOTATION and use a pre-generated CSV directly.
    // Opt-in: params.cicero_per_ct.enabled = true.
    // Per-stratum floor (250 cells) applied in CICERO_TRIPLETS_PER_CT via
    // exit 77 → errorStrategy ignore → clean skip.
    // ========================================================================
    if (params.cicero_per_ct?.enabled == true && atac_completed) {
        def csv_ch
        if (params.cicero_per_ct?.annotation_csv) {
            csv_ch = Channel.value(file(params.cicero_per_ct.annotation_csv))
        } else {
            BUILD_CT_ANNOTATION(
                ch_atac_peak_matrix.first(),
                params.cicero_per_ct?.cell_type_col ?: 'cell_type_broad',
                params.cicero_per_ct?.condition_col  ?: 'condition',
                params.cicero_per_ct?.min_cells       ?: 50,
                params.cicero_per_ct?.min_pct         ?: 0.01
            )
            csv_ch = BUILD_CT_ANNOTATION.out.csv
        }

        def strata_ch = csv_ch
            .splitCsv(header: true)
            .map  { row -> tuple(row.cell_type_v2, row.condition) }
            .filter { ct, cond -> ct != 'EXCLUDED' }
            .unique()

        CICERO_TRIPLETS_PER_CT(
            ch_atac_peak_matrix.first(),
            csv_ch.first(),
            strata_ch
        )
        // Phase 2a: estimate distance parameter per stratum
        CICERO_ESTIMATE_DP_PER_CT(CICERO_TRIPLETS_PER_CT.out.triplets)

        // Phase 2b: fan-out by chromosome within each stratum
        def cicero_chroms_per_ct = (params.species == 'mouse' ? (1..19) : (1..22))
                                       .collect { "chr${it}" } + ['chrX', 'chrY', 'chrM']

        def per_ct_chrom_in = CICERO_ESTIMATE_DP_PER_CT.out.all_out
            .combine(Channel.fromList(cicero_chroms_per_ct))
            .map { ct, cond, dp, cds, gene_ann, ordered_cds, chrom ->
                tuple(ct, cond, chrom, cds, gene_ann, dp)
            }

        CICERO_FULL_CHROM_PER_CT(per_ct_chrom_in)

        // Phase 2c: collect per-stratum chrom connections, join with ordered CDS
        def per_ct_chrom_conns = CICERO_FULL_CHROM_PER_CT.out.chrom_conns
            .groupTuple(by: [0, 1], size: cicero_chroms_per_ct.size())
            .map { ct, cond, _chroms, files -> tuple(ct, cond, files) }

        def per_ct_ordered_cds = CICERO_ESTIMATE_DP_PER_CT.out.all_out
            .map { ct, cond, dp, cds, gene_ann, ordered_cds -> tuple(ct, cond, ordered_cds) }

        def per_ct_join_in = per_ct_chrom_conns.join(per_ct_ordered_cds, by: [0, 1])

        CICERO_JOIN_PER_CT(per_ct_join_in, Channel.value(params.cicero.gtf_full))
    }

    // SHI_FIGURES — Shi et al. 2025 figure equivalents (Tier A always; Tier B
    // gates on existence of differential outputs, so single-condition runs
    // produce 1E/2B/2D-foundation only).
    // 2026-05-04 refactor: SHI now consumes upstream channels (peak_matrix
    // from ATAC; bigwig manifest+dir from EFR), so Tier A fires on first run
    // without requiring a -resume cycle.
    if (params.shi_figures?.enabled == true) {
        def shi_bw_manifest_ch = Channel.empty()
        def shi_bw_dir_ch      = Channel.empty()
        if (params.enhancer_footprinting.run && atac_completed &&
            params.cicero.run && params.chromvar.run && params.scprinter.run) {
            shi_bw_manifest_ch = ENHANCER_FOOTPRINTING_RECIPES.out.shi_bw_manifest
            shi_bw_dir_ch      = ENHANCER_FOOTPRINTING_RECIPES.out.shi_bw_dir
        }
        SHI_FIGURES(ch_atac_peak_matrix, shi_bw_manifest_ch, shi_bw_dir_ch, atac_broad_cell_type_key)
    }
}

// ============================================================================
// VIZ_ONLY — read-only viz entry workflow
// D2: absorbed from SDas_nf. Re-renders post-hoc QC + Cicero target-gene plots
// from already-persisted artifacts under results/. Does NOT re-run upstream
// compute. Invoke via:
//
//   nextflow run main.nf -entry VIZ_ONLY \
//     --viz_only.peak_matrix_h5ad results/atac/final/peak_matrix_annotated.h5ad \
//     --viz_only.cicero_connections results/cicero/.../cicero_connections.tsv.gz \
//     --viz_only.cicero_ccan        results/cicero/.../CCAN_assignments.tsv.gz \
//     --viz_only.cicero_cds         results/cicero/.../input_cds_ordered.rds \
//     --viz_only.target_genes 'Hmgb1,Hspa8,Atf6'
// ============================================================================
workflow VIZ_ONLY {
    main:

    if (params.viz_only?.peak_matrix_h5ad) {
        POST_QC_REPORT(
            file(params.viz_only.peak_matrix_h5ad),
            params.atac.min_tsse        ?: params.atac.initial_min_tsse   ?: 6,
            params.atac.min_counts      ?: params.atac.initial_min_counts ?: 5000,
            params.atac.max_counts      ?: params.atac.initial_max_counts ?: 100000
        )
    } else {
        log.warn "VIZ_ONLY: --viz_only.peak_matrix_h5ad not provided; skipping POST_QC_REPORT."
    }

    if (params.viz_only?.cicero_connections &&
        params.viz_only?.cicero_ccan        &&
        params.viz_only?.cicero_cds) {

        def viz_genes
        if (params.viz_only?.target_genes) {
            viz_genes = (params.viz_only.target_genes instanceof List) ?
                params.viz_only.target_genes :
                params.viz_only.target_genes.toString().split(',').collect { it.trim() }.findAll { it }
        } else {
            viz_genes = params.cicero.target_genes ?: []
        }

        if (!viz_genes) {
            log.warn "VIZ_ONLY: no target genes provided (--viz_only.target_genes or params.cicero.target_genes); skipping CICERO_TARGET_PLOTS."
        } else {
            CICERO_TARGET_PLOTS(
                file(params.viz_only.cicero_connections),
                file(params.viz_only.cicero_ccan),
                file(params.viz_only.cicero_cds),
                params.cicero.gtf_plot,
                viz_genes
            )
        }
    } else {
        log.warn "VIZ_ONLY: cicero_connections / cicero_ccan / cicero_cds not all provided; skipping CICERO_TARGET_PLOTS."
    }
}

// ============================================================================
// SHI_FIGURES — Shi et al. 2025 figure-equivalents from persisted artifacts.
// D3: absorbed from SDas_nf 2026-05-03. Tier A modules (ANNOTATE_PEAK_TYPES,
// MARKER_COVERAGE_TRACKS) run unconditionally; Tier B (NMF_ENHANCER_PROGRAMS,
// DA/TF differential plots, curated networks, locus bars) auto-skip on
// single-condition runs.
// ENG-24 (2026-05-04): NMF moved to Tier B because nmf_enhancer_programs.py
//   requires obs[condition_key] which is only guaranteed under >=2 conditions.
//
// Path defaults assume the FORGE outdir layout (different from SDas_nf):
//   peak_matrix     — ${outdir}/atac/final/peak_matrix.h5ad
//   diff DA peaks   — ${outdir}/differential/DA_peaks_*.csv
//   diff TF csvs    — ${outdir}/differential_tf/tf_differential_*.csv
//   ccan gene links — ${outdir}/enhancer_footprinting/ccan_enhancers/
//   tf adjacency    — ${outdir}/enhancer_footprinting/network/tf_gene_adjacency.tsv
//   bigwig manifest — ${outdir}/enhancer_viz/bigwigs/bigwigs/manifest.json
// Override any with --shi_figures.<key> on the CLI.
// ============================================================================
workflow SHI_FIGURES {

    take:
        peak_matrix_ch        // upstream-emitter channel (ATAC consolidated/onramp).
                              //   Replaces ${outdir}/atac/final/peak_matrix.h5ad
                              //   path probe so Tier A fires deterministically on
                              //   first run instead of waiting for a -resume cycle.
        bw_manifest_ch        // EXPORT_ATAC_BIGWIGS.out.manifest, or Channel.empty()
                              //   when the scenicplus / no-enhancer-viz path fires.
        bw_dir_ch             // EXPORT_ATAC_BIGWIGS.out.bigwigs (queue), or empty.
        cell_type_col         // broad ATAC cell-type obs column on peak_matrix
                              //   (e.g. 'cell_type_broad'); flows in from
                              //   atac_broad_cell_type_key at the call site so
                              //   downstream NMF doesn't read a flat param.

    main:

    // FIX ENG-13/MIS-08: use canonical params (gtf_human_full / gtf_mouse_full)
    // matching the resolution used elsewhere (main.nf:574, 643). Fall back to
    // params.cicero.gtf_full if set. Hard-fail rather than pass null downstream.
    def gtf = (params.species == 'human' ? params.gtf_human_full : params.gtf_mouse_full) ?:
              params.cicero?.gtf_full
    if (!gtf || !file(gtf).exists()) {
        log.warn "SHI_FIGURES: GTF not found (params.gtf_${params.species}_full=${gtf}); skipping all GTF-dependent panels (ANNOTATE_PEAK_TYPES, MARKER_COVERAGE_TRACKS)."
        gtf = null
    }
    def trt_label  = params.shi_figures?.treatment ?: ''
    def ctrl_label = params.shi_figures?.control   ?: ''

    // Tier B file-glob paths still resolve at construction time (ENG-02:
    // intentionally allows consumption of differential CSVs from prior runs).
    def conn_ctrl_path = params.shi_figures?.connections_ctrl ?:
        "${params.outdir}/cicero/stratified/${ctrl_label}/cicero_connections.tsv.gz"
    def conn_trt_path  = params.shi_figures?.connections_trt ?:
        "${params.outdir}/cicero/stratified/${trt_label}/cicero_connections.tsv.gz"
    def gene_links_ctrl_path = params.shi_figures?.gene_links_ctrl ?:
        "${params.outdir}/enhancer_footprinting/ccan_enhancers/ccan_enhancer_gene_links_${ctrl_label}.tsv"
    def gene_links_trt_path  = params.shi_figures?.gene_links_trt ?:
        "${params.outdir}/enhancer_footprinting/ccan_enhancers/ccan_enhancer_gene_links_${trt_label}.tsv"
    def gene_links_global    = params.shi_figures?.gene_links_global ?:
        "${params.outdir}/enhancer_footprinting/ccan_enhancers/ccan_enhancer_gene_links.tsv"
    def enh_bed_path  = params.shi_figures?.enhancer_bed ?:
        "${params.outdir}/enhancer_footprinting/ccan_enhancers/ccan_enhancer_peaks.bed.gz"
    def adjacency_path = params.shi_figures?.adjacency ?:
        "${params.outdir}/enhancer_footprinting/network/tf_gene_adjacency.tsv"
    def have_two_cond = (trt_label && ctrl_label)

    // ----- Tier A — single-condition compatible (channel-driven 2026-05-04) -----

    // Foundation: peak biotype classification. peak_matrix_ch is empty when
    // ATAC is disabled or the upstream emitter never fired, so the process
    // naturally skips (no construction-time path probe).
    if (gtf) {
        ANNOTATE_PEAK_TYPES(peak_matrix_ch, file(gtf))
    } else {
        log.warn "SHI_FIGURES: GTF unresolved; skipping ANNOTATE_PEAK_TYPES."
    }

    // 1E — broad-class bigWigs + marker coverage tracks. bw_manifest_ch is
    // empty when EXPORT_ATAC_BIGWIGS didn't fire (scenicplus pseudobulk path
    // or enhancer_viz disabled), so MARKER_COVERAGE_TRACKS naturally skips.
    if (gtf) {
        MARKER_COVERAGE_TRACKS(bw_manifest_ch, bw_dir_ch.collect(), file(gtf))
    }

    // ----- Tier B — require >=2 conditions; auto-skip otherwise -----

    if (!have_two_cond) {
        log.info "SHI_FIGURES: shi_figures.treatment / .control not both set; skipping " +
                 "Tier B (single-condition mode), including NMF_ENHANCER_PROGRAMS."
        return
    }

    // ENG-19 (2026-05-04): tighten globs by current trt/ctrl labels so stale
    // CSVs from a prior comparison label can't be silently consumed. Filenames
    // emitted by snapatac_diff.nf are DA_peaks_<ct>__<trt>_vs_<ctrl>.csv and
    // by differential_tf_accessibility.py are tf_differential_<ct>_<trt>_vs_<ctrl>.csv.
    def diff_glob_strict    = "${params.outdir}/differential/DA_peaks_*__${trt_label}_vs_${ctrl_label}.csv"
    def diff_tf_glob_strict = "${params.outdir}/differential_tf/tf_differential_*_${trt_label}_vs_${ctrl_label}.csv"
    def have_diff_da_strict = !file(diff_glob_strict).isEmpty()
    def have_diff_tf_strict = !file(diff_tf_glob_strict).isEmpty()

    // ENG-02 (2026-05-04): when params.differential.run=false but Tier B is
    // consuming pre-existing CSVs, log each filename + mtime so the user can
    // see what's being re-rendered against. Does not block — re-rendering from
    // a prior differential run is a legitimate use case.
    if (!params.differential?.run && (have_diff_da_strict || have_diff_tf_strict)) {
        def staleFiles = []
        if (have_diff_da_strict) staleFiles.addAll(file(diff_glob_strict))
        if (have_diff_tf_strict) staleFiles.addAll(file(diff_tf_glob_strict))
        log.warn "SHI_FIGURES: differential.run=false but Tier B is consuming " +
                 "${staleFiles.size()} pre-existing CSV(s) from a prior run:"
        def mtfmt = new java.text.SimpleDateFormat('yyyy-MM-dd HH:mm:ss')
        staleFiles.each { f ->
            log.warn "  - ${f.name} (mtime: ${mtfmt.format(new Date(f.lastModified()))})"
        }
        log.warn "  Verify these reflect your current parameters; otherwise re-run differential."
    }

    Channel.fromPath(diff_tf_glob_strict).collect().set { ch_diff_tf }
    Channel.fromPath(diff_glob_strict).collect().set    { ch_diff }

    // 2B — NMF on enhancer × pseudobulk. peak_matrix_ch upstream-driven; enhancer
    // BED still file-based since it can come from a prior EFR run. cell_type_col
    // flows in via take: so the obs column name is decided once at the top-level
    // call site (atac_broad_cell_type_key) instead of in a per-process params lookup.
    if (file(enh_bed_path).exists()) {
        NMF_ENHANCER_PROGRAMS(peak_matrix_ch, file(enh_bed_path), cell_type_col)
    } else {
        log.warn "SHI_FIGURES: enhancer BED not found at ${enh_bed_path}; skipping NMF_ENHANCER_PROGRAMS (2B)."
    }

    // 2D — differential peak biotype breakdown (depends on DA peaks + peak annotation)
    if (have_diff_da_strict) {
        DA_PEAK_BREAKDOWN(ch_diff, ANNOTATE_PEAK_TYPES.out.tsv)
        // 2E — per-celltype log2FC heatmaps
        if (file(gene_links_global).exists()) {
            DA_LOG2FC_HEATMAPS(ch_diff, ANNOTATE_PEAK_TYPES.out.tsv, file(gene_links_global))
        } else {
            log.warn "SHI_FIGURES: gene links global not found at ${gene_links_global}; skipping DA_LOG2FC_HEATMAPS (2E)."
        }
    } else {
        log.info "SHI_FIGURES: no DA-peak CSVs at ${diff_glob_strict}; skipping 2D/2E."
    }

    // 2C — co-accessibility correlation (requires both stratified Cicero outputs)
    if (file(conn_ctrl_path).exists() && file(conn_trt_path).exists()) {
        COACC_CORRELATION_MATRIX(
            file(conn_ctrl_path),
            file(conn_trt_path),
            ANNOTATE_PEAK_TYPES.out.tsv,
        )
    } else {
        log.info "SHI_FIGURES: stratified Cicero connections not found for both conditions; skipping COACC_CORRELATION_MATRIX (2C)."
    }

    // 4B/5A/5C/5E — TF differential volcano (depends on diff_tf CSVs)
    if (have_diff_tf_strict) {
        TF_DIFFERENTIAL_VOLCANO(ch_diff_tf)
    } else {
        log.info "SHI_FIGURES: no TF-diff CSVs at ${diff_tf_glob_strict}; skipping TF_DIFFERENTIAL_VOLCANO (4B/5A/5C/5E)."
    }

    // Curated panels — gate on Cicero stratified + ccan gene links + adjacency + diff_tf
    def have_curated_inputs = have_diff_tf_strict &&
        file(conn_ctrl_path).exists() && file(conn_trt_path).exists() &&
        file(gene_links_ctrl_path).exists() && file(gene_links_trt_path).exists() &&
        file(adjacency_path).exists()

    if (have_curated_inputs) {
        SELECT_SHI_CANDIDATES(
            ch_diff_tf,
            file(conn_ctrl_path),
            file(conn_trt_path),
            file(gene_links_ctrl_path),
            file(gene_links_trt_path),
            file(adjacency_path),
        )

        // 4C/5B/5D/5F — curated TF→target subgraphs
        CURATED_TF_NETWORKS(
            file(adjacency_path),
            SELECT_SHI_CANDIDATES.out.candidates,
            ch_diff_tf,
        )

        // 4E — locus-restricted TF binding bars
        LOCUS_TF_BINDING(
            file(adjacency_path),
            ch_diff_tf,
            SELECT_SHI_CANDIDATES.out.candidates,
        )
    } else {
        log.info "SHI_FIGURES: candidate-selector inputs incomplete; skipping SELECT_SHI_CANDIDATES + CURATED_TF_NETWORKS + LOCUS_TF_BINDING."
    }
}

// ============================================================================
// WORKFLOW COMPLETION HANDLER
// ============================================================================
workflow.onComplete {
    def status_text = workflow.success ? 'SUCCESS' : 'FAILED'

    log.info """

    PIPELINE EXECUTION SUMMARY

    Status:    ${status_text}
    Duration:  ${workflow.duration}
    Completed: ${workflow.complete}

    Results Directory: ${params.outdir}/

    Output Locations:

    RNA Results:
      - CellBender reports:  ${params.outdir}/cellbender/
      - QC h5ad files:       ${params.outdir}/rna_qc/
      - Integrated h5ad:     ${params.outdir}/integration/
      - CellChat results:    ${params.outdir}/cellchat/
      - hdWGCNA results:     ${params.outdir}/hdwgcna/
      - Differential expr:   ${params.outdir}/rna_differential/

    ATAC Results:
      - Initial QC:          ${params.outdir}/atac/initial_qc/
      - Final QC:            ${params.outdir}/atac/final/
      - Peak matrix:         ${params.outdir}/atac/final/peak_matrix.h5ad
      - Cell-type annot:     ${params.outdir}/atac/final/celltype_annotations.json
      - Differential peaks:  ${params.outdir}/differential/

    Regulatory Analysis:
      - Cicero connections:  ${params.outdir}/cicero/
      - ChromVAR results:    ${params.outdir}/chromvar/
      - scPRINT footprints:  ${params.outdir}/scprinter/

    Multiome Integration:
      - MuData object:       ${params.outdir}/multiome/mudata/
      - MOFA factors:        ${params.outdir}/multiome/mofa/
      - Bootstrap results:   ${params.outdir}/multiome/mofa_bootstrap/
      - MultiVI model:       ${params.outdir}/multiome/multivi/
      - MultiVI plots:       ${params.outdir}/multiome/multivi/visualizations/

    Enhancer Visualization:
      - Track configs:       ${params.outdir}/enhancer_viz/tracks/
      - Composite figures:   ${params.outdir}/enhancer_viz/composites/

    Logs:
      - Execution trace:     logs/nextflow/trace.txt
      - Timeline:            logs/nextflow/timeline.html
      - Report:              logs/nextflow/report.html

    """.stripIndent()

    if (!workflow.success) {
        log.error """

        TROUBLESHOOTING TIPS:

        1. Check error logs: logs/nextflow/trace.txt
        2. Resume failed run: nextflow run main.nf -resume
        3. Check individual process logs in work/ directory
        4. Verify input files exist and are readable

        """
    }
}

// ============================================================================
// WORKFLOW ERROR HANDLER
// ============================================================================
workflow.onError {
    log.error """

    PIPELINE ERROR

    Error message: ${workflow.errorMessage}
    Error report: ${workflow.errorReport}

    Check the logs directory for details:
      ${params.outdir}/logs/nextflow/

    """
}
