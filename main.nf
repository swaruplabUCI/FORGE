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
include { CICERO_FULL_CHROM  } from './modules/cicero/cicero_full_chrom'
include { CICERO_JOIN        } from './modules/cicero/cicero_join'
include { CICERO_TARGET_PLOTS } from './modules/cicero/cicero_target_plots'

// ChromVAR modules
include { GPU_CHROMVAR } from './modules/chromvar/gpu_chromvar'
include { VIS_CHROMVAR } from './modules/chromvar/vischromvar'
include { EXTRACT_CHROMVAR_MOTIFS } from './modules/chromvar/extract_motifs'
include { MAP_TF_TO_TARGET_GENES } from './modules/chromvar/map_tf_targets'

// scPRINT modules
include { SCPRINTER_BARCODES } from './modules/scprint/barcodes'
include { RESOLVE_GENE_COORDINATES } from './modules/scprint/resolve_coordinates'
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
include { PYCISTOPIC_PREPARE } from './modules/multiome/pycistopic_prepare'
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


// ============================================================================
// GLOBAL: Compute cell_type_key ONCE
// ============================================================================
def has_reference = (params.species == 'human' && params.ref_dir_human_integrated) ||
                    (params.species == 'mouse' && params.ref_dir_mouse_integrated)
// Unified canonical cell-type column. BUILD_MUDATA and PLOT_POST_SCANVI write
// 'cell_type' by precedence (scanvi > celltypist > marker) and stamp provenance
// in 'cell_type_source'. Downstream modules use this single key regardless of
// which annotation tool ran.
def cell_type_key = 'cell_type'

// ATAC cell type column: depends on annotation mode
// marker_file → 'cell_type', scATAnno → 'cell_type_prediction', CellTypist → 'celltypist_prediction'
def atac_cell_type_key = params.atac.marker_file ? 'cell_type' :
    (params.atac.annotation_method == 'scatanno' ? 'cell_type_prediction' : 'celltypist_prediction')


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

                // P0-11 (A4): Require condition_group column if differential analysis is enabled
                def cgIdx = header.findIndexOf { it == 'condition_group' }
                if (cgIdx < 0 && (params.differential?.run ?: false)) {
                    errors << "Manifest CSV missing 'condition_group' column but differential analysis is enabled. " +
                              "All samples will default to 'Control', producing zero DE genes."
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
            errors << "Missing container files: ${containerMissing}. " +
                      "Expected in singularity_cache/. Run container build/pull first."
        } else {
            checks_passed << "Containers (${uniqueSifs.size()} SIF files)"
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

        log.info "Step 6b: Running CellTypist annotation..."
        RUN_CELLTYPIST(TRAIN_SCANVI.out.annotated)

        log.info "Step 7: Generating post-integration visualizations..."
        PLOT_POST_SCANVI(RUN_CELLTYPIST.out.annotated_h5ad, cell_type_key)

    } else {
        // ============================================================
        // PATH B: No reference atlas -> CellTypist (direct)
        // ============================================================
        log.info "No reference atlas provided -- using CellTypist (direct) path"

        log.info "Step 5-alt: Running CellTypist annotation..."
        RUN_CELLTYPIST(CONCAT_BATCHES.out.concatenated)

        log.info "Step 7: Generating post-integration visualizations..."
        PLOT_POST_SCANVI(RUN_CELLTYPIST.out.annotated_h5ad, cell_type_key)
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
        log.info "  Step 9d: Tier 1 -- enrichment & network visualization..."
        HDWGCNA_ENRICHMENT(
            HDWGCNA_PER_CELLTYPE.out.results.map { ct, rds -> rds },
            HDWGCNA_PER_CELLTYPE.out.results.map { ct, rds -> ct },
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

            // FIX: .out.results is tuple(val(cell_type), path(rds)) — destructure properly
            HDWGCNA_DIFFERENTIAL(
                HDWGCNA_PER_CELLTYPE.out.results.map { ct, rds -> rds },
                HDWGCNA_PER_CELLTYPE.out.results.map { ct, rds -> ct },
                cell_type_key,
                params.hdwgcna.condition_key,
                params.hdwgcna.control_condition,
                params.hdwgcna.treatment_condition,
                traits_str
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

    SNAPATAC_DIFFERENTIAL(
        peak_matrix,
        metadata,
        ch_tasks.map { tuple(it[0], it[1]) },
        ch_tasks.map { it[2] }
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
    // ================================================================
    if (params.cicero.run) {
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
            params.cicero.gtf_full
        )

        if (!use_chromvar_for_cicero && params.cicero.target_genes && !params.cicero.target_genes.isEmpty()) {
            log.info "Rendering Cicero target plots with static gene list: ${params.cicero.target_genes}"
            CICERO_TARGET_PLOTS(
                CICERO_JOIN.out.connections,
                CICERO_JOIN.out.ccan,
                CICERO_JOIN.out.cds,
                params.cicero.gtf_plot,
                params.cicero.target_genes
            )
        }
    }

    // ================================================================
    // PARALLEL LEG B: GPU ChromVAR TF Motif Enrichment
    // ================================================================
    if (params.chromvar.run) {
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

        if (has_da_peaks) {
            VIS_CHROMVAR(
                GPU_CHROMVAR.out.chromvar_dev
            )
        } else {
            log.info "Skipping VIS_CHROMVAR (requires differential conditions for permutation tests)"
        }

        if (is_discovery_mode && params.scprinter.run) {
            EXTRACT_CHROMVAR_MOTIFS(
                GPU_CHROMVAR.out.chromvar_dev,
                params.chromvar.top_n_per_celltype,
                params.chromvar.min_motif_zscore,
                atac_cell_type_key
            )

            EXTRACT_CHROMVAR_MOTIFS.out.report.view {
                "\n==== PER-CELL-TYPE CHROMVAR MOTIFS ====\n${it.text}\n======================================="
            }
        }
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
            CICERO_JOIN.out.connections,
            CICERO_JOIN.out.ccan,
            CICERO_JOIN.out.cds,
            params.cicero.gtf_plot,
            ch_chromvar_target_genes
        )
    }

    // ================================================================
    // STEP 3: scPRINTER TF Footprinting
    // ================================================================
    if (params.scprinter.run) {
        log.info "Running scPRINT TF footprinting workflow..."

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
                GPU_CHROMVAR.out.chromvar_raw,
                EXTRACT_CHROMVAR_MOTIFS.out.motif_list,
                CICERO_JOIN.out.ccan,
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
                CICERO_JOIN.out.connections.ifEmpty(file('NO_FILE')) :
                Channel.value(file('NO_FILE'))
            def pfm_for_fp = params.scprinter.pfms ?
                Channel.value(file(params.scprinter.pfms)) :
                Channel.value(file('NO_FILE'))

            SCPRINTER_FOOTPRINTING(
                peak_matrix,
                metadata,
                SCPRINTER_BUILD_PRINTER.out.printer,
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

                SCPRINTER_FOOTPRINTING_DIFF(
                    peak_matrix,
                    metadata,
                    SCPRINTER_BUILD_PRINTER.out.printer,
                    ch_fp.map { it[0] },
                    ch_fp.map { it[1] },
                    ch_fp.map { it[2] },
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
                    SCPRINTER_BUILD_PRINTER.out.printer,
                    ch_fp.map { it[0] },
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

            SCPRINTER_FOOTPRINTING(
                peak_matrix,
                metadata,
                SCPRINTER_BUILD_PRINTER.out.printer,
                'targeted',
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
                    SCPRINTER_BUILD_PRINTER.out.printer,
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
                    SCPRINTER_BUILD_PRINTER.out.printer,
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
    // EMIT
    // ================================================================
    emit:
    cicero_connections   = params.cicero.run ? CICERO_JOIN.out.connections : Channel.empty()
    cicero_ccan          = params.cicero.run ? CICERO_JOIN.out.ccan : Channel.empty()
    chromvar_deviations  = params.chromvar.run ? GPU_CHROMVAR.out.chromvar_dev : Channel.empty()
    chromvar_per_ct      = (is_discovery_mode && params.chromvar.run && params.scprinter.run) ?
        EXTRACT_CHROMVAR_MOTIFS.out.motif_list : Channel.empty()
    scprinter_printer    = params.scprinter.run ?
        SCPRINTER_BUILD_PRINTER.out.printer : Channel.empty()
    scprinter_footprints = params.scprinter.run ?
        SCPRINTER_FOOTPRINTING.out.footprints : Channel.empty()
    scprinter_diff       = (has_da_peaks && params.scprinter.run) ?
        SCPRINTER_FOOTPRINTING_DIFF.out.footprints : Channel.empty()
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

    // pycisTopic preparation
    if (params.pycistopic.run) {
        log.info "Running pycisTopic preparation..."

        PYCISTOPIC_PREPARE(
            metadata_csv,
            rna_h5ad,       // Use CellTypist-annotated RNA (not un-annotated MuData export)
            // FIX-P0-8: Auto-derive pycistopic species from params.species
            params.pycistopic.species ?: [human: 'hsapiens', mouse: 'mmusculus'].get(params.species, params.species),
            mudata_stats,
            blacklist_bed,
            file(params.pycistopic.gtf),
            cell_type_key
        )
    }

    // SCENIC+ via Snakemake
    if (params.scenicplus.run && params.pycistopic.run) {
        log.info "Running SCENIC+ Snakemake pipeline..."

        SCENICPLUS_RUN(
            PYCISTOPIC_PREPARE.out.cistopic_obj,
            rna_for_dorc,
            PYCISTOPIC_PREPARE.out.region_sets,
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
    cistopic_obj   = params.pycistopic.run ? PYCISTOPIC_PREPARE.out.cistopic_obj   : Channel.empty()
    region_sets    = params.pycistopic.run ? PYCISTOPIC_PREPARE.out.region_sets    : Channel.empty()
    gene_activity  = (params.pycistopic.run && PYCISTOPIC_PREPARE.out.gene_activity) ?
                     PYCISTOPIC_PREPARE.out.gene_activity : Channel.empty()
    pseudobulk_bigwigs = params.pycistopic.run ? PYCISTOPIC_PREPARE.out.pseudobulk_bigwigs : Channel.empty()

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

    main:

    // ================================================================
    // PHASE 1: ATAC-Only Enhancer Footprinting (Recipe A)
    // ================================================================
    log.info "ENHANCER FOOTPRINTING RECIPES: Phase 1 (ATAC-only)"

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
        cell_type_col
    )

    // ================================================================
    // PHASE 2: Multiome Integration (Recipe B)
    // ================================================================
    def has_scenic = (params.scenicplus.run ?: false) && (params.pycistopic.run ?: false)
    def has_dorc = (params.dorc.run ?: false)

    if (has_scenic) {
        log.info "ENHANCER FOOTPRINTING RECIPES: Phase 2 (Multiome integration)"

        def dorc_sig_file = has_dorc ?
            dorc_sig_ch : Channel.value(file('NO_FILE_dorc'))

        EXTRACT_EREGULON_REGIONS(
            ereg_direct_ch,
            r2g_ch,
            dorc_sig_file
        )

        CROSS_MODAL_VALIDATION(
            ENHANCER_FOOTPRINTING.out.footprints.collect(),
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

        def ereg_regions_ch = has_scenic ?
            EXTRACT_EREGULON_REGIONS.out.region_sets :
            Channel.value(file('NO_FILE_ereg_regions'))
        def ereg_manifest_ch = has_scenic ?
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

        SIGNAL_CHAIN_CORRELATION(
            ENHANCER_FOOTPRINTING.out.footprints.collect(),
            rna_h5ad_ch,
            EXTRACT_SIGNALING_TARGETS.out.metadata
        )
    }

    // ================================================================
    // PHASE 4: Composite Enhancer Visualization (Recipe D)
    // ================================================================
    if (params.enhancer_viz.run) {
        log.info "ENHANCER FOOTPRINTING RECIPES: Phase 4 (Composite Visualization)"

        def bigwig_dir_ch = pseudobulk_bigwigs_ch
            .ifEmpty(file('NO_BIGWIGS'))
            .collect()

        PREPARE_ENHANCER_VIZ_TRACKS(
            cicero_conns_ch,
            EXTRACT_CCAN_ENHANCERS.out.enhancer_peaks,
            bigwig_dir_ch,
            file(params.species == 'human' ? params.scprinter.gtf_human : params.scprinter.gtf_mouse),
            params.enhancer_viz.target_genes,
            MOTIF_SCAN_ENHANCERS.out.manifest
        )

        ch_viz_tasks = MOTIF_SCAN_ENHANCERS.out.manifest
            .flatMap { manifest_file ->
                def data = new groovy.json.JsonSlurper().parseText(manifest_file.text)
                def seen = new HashSet()
                data.region_sets.collect { entry ->
                    def key = "${entry.cell_type}_${entry.tf}"
                    if (!seen.contains(key)) {
                        seen.add(key)
                        tuple(entry.tf, entry.tf)
                    } else {
                        null
                    }
                }.findAll { it != null }
            }

        if (params.enhancer_viz.target_genes && !params.enhancer_viz.target_genes.isEmpty()) {
            ch_viz_tasks = MOTIF_SCAN_ENHANCERS.out.manifest
                .flatMap { manifest_file ->
                    def data = new groovy.json.JsonSlurper().parseText(manifest_file.text)
                    def tfs = data.region_sets.collect { it.tf }.unique()
                    def genes = params.enhancer_viz.target_genes
                    def pairs = []
                    genes.each { gene ->
                        tfs.each { tf ->
                            pairs.add(tuple(gene, tf))
                        }
                    }
                    pairs
                }
        }

        // FIX-P0-2: Use channel output instead of reading from publishDir
        def fp_pngs_ch = ENHANCER_FOOTPRINTING.out.plots.collect().ifEmpty([])

        COMPOSITE_ENHANCER_VIZ(
            PREPARE_ENHANCER_VIZ_TRACKS.out.track_manifest,
            PREPARE_ENHANCER_VIZ_TRACKS.out.track_inis.collect(),
            MOTIF_SCAN_ENHANCERS.out.motif_scan,
            ch_viz_tasks.map { it[0] },
            ch_viz_tasks.map { it[1] },
            fp_pngs_ch
        )
    }

    emit:
    enhancer_peaks     = EXTRACT_CCAN_ENHANCERS.out.enhancer_peaks
    motif_scan         = MOTIF_SCAN_ENHANCERS.out.motif_scan
    region_sets        = MOTIF_SCAN_ENHANCERS.out.region_sets
    enhancer_fps       = ENHANCER_FOOTPRINTING.out.footprints
    cross_modal        = has_scenic ? CROSS_MODAL_VALIDATION.out.validation_table : Channel.empty()
    evidence_tiers     = params.enhancer_recipe_c.run ? SIGNAL_CHAIN_CORRELATION.out.evidence_tiers : Channel.empty()
    enhancer_viz       = (params.enhancer_viz.run ?: false) ? COMPOSITE_ENHANCER_VIZ.out.composite_png : Channel.empty()
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

    On-ramps:
      rna_integrated_h5ad:    ${params.onramp?.rna_integrated_h5ad ?: 'none'}
      atac_peak_matrix_h5ad:  ${params.onramp?.atac_peak_matrix_h5ad ?: 'none'}
      mudata_h5mu:            ${params.onramp?.mudata_h5mu ?: 'none'}
      printer_h5ad:           ${params.onramp?.printer_h5ad ?: 'none'}
      cicero_connections:     ${params.onramp?.cicero_connections ?: 'none'}
      chromvar_deviations:    ${params.onramp?.chromvar_deviations ?: 'none'}

    """

    // ========================================================================
    // RNA PROCESSING (with on-ramp)
    // ========================================================================
    if (params.onramp?.rna_integrated_h5ad) {
        log.info "ON-RAMP: Using pre-integrated RNA: ${params.onramp.rna_integrated_h5ad}"
        ch_integrated_rna = Channel.value(file(params.onramp.rna_integrated_h5ad))
        rna_completed = true
        rna_from_onramp = true
    } else if (params.rna.run) {
        log.info "Starting RNA workflow..."
        RNA()
        ch_integrated_rna = RNA.out.integrated_rna
        rna_completed = true
        rna_from_onramp = false
    } else {
        log.info "Skipping RNA workflow (disabled in config)"
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
    // ATAC PROCESSING (with on-ramp)
    // ========================================================================
    if (params.onramp?.atac_peak_matrix_h5ad) {
        log.info "ON-RAMP: Using pre-computed ATAC peak matrix: ${params.onramp.atac_peak_matrix_h5ad}"
        ch_atac_peak_matrix = Channel.value(file(params.onramp.atac_peak_matrix_h5ad))
        atac_completed = true
        atac_from_onramp = true
    } else if (params.atac.run && params.atac.run_initial_qc) {
        log.info "Starting ATAC workflow..."

        // Stage 1: Initial QC
        ATAC_INITIAL()

        // Stage 2: Final QC with sample-specific thresholds
        ATAC_FINAL(ATAC_INITIAL.out.thresholds)

        ch_atac_peak_matrix = ATAC_FINAL.out.peak_matrix
        atac_completed = true
        atac_from_onramp = false

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

        // REGULATORY ANALYSIS (Cicero, ChromVAR, scPRINT)
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
            ATAC_FINAL.out.peak_matrix,
            ATAC_FINAL.out.individual_samples,
            da_peaks_ch,
            file(params.atac.sample_metadata),
            ch_reg_fragments,
            atac_cell_type_key
        )

    } else {
        log.info "Skipping ATAC workflow (disabled in config)"
        atac_completed = false
        atac_from_onramp = false
    }

    // ========================================================================
    // MULTIOME INTEGRATION (requires both RNA and ATAC)
    // ========================================================================
    if (params.onramp?.mudata_h5mu) {
        log.info "ON-RAMP: Using pre-computed MuData: ${params.onramp.mudata_h5mu}"
        // Skip MULTIOME_INTEGRATION entirely -- downstream GRN workflows
        // should read from the on-ramp mudata directly.
    } else if (params.run_multiome_integration && rna_completed && atac_completed && !atac_from_onramp && !rna_from_onramp) {
        log.info """
        MULTIOME INTEGRATION
        Waiting for RNA and ATAC workflows to complete...
        """

        // Explicit dependency: wait for ATAC_FINAL to complete
        ATAC_FINAL.out.individual_samples
            .collect()
            .subscribe { atac_files ->
                log.info "ATAC processing complete: ${atac_files.size()} samples"
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
        ch_rna_mapped = RNA.out.qc_h5ads
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
        ch_atac_mapped = ATAC_FINAL.out.individual_samples
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
    // ========================================================================
    if (rna_completed && atac_completed && !atac_from_onramp && !rna_from_onramp &&
        params.run_multiome_integration && !params.onramp?.mudata_h5mu &&
        (params.pycistopic.run || params.scenicplus.run || params.dorc.run)) {
        log.info """
        MULTIOME GRN ANALYSIS
        Running pycisTopic, SCENIC+ and/or DORC as requested in config
        """
        MULTIOME_GRN(
            ch_integrated_rna,
            ATAC_FINAL.out.peak_matrix,
            file(params.metadata_file),
            MULTIOME_INTEGRATION.out.rna_for_dorc,
            MULTIOME_INTEGRATION.out.stats,
            file(params.pycistopic.blacklist_bed),
            cell_type_key
        )
    }

    // ========================================================================
    // ENHANCER FOOTPRINTING RECIPES (A/B/C/D)
    // Requires: ATAC completed + REGULATORY_ANALYSIS (Cicero, ChromVAR, scPRINTER)
    // ========================================================================
    if (params.enhancer_footprinting.run && atac_completed && !atac_from_onramp &&
        params.cicero.run && params.chromvar.run && params.scprinter.run) {
        log.info """
        ENHANCER FOOTPRINTING RECIPES
        Phase 1 (ATAC-only):     always
        Phase 2 (Multiome):      ${(params.scenicplus.run ?: false) ? 'enabled' : 'disabled'}
        Phase 3 (CellChat):      ${(params.enhancer_recipe_c.run ?: false) ? 'enabled' : 'disabled'}
        Phase 4 (Visualization): ${(params.enhancer_viz.run ?: false) ? 'enabled' : 'disabled'}
        """

        def has_grn = rna_completed && atac_completed && !atac_from_onramp && !rna_from_onramp &&
                      params.run_multiome_integration && !params.onramp?.mudata_h5mu &&
                      (params.pycistopic.run || params.scenicplus.run || params.dorc.run)

        def ereg_direct_ch = (has_grn && params.scenicplus.run && params.pycistopic.run) ?
            MULTIOME_GRN.out.scplus_ereg_direct : Channel.empty()
        def r2g_ch = (has_grn && params.scenicplus.run && params.pycistopic.run) ?
            MULTIOME_GRN.out.scplus_r2g : Channel.empty()
        def dorc_sig_ch = (has_grn && params.dorc.run) ?
            MULTIOME_GRN.out.dorc_sig : Channel.empty()
        def rna_ch = rna_completed ? ch_integrated_rna : Channel.empty()
        def cellchat_csv_ch = (rna_completed && !rna_from_onramp && params.cellchat.run) ?
            RNA.out.cellchat_csv : Channel.empty()
        def bigwigs_ch = (has_grn && params.pycistopic.run) ?
            MULTIOME_GRN.out.pseudobulk_bigwigs : Channel.empty()

        ENHANCER_FOOTPRINTING_RECIPES(
            ATAC_FINAL.out.peak_matrix,
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
            atac_cell_type_key
        )
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
