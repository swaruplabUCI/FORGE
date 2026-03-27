# Refactor5 Pipeline: Prioritized Breakage Risk List (Round 2)

**Date:** 2026-03-25
**Reviewer:** Senior Bioinformatics Engineering Review
**Branch:** dev (commit 9877678)
**Round 2 expansion:** Merged 25 original Round 1 issues with 37 new input-validation findings from systematic "Useful Idiot" analysis.

## What Changed from Round 1

Round 1 identified 25 issues from code review of the Nextflow orchestration, config system, and key Python scripts. Round 2 added 37 issues from a systematic input-validation walkthrough covering manifest CSV handling, reference file assumptions, ATAC/RNA-specific inputs, configuration parameter validation, and intermediate file coupling. After deduplication (8 overlaps identified), the merged list contains **54 unique issues**, re-ranked from scratch using stricter severity criteria.

**Key re-ranking changes:**
- Several Round 2 "silent corruption" findings (barcode mismatch, GTF genome build mismatch, species parameter mismatch) enter at P0 because they produce wrong biology with no warning.
- Round 1 items that were P2 (e.g., hardcoded DE comparisons) are promoted because Round 2 confirmed they silently produce empty results for any non-default experiment design.
- Round 1's CSV validation issue (formerly rank 21, P2) is now decomposed into specific failure modes (A1-A7) and distributed across tiers by severity.

## Ranking Criteria

1. **Silent corruption** (wrong biology, wrong cell pairing) >> **Silent empty results** (analysis runs but produces nothing) >> **Late crashes** (after expensive GPU compute) >> **Immediate crashes** >> **Cosmetic/resource waste**
2. **Affects any new user** >> **Affects only edge cases**
3. **Simple fix** prioritized as quick wins within each severity band

## Priority Tier Definitions

| Tier | Meaning |
|------|---------|
| P0   | Fix before ANY production run -- data integrity at risk |
| P1   | Fix before sharing with collaborators -- portability/usability blockers |
| P2   | Fix before publication/release -- robustness/maintainability |
| P3   | Nice-to-have cleanup |

---

## Summary Table

| Rank | ID | Title | Severity | Tier | Category | Source |
|------|----|-------|----------|------|----------|--------|
| 1 | R1-2 | Global `errorStrategy = 'ignore'` silently swallows critical failures | Silent corruption | P0 | Nextflow config | Round 1 |
| 2 | R1-18 | Reading from publishDir instead of process output channels | Silent corruption | P0 | Nextflow config | Round 1 |
| 3 | R1-12 / A3 | `ad.concat(index_unique=None)` + duplicate sample_id corrupt cell pairing | Silent corruption | P0 | Data integrity | Round 1 + Round 2 |
| 4 | B1 | GTF genome build mismatch with fragment files | Silent corruption | P0 | Reference files | NEW |
| 5 | C2 | RNA-ATAC barcode format mismatch causes zero paired cells | Silent corruption | P0 | ATAC input | NEW |
| 6 | B2 | Species parameter mismatch with GTF/reference atlas | Silent corruption | P0 | Reference files | NEW |
| 7 | A2 | sample_type case mismatch silently skips entire workflow arms | Silent corruption | P0 | Manifest CSV | NEW |
| 8 | E1 / E2 | species='mouse' with human references / pycistopic.species divergence | Silent corruption | P0 | Config params | NEW |
| 9 | B4 | SCENIC+ feather/rankings wrong species or version | Silent corruption | P0 | Reference files | NEW |
| 10 | F11 | Positional barcode matching for BD silently pairs wrong cells | Silent corruption | P0 | Intermediate files | NEW |
| 11 | R1-1 / A4 | Undefined legacy params + missing condition_group silently defaults | Silent corruption | P0 | Nextflow config | Round 1 + Round 2 |
| 12 | R1-10 | Fragile NO_FILE sentinel pattern across main.nf | Silent corruption | P0 | Nextflow config | Round 1 |
| 13 | B5 | Blacklist BED file for wrong genome build | Silent corruption | P0 | Reference files | NEW |
| 14 | E3 | cellbender.expected_cells wildly wrong causes quality degradation | Silent corruption | P0 | Config params | NEW |
| 15 | R1-5 | ~~User-specific conda path in launch.sh~~ | ~~Immediate crash~~ | ~~P1~~ | ~~Portability~~ | **FIXED** — Removed dead conda source line; pipeline is fully containerized, nextflow is a standalone binary on PATH |
| 16 | R1-24 | ~~SLURM account hardcoded to vswarup_lab~~ | ~~Immediate crash~~ | ~~P1~~ | ~~Portability~~ | **RESOLVED (no change needed)** — Already parameterized via `params.slurm_account`/`params.slurm_account_gpu` in nextflow.config; collaborators override in dataset config. Only hardcoded instance is `#SBATCH -A` in launch.sh for the launcher job itself, which is expected. |
| 17 | R1-4 | ~~Hardcoded /dfs7 absolute paths throughout configs~~ | ~~Immediate crash~~ | ~~P1~~ | ~~Portability~~ | **RESOLVED (no change needed)** — Dataset configs intentionally contain site-specific data/reference paths (collaborators write their own). launch.sh environment block (`PATH`, `SINGULARITY_BINDPATH`) is inherently site-specific infrastructure. |
| 18 | R1-14 | ~~Hardcoded SLURM partition names~~ | ~~Immediate crash~~ | ~~P1~~ | ~~Portability~~ | **RESOLVED (no change needed)** — Default partitions already parameterized in nextflow.config. Per-process overrides (`highmem`, `hugemem`, `maxmem`) in resource tier configs are site-specific HPC tuning that collaborators would re-profile for their own cluster. |
| 19 | R1-6 | Missing medium resource tier | Silent degradation | P1 | Config system | Round 1 |
| 20 | R1-23 | pycistopic.gtf = null but enabled by default | Immediate crash | P1 | Config params | Round 1 |
| 21 | E4 | mt_threshold=0 filters all cells | Silent empty results | P1 | Config params | NEW |
| 22 | E5 / R1-15 | Hardcoded Group1-4 DE comparisons vs. user group names | Silent empty results | P1 | Config params | Round 1 + Round 2 |
| 23 | D1 | BD MEX directory with missing matrix files | Immediate crash | P1 | RNA input | NEW |
| 24 | A1 | No schema validation of CSV column names | Silent corruption | P1 | Manifest CSV | NEW |
| 25 | C4 | annotation_resolution referencing non-existent clustering resolution | Silent empty results | P1 | ATAC input | NEW |
| 26 | R1-7 | Container name mismatch between launch.sh check and config | Delayed crash | P1 | Config system | Round 1 |
| 27 | R1-16 | GPU profile applies to ALL processes including CPU-only | Resource waste | P1 | Nextflow config | Round 1 |
| 28 | R1-22 | SINGULARITY_BINDPATH has no effect with --contain | Portability trap | P1 | Nextflow config | Round 1 |
| 29 | F1 | cell_type_key mismatch between pipeline stages | Silent empty results | P1 | Intermediate files | NEW |
| 30 | B3 | Reference atlas directory with no h5ad files | Late crash | P1 | Reference files | NEW |
| 31 | B6 | CellTypist model inappropriate for tissue type | Silent wrong annotations | P2 | Reference files | NEW |
| 32 | C1 | Fragment file sort order mismatch for scPRINTER | Silent noisy results | P2 | ATAC input | NEW |
| 33 | R1-9 | Only human MT genes detected; mouse mt- ignored | Silent QC bypass | P2 | Species handling | Round 1 |
| 34 | R1-3 / D3 / F2 | Hardcoded batch detection (june/july/nov) + extract_sample_id BD naming | Silent wrong batch | P2 | Batch handling | Round 1 + Round 2 |
| 35 | D4 | Scrublet on very small samples silently filters all cells | Silent empty results | P2 | RNA input | NEW |
| 36 | E6 | chromvar.min_motif_zscore too high filters all TFs | Silent empty results | P2 | Config params | NEW |
| 37 | E7 | differential.cell_types empty while differential.run=true | Silent empty results | P2 | Config params | NEW |
| 38 | F4 | run_snapatac_diff.py hardcodes cell_type/condition_group columns | Silent empty results | P2 | Intermediate files | NEW |
| 39 | F5 | make_atac_thresholds_from_csv.py assumes specific column names | Silent empty results | P2 | Intermediate files | NEW |
| 40 | A5 | Whitespace in CSV fields silently breaks file path resolution | Silent corruption | P2 | Manifest CSV | NEW |
| 41 | C3 | ATAC sample_metadata diverging from main manifest | Silent inconsistency | P2 | ATAC input | NEW |
| 42 | C5 | Non-standard fragment file column format | Silent misinterpretation | P2 | ATAC input | NEW |
| 43 | R1-13 | /tmp as TMPDIR on shared HPC (disk exhaustion) | Late crash | P2 | HPC config | Round 1 |
| 44 | R1-11 | MOFA gets GPU container but CPU partition in small tier | Config mismatch | P2 | Nextflow config | Round 1 |
| 45 | R1-17 / F9 | pathway_to_tfs.json not validated before use | Silent empty output | P2 | Intermediate files | Round 1 + Round 2 |
| 46 | D2 | RNA h5 file with non-standard structure | Cryptic error | P2 | RNA input | NEW |
| 47 | E8 | mofa.mode invalid string | Late crash | P2 | Config params | NEW |
| 48 | F3 | concat_batches.py outer join creates NaN genes across batches | Silent data quality | P2 | Intermediate files | NEW |
| 49 | A6 | Empty rows/trailing newlines produce phantom samples | Silent corruption | P2 | Manifest CSV | NEW |
| 50 | A7 | File paths not checked in dry-run/preview mode | UX gap | P2 | Manifest CSV | NEW |
| 51 | R1-8 | eval EXTRA_ARGS shell injection risk | Security | P2 | launch.sh | Round 1 |
| 52 | F12 | params.scprinter.genome not validated against ATAC genome | Silent mismatch | P2 | Intermediate files | NEW |
| 53 | F6 | Barcode delimiter inconsistency across pipeline stages | Silent mismatch | P2 | Intermediate files | NEW |
| 54 | F7 | scenicplus.bc_transform_func is an eval'd lambda string | Security/fragility | P2 | Intermediate files | NEW |
| 55 | R1-23 / F8 | hdWGCNA filename parsing for cell type is fragile | Fragile parsing | P3 | Intermediate files | Round 1 + Round 2 |
| 56 | R1-19 | builtins.long = int monkey-patch fragility | Fragile workaround | P3 | Python compat | Round 1 |
| 57 | R1-20 | Duplicate process definitions across config files | Maintenance burden | P3 | Nextflow config | Round 1 |
| 58 | F10 | launch.sh resource tier auto-detection counts blank lines | Minor bug | P3 | launch.sh | NEW |

---

## Detailed Issues by Tier

### P0: Fix Before ANY Production Run (14 issues)

**1. R1-2: Global `errorStrategy = 'ignore'` silently swallows critical failures** (Round 1)
The cluster profile sets the default error strategy to `ignore`, meaning any core process that fails is silently skipped. Downstream steps receive empty channels and produce incomplete results with no error. This is the single worst failure mode for a scientific pipeline -- it makes every other bug harder to detect.
*Fix: Change default to `'retry'` with `maxRetries = 2`, use `'ignore'` only on explicitly optional processes.*

**2. R1-18: Reading from publishDir instead of process output channels** (Round 1)
If any process reads files from publishDir rather than Nextflow output channels, results depend on filesystem timing. A slow NFS can produce silent data loss or stale inputs. This breaks Nextflow's DAG guarantees and is impossible to debug.
*Fix: Audit all process inputs to ensure they come from channel emissions, never from publishDir paths.*

**3. R1-12 / A3: `ad.concat(index_unique=None)` + duplicate sample_id corrupt cell pairing** (Round 1 + Round 2)
Found in 3+ scripts including concat_batches.py. Duplicate barcodes across samples are silently collapsed, corrupting all downstream per-cell analyses. Round 2 adds that duplicate sample_id in the manifest CSV causes wrong RNA-ATAC pairing in channel joins, compounding the barcode collision.
*Fix: Add `index_unique='-'` or equivalent, and validate sample_id uniqueness at manifest parse time.*

**4. B1: GTF genome build mismatch with fragment files** (NEW)
If the user provides an hg38 GTF but fragment files were generated against hg19 (or vice versa), all coordinate-based operations (peak calling, gene activity scoring, pycisTopic) silently use wrong genomic coordinates. The pipeline produces plausible-looking but biologically meaningless results.
*Fix: Extract genome build from GTF header and fragment file header; assert they match at pipeline start.*

**5. C2: RNA-ATAC barcode format mismatch causes zero paired cells** (NEW)
10x multiome and BD Rhapsody use different barcode formats. If RNA barcodes have a `-1` suffix but ATAC barcodes do not (or vice versa), the pairing join produces zero matches. With `errorStrategy = 'ignore'`, this silently drops all multiome integration.
*Fix: Normalize barcode format (strip/add suffix) before any RNA-ATAC join; warn if join yield < 10%.*

**6. B2: Species parameter mismatch with GTF/reference atlas** (NEW)
Setting `params.species = 'human'` but providing a mouse GTF or mouse reference atlas produces zero marker gene matches, zero cell type annotations, and zero regulatory networks -- all silently, since empty results are not treated as errors.
*Fix: Validate species consistency across GTF gene name prefixes, reference atlas, and params.species at startup.*

**7. A2: sample_type case mismatch silently skips entire workflow arms** (NEW)
The manifest CSV `sample_type` column is matched with exact string comparison. A value of `rna` instead of `RNA` (or `Rna`, `rna `) causes the sample to be excluded from the RNA workflow arm entirely. With `errorStrategy = 'ignore'`, the pipeline completes successfully with missing samples.
*Fix: Normalize sample_type to uppercase and strip whitespace at manifest parse time.*

**8. E1 / E2: species='mouse' with human references / pycistopic.species divergence** (NEW)
`params.species` controls MT gene prefix detection and marker gene lists, but `params.pycistopic.species` is a separate parameter. If they diverge, pycisTopic uses wrong genome annotations while the rest of the pipeline uses another. Additionally, setting species='mouse' but leaving all reference paths at their human defaults produces silent wrong results everywhere.
*Fix: Derive pycistopic.species from params.species; validate all reference paths contain species-consistent content.*

**9. B4: SCENIC+ feather/rankings files wrong species or version** (NEW)
SCENIC+ cisTarget databases are species- and genome-build-specific. Using human rankings with mouse data (or hg19 rankings with hg38 peaks) produces eRegulon predictions that look plausible but are biologically wrong. No runtime error occurs.
*Fix: Validate that cisTarget database filenames contain the expected species/build string; warn on mismatch.*

**10. F11: Positional barcode matching for BD silently pairs wrong cells** (NEW)
BD Rhapsody barcode pairing relies on positional matching (same row index in RNA and ATAC barcode lists). If the lists are sorted differently or have different filtering, position N in RNA maps to a completely different cell than position N in ATAC -- producing silent chimeric cells.
*Fix: Use explicit barcode sequence matching rather than positional indexing; add a join-yield sanity check.*

**11. R1-1 / A4: Undefined legacy params + missing condition_group defaults to 'Control'** (Round 1 + Round 2)
Legacy params (`june_metadata`, `july_souporcell_dir`) are still referenced in main.nf. For any non-BD-brain dataset they resolve to NO_FILE sentinels via a fragile implicit contract. Round 2 adds that a missing `condition_group` column in the manifest silently defaults all samples to 'Control', making all differential analyses produce zero DE genes (no contrast).
*Fix: Remove legacy params; require condition_group in manifest schema; fail early if missing.*

**12. R1-10: Fragile NO_FILE sentinel pattern across main.nf** (Round 1)
Used in 15+ locations with inconsistent sentinel names (`NO_FILE`, `NO_FILE_JULY`, `NO_FILE_dorc`). Any process that does not explicitly check for these strings attempts to open a nonexistent file, and `errorStrategy = 'ignore'` hides the crash.
*Fix: Use Nextflow's native optional input pattern; remove string-based sentinel checking.*

**13. B5: Blacklist BED file for wrong genome build** (NEW)
ENCODE blacklist regions are genome-build-specific. Using an hg19 blacklist with hg38 data (or vice versa) removes the wrong genomic regions from ATAC peak calling, introducing noise or removing real signal. No error is raised.
*Fix: Validate blacklist BED chromosome names against fragment file chromosomes at startup.*

**14. E3: cellbender.expected_cells wildly wrong causes quality degradation** (NEW)
CellBender's `expected_cells` parameter strongly influences ambient RNA removal. Setting it 10x too high or too low produces either over-correction (removing real signal) or under-correction (leaving ambient RNA). Both are silent and affect all downstream RNA analyses.
*Fix: Add a sanity check that expected_cells is within 2x of the knee-point estimate from the raw matrix; warn if not.*

---

### P1: Fix Before Sharing with Collaborators (16 issues)

**15. R1-5: User-specific conda path in launch.sh** (Round 1) — **FIXED**
~~Any collaborator who clones the repo gets an immediate launch failure due to `/pub/lesolano/miniconda3`.~~
*Resolution: Removed dead `source /pub/lesolano/miniconda3/etc/profile.d/conda.sh` line. Pipeline is fully containerized; Nextflow is a standalone binary available via PATH. No conda environment is activated or needed.*

**16. R1-24: SLURM account hardcoded to vswarup_lab** (Round 1) — **RESOLVED (no change needed)**
~~All clusterOptions use `params.slurm_account` defaulting to `vswarup_lab`. Any external collaborator's jobs are rejected by the scheduler immediately.~~
*Resolution: Upon review, SLURM account is already properly parameterized in nextflow.config (`params.slurm_account`, `params.slurm_account_gpu`). All process configs reference the parameter, not the raw string. Collaborators override via their dataset config. The only hardcoded instance (`#SBATCH -A vswarup_lab` in launch.sh) is for the lightweight launcher job and is expected to be edited alongside other SBATCH directives.*

**17. R1-4: Hardcoded /dfs7 absolute paths throughout configs** (Round 1) — **RESOLVED (no change needed)**
~~Binds the pipeline to UCI HPC3. Affects launch.sh, singularity runOptions, and tool paths.~~
*Resolution: The `/dfs7` paths fall into two categories: (1) Dataset configs contain site-specific data/reference paths by design — collaborators write their own dataset config with their own paths. (2) launch.sh environment block (PATH, SINGULARITY_BINDPATH) is inherently site-specific infrastructure that any collaborator would customize for their HPC. No abstraction adds value here.*

**18. R1-14: Hardcoded SLURM partition names** (Round 1) — **RESOLVED (no change needed)**
~~20+ clusterOptions lines reference UCI-specific partitions (`highmem`, `hugemem`, `maxmem`). Job submission fails at any other site.~~
*Resolution: Default partitions (`standard`, `gpu`, `gpu-hugemem`) are already parameterized in nextflow.config. Per-process overrides in resource tier configs (e.g., `highmem` for ANNOTATE_ATAC, `maxmem` for BUILD_MUDATA) are site-specific HPC tuning based on memory profiling. A collaborator would need to re-profile resource requirements for their own cluster regardless — parameterizing partition names alone doesn't solve the portability problem.*

**19. R1-6: Missing medium resource tier** (Round 1)
launch.sh assigns `medium` for 6-50 samples but only `small` and `large` configs exist. Falls back to `small`, causing OOM kills on moderately sized datasets.
*Fix: Create medium.config or adjust tier boundaries.*

**20. R1-23: pycistopic.gtf = null but enabled by default** (Round 1)
pycisTopic is enabled (`run = true`) but its required GTF defaults to null. Runtime fails with an unhelpful null-pointer error.
*Fix: Validate that all required params for enabled modules are non-null at startup.*

**21. E4: mt_threshold=0 filters all cells** (NEW)
Setting `mt_threshold = 0` (intending "no filter") actually filters every cell since all cells have >= 0% mitochondrial reads. The parameter semantics are inverted from user expectation.
*Fix: Use null or -1 as "disabled" sentinel; document that 0 means "remove all cells with any mito reads."*

**22. E5 / R1-15: Hardcoded Group1-4 DE comparisons vs. user group names** (Round 1 + Round 2)
RNA_DIFFERENTIAL hardcodes 6 pairwise comparisons between Group1-Group4. Any experiment with different condition names produces silent empty DE results (zero genes, no error). Round 2 confirms this affects both RNA and ATAC differential modules.
*Fix: Read comparisons from params or auto-generate from unique condition_group values in the manifest.*

**23. D1: BD MEX directory with missing matrix files** (NEW)
BD Rhapsody MEX output requires matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz. A missing file causes a crash, but the error message is cryptic (scanpy read error, not "missing barcodes.tsv.gz").
*Fix: Check for all three required files at manifest validation time; report which file is missing.*

**24. A1: No schema validation of CSV column names** (NEW)
Misspelled column names (e.g., `sampl_type` instead of `sample_type`) return null in Groovy channel operations. With `errorStrategy = 'ignore'`, the entire pipeline can run with all samples silently excluded.
*Fix: Validate required columns (sample_id, sample_type, data_dir, condition_group, fragment_file) at manifest parse; fail with column name suggestions on mismatch.*

**25. C4: annotation_resolution referencing non-existent clustering resolution** (NEW)
If `params.annotation_resolution` references a resolution not computed in the clustering step, downstream annotation and differential analysis silently produce empty results.
*Fix: Validate that annotation_resolution is in the list of computed resolutions before proceeding.*

**26. R1-7: Container name mismatch between launch.sh check and config** (Round 1)
launch.sh checks for `cicero.sif` but config defines `r_cicero` pointing to `cicero.sif`. If naming drifts, pre-flight passes but process fails.
*Fix: Generate container check list from the same config source used at runtime.*

**27. R1-16: GPU profile applies to ALL processes including CPU-only** (Round 1)
The `gpu` profile sets `--gres=gpu` globally. CPU-only processes request and waste GPU allocations, increasing queue wait and burning credits.
*Fix: Apply GPU clusterOptions only to GPU-accelerated process labels.*

**28. R1-22: SINGULARITY_BINDPATH has no effect with --contain** (Round 1)
The singularity profile uses `--contain` which ignores the environment variable. The actual bind mounts come from runOptions. Misleading for anyone adding new bind paths.
*Fix: Remove the SINGULARITY_BINDPATH export; document that bind paths must go in runOptions.*

**29. F1: cell_type_key mismatch between pipeline stages** (NEW)
Different modules expect the cell type annotation in different `.obs` column names (`cell_type`, `celltype`, `CellType`, `cell_type_key`). If a module writes `cell_type` but the next reads `celltype`, all downstream grouping is wrong or empty.
*Fix: Standardize on a single cell_type_key parameter used across all modules; validate at handoff points.*

**30. B3: Reference atlas directory with no h5ad files** (NEW)
If the reference atlas directory exists but contains no .h5ad files (wrong path, wrong extension), the pipeline proceeds through expensive GPU-based CellBender and scVI steps before crashing at the annotation step.
*Fix: Validate reference atlas file existence at pipeline startup, before any compute.*

---

### P2: Fix Before Publication/Release (24 issues)

**31. B6: CellTypist model inappropriate for tissue type** (NEW)
Using a brain-trained CellTypist model on kidney data produces plausible but wrong cell type labels. No error occurs. This is fundamentally a user-responsibility issue but the pipeline can help.
*Fix: Log which CellTypist model is being used; optionally validate model tissue type against a user-declared tissue parameter.*

**32. C1: Fragment file sort order mismatch for scPRINTER** (NEW)
scPRINTER requires position-sorted fragment files. If the input is name-sorted, footprinting results are noisy/wrong with no error.
*Fix: Check sort order of first 1000 lines at process start; re-sort if needed or fail with a clear message.*

**33. R1-9: Only human MT genes detected; mouse mt- ignored** (Round 1)
concat_batches.py only checks `MT-` prefix. Mouse datasets have 0% mito genes, disabling QC filtering and passing low-quality cells through.
*Fix: Use species parameter to select MT-/mt- prefix; or detect both.*

**34. R1-3 / D3 / F2: Hardcoded batch detection + extract_sample_id BD naming** (Round 1 + Round 2)
concat_batches.py batch inference is hardcoded for BD naming conventions (`june`/`july`/`nov`). extract_sample_id() filename parsing also assumes BD directory structure. Non-BD datasets get wrong batch assignments silently.
*Fix: Read batch from manifest CSV rather than inferring from filenames.*

**35. D4: Scrublet on very small samples silently filters all cells** (NEW)
Scrublet's doublet detection on samples with fewer than ~100 cells can flag nearly all cells as doublets. The pipeline silently removes them, producing an empty AnnData that propagates downstream.
*Fix: Skip Scrublet for samples below a minimum cell count; log a warning.*

**36. E6: chromvar.min_motif_zscore too high filters all TFs** (NEW)
If `min_motif_zscore` is set too aggressively (e.g., > 10), all transcription factors are filtered out and the regulatory analysis produces empty output silently.
*Fix: Warn if zero TFs pass the threshold; suggest a lower value.*

**37. E7: differential.cell_types empty while differential.run=true** (NEW)
If `differential.run = true` but `differential.cell_types` is an empty list, the differential module runs with no cell types and produces no output -- silently.
*Fix: If cell_types is empty and run is true, default to all cell types found in the data; log the decision.*

**38. F4: run_snapatac_diff.py hardcodes cell_type/condition_group column names** (NEW)
The SnapATAC differential script hardcodes `.obs` column names rather than reading them from parameters. Datasets with different column naming conventions produce empty or wrong groupings.
*Fix: Accept column names as script arguments; pass from Nextflow params.*

**39. F5: make_atac_thresholds_from_csv.py assumes specific column names** (NEW)
Similar to F4 -- the threshold CSV parser assumes specific headers. Non-standard headers produce silent wrong thresholds.
*Fix: Validate expected columns at parse time; fail with a clear message listing expected vs. found columns.*

**40. A5: Whitespace in CSV fields silently breaks file path resolution** (NEW)
Trailing spaces in manifest CSV fields (e.g., `data_dir = "/path/to/data "`) cause file-not-found errors that are hidden by `errorStrategy = 'ignore'`.
*Fix: Strip whitespace from all CSV fields at parse time.*

**41. C3: ATAC sample_metadata diverging from main manifest** (NEW)
Some ATAC processes use a separate sample_metadata file. If it diverges from the main manifest (different sample IDs, different condition groups), RNA and ATAC arms of the pipeline analyze different sample sets.
*Fix: Generate ATAC sample_metadata from the main manifest programmatically; do not accept a separate file.*

**42. C5: Non-standard fragment file column format** (NEW)
Fragment files are expected to be BED-like with columns: chr, start, end, barcode, count. Files with extra columns or different ordering are silently misinterpreted.
*Fix: Validate fragment file header/first line format at startup.*

**43. R1-13: /tmp as TMPDIR on shared HPC** (Round 1)
Large intermediate files can fill the shared /tmp partition, crashing other users' jobs.
*Fix: Set TMPDIR to a job-specific directory under the work directory.*

**44. R1-11: MOFA gets GPU container but CPU partition in small tier** (Round 1)
MOFA_INTEGRATE is in the GPU accelerator block but small.config routes it to CPU partition. Job either fails or wastes a GPU slot.
*Fix: Ensure GPU processes always get GPU partitions regardless of resource tier.*

**45. R1-17 / F9: pathway_to_tfs.json not validated before use** (Round 1 + Round 2)
The CellChat-to-TF hypothesis script loads this JSON without schema checks. Malformed or empty file produces silent empty output.
*Fix: Validate JSON schema at load time; fail if required keys are missing.*

**46. D2: RNA h5 file with non-standard structure** (NEW)
Non-10x h5 files (e.g., from BD Rhapsody or custom pipelines) may have different internal group names. scanpy.read_10x_h5() fails with a cryptic HDF5 key error.
*Fix: Detect h5 format and route to appropriate reader; provide a clear error message for unsupported formats.*

**47. E8: mofa.mode invalid string** (NEW)
An invalid `mofa.mode` string (e.g., `'muliti'` instead of `'multi'`) is not caught until MOFA training starts, after potentially hours of upstream compute.
*Fix: Validate mofa.mode against allowed values at pipeline startup.*

**48. F3: concat_batches.py outer join creates NaN genes across batches** (NEW)
Using outer join in concatenation introduces NaN values for genes present in some batches but not others. Downstream analyses may silently treat NaN as zero or fail on NaN arithmetic.
*Fix: Use inner join by default or impute zeros explicitly; log the number of genes lost/gained.*

**49. A6: Empty rows/trailing newlines produce phantom samples** (NEW)
Trailing newlines or empty rows in the manifest CSV create phantom sample entries with all-null fields, causing cryptic downstream errors.
*Fix: Filter empty rows at manifest parse time.*

**50. A7: File paths not checked in dry-run/preview mode** (NEW)
Nextflow's `-preview` mode does not trigger file existence checks. Users cannot validate their manifest paths without starting actual compute.
*Fix: Add a validate_manifest process at the start of the pipeline that checks all paths regardless of execution mode.*

**51. R1-8: eval EXTRA_ARGS shell injection risk** (Round 1)
EXTRA_ARGS from positional arguments is eval'd. Low risk on authenticated HPC but violates defense-in-depth.
*Fix: Use an array and proper quoting instead of eval.*

**52. F12: params.scprinter.genome not validated against ATAC genome** (NEW)
scPRINTER has its own genome parameter that can diverge from the ATAC genome used for peak calling. Mismatched genomes produce wrong footprinting coordinates silently.
*Fix: Derive scprinter.genome from the pipeline-wide genome parameter.*

**53. F6: Barcode delimiter inconsistency across pipeline stages** (NEW)
Some stages use `_` as barcode delimiter, others use `-`, and auto-detection heuristics can fail on edge cases. This causes silent barcode mismatches at integration points.
*Fix: Standardize on a single delimiter; normalize at each stage boundary.*

**54. F7: scenicplus.bc_transform_func is an eval'd lambda string** (NEW)
The barcode transform function is stored as a string and eval'd at runtime. This is both a security concern and a fragility issue (syntax errors produce cryptic failures).
*Fix: Replace with a named transform option (e.g., 'strip_suffix', 'add_prefix') rather than arbitrary code execution.*

---

### P3: Nice-to-Have Cleanup (4 issues)

**55. R1-23 / F8: hdWGCNA filename parsing for cell type is fragile** (Round 1 + Round 2)
The `.replaceAll(/^hdwgcna_/, '')` pattern extracts cell type from filename. Cell types with underscores or special characters are mangled. Round 2 confirms this is a general pattern across several scripts.
*Fix: Pass cell type as an explicit process argument rather than encoding in filenames.*

**56. R1-19: builtins.long = int monkey-patch fragility** (Round 1)
Present in 3 scPRINTER scripts. Works around a Python 3 incompatibility in an upstream library but will break if the library updates.
*Fix: Pin the library version or submit an upstream fix; remove the monkey-patch when possible.*

**57. R1-20: Duplicate process definitions across config files** (Round 1)
Container and resource assignments are duplicated between nextflow.config, small.config, and large.config. Known Nextflow limitation.
*Fix: Extract shared definitions into a base config included by all tiers.*

**58. F10: launch.sh resource tier auto-detection counts blank lines** (NEW)
The sample count logic in launch.sh counts blank lines in the manifest CSV, potentially assigning a larger resource tier than needed.
*Fix: Filter blank lines before counting: `tail -n +2 "$csv" | grep -c '[^[:space:]]'`.*

---

## Summary by Tier

| Tier | Count | Action |
|------|-------|--------|
| P0 | 14 | Fix immediately -- silent data corruption or masked failures |
| P1 | 16 | Fix before any collaborator runs the pipeline -- portability and usability blockers |
| P2 | 24 | Fix before publication -- correctness edge cases, hardening, and robustness |
| P3 | 4 | Cleanup when convenient -- tech debt, minor UX |
| **Total** | **58** | |

## Deduplication Notes

The following Round 1 / Round 2 items were merged (8 overlaps):

| Merged ID | Round 1 Source | Round 2 Source | Rationale |
|-----------|---------------|----------------|-----------|
| Rank 3 | R1-12 (ad.concat barcode merge) | A3 (duplicate sample_id) | Both cause wrong cell pairing via barcode collision |
| Rank 8 | -- | E1 + E2 | Both are species parameter consistency issues |
| Rank 11 | R1-1 (legacy params) | A4 (missing condition_group) | Both cause silent wrong defaults from missing manifest data |
| Rank 22 | R1-15 (hardcoded DE comparisons) | E5 (Group1-4 vs user names) | Same issue described from config and user perspectives |
| Rank 34 | R1-3 (batch detection) | D3 + F2 (BD naming hardcoded) | All three are manifestations of BD-specific filename assumptions |
| Rank 45 | R1-17 (pathway_to_tfs.json) | F9 (same file not validated) | Identical issue found independently |
| Rank 55 | R1-23 (hdWGCNA parsing) | F8 (same fragile parsing) | Identical issue found independently |
| Rank 24 | R1-21 (CSV validation) | A1 (column name validation) | R1-21 was general; A1 is the specific actionable version |
