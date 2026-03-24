// modules/chromvar/gpu_chromvar.nf
// GPU-accelerated ChromVAR motif enrichment using scprinter + cupy
//
// FIX-23b: Updated NO_FILE sentinel checks to use .startsWith('NO_FILE')
//          instead of exact == 'NO_FILE' match.  main.nf now emits unique
//          sentinel names (NO_FILE_METADATA, NO_FILE_DA_PEAKS, etc.) to
//          avoid Nextflow file-name collisions when multiple path inputs
//          share the same staged filename.
//
// FIX-24:  Removed Cicero inputs (cicero_conns, ccan_assignments).
//          GPU_CHROMVAR now runs independently of Cicero, enabling parallel
//          execution.  Cicero co-accessibility and ChromVAR motif enrichment
//          answer complementary questions; neither requires the other's output.
//          ChromVAR-discovered TFs now drive BOTH Cicero target plots AND
//          scPRINTER footprinting in discovery mode.
//
// FIX-31:  Added --chunk-size passthrough from params.chromvar.chunk_size.
//          A100 (80GB VRAM) can handle 30-40K cells per chunk vs 10K on A30,
//          reducing total chunks from 12 to 4 and eliminating most RMM spilling.

process GPU_CHROMVAR {
    label 'process_gpu'
    publishDir "${params.outdir}/chromvar", mode: 'copy'

    input:
    path peak_matrix
    path metadata
    val cache_dir
    val pfms
    val genome
    path da_peaks

    output:
    path "chromvar_matrix_deviations.h5ad", emit: chromvar_dev
    path "chromvar_matrix.h5ad", emit: chromvar_raw

    script:
    // FIX-23b: sentinel names are now NO_FILE_DA_PEAKS, NO_FILE_METADATA, etc.
    // Use .startsWith('NO_FILE') so any NO_FILE* variant is treated as absent.
    def da_arg    = !da_peaks.name.startsWith('NO_FILE')          ? "--da-peaks ${da_peaks}"           : "--da-peaks ''"
    def meta_arg  = !metadata.name.startsWith('NO_FILE')          ? "--metadata ${metadata}"           : ""
    def pfms_arg  = pfms                                         ? "--pfms '${pfms}'"                 : ""
    def chunk_sz  = params.chromvar.chunk_size ?: 10000

    """
    gpu_chromvar_nf.py \\
      --peak-matrix ${peak_matrix} \\
      ${meta_arg} \\
      --cache-dir '${cache_dir}' \\
      ${pfms_arg} \\
      --genome '${genome}' \\
      --out-prefix chromvar_matrix \\
      --chunk-size ${chunk_sz} \\
      ${da_arg}
    """
}
