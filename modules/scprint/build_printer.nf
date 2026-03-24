// modules/scprint/build_printer.nf
//
// Build the scPrinter printer object from fragment files + barcode whitelists.
//
// FIX-37d: Accept fragment files as an explicit input channel instead of
//          hardcoding BD-specific demux directory paths (june/july/nov).
//          This makes the module work for any dataset (BD, 10x PBMC, etc.).
process SCPRINTER_BUILD_PRINTER {
    label 'process_high'
    publishDir "${params.outdir}/scprinter", mode: 'copy'

    input:
    // All barcode files produced by SCPRINTER_BARCODES
    path barcode_files
    // Fragment files resolved from the manifest (same as ATAC workflows)
    path fragment_files_input

    output:
    // Printer file written in the working directory with this name
    path "human_multiome_scprinter.h5ad", emit: printer

    script:
    // Extract sample names from barcode file names
    // Remove _barcodes.txt suffix to get sample names
    def sample_names = barcode_files.collect { file ->
        def name = file.name.replace('_barcodes.txt', '')
        // Skip non-sample barcode files
        if (name.contains('gene_matrix') || name.contains('peak_matrix')) {
            return null
        }
        return name
    }.findAll { it != null }

    // Fragment files are now passed as input — convert to list
    def frag_files_list
    if (fragment_files_input instanceof List) {
        frag_files_list = fragment_files_input
    } else {
        frag_files_list = [fragment_files_input]
    }

    // Filter barcode files to only sample-specific ones
    def filtered_barcode_files = barcode_files.findAll { file ->
        !file.name.contains('gene_matrix') && !file.name.contains('peak_matrix')
    }

    // Create lists for command
    def frag_list = frag_files_list.collect { "'${it}'" }.join(' ')
    def sample_list = sample_names.collect { "'${it}'" }.join(' ')
    def barcode_list = filtered_barcode_files.collect { "'${it}'" }.join(' ')

    // Debug output
    println "Fragment files: ${frag_files_list}"
    println "Sample names: ${sample_names}"
    println "Barcode files: ${filtered_barcode_files}"

    """
    # FIX-35c: Set SCPRINTER_DATA so pooch caches model files to the
    # persistent ref dir on DFS7 instead of ephemeral /tmp under --contain.
    export SCPRINTER_DATA='${params.scprinter.cache_dir}'

    # Debug: Check if fragment files exist
    for frag in ${frag_list}; do
        if [[ -f "\$frag" ]]; then
            echo "Found fragment file: \$frag"
        else
            echo "WARNING: Fragment file not found: \$frag"
        fi
    done

    build_scprinter_printer.py \\
      --cache-dir '${params.scprinter.cache_dir}' \\
      --output 'human_multiome_scprinter.h5ad' \\
      --genome '${params.scprinter.genome}' \\
      --fragment-files ${frag_list} \\
      --barcode-lists ${barcode_list} \\
      --sample-names ${sample_list} \\
      --n-jobs ${task.cpus}
    """
}
