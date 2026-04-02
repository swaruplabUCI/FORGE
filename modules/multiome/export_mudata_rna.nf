process EXPORT_MUDATA_RNA {
    label 'process_low'
    tag "export_mudata_rna"

    input:
    path mudata_file

    output:
    path "integrated_rna_for_dorc.h5ad"

    script:
    """
    export_mudata_rna.py \
        --mudata ${mudata_file} \
        --out integrated_rna_for_dorc.h5ad
    """
}
