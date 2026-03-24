// modules/multiome/scenicplus_run.nf

nextflow.enable.dsl=2

process SCENICPLUS_RUN {

    tag { "scenicplus_snakemake" }

    publishDir "${params.outdir}/scenicplus", mode: 'copy', overwrite: true

    // Conda is handled in nextflow.config

    input:
    // pycisTopic object with regions/cell metadata
    path cistopic_obj
    // RNA AnnData (integrated & QC'd)
    path rna_anndata
    // Region set directory (with subfolders Topics_otsu, Topics_top_3k, DARs_cell_type)
    path region_sets_dir
    // cisTarget ranking database (.rankings.feather)
    val ctx_db_rankings
    // cisTarget scores database (.scores.feather)
    val ctx_db_scores
    // Motif annotation table (.tbl)
    val motif_annotations_tbl
    // Barcode transform lambda string
    val bc_transform

    output:
    path "outs/ACC_GEX.h5mu", emit: acc_gex_mudata
    path "outs/ctx_results.hdf5", emit: ctx_h5
    path "outs/ctx_results.html", emit: ctx_html
    path "outs/dem_results.hdf5", emit: dem_h5
    path "outs/dem_results.html", emit: dem_html
    path "outs/cistromes_direct.h5ad", emit: cistromes_direct
    path "outs/cistromes_extended.h5ad", emit: cistromes_extended
    path "outs/tf_names.txt", emit: tf_names
    path "outs/genome_annotation.tsv", emit: genome_annotation
    path "outs/chromsizes.tsv", emit: chromsizes
    path "outs/search_space.tsv", emit: search_space
    path "outs/tf_to_gene_adj.tsv", emit: tf2g
    path "outs/region_to_gene_adj.tsv", emit: r2g
    path "outs/eRegulon_direct.tsv", emit: ereg_direct
    path "outs/eRegulons_extended.tsv", emit: ereg_extended
    path "outs/AUCell_direct.h5mu", emit: aucell_direct
    path "outs/AUCell_extended.h5mu", emit: aucell_extended
    path "outs/scplusmdata.h5mu", emit: scplus_mudata
    path "scplus_pipeline/Snakemake/.snakemake/log/*", optional: true, emit: snakemake_logs

    script:
    """
    export PYTHONPATH="${projectDir}/tools/scenicplus-main/src:\${PYTHONPATH:-}"

    mkdir -p scplus_pipeline

    python ${projectDir}/bin/run_scenicplus_snakemake.py \\
        --cistopic ${cistopic_obj} \\
        --rna ${rna_anndata} \\
        --region-sets ${region_sets_dir} \\
        --ctx-rankings ${ctx_db_rankings} \\
        --ctx-scores ${ctx_db_scores} \\
        --motif-annotations ${motif_annotations_tbl} \\
        --bc-transform "${bc_transform}" \\
        --gtf ${params.scenicplus.gtf} \\
        --fai ${params.scenicplus.fai} \\
        --outdir .

    """
}
