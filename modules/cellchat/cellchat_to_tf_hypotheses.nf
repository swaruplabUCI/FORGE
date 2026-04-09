// modules/cellchat/cellchat_to_tf_hypotheses.nf
//
// CELLCHAT_TO_TF_HYPOTHESES — Recipe Steps C2+C3
// Map CellChat signaling pathways to downstream TFs, validate with
// chromVAR deviations in receiver cells (Mann-Whitney U test).

process CELLCHAT_TO_TF_HYPOTHESES {

    tag "cellchat_tf_hypotheses"
    label 'process_medium'
    publishDir "${params.outdir}/cellchat/tf_hypotheses", mode: 'copy'

    input:
    path cellchat_csv         // CellChat results CSV from RUN_CELLCHAT
    path chromvar_dev         // chromvar_matrix_deviations.h5ad from GPU_CHROMVAR
    path pathway_to_tfs_json  // data/pathway_to_tfs.json

    output:
    path "receiver_tf_hypotheses.tsv",       emit: hypotheses
    path "validated_receiver_tf_pairs.tsv",  emit: validated_pairs
    path "unmapped_pathways.txt",            emit: unmapped
    path "hypothesis_summary.txt",           emit: summary

    script:
    def annotation_key = params.atac.marker_file ? 'cell_type' : (params.atac.annotation_method == 'scatanno' ? 'cell_type_prediction' : 'celltypist_prediction')
    """
    python ${projectDir}/bin/cellchat_to_tf_hypotheses.py \\
        --cellchat-csv '${cellchat_csv}' \\
        --chromvar-dev '${chromvar_dev}' \\
        --pathway-to-tfs '${pathway_to_tfs_json}' \\
        --cell-type-key ${annotation_key} \\
        --pval-threshold 0.05 \\
        --outdir .
    """
}
