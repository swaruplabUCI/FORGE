// modules/cellchat/signal_chain_correlation.nf
//
// SIGNAL_CHAIN_CORRELATION — Recipe Step C6
// Full-chain correlation: ligand (sender) ~ footprint (receiver) ~
// target gene (receiver).  Classifies interactions into evidence tiers.

process SIGNAL_CHAIN_CORRELATION {

    tag "signal_chain_correlation"
    label 'process_medium'
    publishDir "${params.outdir}/cellchat/signal_chain", mode: 'copy'

    input:
    path footprint_h5ads      // Collected enhancer footprint h5ads
    path rna_h5ad             // Integrated RNA h5ad
    path signaling_metadata   // signaling_metadata.tsv

    output:
    path "signal_chain_correlations.tsv",   emit: correlations
    path "evidence_tiers.tsv",              emit: evidence_tiers
    path "chain_plots/",                    emit: plots, optional: true
    path "chain_summary.txt",              emit: summary

    script:
    """
    python ${projectDir}/bin/signal_chain_correlation.py \\
        --footprint-h5ads ${footprint_h5ads} \\
        --rna-h5ad '${rna_h5ad}' \\
        --signaling-metadata '${signaling_metadata}' \\
        --outdir .
    """
}
