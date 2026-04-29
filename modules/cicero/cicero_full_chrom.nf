// modules/cicero/cicero_full_chrom.nf
//
// Cicero Phase 3 — step 2 of 3: per-chromosome generate_cicero_models +
// assemble_connections using the pinned global distance_parameter from step 1.
// Fan-out across chromosomes provides the ~4x speedup validated by FORGE
// (rho=1.0 vs monolithic baseline).
// Ported from src_FORGE 2026-04-23; used verbatim for plain + stratified legs
// via include-as aliasing in main.nf.

process CICERO_FULL_CHROM {
    tag "${chrom}"
    label 'process_medium'
    errorStrategy 'terminate'

    input:
    tuple val(chrom), path(cicero_cds_rds), path(gene_ann_rds), path(dp_file)

    output:
    tuple val(chrom), path("conns_${chrom}.tsv.gz"), emit: chrom_conns

    script:
    """
    export HOME=/tmp/container_home
    mkdir -p /tmp/container_home
    export R_LIBS_USER=""
    DP=\$(cat ${dp_file})
    cicero_full_chrom.R \\
      --chrom "${chrom}" \\
      --cicero_cds "${cicero_cds_rds}" \\
      --gene_annotation "${gene_ann_rds}" \\
      --dp "\${DP}" \\
      --window_bp 500000 \\
      --outdir "."
    """
}
