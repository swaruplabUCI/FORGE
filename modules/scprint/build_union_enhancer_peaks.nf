// modules/scprint/build_union_enhancer_peaks.nf
//
// BUILD_UNION_ENHANCER_PEAKS — methodological prerequisite for cis-rewiring.
// Concatenates per-condition EXTRACT_CCAN_ENHANCERS BEDs, sorts, and merges
// overlapping intervals (awk-based; no bedtools dependency) → a single union
// peakset covering BOTH conditions. Subsequent motif scan against this BED
// captures motifs at treatment-only enhancers that a control-only scan would
// miss (the caveat that drove the SREBF1/Cers2 0/12 conservative call).

process BUILD_UNION_ENHANCER_PEAKS {

    tag "build_union_enhancer_peaks"
    label 'process_low'
    publishDir "${params.outdir}/enhancer_footprinting/ccan_enhancers",
               mode: 'copy'
    errorStrategy 'terminate'

    input:
    path peaks_ctrl       // ccan_enhancer_peaks_<ctrl>.bed.gz
    path peaks_trt        // ccan_enhancer_peaks_<trt>.bed.gz

    output:
    path "ccan_enhancer_peaks_union.bed.gz",   emit: peaks_union
    path "ccan_enhancer_peaks_union_summary.txt", emit: summary

    script:
    """
    set -eo pipefail
    zcat '${peaks_ctrl}' '${peaks_trt}' \\
        | awk 'BEGIN{OFS="\\t"} \$1!~/^#/ && NF>=3 {print \$1,\$2,\$3}' \\
        | sort -k1,1 -k2,2n \\
        | awk 'BEGIN{OFS="\\t"}
            NR==1 {c=\$1; s=\$2; e=\$3; next}
            \$1==c && \$2<=e { if (\$3>e) e=\$3; next }
            { print c, s, e; c=\$1; s=\$2; e=\$3 }
            END { if (NR>0) print c, s, e }' \\
        | gzip -c > ccan_enhancer_peaks_union.bed.gz

    n_ctrl=\$(zcat '${peaks_ctrl}' | awk 'NF>=3 && \$1!~/^#/' | wc -l)
    n_trt=\$(zcat  '${peaks_trt}'  | awk 'NF>=3 && \$1!~/^#/' | wc -l)
    n_union=\$(zcat ccan_enhancer_peaks_union.bed.gz | wc -l)
    {
        echo "ctrl input intervals : \$n_ctrl"
        echo "trt  input intervals : \$n_trt"
        echo "merged union intervals: \$n_union"
    } > ccan_enhancer_peaks_union_summary.txt
    """

    stub:
    """
    echo -e "chr1\\t100\\t200" | gzip -c > ccan_enhancer_peaks_union.bed.gz
    echo "stub" > ccan_enhancer_peaks_union_summary.txt
    """
}
