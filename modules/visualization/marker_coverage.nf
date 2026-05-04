// modules/visualization/marker_coverage.nf
// Render pyGenomeTracks coverage panels at canonical marker promoters
// (Shi Fig 1E equivalent) using the broad-class bigWigs.

nextflow.enable.dsl=2

process MARKER_COVERAGE_TRACKS {

    tag "marker_coverage_tracks"
    label 'process_medium'
    errorStrategy 'terminate'
    publishDir "${params.outdir}/shi_equivalent/1E_marker_tracks", mode: 'copy'

    input:
    path bigwig_manifest      // manifest.json from EXPORT_ATAC_BIGWIGS (broad ct)
    path "bigwigs/*"          // collected .bw files staged into bigwigs/
    path gtf

    output:
    path "*_coverage.{pdf,png}",            emit: tracks
    path "marker_tracks.ini",               emit: ini
    path "marker_tracks_manifest.json",     emit: manifest

    script:
    def markers = params.shi_figures?.markers_json ?: ''
    def window  = params.shi_figures?.marker_window ?: 5000
    def width   = params.shi_figures?.marker_width  ?: 18
    def markers_arg = markers ?
        "--markers '${markers}'" :
        "--markers 'Astrocyte=Gfap,Excitatory=Slc17a7,Inhibitory=Gad2,Microglia=Csf1r,Oligodendrocyte=Mobp,OPC=Pdgfra,Vascular=Cldn5'"
    def keep_csv = (params.shi_figures?.marker_lineages instanceof List) ?
        params.shi_figures.marker_lineages.join(',') :
        (params.shi_figures?.marker_lineages ?:
            'Astro-Epen,IT-ET Glut,CTX-MGE GABA,Immune,OPC-Oligo,Vascular')
    """
    # Rewrite the staged manifest:
    #   - point at staged work-dir bigwig basenames
    #   - subset to the lineages requested via marker_lineages so the figure stays legible.
    python <<'PYEOF'
import json, os
m = json.load(open("${bigwig_manifest}"))
bws = m.get("bigwigs", m) if isinstance(m, dict) else m
keep = {s.strip() for s in "${keep_csv}".split(",") if s.strip()}
fixed = {}
for k, v in bws.items():
    if keep and k not in keep:
        continue
    fixed[k] = os.path.join("bigwigs", os.path.basename(v))
json.dump(fixed, open("bigwigs_local.json", "w"), indent=2)
print("[marker_coverage] using lineages:", list(fixed.keys()))
PYEOF

    python ${projectDir}/bin/plot_marker_coverage_tracks.py \\
        --bigwig-manifest bigwigs_local.json \\
        ${markers_arg} \\
        --gtf '${gtf}' \\
        --window ${window} \\
        --width ${width} \\
        --outdir .
    """
}
