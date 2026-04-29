// modules/scprint/manifest.nf
//
// Build a scPrinter manifest TSV mapping:
//   sample_id  fragment_file  barcode_file
//
// The tricky part: barcode sample IDs may differ from metadata sample_id.
// We support a regex mapping:
//   params.scprinter.sample_id_regex
// which must contain a named capture group "id" that returns the metadata sample_id.
//
// Example:
//   barcode file: sample1_batch1_barcodes.txt
//   regex: '^(?<id>.+?)_batch\\d+$'
//   => id = sample1
//
process SCPRINTER_BUILD_MANIFEST {
    label 'process_small'
    publishDir "${params.outdir}/scprinter", mode: 'copy'

    input:
    path metadata_csv
    path barcode_files
    val fragments_dir
    val sample_id_regex

    output:
    path "scprinter_manifest.tsv", emit: manifest

    script:
    """
    python3 - << 'PY'
    import csv
    import re
    from pathlib import Path
    import sys

    metadata_csv = Path("${metadata_csv}")
    fragments_dir = Path("${fragments_dir}")
    sample_id_regex = "${sample_id_regex}".strip()

    # Read metadata -> sample_id => fragment_file (resolved to absolute if needed)
    meta = {}
    with metadata_csv.open() as f:
        reader = csv.DictReader(f)
        if 'sample_id' not in reader.fieldnames or 'fragment_file' not in reader.fieldnames:
            raise SystemExit(
                f"ERROR: metadata must contain columns sample_id and fragment_file; found: {reader.fieldnames}"
            )
        for row in reader:
            sid = (row.get('sample_id') or '').strip()
            frag = (row.get('fragment_file') or '').strip()
            if not sid or not frag:
                continue
            frag_path = Path(frag)
            if not frag_path.is_absolute():
                frag_path = fragments_dir / frag_path
            meta[sid] = str(frag_path)

    # In work dir, Nextflow stages barcode_files. Discover them here.
    barcode_paths = sorted(Path('.').glob('*_barcodes.txt'))
    barcode_paths = [p for p in barcode_paths if not (p.name.startswith('gene_matrix') or p.name.startswith('peak_matrix'))]

    if not barcode_paths:
        raise SystemExit("ERROR: No *_barcodes.txt files found in work directory")

    # Optional regex mapping from barcode sid -> metadata sid
    rx = None
    if sample_id_regex:
        rx = re.compile(sample_id_regex)

    def map_barcode_sid_to_meta_sid(barcode_sid: str) -> str:
        if rx is None:
            return barcode_sid
        m = rx.match(barcode_sid)
        if not m:
            raise SystemExit(
                f"ERROR: sample_id_regex did not match barcode sample id '{barcode_sid}'. "
                f"Regex: {sample_id_regex}"
            )
        if 'id' not in m.groupdict():
            raise SystemExit(
                "ERROR: sample_id_regex must include a named capture group 'id'. "
                f"Regex: {sample_id_regex}"
            )
        return m.group('id')

    missing = []
    rows = []
    for bc in barcode_paths:
        barcode_sid = bc.name.replace('_barcodes.txt', '')
        meta_sid = map_barcode_sid_to_meta_sid(barcode_sid)

        if meta_sid not in meta:
            missing.append((barcode_sid, meta_sid))
            continue

        rows.append((meta_sid, meta[meta_sid], str(bc.resolve())))

    if missing:
        msg = "ERROR: Barcode sample IDs could not be mapped to metadata sample_id:\\n"
        for barcode_sid, meta_sid in missing:
            msg += f"  barcode='{barcode_sid}' -> meta='{meta_sid}' (NOT FOUND IN METADATA)\\n"
        raise SystemExit(msg)

    # Write manifest TSV with header
    out = Path("scprinter_manifest.tsv")
    with out.open("w") as f:
        f.write("sample_id\\tfragment_file\\tbarcode_file\\n")
        for sid, frag, bc in rows:
            f.write(f"{sid}\\t{frag}\\t{bc}\\n")

    print("Wrote scprinter_manifest.tsv with", len(rows), "samples")
    PY

    echo "=== scprinter_manifest.tsv (head) ==="
    head -n 10 scprinter_manifest.tsv
    echo "====================================="
    """
}
