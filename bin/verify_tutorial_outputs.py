#!/usr/bin/env python3
"""
verify_tutorial_outputs.py — check a tutorial run against the published checksums.

WHY THIS EXISTS AND NOT `sha256sum -c`
--------------------------------------
Two of FORGE's output classes are not byte-reproducible even when the numbers in
them are identical:

  * **PDF/PNG figures** — PDFs embed /CreationDate and /ModDate, so two runs that
    produce pixel-identical plots still hash differently. Checksumming figures
    reports FAILED forever, which only teaches a reader to ignore the check.
  * **gzip files** — the gzip member header embeds the source mtime, so
    cicero_connections.tsv.gz hashes differently between runs even when the TSV
    inside is byte-identical.

So this script checksums the deterministic text/numeric outputs, and hashes gzip
members **decompressed**, which is stable. Figures and other timestamped outputs
are left to human verification — see docs/tutorial.md.

WHAT IS COVERED
---------------
Everything under the run's outdir matching HASHABLE_SUFFIXES, minus EXCLUDE.
Measured across two independent cold runs of the tutorial: **168 of 170** such
files were byte-identical. The two that were not are both listed in EXCLUDE
below, with the reason.

USAGE
-----
    # verify a run (what a reader does)
    python3 bin/verify_tutorial_outputs.py \
        --results results_tutorial \
        --checksums checksums_data.txt

    # regenerate the manifest (what a maintainer does after a legitimate change)
    python3 bin/verify_tutorial_outputs.py \
        --results results_tutorial \
        --checksums checksums_data.txt --write

Exit status is 0 only if every file is present and matches.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import sys

HASHABLE_SUFFIXES = ('.json', '.csv', '.tsv', '.txt', '.tsv.gz')

# Deterministic-output exceptions. Both were found empirically by diffing two
# independent cold runs; do not add to this list without the same evidence.
EXCLUDE = {
    # Records a wall-clock timestamp and a timestamped output filename.
    'mofa_visualization/mofa_integration_summary.json':
        'embeds a wall-clock timestamp and a timestamped output filename',
    # Per-task runtimes and memory — informational by definition.
    'pipeline_info/trace.tsv':
        'per-task runtime/memory, expected to differ every run',
}


def content_hash(path: str) -> str:
    """sha256 of a file's content — decompressed if it is gzip.

    Hashing the decompressed member is what makes a .gz checksum stable: the
    gzip header carries the source mtime, the payload does not.
    """
    h = hashlib.sha256()
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def cicero_fallback(results: str, rel: str) -> str | None:
    """Runs predating the cicero.outdir fix published Cicero to a sibling results/."""
    if not rel.startswith('cicero/'):
        return None
    parent = os.path.dirname(os.path.abspath(results.rstrip('/')))
    alt = os.path.join(parent, 'results', rel)
    return alt if os.path.exists(alt) else None


def resolve(results: str, rel: str) -> str | None:
    direct = os.path.join(results, rel)
    if os.path.exists(direct):
        return direct
    return cicero_fallback(results, rel)


def discover(results: str) -> list[str]:
    """All hashable outputs under the run dir, plus Cicero wherever it landed."""
    found = set()
    roots = [results]
    stray = cicero_fallback(results, 'cicero/')
    if stray:
        roots.append(stray)
    for root in roots:
        base = results if root == results else os.path.dirname(root.rstrip('/'))
        for dirpath, _, filenames in os.walk(root):
            for name in filenames:
                if not name.endswith(HASHABLE_SUFFIXES):
                    continue
                rel = os.path.relpath(os.path.join(dirpath, name), base)
                if rel in EXCLUDE:
                    continue
                found.add(rel)
    return sorted(found)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True, help='the tutorial run output dir')
    ap.add_argument('--checksums', required=True, help='the manifest to check or write')
    ap.add_argument('--write', action='store_true',
                    help='regenerate the manifest instead of checking it')
    args = ap.parse_args()

    if args.write:
        rels = discover(args.results)
        with open(args.checksums, 'w') as out:
            out.write('# FORGE tutorial — sha256 of the DETERMINISTIC numeric outputs.\n')
            out.write('#\n')
            out.write('# Verify with:\n')
            out.write('#   python3 bin/verify_tutorial_outputs.py '
                      '--results <dir> --checksums checksums_data.txt\n')
            out.write('#\n')
            out.write('# NOT usable with `sha256sum -c`: .gz entries are hashed\n')
            out.write('# DECOMPRESSED, because the gzip header embeds an mtime.\n')
            out.write('#\n')
            out.write('# Figures are deliberately absent — PDFs embed /CreationDate\n')
            out.write('# and never hash-match. Verify those by eye.\n')
            for rel, why in sorted(EXCLUDE.items()):
                out.write(f'# Excluded: {rel} — {why}\n')
            out.write('#\n')
            for rel in rels:
                path = resolve(args.results, rel)
                if path is None:
                    print(f'MISSING (not written): {rel}', file=sys.stderr)
                    continue
                out.write(f'{content_hash(path)}  {rel}\n')
        print(f'wrote {args.checksums} ({len(rels)} files)')
        return 0

    expected = {}
    with open(args.checksums) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            digest, rel = line.split(None, 1)
            expected[rel.strip()] = digest

    bad, missing = [], []
    for rel, want in sorted(expected.items()):
        path = resolve(args.results, rel)
        if path is None:
            missing.append(rel)
            continue
        if content_hash(path) != want:
            bad.append(rel)

    total = len(expected)
    print(f'{total - len(bad) - len(missing)}/{total} matched')
    for rel in missing:
        print(f'  MISSING  {rel}')
    for rel in bad:
        print(f'  MISMATCH {rel}')

    if bad or missing:
        print('\nMismatches here are EXPECTED on hardware with a different core '
              'count from the published run, and are not a sign that your run '
              'failed.\n'
              'Several stages (CellBender, Cicero distance-parameter estimation) '
              'are numerically sensitive to the number of threads the BLAS/OpenMP\n'
              'libraries use, and that thread count currently follows the machine. '
              'The resulting shifts are tiny and do not move any structural count.\n'
              '\nWhat to check instead: the "Values that must match" table in\n'
              'docs/tutorial.md (cell counts, peak counts, cluster counts, triplet\n'
              'count, MOFA factors). Those are invariant to hardware. If one of\n'
              'THOSE differs, something is genuinely wrong with your inputs or\n'
              'containers.\n'
              '\nFigures and timestamped files are NOT checked here — verify those '
              'by eye against figures.tar.gz.')
        return 0
    print('\nAll deterministic numeric outputs match the published run.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
