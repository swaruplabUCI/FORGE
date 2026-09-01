# Building mm39 cisTarget databases by liftOver

Used by: `SCENICPLUS_RUN` (and `PYCISTARGET_*`) for any mouse dataset aligned to **GRCm39 / mm39**.

!!! warning "There is no upstream mm39 cisTarget database"
    The [aertslab cisTarget resources](https://resources.aertslab.org/cistarget/databases/)
    publish `v10nr_clust` databases for **hg38** and **mm10** only. If your mouse data is
    aligned to GRCm39 (as the BD Rhapsody datasets in this project are), the mm10 database
    coordinates will not match your consensus peaks and SCENIC+ will silently return
    near-empty motif enrichment.

    The two options are: (a) realign/liftOver your **peaks** to mm10, or (b) liftOver the
    **database region coordinates** to mm39, which is what this page describes. Option (b)
    is preferable — it is done once per database, not once per dataset, and it leaves your
    analysis coordinates untouched.

## What this produces

| Output | Size |
|--------|------|
| `mm39_region_based.scores.feather` | 10.2 GiB |
| `mm39_region_based.rankings.feather` | 21.0 GiB |

Reference run statistics (mm10 `v10nr_clust` region-based database → mm39):

| Metric | Value |
|--------|-------|
| Input regions (feather columns) | 1,110,655 |
| Lifted 1:1 to mm39 | 1,110,637 (99.998%) |
| Unmapped (partially deleted / split in new) | 18 |
| Split regions dropped | 0 |
| Duplicate mm39 targets dropped | 0 |

Total cost is roughly 5 CPU-hours, dominated by step 4. Step 4 needs **256 GB RAM** — the
rankings feather is rewritten column-wise in memory.

## Why the coordinates live in the column names

A cisTarget region-based database is a feather table whose **columns are genomic regions**
(`chr1:3050000-3050300`) and whose rows are motifs, plus a trailing `motifs` column. So
lifting the database over means rewriting ~1.1 million column *names* and dropping the
columns that fail to lift — the matrix values themselves are unchanged. That is why this
works at all, and why it is cheap relative to rebuilding a database from scratch with
`create_cisTarget_databases`.

Regions that lift to **multiple** mm39 intervals are dropped rather than duplicated: a
motif score is not meaningfully assignable to two disjoint targets.

---

## Prerequisites

- The mm10 region-based database and scores, already downloaded (see
  [Reference files](references.md#scenic-cistarget-rankings-and-scores)):
  `mm10_region_based.rankings.feather`, `mm10_region_based.scores.feather`
- ~65 GB free scratch (inputs stay in place; outputs are written alongside)
- A conda env with `pyarrow>=14` and `pandas>=2`
- The UCSC `liftOver` binary and the `mm10ToMm39.over.chain.gz` chain file

Set a working directory once:

```bash
export LO=/path/to/UCSC_Liftover
export RESOURCE_DIR=/path/to/ref/scenic_plus_resources/mouse
mkdir -p "$LO"/{bin,logs}
```

---

## Step 1 — environment, liftOver binary, chain file

```bash
eval "$(conda shell.bash hook)"
conda create -y -p "$LO/env" -c conda-forge "pyarrow>=14" "pandas>=2"

wget -q -O "$LO/bin/liftOver" \
  https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x "$LO/bin/liftOver"

wget -q -P "$LO/" \
  https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz
```

!!! note "Use conda, not mamba, on DFS filesystems"
    Solver I/O on a parallel filesystem is slow enough that mamba's speed advantage
    disappears and its failure modes get worse. Give the job generous walltime (6 h).

The `liftOver` binary is unversioned — record the **retrieval date** for your methods
section, since UCSC replaces the binary in place.

## Step 2 — extract region coordinates to BED4

Reads only the feather **schema**, so memory use is trivial. The 4th BED column carries
the original mm10 region string, which is what lets step 4 map old name → new name.

```python title="extract_bed.py"
#!/usr/bin/env python3
"""Extract genomic region column names from a cisTarget feather file to BED4."""
import argparse
import pyarrow.ipc as ipc

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--feather", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    schema = ipc.open_file(args.feather).schema
    n_regions = n_skipped = 0

    with open(args.output, "w") as f:
        for i in range(len(schema)):
            name = schema.field(i).name
            if not name.startswith("chr"):
                n_skipped += 1          # the trailing 'motifs' column
                continue
            try:
                chrom, coords = name.split(":")
                start, end = coords.split("-")
                f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
                n_regions += 1
            except ValueError:
                print(f"  WARNING: could not parse column '{name}', skipping")
                n_skipped += 1

    print(f"Extracted {n_regions} regions ({n_skipped} non-region columns skipped)")

if __name__ == "__main__":
    main()
```

```bash
conda activate "$LO/env"
python3 "$LO/extract_bed.py" \
  --feather "$RESOURCE_DIR/mm10_region_based.scores.feather" \
  --output  "$LO/mm10_regions.bed"
wc -l "$LO/mm10_regions.bed"      # expect 1110655
```

Extract from the **scores** file (smaller) — scores and rankings share an identical
column set, so one BED serves both.

## Step 3 — liftOver mm10 → mm39

```bash
"$LO/bin/liftOver" \
  "$LO/mm10_regions.bed" \
  "$LO/mm10ToMm39.over.chain.gz" \
  "$LO/mm39_regions.bed" \
  "$LO/mm10_unmapped.bed" \
  -minMatch=0.95
```

`-minMatch=0.95` requires 95% of the interval to map. This is deliberately strict: a
cisTarget region is a ~300 bp regulatory window, and a partially lifted window no longer
corresponds to the sequence the motif scores were computed on.

Check the rate before continuing — anything below ~99% means the wrong chain file:

```bash
MAPPED=$(wc -l < "$LO/mm39_regions.bed")
UNMAPPED=$(grep -c '^[^#]' "$LO/mm10_unmapped.bed" || echo 0)
echo "mapped=$MAPPED unmapped=$UNMAPPED"
head "$LO/mm10_unmapped.bed"      # reasons are '#'-prefixed
```

Expected failure reasons are `#Partially deleted in new` and `#Split in new`, 18 regions total.

## Step 4 — rebuild the feather files (256 GB)

Builds the mm10→mm39 name mapping, then selects, renames, coordinate-sorts and rewrites
both feathers. Regions that split, failed, or collide on the same mm39 target are dropped.

```python title="rebuild_feather.py"
#!/usr/bin/env python3
"""Rebuild SCENIC+ cisTarget feather databases with mm39 coordinates."""
import argparse
import os
import pyarrow.feather as pf
import pyarrow.ipc as ipc


def load_liftover_mapping(bed_mapped_path, bed_unmapped_path):
    """Build mm10_region -> mm39_region mapping, excluding splits and collisions."""
    name_counts = {}
    with open(bed_mapped_path) as f:
        for line in f:
            mm10_name = line.strip().split("\t")[3]
            name_counts[mm10_name] = name_counts.get(mm10_name, 0) + 1

    mapping, mm39_seen, skipped_split = {}, {}, 0
    with open(bed_mapped_path) as f:
        for line in f:
            chrom, start, end, mm10_name = line.strip().split("\t")[:4]
            if name_counts[mm10_name] > 1:        # lifted to >1 interval
                skipped_split += 1
                continue
            mm39_name = f"{chrom}:{start}-{end}"
            if mm39_name in mm39_seen:            # two mm10 regions collapsed onto one
                print(f"  WARNING: duplicate mm39 target {mm39_name} "
                      f"(from {mm39_seen[mm39_name]} and {mm10_name}), keeping first")
                continue
            mm39_seen[mm39_name] = mm10_name
            mapping[mm10_name] = mm39_name

    unmapped = sum(1 for line in open(bed_unmapped_path) if not line.startswith("#"))
    unique_split = len([k for k, v in name_counts.items() if v > 1])

    print("Liftover mapping summary:")
    print(f"  Successfully mapped (1:1): {len(mapping)}")
    print(f"  Split regions (dropped):   {unique_split}")
    print(f"  Unmapped regions:          {unmapped}")
    return mapping


def sort_key(region_str):
    chrom, coords = region_str.split(":")
    start, end = coords.split("-")
    c = chrom.replace("chr", "")
    try:
        c = int(c)
    except ValueError:
        c = {"X": 100, "Y": 101, "M": 102, "MT": 102}.get(c, 200)
    return (c, int(start), int(end))


def rebuild_feather(input_path, output_path, mapping):
    print(f"\nProcessing: {input_path}\n  Output:   {output_path}")
    schema = ipc.open_file(input_path).schema

    keep, new_names, region_new = [], [], []
    for i in range(len(schema)):
        col = schema.field(i).name
        if col.startswith("chr"):
            if col in mapping:
                keep.append(col)
                new_names.append(mapping[col])
                region_new.append(mapping[col])
        else:
            keep.append(col)
            new_names.append(col)

    print(f"  Original columns: {len(schema)}")
    print(f"  Kept columns:     {len(keep)}")
    print(f"  Dropped columns:  {len(schema) - len(keep)}")

    table = pf.read_table(input_path, columns=keep)
    table = table.rename_columns(new_names)

    non_region = [n for n in new_names if not n.startswith("chr")]
    region_new.sort(key=sort_key)
    table = table.select(region_new + non_region)

    pf.write_feather(table, output_path, compression="lz4")
    print(f"  Done. {os.path.getsize(output_path) / 1024**3:.1f} GiB, "
          f"{len(region_new) + len(non_region)} columns")
    del table


def main():
    p = argparse.ArgumentParser()
    for a in ("--bed-mapped", "--bed-unmapped", "--scores-in",
              "--rankings-in", "--scores-out", "--rankings-out"):
        p.add_argument(a, required=True)
    args = p.parse_args()

    mapping = load_liftover_mapping(args.bed_mapped, args.bed_unmapped)
    rebuild_feather(args.scores_in, args.scores_out, mapping)      # smaller first
    rebuild_feather(args.rankings_in, args.rankings_out, mapping)


if __name__ == "__main__":
    main()
```

```bash
python3 "$LO/rebuild_feather.py" \
  --bed-mapped   "$LO/mm39_regions.bed" \
  --bed-unmapped "$LO/mm10_unmapped.bed" \
  --scores-in    "$RESOURCE_DIR/mm10_region_based.scores.feather" \
  --rankings-in  "$RESOURCE_DIR/mm10_region_based.rankings.feather" \
  --scores-out   "$RESOURCE_DIR/mm39_region_based.scores.feather" \
  --rankings-out "$RESOURCE_DIR/mm39_region_based.rankings.feather"
```

Scores are processed first on purpose: it is the smaller file, so a schema or memory
problem surfaces in minutes rather than after an hour on the rankings.

## Step 5 — verify

Do not skip this. A database whose columns are silently still mm10 produces a SCENIC+ run
that completes successfully and returns meaningless enrichment.

```bash
python3 - <<'PY'
import pyarrow.ipc as ipc, pyarrow.feather as pf
RESOURCE_DIR = "/path/to/ref/scenic_plus_resources/mouse"

for suffix in ("scores", "rankings"):
    mm10 = f"{RESOURCE_DIR}/mm10_region_based.{suffix}.feather"
    mm39 = f"{RESOURCE_DIR}/mm39_region_based.{suffix}.feather"
    r10, r39 = ipc.open_file(mm10), ipc.open_file(mm39)
    print(f"\n=== {suffix} ===")
    print(f"  mm10 columns {len(r10.schema)} -> mm39 columns {len(r39.schema)} "
          f"(dropped {len(r10.schema) - len(r39.schema)})")

    names = [r39.schema.field(i).name for i in range(len(r39.schema))]
    regions = [n for n in names if n.startswith("chr")]
    print(f"  region columns: {len(regions)}  non-region: {len(names) - len(regions)}")
    print(f"  chromosomes: {sorted({n.split(':')[0] for n in regions})}")
    print(f"  motifs column last: {names[-1] == 'motifs'}")
    print(f"  sorted from chr1:  {names[0].startswith('chr1:')}")

    t = pf.read_table(mm39, columns=names[:3])
    print(f"  data slice rows: {t.num_rows}; sample: {t.column(0).to_pylist()[:3]}")
PY
```

Expect: 18 columns dropped, `motifs` still the final column, first column on `chr1`, and a
non-empty data slice.

### Spot-check against a known locus

Schema checks cannot tell you the coordinates are *correct*, only that they changed. Confirm
one known mm39 promoter appears in the new column set — pick a gene from your own analysis,
take its mm39 TSS from the GENCODE vM37 GTF, and verify a database region overlaps it:

```bash
python3 - <<'PY'
import pyarrow.ipc as ipc
names = [f.name for f in ipc.open_file(
    "/path/to/mm39_region_based.scores.feather").schema]
CHROM, TSS = "chr11", 69_588_000        # replace with your locus, mm39 coordinates
hits = [n for n in names if n.startswith(CHROM + ":")
        and int(n.split(":")[1].split("-")[0]) <= TSS <= int(n.split("-")[1])]
print(hits or "NO OVERLAP — check the chain file direction")
PY
```

---

## Wiring it into FORGE

Point the dataset config at the new files:

```groovy
params {
    scenicplus {
        rankings_db = "/path/to/ref/scenic_plus_resources/mouse/mm39_region_based.rankings.feather"
        scores_db   = "/path/to/ref/scenic_plus_resources/mouse/mm39_region_based.scores.feather"
        motif_annot = "/path/to/ref/scenic_plus_resources/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
    }
}
```

The motif annotation table (`motifs-v10nr_clust-nr.mgi-*.tbl`) is **not** coordinate-based
and needs no liftOver — use the mm10/MGI table unchanged.

!!! tip "Mouse TF symbol casing"
    The MGI motif annotation table is the right choice for mouse. If you substitute the
    HGNC table, TF names arrive in human uppercase and cistrome construction silently
    fails to match mouse gene symbols.

## Reporting this in a manuscript

Record: the chain file (`mm10ToMm39.over.chain.gz`, UCSC goldenPath), the `liftOver`
binary retrieval date, `-minMatch=0.95`, the source database version (`v10nr_clust`,
region-based), and the retained-region count (1,110,637 of 1,110,655). The derived
databases are a resource in their own right — deposit them, or cite this page as the
build recipe.
