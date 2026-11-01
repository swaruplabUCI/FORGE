#!/usr/bin/env python3
"""Regression test for BUILD_MUDATA's cell-type column precedence.

The canonical 'cell_type' column was missing from _CT_CANDIDATES, so an
on-ramped integrated h5ad -- which stamps cell_type from the user's own column
name (subclass_label on Allen and SEA-AD) rather than from a tool-specific one
-- resolved to nothing and the MuData's RNA labels ended as NaN. The native
path was unaffected only because plot_postscanvi.py leaves the tool-specific
column in place next to the canonical one.

Run: python3 tests/onramp_rds/test_build_mudata_celltype.py
"""
import os, re, sys, pandas as pd
_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "bin", "build_mudata_batched.py")
src = open(sys.argv[1] if len(sys.argv) > 1 else _DEFAULT).read()

# lift the real _CT_CANDIDATES and _resolve_ct_source out of the module,
# verbatim, so this tests the shipped code rather than a paraphrase of it.
start = src.index("    _CT_CANDIDATES = [")
end   = src.index("    # Lookup keyed by", start)
block = "\n".join(l[4:] if l.startswith("    ") else l
                  for l in src[start:end].splitlines())
ns = {}
exec(block, ns)
_CT_CANDIDATES = ns["_CT_CANDIDATES"]
_resolve = ns["_resolve_ct_source"]

def select(obs):
    """The selection loop as it appears at the scanvi_file site."""
    for candidate, source in _CT_CANDIDATES:
        if candidate in obs.columns:
            vals = obs[candidate].dropna().unique()
            if len(vals) > 1 or (len(vals) == 1 and vals[0] != 'Unknown'):
                return candidate, _resolve(obs, source)
    return None, None

ok = True
def chk(label, got, want):
    global ok
    good = got == want
    print(f"  {label:<56} {'PASS' if good else 'FAIL'}   got={got}")
    if not good: ok = False

print("=== native integrated h5ad (plot_postscanvi output) ===")
native = pd.DataFrame({
    'scanvi_prediction': ['Astro', 'Oligo'],
    'cell_type':         ['Astro', 'Oligo'],
    'cell_type_source':  ['scanvi', 'scanvi'],
})
chk("picks cell_type, provenance read from obs", select(native), ('cell_type', 'scanvi'))
chk("same LABELS as the old scanvi_prediction pick",
    list(native['cell_type']), list(native['scanvi_prediction']))

print("=== seurat_rds on-ramp h5ad (the bug) ===")
onramp = pd.DataFrame({
    'subclass_label':   ['L6 IT', 'Sst'],
    'cell_type':        ['L6 IT', 'Sst'],
    'cell_type_source': ['seurat_onramp', 'seurat_onramp'],
})
chk("now resolves (was None -> NaN labels)", select(onramp), ('cell_type', 'seurat_onramp'))

print("=== on-ramp carrying a STALE tool column ===")
stale = pd.DataFrame({
    'subclass_label':    ['L6 IT', 'Sst'],
    'scanvi_prediction': ['WRONG_A', 'WRONG_B'],   # left over from an older run
    'cell_type':         ['L6 IT', 'Sst'],
    'cell_type_source':  ['seurat_onramp', 'seurat_onramp'],
})
chk("canonical wins over the stale column", select(stale), ('cell_type', 'seurat_onramp'))

print("=== legacy h5ad with no canonical column ===")
legacy = pd.DataFrame({'celltypist_prediction': ['A', 'B']})
chk("still falls back to tool columns", select(legacy), ('celltypist_prediction', 'celltypist'))

print("=== canonical present but no provenance stamped ===")
noprov = pd.DataFrame({'cell_type': ['A', 'B']})
chk("provenance defaults to 'canonical'", select(noprov), ('cell_type', 'canonical'))

print("=== degenerate: cell_type all 'Unknown' ===")
unk = pd.DataFrame({'cell_type': ['Unknown', 'Unknown'],
                    'scanvi_prediction': ['A', 'B']})
chk("skips unusable canonical, uses scanvi", select(unk), ('scanvi_prediction', 'scanvi'))

print("=== nothing usable ===")
chk("returns None (now warns at the call site)",
    select(pd.DataFrame({'foo': [1, 2]})), (None, None))

print("\nALL PASS" if ok else "\nFAILURES ABOVE")
sys.exit(0 if ok else 1)
