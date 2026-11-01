# Seurat `.rds` on-ramp — verification gates

These gates exist because every defect found in the `.rds` on-ramp so far was
**silent**. Nothing crashed, no warning was printed, and the numbers that came
out were wrong in ways only a comparison could reveal. Three of them shipped
past a gate that could not see them:

| Defect | What it did | Caught by |
|---|---|---|
| `logcounts` dropped on the h5ad → Seurat leg | hdWGCNA metacells composed from a raw `data` layer (991 vs 1005 metacells) | the **equivalence** gate |
| CellChat reading `X` as normalised | raw counts fed straight into L–R inference, no metacell step to launder them | the **equivalence** gate |
| no normalisation carried for an unnormalised object | `ScaleData`/PCA/kNN on raw values | the **equivalence** gate, on an unnormalised subject |

The **fidelity** gate saw none of them: the matrix was bit-exact in every case.
That is the whole argument for the analysis-level gates existing.

---

## Naming

`docs/verification.md` already uses **Tier 1 / 2 / 3** for the whole pipeline
(pre-flight → tiny dataset → full datasets). These gates are a different axis
and are named for what they prove:

| Name | Asks | Scripts |
|---|---|---|
| **fidelity** | did the matrix survive conversion? | `summarize_seurat.R`, `summarize_h5ad.py`, `compare.py` |
| **equivalence** | does an *analysis* on that matrix land in the same place? | `roundtrip_hdwgcna.sbatch`, `roundtrip_cellchat.sbatch`, `compare_*.R` |
| **scale** | does it hold at atlas size, and what does it cost? | the same equivalence harness on a large object |
| **determinism control** | is the equivalence gate even readable? | `determinism_control.sbatch` |
| **unit** | the one gate that needs no cluster | `test_build_mudata_celltype.py` |

Commit messages and session logs written before 2026-10-30 call these
Tier 1 / Tier 2 / Tier 3. Same things.

---

## Run the unit test (seconds, no cluster)

```bash
python3 tests/onramp_rds/test_build_mudata_celltype.py
```

Covers `BUILD_MUDATA`'s cell-type column precedence. It `exec`s the real
`_CT_CANDIDATES` and `_resolve_ct_source` out of `bin/build_mudata_batched.py`
rather than restating them, so it cannot drift from the shipped code. Needs
only `pandas`. This one runs in CI.

---

## Fidelity gate

Three scripts, deliberately split so the two sides never share a runtime — not
even a container. Everything crosses the boundary as CSV/JSON.

```bash
# reference side, from the ORIGINAL object
Rscript tests/onramp_rds/summarize_seurat.R --rds in.rds --outdir ref/
# converted side, from the h5ad the converter produced
python3 tests/onramp_rds/summarize_h5ad.py --h5ad out.h5ad --outdir test/
# diff them
python3 tests/onramp_rds/compare.py --ref ref/ --test test/
```

Exit 0 = every check passed. Tolerances are tight on purpose: the converter
reuses the sparse slots verbatim, so the only expected difference is float64
round-tripping through HDF5, which is exact.

`plot_compare.py` renders the two artifact sets side by side when a number
disagrees and you want to see where.

---

## Equivalence gates (SLURM, hours)

Both run two arms on the same cell types and compare the *analysis*, not the
matrix. The native arm is pre-conditioned so the only variable under test is
the converter round trip — see the step 0b comment block in
`roundtrip_hdwgcna.sbatch` for which normalisations are mirrored and why.

```bash
sbatch --export=ALL,\
RDS=/path/object.rds,CT_KEY=subclass_label,N_CT=3,OUT=/path/outdir \
  tests/onramp_rds/roundtrip_hdwgcna.sbatch

sbatch --export=ALL,\
RDS=/path/object.rds,CT_KEY=subclass_label,N_CT=6,OUT=/path/outdir \
  tests/onramp_rds/roundtrip_cellchat.sbatch
```

`N_CT` picks the N most abundant cell types rather than hardcoding names, so
these run on any object and never build networks from too few cells.
Per-step results are **cached on their outputs**: a resubmit after a harness
fix skips what already succeeded. That also means *a stale result will be
silently reused* — delete `$OUT/native/` (or the whole `$OUT`) when you change
anything that should invalidate it.

Use a **fresh `$OUT` per question.** Reusing a green run's directory skips
every step and proves nothing.

### Reading a failure

`State=FAILED` means **either** a metric drifted **or** an arm crashed. Read
`$OUT/logs/04_compare_*.log` to tell which. Do not conflate them — three of the
first four failures here were harness bugs, not product defects.

**Do not loosen a threshold to make the gate green.** Thresholds live in one
`THRESH` block at the top of `compare_hdwgcna.R` with the rationale. If a
metric fails, the first question is whether a package or container version
changed, not whether the number is too strict.

---

## Determinism control — run this first, and after any upgrade

```bash
sbatch --export=ALL,\
RDS=/path/object.rds,CT=<cell type>,CT_KEY=<column>,OUT=/path/outdir \
  tests/onramp_rds/determinism_control.sbatch
```

Two replicates of the **same** arm on the **same** input. It answers a question
the equivalence gates cannot answer about themselves: *does hdWGCNA have a
noise floor?* If it did, the thresholds would be scoring scatter and no
converter work could ever turn the gate green.

Measured floor: **zero** — all nine metrics at exactly 1.0000 (job 58561738).
That is what licenses the 0.99 thresholds and what makes any deviation from
1.0 converter signal rather than noise.

**A zero noise floor is a property of the installed hdWGCNA, WGCNA, Seurat and
BLAS — not a law.** Re-run this control after any of those change, before
trusting an equivalence result or editing a threshold.

---

## What is gated, and what is not

Gated: the converted matrix; hdWGCNA end to end through both converters;
CellChat L–R inference; `BUILD_MUDATA`'s label resolution (unit).

**Not gated:** MAST differential expression, which goes through the same
`convert_h5ad_to_seurat.R` as hdWGCNA but has no equivalence test of its own;
and the multiome / GRN consumers — MuData assembly end to end, DORC, SCENIC+,
pycisTopic. Check those against a native run before publishing off them.
