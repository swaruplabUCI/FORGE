# Wiring FORGE with an LLM

FORGE asks you to write two files: a [manifest CSV](core/manifest.md) and a
[dataset config](core/config.md). This can be tedious and has a steep learning curve if you are unfamiliar with Nextflow. 

Turns out AI coding assistants do this quite well, provided proper context. This page gives you that prompt, and demonstrates a worked
example.

!!! info "This is a **setup aid**, not an autopilot."
     The prompt produces a manifest and
    a config, proves they pass FORGE's pre-flight checklist, and reports what it
    assumed. You still read the result before you launch. The
    [review checklist](#before-you-launch) at the bottom is the part a human should do.

---

## What your AI/LLM/agentic coding assistant needs

!!! tip "Which model"
    We tested this prompt on two tiers of model (Opus vs Haiku) against the same PBMC setup
    below. **Both produced a launchable config that passed pre-flight on the
    first try.** The smaller, faster model got every gate and every required
    parameter right; where it fell short was in the accuracy of its summary of its own
    work. 

    So: a lightweight model is enough to get a first run going. A stronger one might be worth it if you intend to build out and expand the pipeline later. **Either way, work the
    [review checklist](#before-you-launch) yourself**. Trust the config file and the pre-flight output over the assistant's description of them.

Give it a clone of the repository and the ability to run shell commands. The
prompt tells it to read `docs/core/`, `nextflow.config`, and
`validateStartupParams()` in `main.nf` — so it works from the repository as the
source of truth rather than from whatever it remembers about Nextflow pipelines.

Then paste the prompt below. Set `$FORGE` to your clone first, or edit the paths
in the verification block.

??? example "The prompt — click to expand, then copy"

    ```text
    You are helping me configure FORGE, a Nextflow pipeline for end-to-end
    single-cell / single-nucleus multiome (RNA + ATAC) analysis. You have a clone
    of the repository on disk and can read files and run shell commands.

    Your job is to produce a manifest CSV and a dataset config that correctly
    describe my experiment, and to prove they are valid before I spend compute.

    ## Hard constraints

    1. Write only my manifest CSV and one dataset config. Do not edit main.nf,
       nextflow.config, modules/, bin/, or configs/resource_tiers/.
    2. Most of what you need is a parameter. Some things are not — SLURM accounts
       and partitions are baked in by the cluster profile and the resource tiers
       at parse time, and cannot be changed from a dataset config. If you hit one
       of those, stop and tell me. Do not invent a workaround, and do not leave a
       directive in place that does not fit my cluster.
    3. Every parameter you write must literally exist in nextflow.config. Confirm
       each with grep -n. Nextflow silently ignores undeclared keys, so an
       invented name is not an error — it is a stage that quietly does nothing.
    4. Do not guess paths. ls every path you write.
    5. Do not fabricate values I have to supply. If you do not know, ask me.

    ## Read these first, in this order

    - docs/core/manifest.md — the input contract
    - docs/core/config.md — layering rules and block reference
    - docs/core/architecture.md — which sub-workflow each gate turns on
    - nextflow.config — the authority. Where a doc and this file disagree, this
      file wins.
    - configs/datasets/ — working examples. Treat them as STRUCTURE ONLY. Their
      paths are specific to the authors' cluster, several ship expensive stages
      already enabled, and some carry stale parameter values that main.nf no
      longer reads.

    Read validateStartupParams() in main.nf before writing anything. It is the
    pre-flight checklist in source form. Reading it up front is far cheaper than
    discovering its rules one failed preview at a time.

    ## Interview me before writing anything

    Ask, then wait for my answers:

    1. Species and genome build.
    2. Platform (10x Multiome / BD Rhapsody) and number of samples.
    3. Tissue, and whether there is more than one experimental condition. If so,
       what the groups are called. If not, what single label I want in
       condition_group — it is mandatory even with one condition.
    4. Where the raw files live and what they are named.
    5. What sample_id I want for each sample. These become the prefix on nearly
       every output file and the key that joins RNA to ATAC; changing one later
       invalidates the resume cache.
    6. Which references I have on disk: GTF, ENCODE blacklist, the scATAnno
       reference atlas (.h5ad) for my tissue, cisTarget rankings/scores + motif
       annotations, motif PFMs.
    7. Whether I have an integrated RNA reference atlas for scANVI
       (ref_dir_human_integrated / ref_dir_mouse_integrated). This decides
       whether the RNA arm trains scVI/scANVI or runs CellTypist directly, and it
       is the biggest fork in that arm.
    8. How many cells the chemistry actually recovered, for cellbender.expected_cells.
    9. What I want out of the run — cell types only, or also co-accessibility /
       GRNs / TF footprinting.
    10. Cluster and scheduler: executor, account, partition, GPU availability.
    11. Have the containers been built? ls singularity_cache/*.sif. A missing
        .sif is only a WARNING under -preview but a FATAL error on a real run, so
        a config can pass every check I ask of you and still die on launch.
    12. Where should outputs and the Nextflow work directory go? On a quota'd
        shared filesystem this matters, and the work dir is a CLI flag (-w), not
        a config key.

    If I am vague on 6 or 9, tell me which references are required for what I
    asked for, and which stages you will leave off until I have them.

    ## Rules you cannot derive from reading the config

    1. A dataset config cannot select the resource tier. params.resource_tier is
       read while nextflow.config is parsed, before any -c file merges. Writing
       it in a dataset config is decorative. Pass --resource_tier <tier> on the
       nextflow run command line, or use a profile.

    2. Parse-time interpolation strands every derived default. Any value in the
       params block computed from another parameter is frozen when
       nextflow.config is parsed, before your -c file merges. Do this — do not
       just trust the examples:

           grep -n '\${params\.\|params\.[a-z_]* ==' nextflow.config

       Re-declare, in your own config, every hit whose input you override. Put
       the re-declaration AFTER your own assignment of the input it depends on,
       in the same file — then "${params.outdir}/cicero" resolves correctly,
       because your file's own outdir assignment has already happened. That is
       the convention every shipped dataset config uses. A literal path works
       too. Known cases: cicero.outdir (otherwise
       publishes to results/cicero while everything else goes to your outdir —
       the stage still succeeds, which is what makes it easy to miss) and
       scprinter.gtf_human / gtf_mouse (otherwise resolve to the literal string
       'null'). pipeline_info/ is stranded the same way and cannot be fixed from
       a dataset config or the CLI; expect a stray results/ directory.

    3. atac.sample_metadata is mandatory whenever atac.run = true. main.nf calls
       file() on it unguarded, so null aborts the run. Normally the same path as
       metadata_file.

    4. Per-batch directories need three parameters, not one. docs/core/manifest.md
       documents batch_dirs as the alternative to a per-row data_dir. That covers
       RNA only. ATAC resolves through atac_batch_dirs and atac_coord_batch_dirs,
       neither of which is declared in nextflow.config — they exist only in
       dataset configs. If you omit data_dir, set all three, and see
       configs/datasets/tutorial_pbmc.config for the pattern. The simpler and
       safer choice for a small study is a per-row data_dir.

    5. ATAC annotation means scATAnno. atac.annotation_method = 'scatanno' is the
       only supported value — 'celltypist' aborts the run, because ATAC
       CellTypist-on-gene-activity was removed. scATAnno requires
       scatanno.reference_atlas, a .h5ad atlas for my tissue, and pre-flight
       enforces it. If I do not have one, say so and stop; getting the right
       atlas is the fix. (atac.marker_file exists and bypasses the atlas check,
       but it is a "super-user" escape the shipped tutorial uses to avoid a
       2.76 GB download. Do not reach for it to route around a missing atlas.)
       params.celltypist governs the RNA arm only. RNA and ATAC cell-type labels
       are deliberately separate keys produced by different tools — do not try to
       make them agree.

    6. Expensive blocks default off in nextflow.config but ship ON in the dataset
       configs you were told to read. These all default run = true: chromvar,
       scprinter, pycistopic, scenicplus, enhancer_footprinting,
       enhancer_recipe_c, enhancer_viz, cellchat, hdwgcna. Enhancer footprinting
       alone was 54% of all compute across the four published datasets, and
       cellchat / hdwgcna fan out per cell type. For a first run, turn off
       everything I did not ask for, get a green run, then enable one block at a
       time.

    7. Nested gates: the outer one wins. msfp_strip.enabled does nothing unless
       enhancer_footprinting.msfp_enabled is also true.

    8. Manifest literalism. sample_type must be exactly 'lane', lowercase.
       rna_file / fragment_file are FILENAMES, not paths — the directory comes
       from the row's data_dir or from the batch_dirs family. condition_group
       must be non-empty even for a single-condition study.

    9. celltypist.model defaults to Immune_All_Low.pkl — right for PBMC, wrong
       for every other tissue.

    10. cellbender.total_droplets must be strictly less than the barcode count in
        the raw matrix (equality is an IndexError), and expected_cells should
        match the kit's real recovery. To read the barcode count:

            singularity exec $FORGE/singularity_cache/scgpu_extended.sif python -c \
              "import h5py,sys; print(h5py.File(sys.argv[1])['matrix/shape'][:])" <raw.h5>

    11. Shallow or subset ATAC needs explicit QC thresholds. A null
        atac.min_counts does NOT fall back to atac.initial_min_counts; it falls
        through to the underlying script's own default of 5000, which is right
        for whole-genome 10x and wrong for anything shallower.

    ## GPU

    Assume a GPU is available. CELLBENDER and the MULTIVI_* processes are pinned
    to GPU partitions by the production resource tiers (small / medium / large),
    so those tiers need one. Note that MULTIVI_INTEGRATE and MULTIVI_VISUALIZE
    request a V100 by literal string, not via params.slurm_gpu_type — if I have
    only one card type, check that before launching.

    If I tell you I have no GPU, do not rewrite process directives. Tell me, and
    give me the two supported options:

    - run the CPU-only tutorial tier / profile, which strips GPU directives by
      design, or
    - do my own ambient-RNA correction and QC in scanpy and skip the RNA arm's
      upstream via onramp.rna_integrated_h5ad — pairing it with
      onramp.rna_per_sample_h5ads_dir, which main.nf requires whenever
      run_multiome_integration is true.

    ## Then prove it

    Work from a scratch directory, not the repository. Set $FORGE to the clone.

        # 1. Confirm your keys landed in the merged config.
        nextflow -c <your>.config config $FORGE -profile cluster,singularity | grep -n '<key>'

        # 2. Confirm where work will actually be submitted. A green preview says
        #    nothing about whether my cluster can schedule it.
        nextflow -c <your>.config config $FORGE -profile cluster,singularity | grep clusterOptions

        # 3. Pre-flight + full DAG construction. No compute, ~15 s.
        nextflow run $FORGE/main.nf -profile cluster,singularity -preview \
            -c <your>.config --resource_tier <tier>

    Note that `nextflow config` rejects --resource_tier; only `nextflow run`
    accepts it. So step 3's "[OK] Resource tier (...)" line is the only place the
    tier is actually confirmed — do not conclude from steps 1-2 that it is right.

    Iterate until PRE-FLIGHT CHECKLIST PASSED. The checklist reports every
    problem at once, so fix them in batches. A missing-container warning under
    -preview is expected; it becomes a fatal error on a real run.

    -preview appends a session to .nextflow/history, and a later bare -resume
    picks the most recent session — the empty preview — and re-runs everything.
    Preview from a scratch directory, or pass -resume <prior-session-uuid>
    explicitly.

    ## Report back

    - The two files you wrote, and the full launch command.
    - Which sub-workflows will materialize, and how many process nodes -preview
      constructed. NODES ARE NOT TASKS — per-chromosome and per-cell-type stages
      fan out at runtime, so tell me the real job count you expect.
    - Every parameter you set that is NOT a repository default, one line each on
      why.
    - Anything you assumed, and what you left off and what would turn it on.
    ```

---

## A worked example: 10k human PBMC

The public 10x Genomics *10k Human PBMC, Multiome, nextGEM Chromium X* sample,
wired from scratch. The answers given to the interview were:

| Question | Answer |
|---|---|
| Species / build | Human, GRCh38 / hg38 |
| Platform | 10x Multiome, one sample |
| Tissue / conditions | PBMC, single condition (`Healthy`) |
| References on disk | Gencode v38 GTF, hg38 blacklist, scATAnno PBMC atlas |
| References *not* on disk | cisTarget rankings/scores, motif annotations, PFMs |
| Goal | Cell types on both arms, multiome integration, Cicero. **No** GRN, **no** footprinting |
| Cluster | SLURM, GPU available, `small` tier |

### What came back

The manifest is one row:

```csv
sample_id,batch,sample_type,original_lane_id,rna_file,fragment_file,condition_group,data_dir
10k_PBMC,PBMC,lane,L1,10k_PBMC_..._raw_feature_bc_matrix.h5,10k_PBMC_..._atac_fragments.tsv.gz,Healthy,/data/pbmc
```

The config, with the interesting decisions kept and boilerplate trimmed:

```groovy
params {
    species       = 'human'
    metadata_file = '/study/pbmc/pbmc_10k_manifest.csv'
    outdir        = '/study/pbmc/results_pbmc_firstpass'

    gtf_human_full = '/refs/gencode.v38.annotation.gtf'
    blacklist_bed  = '/refs/hg38-blacklist.v2.bed.gz'

    // No integrated atlas → RNA Path B (CellTypist direct, no scVI/scANVI)
    ref_dir_human_integrated = null

    rna        { run = true; annotation_method = 'celltypist' }
    celltypist { enabled = true; model = 'Immune_All_Low.pkl' }
    cellbender { total_droplets = 20000; expected_cells = 10000 }

    atac {
        run               = true
        // Mandatory: main.nf calls file() on this unguarded.
        sample_metadata   = '/study/pbmc/pbmc_10k_manifest.csv'
        annotation_method = 'scatanno'
        marker_file       = null          // do NOT bypass the atlas
        tissue_type       = 'pbmc'
    }
    scatanno { reference_atlas = '/refs/PBMC_reference_atlas_final.h5ad'; atlas_name = 'pbmc' }

    run_multiome_integration = true
    mofa { run = true }; multivi { run = true }

    cicero {
        run = true
        // Re-declared as a LITERAL: the default interpolates params.outdir at
        // parse time, before this file merges, and would publish to results/cicero.
        outdir               = '/study/pbmc/results_pbmc_firstpass/cicero'
        gtf_full             = '/refs/gencode.v38.annotation.gtf'
        gtf_plot             = '/refs/gencode.v38.annotation.gtf'
        use_chromvar_targets = false      // ChromVAR is off this pass
    }

    // Nine blocks that default ON, switched off for a first pass
    chromvar { run = false }
    scprinter {
        run = false
        // Re-declared: the default is "${params.gtf_human_full}", frozen to the
        // literal string 'null'. Fixed now so enabling this later just works.
        gtf_human = '/refs/gencode.v38.annotation.gtf'
    }
    enhancer_footprinting { run = false; msfp_enabled = false }
    enhancer_recipe_c { run = false }
    enhancer_viz      { run = false }
    pycistopic { run = false }            // no cisTarget references on disk
    scenicplus { run = false }
    cellchat   { run = false }            // fans out per cell type
    hdwgcna    { run = false }            // fans out per cell type

    differential { run = false }          // single condition
    differential_rna { run = false }
}
```

!!! note "Paths above are tidied for readability"
    The structure, comments, and every parameter decision are as produced. Only
    the long working-directory paths were shortened.

### The proof

```text
PRE-FLIGHT CHECKLIST PASSED (8 checks):
    [OK] Manifest schema (1 rows)
    [OK] Species/genome consistency
    [OK] GTF files (1 paths validated)
    [OK] scATAnno reference atlas (PBMC_reference_atlas_final.h5ad)
    [OK] MOFA mode (high_memory)
    [OK] CellTypist model: Immune_All_Low.pkl
    [OK] Containers (5 SIF files)
    [OK] Resource tier (small)
  No warnings.
```

One iteration, 10.8 s, exit 0. It reported **21 process nodes** and — correctly
distinguishing nodes from tasks — **≈45 real SLURM jobs**, because
`CICERO_FULL_CHROM` fans out across `chr1–22 + X/Y/M` for human while every other
stage is a single task at one sample.

Four sub-workflows materialize: `RNA` (Path B), `ATAC_INITIAL` → `ATAC_FINAL`
(scATAnno branch), `REGULATORY_ANALYSIS` (Cicero leg only), and
`MULTIOME_INTEGRATION`. `MULTIOME_GRN`, `ENHANCER_FOOTPRINTING_RECIPES`.

---

## The same reasoning, at tutorial scale

The [tutorial dataset](tutorial.md) is a ~1,000-cell subset of this same PBMC
sample with ATAC restricted to chr21 + chr22. It is the cheapest place to watch
the reasoning play out, because you can run the result end to end on 8 CPUs with
no GPU and no reference downloads.

---

## Before you launch

The prompt makes the assistant prove the configuration parses. It cannot prove
the configuration is what you *meant*. Read these yourself:

- [ ] **The config file is the ground truth.** Read the file,
      not the summary. In testing, a model's prose omitted two blocks it had in
      fact disabled correctly.
- [ ] **`sample_id` values are what you want to live with.** They prefix nearly
      every output and join RNA to ATAC. Renaming one later invalidates the cache.
- [ ] **`condition_group` labels are real.** A single-condition study still needs
      a label, so the assistant invented one if you did not supply it.
- [ ] **`cellbender.expected_cells` matches your actual recovery,** not the kit's
      nominal number.
- [ ] **The CellTypist model fits your tissue.** The default is an immune model.
- [ ] **Every stage you wanted is actually on,** and every stage you did not ask
      for is off. Check the `-preview` process list, not the config.
- [ ] **`clusterOptions` point at partitions you can use.** A green pre-flight
      says nothing about schedulability — `-preview` submits nothing. See
      [Adapting to your cluster](setup/cluster.md).
- [ ] **The containers exist.** Missing `.sif` files are a warning under
      `-preview` and fatal on a real run.
- [ ] **You are not resuming into a preview session.** Check
      `tail -3 .nextflow/history` before a bare `-resume`.

---

## Where this helps least

Be realistic about the boundary:

- **It may not correctly select your biology.** Which atlas, which CellTypist model, which QC thresholds suit your tissue — AI suggestions here should be carefully evaluated.
- **It cannot validate against data it has not seen.** Pre-flight checks that
  files exist and parameters cohere, not that your fragments are any good.
- **Site adaptation is out of scope.** SLURM accounts, partitions and QOS are
  baked in by the cluster profile and resource tiers at parse time. Moving FORGE
  to another cluster means editing those files — see
  [Adapting to your cluster](setup/cluster.md) — and the prompt deliberately
  tells the assistant to stop and say so rather than improvise.
- **Trust the pre-flight, not the prose.** If the assistant's explanation and the
  checklist disagree, manually read and understand the checklist.

---

## Where to go next

- [The three core files](core/index.md) — what the assistant is actually writing
- [Verifying FORGE works](verification.md) — the three verification tiers
- [Tutorial](tutorial.md) — run the result end to end
- [On-ramps & resuming](onramps.md) — skipping stages you have already computed
