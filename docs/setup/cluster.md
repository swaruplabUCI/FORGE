# Adapting to your cluster

FORGE was developed on UCI's HPC3, and the shipped profiles and resource tiers
encode that site's SLURM partitions, accounts, and QOS names. Running elsewhere
means changing those in three places. Nothing about the pipeline logic is
site-specific — only the scheduler details.

---

## What is site-specific

| File | Contains |
|---|---|
| `configs/profiles/hpc3_cluster.config` | Executor, queue names, account, job submission rate limits |
| `configs/resource_tiers/{small,medium,large}.config` | Per-process `cpus`, `memory`, `time`, and `clusterOptions` |
| `nextflow.config` (`profiles` block) | GPU allocation, container engine flags |

Grep for the parameters that name your site:

```bash
grep -rn "slurm_account\|slurm_partition\|slurm_qos\|gres=gpu" \
    nextflow.config configs/
```

---

## Step 1 — Set your scheduler identifiers

The tiers reference SLURM settings through parameters rather than hardcoding them
everywhere, so most of the work is redefining those in your dataset config:

```groovy
params {
    slurm_account                = 'my_lab'
    slurm_partition_cpu          = 'compute'
    slurm_account_gpu            = 'my_lab_gpu'
    slurm_partition_gpu          = 'gpu'
    slurm_partition_gpu_hugemem  = 'gpu-bigmem'
    slurm_qos_gpu_hugemem        = 'normal'
    slurm_gpu_type               = 'a100'   // must match your --gres names
    slurm_gpu_count              = 1
}
```

Confirm the values resolve before running anything:

```bash
nextflow config -profile cluster,singularity -c my_study.config | grep slurm
```

## Step 2 — Check the GPU request syntax

FORGE requests GPUs with an explicit `--gres` string. If your site uses different
GPU type names — or no type qualifier at all — the submission will be rejected.
The relevant lines look like:

```groovy
clusterOptions = "-A ${params.slurm_account_gpu} -p ${params.slurm_partition_gpu} \
                  --gres=gpu:${params.slurm_gpu_type}:${params.slurm_gpu_count}"
```

A few processes pin a specific device (`--gres=gpu:V100:1`, `--gres=gpu:A30:1`)
because they were tuned against known memory ceilings. Search the tier file for
`gres=gpu:` and adjust these to whatever your site provides.

## Step 3 — Review memory and walltime

The shipped values are deliberately generous. On HPC3, memory is requested via
cores and unused walltime is refunded, so over-asking is cheap while a walltime
kill is expensive. **If your site bills reserved rather than used resources, these
defaults will overspend.**

Right-size from the measured evidence rather than guessing: after any run,
`logs/nextflow/report.html` and `logs/nextflow/trace.txt` give per-process peak
RSS and runtime. Adjust the tier to fit what your data actually needed.

---

## Non-SLURM schedulers

FORGE only ships a SLURM profile, but Nextflow supports many executors. Add a
profile for yours:

```groovy
profiles {
    my_cluster {
        process.executor = 'sge'        // or 'lsf', 'pbs', 'awsbatch', 'k8s'
        process.queue    = 'all.q'
    }
}
```

Then run with `-profile my_cluster,singularity`. The catch is that
`clusterOptions` strings throughout the resource tiers are SLURM syntax and will
not translate — you will need to replace them with your scheduler's equivalents,
or drop them and rely on `cpus`/`memory`/`time` alone.

---

## Running without a scheduler

For small tests, run everything locally:

```bash
nextflow run main.nf -profile standard,singularity -c my_study.config
```

This ignores partitions and QOS entirely. It is appropriate for a small dataset
on a large workstation, and for the
[pre-flight validation tier](../verification.md#tier-1-pre-flight-validation) — which
needs no containers or GPU at all.

!!! warning "GPU-dependent stages"
    CellBender, scVI/scANVI, ChromVAR, MOFA+, and MultiVI expect a GPU. Several
    also carry large host-memory requests. On a machine without an NVIDIA GPU,
    disable those stages rather than expecting a CPU fallback.

---

## Container paths

If your `.sif` files live outside the repository — a shared read-only location, or
a scratch filesystem — override the map rather than moving files:

```groovy
params.containers = [
    scgpu:      '/shared/containers/scgpu_extended.sif',
    snapatac:   '/shared/containers/snapatac_extended.sif',
    scenicplus: '/shared/containers/scenicplus.sif',
    r_cicero:   '/shared/containers/cicero.sif',
    r_seurat:   '/shared/containers/seurat_extended.sif',
    r_cellchat: '/shared/containers/seurat_extended.sif',
]
```

You may also need to bind additional filesystems so containers can see your data
and references:

```groovy
singularity.runOptions = '--nv -B /data -B /refs -B /scratch'
```

See [Containers](containers.md) for build and bind-path details.

---

## A checklist

1. `nextflow config -profile cluster,singularity -c my_study.config` — do the
   SLURM parameters resolve to your site's values?
2. `nextflow run main.nf -preview -c my_study.config` — does the pre-flight
   checklist pass?
3. Submit one cheap stage (RNA QC only) and confirm jobs are accepted by the
   scheduler.
4. Check `logs/nextflow/trace.txt` for actual versus requested resources, and
   adjust the tier.
5. Scale up one block at a time.
