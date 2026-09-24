# metaFlye

No finished metaFlye assembly is on disk yet. `badread` 0.4.2 is installed at `/mnt/tank/scratch/dsmutin/envs/badread`. It is simulating 8× ONT from the ten `half_strains` genomes (`--seed 1`) into `examples/half/work/ont/reads.fastq`, and metaFlye will write `examples/half/work/flye`. The queued SLURM job 892944 was cancelled so it would not write that same directory. The notes below record the earlier failed install.

## Reads

`half_strains` Illumina reads are `sample_full_R1.fastq` and `sample_full_R2.fastq` only. A name-limited search (`*ont*`, `*nano*`, `*nanopore*`, `*.fastq` / `*.fq`, maxdepth 7) under

- `/mnt/tank/scratch/partition-metagenomics/smuteam/vaegbin_improved/examples/`
- `/mnt/tank/scratch/dsmutin/tools/my/samovar/`

found no ONT or nanopore FASTQ for `half_strains` or `low75`. Flye was not run on the Illumina reads.

## Simulator install

`badread` is not on `PATH` in `/nfs/home/dsmutin/miniconda3/envs/bubbleblower-asm`. A retry with `http_proxy`, `https_proxy`, `HTTP_PROXY`, `HTTPS_PROXY`, `ALL_PROXY`, and `all_proxy` unset failed while conda was fetching repodata, before any package was linked:

```
OSError: [Errno 122] Disk quota exceeded: '/nfs/home/dsmutin/miniconda3/pkgs/cache/497deca9.7cfa.tmp'
```

Command:

```
conda install -y -n bubbleblower-asm -c conda-forge -c bioconda badread
```

Conda 26.1.1, env `bubbleblower-asm`. NanoSim was not installed. The failure is the home package cache quota, so a second conda install would hit the same path. No reads were simulated.

## Flye

The binary that would be used is

`/mnt/tank/scratch/dsmutin/tools/my/samovar/samovar/.cache/samovar/envs/flye/bin/flye`

`flye --version` reports `2.9.6-b1802`. It was not executed. There is no `examples/half/work/flye` graph, no bubble count, and no minimap2 metric table.

## SLURM scripts that were not submitted

Reference FASTA would be the ten genomes already in

`/mnt/tank/scratch/partition-metagenomics/smuteam/vaegbin_improved/examples/half_strains/work/sim/`

(`GCF_000730425`, `GCF_000743015`, `GCF_000816985`, `GCF_000952955`, `GCF_001518855`, `GCF_001549955`, `GCF_001558215`, `GCF_002741615`, `GCF_002900365`, `GCF_003363755`). Coverage would be `8x` of that set. Job ids: none.

Simulator (`examples/half/work/ont/badread.sbatch`):

```bash
#!/bin/bash
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=bb-badread
#SBATCH --output=examples/half/work/ont/badread-%j.out
#SBATCH --error=examples/half/work/ont/badread-%j.err

set -euo pipefail
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY all_proxy

ROOT=/mnt/tank/scratch/dsmutin/tools/my/BubbleBlower
SIM=/mnt/tank/scratch/partition-metagenomics/smuteam/vaegbin_improved/examples/half_strains/work/sim
OUT=$ROOT/examples/half/work/ont
mkdir -p "$OUT"

cat "$SIM"/GCF_000730425.fna \
    "$SIM"/GCF_000743015.fna \
    "$SIM"/GCF_000816985.fna \
    "$SIM"/GCF_000952955.fna \
    "$SIM"/GCF_001518855.fna \
    "$SIM"/GCF_001549955.fna \
    "$SIM"/GCF_001558215.fna \
    "$SIM"/GCF_002741615.fna \
    "$SIM"/GCF_002900365.fna \
    "$SIM"/GCF_003363755.fna \
    > "$OUT/half_strains_sim.fna"

source /nfs/home/dsmutin/miniconda3/etc/profile.d/conda.sh
conda activate bubbleblower-asm
badread simulate --reference "$OUT/half_strains_sim.fna" --quantity 8x --seed 1 \
    > "$OUT/reads.fastq"
```

Flye, submitted only after the simulator job succeeds (`--dependency=afterok:<badread_jobid>`):

```bash
#!/bin/bash
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=bb-flye
#SBATCH --output=examples/half/work/flye/flye-%j.out
#SBATCH --error=examples/half/work/flye/flye-%j.err

set -euo pipefail
ROOT=/mnt/tank/scratch/dsmutin/tools/my/BubbleBlower
FLYE=/mnt/tank/scratch/dsmutin/tools/my/samovar/samovar/.cache/samovar/envs/flye/bin/flye
"$FLYE" --meta --nano-raw "$ROOT/examples/half/work/ont/reads.fastq" \
    --out-dir "$ROOT/examples/half/work/flye" --threads 16
```

BubbleBlower was not compared with Flye contigs. Misassembly counts are unknown.

## Submitted after the quota failure

The conda package cache was pointed at `/mnt/tank/scratch/dsmutin/conda-pkgs` and the env prefix at `/mnt/tank/scratch/dsmutin/envs/badread`, so the install does not write the home `pkgs/cache`. One SLURM job does both steps: install badread if needed, simulate 8× on the ten `half_strains` genomes, then `flye --meta`.

Job id: `892944` (`bb-flye`, partition `main`, 8 CPUs, 8G). Script: `examples/half/work/ont/badread_flye.sbatch`.

Eight metaSPAdes jobs for `low75`, `low75half`, `high100`, and `half_strains` (`892936`–`892943`) are pending. A 510G `cami3_illumina` job is running, and the extra assemblies sit in `QOSMaxMemoryPerUser` until that memory is free. They were resubmitted at 8G instead of 32G so they fit once the large job ends.
