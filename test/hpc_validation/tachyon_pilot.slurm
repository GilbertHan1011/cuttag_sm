#!/bin/bash
#SBATCH --job-name=cuttag-tachyon-pilot
#SBATCH --partition=intel-sc3
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend/logs/cuttag-%j.out
#SBATCH --error=/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend/logs/cuttag-%j.err

set -euo pipefail

root=/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend
repo=$root/code/cuttag_sm
snake=/storage/zhangkaiLab/hanlitian/software/micromamba/micromamba/envs/py311/bin/snakemake
sample=B51-1_CUT22
export CONDA_PKGS_DIRS=$root/conda-pkgs
export TMPDIR=$root/tmp
export XDG_CACHE_HOME=$root/xdg-cache

run_arm() {
  local arm=$1 config=$2
  local out=$root/results/cuttag/$arm
  "$snake" \
    "$out/Important_processed/Bam/$sample.sorted.markd.bam" \
    "$out/Important_processed/Bam/$sample.sorted.markd.fraglen.tsv" \
    "$out/Important_processed/Track/tracks/$sample.bw" \
    "$out/Important_processed/Peaks/callpeaks/${sample}_peaks.bed" \
    "$out/Report/dtools/fingerprint_$sample.tsv" \
    "$out/Report/preseq/lcextrap_$sample.txt" \
    "$out/Report/plotEnrichment/frip_$sample.tsv" \
    --directory "$repo" --snakefile "$repo/workflow/Snakefile" \
    --workflow-profile none --executor local --cores "$SLURM_CPUS_PER_TASK" \
    --resources mem_mb=120000 --use-conda --conda-prefix "$root/conda" \
    --rerun-incomplete --configfile "$repo/test/hpc_validation/$config"
}

run_arm legacy config_legacy.yml
run_arm tachyon config_tachyon.yml
