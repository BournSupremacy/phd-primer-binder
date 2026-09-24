#!/bin/bash
#SBATCH --job-name=binder-align
#SBATCH --output=logs/02_align_%j.out
#SBATCH --error=logs/02_align_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --partition=htc-el8

# bwa mem against 88 short references is fast per sample, but 12 samples
# run sequentially here - still light, but worth batching rather than
# occupying an interactive shell. No --mem set (standard per-cpu default).
#
# Usage: sbatch slurm/02_align_and_count.sbatch [--ref data/binders.fasta]

set -euo pipefail
mkdir -p logs
scripts/02_align_and_count.sh "$@"
