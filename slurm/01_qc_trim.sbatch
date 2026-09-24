#!/bin/bash
#SBATCH --job-name=binder-qc
#SBATCH --output=logs/01_qc_%j.out
#SBATCH --error=logs/01_qc_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# Add your cluster's partition/account if required, e.g.:
# #SBATCH --partition=htc
# #SBATCH --account=dorrity

# fastp over 12 samples run sequentially - light on resources, but still
# minutes rather than seconds, so batch it rather than tying up a shell.
#
# Usage: sbatch slurm/01_qc_trim.sbatch

set -euo pipefail
mkdir -p logs
scripts/01_qc_trim.sh "$@"
