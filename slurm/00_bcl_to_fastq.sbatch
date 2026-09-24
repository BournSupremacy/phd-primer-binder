#!/bin/bash
#SBATCH --job-name=binder-demux
#SBATCH --output=logs/00_demux_%j.out
#SBATCH --error=logs/00_demux_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
# Add your cluster's partition/account if required, e.g.:
# #SBATCH --partition=htc
# #SBATCH --account=dorrity

# bcl-convert is the heaviest, most multi-threaded step in this pipeline -
# this is the one you most want off the login node. Adjust
# cpus-per-task/mem/time above to your run size; these defaults are
# generous for a MiSeq i100 5M run.
#
# Usage: sbatch slurm/00_bcl_to_fastq.sbatch --run-dir /path/to/MiSeqRun [--force]

set -euo pipefail
mkdir -p logs
scripts/00_bcl_to_fastq.sh "$@"
