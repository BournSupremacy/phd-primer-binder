#!/bin/bash
#SBATCH --job-name=binder-demux
#SBATCH --output=logs/00_demux_%j.out
#SBATCH --error=logs/00_demux_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --partition=htc-el8

# bcl-convert is the heaviest, most multi-threaded step in this pipeline -
# this is the one you most want off the login node. No --mem set (standard
# per-cpu default), and 1h is comfortable headroom for a 5M-read MiSeq i100
# run - bump cpus-per-task/time up if you outgrow this dataset.
#
# Usage: sbatch slurm/00_bcl_to_fastq.sbatch --run-dir /path/to/MiSeqRun [--force]

set -euo pipefail
mkdir -p logs
scripts/00_bcl_to_fastq.sh "$@"
