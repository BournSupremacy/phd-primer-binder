#!/usr/bin/env bash
# Submit the three pipeline steps (scripts/00_bcl_to_fastq.sh,
# scripts/01_qc_trim.sh, scripts/02_align_and_count.sh - each carries its
# own #SBATCH directives, so sbatch can run them directly) as SLURM jobs,
# chained with --dependency=afterok so each only starts once the previous
# one has finished successfully. This is the sbatch equivalent of
# run_pipeline.sh - see README.md "Running on a SLURM cluster" for when to
# use this vs. submitting/running a single step yourself.
#
# Usage:
#   slurm/submit_pipeline.sh --run-dir /path/to/MiSeqRun [--force] [--ref data/binders.fasta]
#
# Run from the repo root (paths in the scripts are relative to it).

set -euo pipefail

RUN_DIR=""
REF="data/binders.fasta"
SAMPLE_SHEET="config/IlluminaSampleSheet.csv"
FORCE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir) RUN_DIR="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --sample-sheet) SAMPLE_SHEET="$2"; shift 2 ;;
    --force) FORCE=true; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$RUN_DIR" ]]; then
  echo "Error: --run-dir is required (path to the raw sequencer run folder)" >&2
  exit 1
fi

mkdir -p logs

demux_args=(--run-dir "$RUN_DIR" --sample-sheet "$SAMPLE_SHEET")
$FORCE && demux_args+=(--force)

demux_id=$(sbatch --parsable scripts/00_bcl_to_fastq.sh "${demux_args[@]}")
echo "Submitted demux job $demux_id"

qc_id=$(sbatch --parsable --dependency=afterok:"$demux_id" scripts/01_qc_trim.sh)
echo "Submitted QC/trim job $qc_id (after $demux_id)"

align_id=$(sbatch --parsable --dependency=afterok:"$qc_id" scripts/02_align_and_count.sh --ref "$REF")
echo "Submitted align/count job $align_id (after $qc_id)"

echo ""
echo "Track progress with: squeue -u \$USER"
echo "Logs land in logs/00_demux_${demux_id}.out, logs/01_qc_${qc_id}.out, logs/02_align_${align_id}.out (and matching .err files)"
echo "Once align_id ($align_id) finishes, knit analysis/DMS_analysis.Rmd"
