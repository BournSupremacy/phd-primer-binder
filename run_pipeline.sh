#!/usr/bin/env bash
# Top-level orchestrator: runs the demux -> QC/trim -> align/count steps in
# order. Each step is also a standalone script you can run and inspect on
# its own (see scripts/) - this is just a convenience wrapper.
#
# Usage:
#   ./run_pipeline.sh --run-dir /path/to/MiSeqRun --ref data/binders.fasta
#
# To dry-run the whole pipeline on synthetic data instead of a real
# sequencer run, skip --run-dir and pass --simulate:
#   ./run_pipeline.sh --simulate --ref data/binders.fasta
#
# Pass --force to let a rerun overwrite a previous demux output directory
# (bcl-convert refuses to run into one that already exists otherwise).

set -euo pipefail

RUN_DIR=""
REF="data/binders.fasta"
SAMPLE_SHEET="config/IlluminaSampleSheet.csv"
SIMULATE=false
FORCE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir) RUN_DIR="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --sample-sheet) SAMPLE_SHEET="$2"; shift 2 ;;
    --simulate) SIMULATE=true; shift ;;
    --force) FORCE=true; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

demux_args=()
$FORCE && demux_args+=(--force)

if $SIMULATE; then
  echo "=== Simulating test data (no real sequencer run needed) ==="
  python3 scripts/simulate_test_data.py --ref "$REF"
elif [[ -n "$RUN_DIR" ]]; then
  echo "=== Step 0: demultiplexing ==="
  scripts/00_bcl_to_fastq.sh --run-dir "$RUN_DIR" --sample-sheet "$SAMPLE_SHEET" "${demux_args[@]}"
else
  echo "Error: pass either --run-dir <path> (real data) or --simulate (dry run)" >&2
  exit 1
fi

echo "=== Step 1: QC / adapter trimming ==="
scripts/01_qc_trim.sh

echo "=== Step 2: alignment + per-binder counting ==="
scripts/02_align_and_count.sh --ref "$REF"

echo ""
echo "All done. Next: open analysis/DMS_analysis.Rmd in RStudio and knit it,"
echo "or run: Rscript -e \"rmarkdown::render('analysis/DMS_analysis.Rmd')\""
