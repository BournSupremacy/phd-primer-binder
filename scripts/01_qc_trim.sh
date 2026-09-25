#!/usr/bin/env bash
#SBATCH --job-name=binder-qc
#SBATCH --output=logs/01_qc_%j.out
#SBATCH --error=logs/01_qc_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --partition=htc-el8
#
# Adapter/quality-trim each sample's paired-end reads with fastp and write a
# per-sample HTML/JSON QC report. This is a light-touch step - we are not
# trying to merge/stitch R1+R2 (see README for why that's unnecessary here),
# just removing adapter read-through and low-quality tails before alignment.
#
# Runs standalone (bash scripts/01_qc_trim.sh) or as a SLURM job
# (sbatch scripts/01_qc_trim.sh) - the #SBATCH lines above are ordinary
# comments to bash and are only read by sbatch.
#
# Usage: scripts/01_qc_trim.sh [--samples config/samples.tsv] [--out-dir data/qc_trimmed]

set -euo pipefail
mkdir -p logs
module load fastp

SAMPLES_TSV="config/samples.tsv"
OUT_DIR="data/qc_trimmed"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES_TSV="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$OUT_DIR"

# mapfile (not `while read < file`/`| while read`) so a samples.tsv missing
# a trailing newline on its last line doesn't silently drop that last
# sample - `read` reports failure on a final line with no newline, which
# would otherwise make the while-loop condition false before the body ever
# runs for it.
mapfile -t sample_lines < <(tail -n +2 "$SAMPLES_TSV")
for line in "${sample_lines[@]}"; do
  IFS=$'\t' read -r sample_id _condition _rep _ara _iptg r1 r2 <<< "$line"
  [[ -z "$sample_id" ]] && continue
  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "WARNING: missing fastq for $sample_id ($r1 / $r2) - skipping. Run 00_bcl_to_fastq.sh first." >&2
    continue
  fi
  sample_out="$OUT_DIR/$sample_id"
  mkdir -p "$sample_out"
  echo "Trimming $sample_id..."
  fastp \
    -i "$r1" -I "$r2" \
    -o "$sample_out/${sample_id}_R1.trimmed.fastq.gz" \
    -O "$sample_out/${sample_id}_R2.trimmed.fastq.gz" \
    --detect_adapter_for_pe \
    --json "$sample_out/${sample_id}_fastp.json" \
    --html "$sample_out/${sample_id}_fastp.html" \
    --thread 4
done

echo "QC/trimming done. Trimmed reads and reports are under $OUT_DIR/"
echo "Tip: open the fastp .html reports to sanity-check read quality and adapter content before aligning."
