#!/usr/bin/env bash
# Align each sample's trimmed paired-end reads to the panel of 88 known binder
# sequences and produce a per-sample, per-binder read count.
#
# Why plain paired-end alignment and not read-stitching:
#   bwa mem (like any standard PE aligner) places R1 and R2 independently onto
#   the same reference contig in FR orientation - it does NOT require the two
#   mates to overlap. So even though our binders (>300bp) are longer than the
#   R1+R2 combined footprint (300bp) and the mates never touch in the middle,
#   bwa mem still correctly identifies which of the 88 binder sequences a read
#   pair came from. Stitching (DiMSum/PEAR-style) is only needed when you must
#   reconstruct one full-length consensus read per molecule (e.g. to call
#   mutations across the whole amplicon) - we don't need that here since we're
#   just counting abundance against a small, known reference panel.
#
# Usage: scripts/02_align_and_count.sh [--samples config/samples.tsv]
#          [--ref data/binders.fasta] [--trimmed-dir data/qc_trimmed]
#          [--out-dir results] [--min-mapq 20]

set -euo pipefail

# This step needs bwa, samtools, and python3 with pysam on PATH. Uncomment
# ONE of the following depending on what your cluster provides - see
# README.md "Environment setup".
#
# Option A - cluster environment modules (check exact names/case first with
# `module avail bwa samtools python`):
# module load BWA SAMtools Python
#
# Option B - conda/mamba env built from environment.yml:
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate rcaT-binder-screen

SAMPLES_TSV="config/samples.tsv"
REF="data/binders.fasta"
TRIMMED_DIR="data/qc_trimmed"
OUT_DIR="results"
MIN_MAPQ=20

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES_TSV="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --trimmed-dir) TRIMMED_DIR="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ ! -f "$REF" ]]; then
  echo "Error: reference fasta not found at $REF" >&2
  exit 1
fi

mkdir -p "$OUT_DIR/bam"

if [[ ! -f "${REF}.bwt" ]]; then
  echo "Indexing reference $REF with bwa..."
  bwa index "$REF"
fi
samtools faidx "$REF"

tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id _condition _rep _ara _iptg _r1 _r2; do
  [[ -z "$sample_id" ]] && continue
  r1_trim="$TRIMMED_DIR/$sample_id/${sample_id}_R1.trimmed.fastq.gz"
  r2_trim="$TRIMMED_DIR/$sample_id/${sample_id}_R2.trimmed.fastq.gz"
  if [[ ! -f "$r1_trim" || ! -f "$r2_trim" ]]; then
    echo "WARNING: missing trimmed fastq for $sample_id - skipping. Run 01_qc_trim.sh first." >&2
    continue
  fi
  bam="$OUT_DIR/bam/${sample_id}.bam"
  echo "Aligning $sample_id -> $bam"
  bwa mem -t 4 "$REF" "$r1_trim" "$r2_trim" \
    | samtools sort -@ 2 -o "$bam" -
  samtools index "$bam"
done

echo "Alignment done. Tallying per-binder read counts (MAPQ >= $MIN_MAPQ, properly paired only)..."
python3 scripts/count_reads_per_binder.py \
  --samples "$SAMPLES_TSV" \
  --bam-dir "$OUT_DIR/bam" \
  --min-mapq "$MIN_MAPQ" \
  --out "$OUT_DIR/binder_counts.tsv"

echo "Done. Combined counts table written to $OUT_DIR/binder_counts.tsv - load this into analysis/DMS_analysis.Rmd"
