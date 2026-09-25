#!/usr/bin/env bash
#SBATCH --job-name=binder-demux
#SBATCH --output=logs/00_demux_%j.out
#SBATCH --error=logs/00_demux_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --partition=htc-el8
#
# Demultiplex the MiSeq i100 run into per-sample paired-end fastq.gz files.
#
# The MiSeq i100 series writes base calls in a newer format than classic
# bcl2fastq2 was built for; Illumina's supported demultiplexer for it is
# BCL Convert, which we use by default. If you specifically need legacy
# bcl2fastq (e.g. an older run), pass --demux bcl2fastq - but check it
# actually supports your RunInfo/RTA version first.
#
# Runs standalone (bash scripts/00_bcl_to_fastq.sh ...) or as a SLURM job
# (sbatch scripts/00_bcl_to_fastq.sh ...) - the #SBATCH lines above are
# ordinary comments to bash and are only read by sbatch. This is the
# heaviest, most multi-threaded step in the pipeline - the one you most
# want off a login node. No --mem set (standard per-cpu default); bump
# cpus-per-task/time up if you outgrow a 5M-read MiSeq i100 run.
#
# Usage:
#   scripts/00_bcl_to_fastq.sh --run-dir /path/to/MiSeqRun \
#     [--sample-sheet config/IlluminaSampleSheet.csv] \
#     [--out-dir data/demux] [--demux bcl-convert|bcl2fastq] [--force]
#
# --force: bcl-convert refuses to write into an output folder that already
# exists (a safety check against accidentally overwriting someone else's
# results). Pass --force to let it overwrite - safe here since demux output
# is fully reproducible from the raw run + sample sheet, so re-running (e.g.
# after fixing an index in the sample sheet) is expected.

set -euo pipefail
mkdir -p logs
module load bcl-convert

RUN_DIR=""
SAMPLE_SHEET="config/IlluminaSampleSheet.csv"
OUT_DIR="data/demux"
DEMUX_TOOL="bcl-convert"
FORCE=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir) RUN_DIR="$2"; shift 2 ;;
    --sample-sheet) SAMPLE_SHEET="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --demux) DEMUX_TOOL="$2"; shift 2 ;;
    --force) FORCE=true; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$RUN_DIR" ]]; then
  echo "Error: --run-dir is required (path to the raw sequencer run folder)" >&2
  exit 1
fi
if [[ ! -f "$SAMPLE_SHEET" ]]; then
  echo "Error: sample sheet not found at $SAMPLE_SHEET (copy/edit config/IlluminaSampleSheet_template.csv)" >&2
  exit 1
fi

# Deliberately NOT mkdir-ing $OUT_DIR here - bcl-convert creates it itself
# and errors out if it already exists (see --force above).

bcl_convert_args=()
$FORCE && bcl_convert_args+=(--force)

if [[ "$DEMUX_TOOL" == "bcl-convert" ]]; then
  echo "Running BCL Convert..."
  bcl-convert \
    --bcl-input-directory "$RUN_DIR" \
    --output-directory "$OUT_DIR" \
    --sample-sheet "$SAMPLE_SHEET" \
    "${bcl_convert_args[@]}"
elif [[ "$DEMUX_TOOL" == "bcl2fastq" ]]; then
  echo "Running legacy bcl2fastq (make sure your sample sheet uses the [Data] format, not [BCLConvert_Data])..."
  bcl2fastq \
    --runfolder-dir "$RUN_DIR" \
    --output-dir "$OUT_DIR" \
    --sample-sheet "$SAMPLE_SHEET" \
    --no-lane-splitting
else
  echo "Error: --demux must be bcl-convert or bcl2fastq (got '$DEMUX_TOOL')" >&2
  exit 1
fi

# Both tools name output files <Sample_ID>_S<N>_R{1,2}_001.fastq.gz, all
# together in $OUT_DIR (no per-sample subdirectories). Rename into the flat
# data/demux/<sample_id>_R{1,2}.fastq.gz layout that config/samples.tsv and
# the rest of the pipeline expect.
#
# find_with_retry: on shared/network filesystems, a file bcl-convert just
# finished writing can take a moment to become visible to a fresh `find` -
# retry a few times with a short sleep before giving up, rather than
# leaving the last-written sample's fastq.gz stranded at the top level.
find_with_retry() {
  local pattern="$1"
  local hit=""
  for _ in 1 2 3 4 5; do
    hit=$(find "$OUT_DIR" -maxdepth 3 -name "$pattern" | head -n1 || true)
    [[ -n "$hit" ]] && break
    sleep 2
  done
  echo "$hit"
}

echo "Renaming fastq files in $OUT_DIR/..."
while IFS=$'\t' read -r sample_id _condition _rep _ara _iptg _r1 _r2; do
  [[ "$sample_id" == "sample_id" ]] && continue
  found_r1=$(find_with_retry "${sample_id}_S*_R1_001.fastq.gz")
  found_r2=$(find_with_retry "${sample_id}_S*_R2_001.fastq.gz")
  if [[ -z "$found_r1" || -z "$found_r2" ]]; then
    echo "  WARNING: could not find demuxed fastq for sample '$sample_id' - check the sample sheet Sample_ID matches config/samples.tsv" >&2
    continue
  fi
  mv "$found_r1" "$OUT_DIR/${sample_id}_R1.fastq.gz"
  mv "$found_r2" "$OUT_DIR/${sample_id}_R2.fastq.gz"
  echo "  $sample_id -> $OUT_DIR/${sample_id}_R{1,2}.fastq.gz"
done < config/samples.tsv

echo "Demultiplexing done. Fastq files are under $OUT_DIR/"
