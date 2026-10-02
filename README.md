# RcaT binder pooled screen - a minimal DMS-style pipeline

A small, from-scratch pipeline for the EMBL PhD Primer project *"Designing
functional protein interactions: computational protein design and
high-throughput biochemistry"* (Dorrity/Mahamid/Eustermann groups, MSB Unit).

## The experiment

88 computationally-designed binders against the RcaT toxin are pooled and
expressed in *E. coli* carrying the RcaT toxin-antitoxin system. RcaT kills
the cell unless something inhibits it, so **a binder's relative abundance
after selection is a readout of how well it protects the cell**: the better
a binder blocks RcaT, the more its carrier cells survive, and the more reads
it contributes to the pool.

12 samples: 6 conditions x 2 biological replicates (A/B).

| Condition | Ara (toxin inducer) | IPTG (binder inducer) |
|---|---|---|
| no_inducer | - | - |
| toxin_only | 0.2% | - |
| library_only_0.1mM | - | 0.1 mM |
| library_only_1mM | - | 1 mM |
| toxin_plus_library_0.1mM | 0.2% | 0.1 mM |
| toxin_plus_library_1mM | 0.2% | 1 mM |

Sequencing: Illumina MiSeq i100 5M, 150bp paired-end, dual PCR barcoding for
demultiplexing.

## Repository layout

```
config/
  samples.tsv                     sample sheet used by every downstream step
  IlluminaSampleSheet_template.csv template BCL Convert sample sheet for demux (fill in real indices)
data/
  binders.fasta                   the 88 designed RcaT binder sequences - the alignment reference used throughout
  plate_assay_reference_hits.tsv  collaborator's published plate-based screen results, for the comparison in DMS_analysis.Rmd
  demux/                          demultiplexed fastq.gz land here, flat, as <sample_id>_R{1,2}.fastq.gz
scripts/
  00_bcl_to_fastq.sh               demultiplex raw sequencer output -> data/demux/ (also an sbatch script - see below)
  01_qc_trim.sh                    fastp adapter/quality trimming (also an sbatch script)
  02_align_and_count.sh            bwa mem alignment to data/binders.fasta + counting (also an sbatch script)
  count_reads_per_binder.py        per-binder read-pair tallying used by 02_align_and_count.sh
  simulate_test_data.py            generate synthetic reads so you can dry-run the whole thing now
run_pipeline.sh                    convenience wrapper that runs 00->01->02 in order
slurm/submit_pipeline.sh           submits the same 3 scripts via sbatch, chained with job dependencies - see "Running on a SLURM cluster"
environment.yml                    conda env for the demux/QC/alignment tools (fallback if your cluster has no modules for them)
analysis/
  DMS_analysis.Rmd                 the actual DMS-style analysis (normalisation, log2FC, replicate QC, hit table)
  DMS_analysis_expression_batches.Rmd  variant that splits the panel into high-/low-expressing batches before
                                        normalising, so a handful of dominant binders can't compositionally
                                        drown out the comparison for everyone else - see its own Section 0
results/                           binder_counts.tsv and everything the R notebook produces (gitignored except final CSV)
```

## Running it

### 1. With real sequencing data

`data/binders.fasta` already contains the 88 designed binder sequences
(headers = their design names) and is used as the alignment reference
throughout.

1. Copy `config/IlluminaSampleSheet_template.csv` to
   `config/IlluminaSampleSheet.csv` and fill in the real i7/i5 index
   sequences for each of the 12 PCR-barcoded samples. A few notes on this
   file:
   - Keep the `Sample_ID` values identical to `sample_id` in
     `config/samples.tsv` so downstream steps can find the right fastq
     files automatically.
   - Don't add comments or any other free text to this file - BCL Convert's
     sample sheet parser doesn't support a comment syntax, and stray text
     before `[Header]` can stop it from recognizing the sheet as v2 format
     at all, silently falling back to legacy v1 parsing (which expects a
     `[Data]` section instead of `[BCLConvert_Data]`) and failing with
     "File has no valid [Data] section". Keep it as plain, valid CSV.
   - If you're using legacy `bcl2fastq` instead of BCL Convert (not
     recommended for MiSeq i100 output - see the comment in
     `scripts/00_bcl_to_fastq.sh`), it expects a `[Data]` section (not
     `[BCLConvert_Data]`) with columns
     `Sample_ID,Sample_Name,index,index2,Sample_Project` and no
     `[BCLConvert_Settings]` section - a different file, not a variant of
     this template.
   - Type both `index` (i7) and `index2` (i5) as the plain 5'->3' primer
     sequence exactly as ordered/synthesized - do NOT manually
     reverse-complement either one. On many older 2-channel Illumina
     instruments (NextSeq, NovaSeq, iSeq, MiniSeq) you historically had to
     manually enter the i5 as its reverse complement, which trips a lot of
     people up; on the MiSeq i100 with a standard BCL Convert sample sheet
     like this one, that flip is handled automatically from a flag in the
     run's `RunInfo.xml`, so plain forward orientation for both indices is
     correct. Still worth a 30-second confirmation with GeneCore before a
     real run, since a flipped i5 sends most of your reads to
     `Undetermined` instead of a sample.
2. Check `config/samples.tsv` - the `sample_id` column must match the
   `Sample_ID`s in your Illumina sample sheet.
3. Run:
   ```bash
   ./run_pipeline.sh --run-dir /path/to/raw/MiSeqRun --ref data/binders.fasta
   ```
   or the individual `scripts/00_*`, `01_*`, `02_*` steps one at a time -
   each is a standalone script with `--help`-able arguments, useful for
   understanding/debugging one step at a time.
4. Open `analysis/DMS_analysis.Rmd` and knit it.

## Running on a SLURM cluster

`bcl-convert` in particular is a heavy, multi-threaded job (tens of minutes
on real sequencer output) - don't run it on a login/interactive node. Each
of `scripts/00_bcl_to_fastq.sh`, `scripts/01_qc_trim.sh`, and
`scripts/02_align_and_count.sh` carries its own `#SBATCH` directives right
after the shebang, so the same file works two ways:

```bash
bash scripts/00_bcl_to_fastq.sh --run-dir /path/to/raw/MiSeqRun    # runs directly
sbatch scripts/00_bcl_to_fastq.sh --run-dir /path/to/raw/MiSeqRun  # submits as a SLURM job
```

(`#SBATCH` lines are just `#`-comments to bash, so running a script directly
ignores them - only `sbatch` reads them.)

To run all three chained together with `--dependency=afterok` (so step 2
only starts once step 1 has actually finished successfully), use:

```bash
slurm/submit_pipeline.sh --run-dir /path/to/raw/MiSeqRun
```

This submits all three jobs at once (they just wait on each other in the
queue) and prints the job IDs and log file paths
(`logs/00_demux_<jobid>.out`, etc.) so you can `squeue -u $USER` and tail
the logs. Add `--force` if you're re-running over a previous demux output.
Resource requests are sized for this dataset (88 binders, 12 samples, a 5M
MiSeq i100 run) on the `htc-el8` partition - bump `--cpus-per-task`/`--mem`/
`--time` up in the relevant script's `#SBATCH` lines if you outgrow it, and
add an `--account` directive if your cluster requires one.

An interactive node (`salloc`/`srun --pty bash`) is fine for quick
debugging or running `scripts/simulate_test_data.py`, but isn't needed to
run the real pipeline steps themselves.

## Environment setup

Each script loads what it needs via cluster environment modules, right
after `set -euo pipefail` - no conda required on a cluster that provides
these:

- `scripts/00_bcl_to_fastq.sh`: `module load bcl-convert`
- `scripts/01_qc_trim.sh`: `module load fastp`
- `scripts/02_align_and_count.sh`: `module load BWA SAMtools Python Pysam`
  (`Python`/`Pysam` are for `count_reads_per_binder.py`, called at the end
  of that script - BWA/SAMtools alone don't cover it)

- **Analysis** (`analysis/DMS_analysis.Rmd`): R with `tidyverse`, `scales`,
  and `DESeq2` (Bioconductor, used for the count-based hit-calling in
  Section 5):
  ```r
  install.packages(c("tidyverse", "scales", "rmarkdown"))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
  ```

## The analysis logic (see `analysis/DMS_analysis.Rmd` for the full walkthrough)

This mirrors the lab's existing plate-based version of the same assay: the
same core comparison (toxin+library vs. toxin-only), with `no_inducer` and
`library_only` used as QC/context steps rather than as normalisation inputs.

1. Load `results/binder_counts.tsv` (long format: sample, binder, count) +
   `config/samples.tsv` metadata; reshape into a binder x sample table.
2. QC: read depth and number of binders detected per sample, a **library
   complexity check** - since this is a competitive pooled assay, a handful
   of binders dominating the pool squeezes everyone else into very few
   reads regardless of biology, so check for that before trusting anything
   downstream - and a **viability filter** - drop any binder whose mean read
   count under `no_inducer` (antibiotic only) falls below a threshold, since
   a plasmid that doesn't even survive plain selection isn't a usable data
   point for anything downstream.
3. Normalise: each binder's count -> relative frequency within its sample
   (total-read-depth normalisation), with a pseudocount so zero-count
   binders don't blow up on the log scale.
4. Innate toxicity (context only): log2FC of `library_only` vs. `no_inducer`
   at each IPTG level - how costly is overexpressing this binder on its own,
   with no toxin present? Reported for interpretation, not used to correct
   the main result.
5. **The actual answer**: a quick simple frequency-ratio log2FC of
   toxin+library vs. toxin-only first (no GFP/neutral-overexpression control
   is available in this sequencing-based run, so unlike the plate-based
   assay this ratio isn't further normalised against one - see the caveat in
   the notebook), then a proper count-based test (**DESeq2**) on the same
   comparison - it models counts with a negative binomial distribution and
   gives each binder an honestly wide or narrow uncertainty estimate
   depending on how much read depth actually backs it, rather than treating
   every point estimate as equally trustworthy the way the simple ratio
   does.
6. Replicate reproducibility: correlate normalised frequency between
   replicate A and B (DESeq2 already folds replicate variation into its
   model in step 5, so there's no separate manual consistency check needed
   for its result).
7. Final ranked hit table: the DESeq2 result ranked by q-value (or the
   simple ratio, marked provisional, if only one replicate exists yet),
   with binders averaging too few reads flagged (not excluded) as
   `low_count`, written to `results/final_binder_ranking.csv`.
8. Comparison to the collaborator's published plate-based screen of the same
   binders (`data/plate_assay_reference_hits.tsv`, filtered to the Sen2
   retron) - a rough sanity check on direction/ranking, not a strict
   validation (their effect size uses a different baseline), written to
   `results/plate_assay_comparison.csv`.
