# RcaT binder pooled screen - a minimal DMS-style pipeline

A small, from-scratch pipeline for the EMBL PhD Primer project *"Designing
functional protein interactions: computational protein design and
high-throughput biochemistry"* (Dorrity/Mahamid/Eustermann/Typas groups, MSB Unit).

## The experiment

88 computationally-designed binders against the RcaT toxin are pooled and
expressed in *E. coli* carrying the RcaT toxin-antitoxin system. RcaT kills
the cell unless something inhibits it, so **a binder's relative abundance
after selection is a readout of how well it protects the cell**: the better
a binder blocks RcaT, the more its carrier cells survive, and the more reads
it contributes to the pool.

5 samples, one per condition:

| Condition | Ara (toxin inducer) | IPTG (binder inducer) |
|---|---|---|
| plasmid_library | - | - |
| no_inducer | - | - |
| toxin_only | 0.2% | - |
| library_only_0.1mM | - | 0.1 mM |
| toxin_plus_library_0.1mM | 0.2% | 0.1 mM |

`plasmid_library` is the pooled plasmid prep itself, sequenced directly -
no transformation or growth at all. It's there to check whether the pool is
already unevenly represented before a single cell is grown, or whether any
skew only shows up afterwards (see `analysis/DMS_analysis.Rmd` Section 2).

Sequencing: Illumina MiSeq i100 5M, 50bp paired-end, dual PCR barcoding for
demultiplexing.

## Repository layout

```
config/
  samples_template.tsv             sample sheet used by every downstream step
  IlluminaSampleSheet_template.csv template BCL Convert sample sheet for demux (fill in real indices)
  barcodes.txt                     reference table of the PCR barcode sequences used to fill in the sheet above
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
environment.yml                    conda env for the demux/QC/alignment tools (fallback if your cluster has no modules for them)
analysis/
  DMS_analysis.Rmd                 the actual DMS-style analysis (normalisation, fitness, hit table)
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
   sequences for each of the 5 PCR-barcoded samples. A few notes on this
   file:
   - Keep the `Sample_ID` values identical to `sample_id` in
     `config/samples_template.tsv` so downstream steps can find the right fastq
     files automatically.
   - Don't add comments or any other free text to this file - BCL Convert's
     sample sheet parser doesn't support a comment syntax, and stray text
     before `[Header]` can stop it from recognizing the sheet as v2 format
     at all, silently falling back to legacy v1 parsing (which expects a
     `[Data]` section instead of `[BCLConvert_Data]`) and failing with
     "File has no valid [Data] section". Keep it as plain, valid CSV.
   - Type both `index` (i7) as the plain 5'->3' primer sequence exactly
     as ordered/synthesised, and `index2` (i5) as the reverse complement
     sequence of what was ordered. This is what this version of BCL Convert
     is expecting, but always make sure to check based on the sequencing
     platform and program you're using!
2. Check `config/samples_template.tsv` - the `sample_id` column must match the
   `Sample_ID`s in your Illumina sample sheet.
3. Submit each step yourself, one at a time, checking the output before
   moving to the next. This can all be done on the login node, and make
   sure you're in the main directory of the GitHub repo:
   ```bash
   sbatch scripts/00_bcl_to_fastq.sh --run-dir /path/to/raw/MiSeqRun --sample-sheet config/IlluminaSampleSheet.csv --out-dir data/demux/
   # wait for it to finish, check data/demux/, then:
   sbatch scripts/01_qc_trim.sh --samples config/samples_template.tsv --out-dir data/qc_trimmed
   # wait for it to finish, check the fastp .html reports, then:
   sbatch scripts/02_align_and_count.sh --samples config/samples.tsv --trimmed-dir data/qc_trimmed --out-dir results --min-mapq 20--ref data/binders.fasta
   ```
5. Open `analysis/DMS_analysis.Rmd` and knit it.

## Running on a SLURM cluster

Submit each of the three steps yourself, one at a time, once the previous
one has actually finished - there's no chained/all-in-one submission script
here by design, since it's worth looking at each step's output before
moving on. Track progress with `squeue -u $USER`; logs land in
`logs/00_demux_<jobid>.out` etc. (and matching `.err` files). Add `--force`
to `00_bcl_to_fastq.sh` if you're re-running over a previous demux output.
Resource requests are sized for this dataset (88 binders, 5 samples, a 25M
MiSeq i100 run) on the `htc-el8` partition - bump `--cpus-per-task`/`--mem`/
`--time` up in the relevant script's `#SBATCH` lines if you outgrow it, and
add an `--account` directive if your cluster requires one.

## The analysis logic (see `analysis/DMS_analysis.Rmd` for the full walkthrough)

This mirrors the lab's existing plate-based version of the same assay.
There's a single sample per condition, so every fitness value below is a
one-off point estimate.

1. Load `results/binder_counts.tsv` (long format: sample, binder, count) +
   `config/samples_template.tsv` metadata; reshape into a binder x sample table.
2. **Library representation**: compare `plasmid_library` (the raw plasmid
   prep, no cells involved) against `no_inducer` (grown, unselected) - does
   any skew in the pool already exist before a single cell was grown, or
   does it only show up after growth?
3. QC: read depth and number of binders detected per sample, a **library
   complexity check** - since this is a competitive pooled assay, a handful
   of binders dominating the pool squeezes everyone else into very few
   reads regardless of biology, so check for that before trusting anything
   downstream - a **viability filter** - drop any binder whose read count
   under `no_inducer` (antibiotic only) falls below a threshold, since a
   plasmid that doesn't even survive plain selection isn't a usable data
   point for anything downstream - and a **low-count filter** - drop any
   binder whose mean read count across the remaining samples falls below a
   threshold (5 by default), since the viability filter alone doesn't catch
   a binder that's merely thin everywhere else.
4. Normalise: each binder's count -> relative frequency within its sample
   (total-read-depth normalisation), with a pseudocount so zero-count
   binders don't blow up on the log scale.
5. Innate toxicity (context only): log2FC of `library_only_0.1mM` vs. `no_inducer`
   - how costly is overexpressing this binder on its own, with no toxin
   present? Reported for interpretation, not used to correct the main
   result.
6. **The actual answer**: a DiMSum-style fitness score for toxin+library vs.
   **library-only** - a ratio-of-ratios relative to one reference binder
   (chosen by inspecting the QC plots; see the notebook) rather than a
   plain frequency ratio.
7. Final ranked hit table, written to `results/final_binder_ranking.csv`.
8. Comparison to the collaborator's published plate-based screen of the same
   binders (`data/plate_assay_reference_hits.tsv`, filtered to the Sen2
   retron) - a rough sanity check on direction/ranking, not a strict
   validation (their effect size uses a different baseline), written to
   `results/plate_assay_comparison.csv`.
