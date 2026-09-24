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
  demux/                          per-sample demultiplexed fastq.gz land here
scripts/
  00_bcl_to_fastq.sh               demultiplex raw sequencer output -> data/demux/<sample_id>/
  01_qc_trim.sh                    fastp adapter/quality trimming
  02_align_and_count.sh            bwa mem alignment to data/binders.fasta + counting
  count_reads_per_binder.py        per-binder read-pair tallying used by 02_align_and_count.sh
  simulate_test_data.py            generate synthetic reads so you can dry-run the whole thing now
run_pipeline.sh                    convenience wrapper that runs 00->01->02 in order
environment.yml                    conda env for the demux/QC/alignment tools
analysis/
  DMS_analysis.Rmd                 the actual DMS-style analysis (normalisation, log2FC, replicate QC, hit table)
results/                           binder_counts.tsv and everything the R notebook produces (gitignored except final CSV)
```

## Running it

### 1. With real sequencing data

`data/binders.fasta` already contains the 88 designed binder sequences
(headers = their design names) and is used as the alignment reference
throughout.

1. Copy `config/IlluminaSampleSheet_template.csv` to
   `config/IlluminaSampleSheet.csv` and fill in the real i7/i5 index
   sequences for each of the 12 PCR-barcoded samples.
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

## Tool requirements

- **Demux/QC/alignment** (`scripts/`): see `environment.yml`
  (`conda env create -f environment.yml`). `bcl-convert` isn't on conda -
  install it separately from Illumina and make sure it's on your `PATH`.
- **Analysis** (`analysis/DMS_analysis.Rmd`): R with `tidyverse` and
  `scales`:
  ```r
  install.packages(c("tidyverse", "scales", "rmarkdown"))
  ```

## The analysis logic (see `analysis/DMS_analysis.Rmd` for the full walkthrough)

1. Load `results/binder_counts.tsv` (long format: sample, binder, count) +
   `config/samples.tsv` metadata; reshape into a binder x sample table.
2. QC: read depth and number of binders detected per sample.
3. Normalise: each binder's count -> relative frequency within its sample
   (total-read-depth normalisation), with a pseudocount so zero-count
   binders don't blow up on the log scale.
4. log2 fold-changes:
   - each induced condition vs. no-inducer, as a sanity check that
     induction is doing something at all;
   - **the actual answer**: toxin+library vs. toxin-only (at each IPTG
     level) - this isolates the protective effect of inducing the binder
     while the toxin is present, i.e. which binders rescue growth.
5. Replicate consistency: correlate A vs. B, both on normalised frequency
   and on the derived log2FC, and flag any binder where the two replicates
   disagree before trusting it as a hit.
6. Final ranked hit table + plot of the top consistently-enriched binders,
   written to `results/final_binder_ranking.csv`.
