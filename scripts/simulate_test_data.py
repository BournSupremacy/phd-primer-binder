#!/usr/bin/env python3
"""Generate a small synthetic dataset so the pipeline can be dry-run end to end
before real sequencing data exists (or while waiting for GeneCore).

It creates:
  - data/binders.fasta                (if it doesn't already exist): N synthetic
    "binder" sequences of realistic length (300-450bp), a handful of which are
    designated as true strong RcaT binders.
  - data/demux/<sample_id>_R{1,2}.fastq.gz for every sample in
    config/samples.tsv, with per-binder relative abundances that reflect the
    expected biology:
      * no_inducer / library_only (any IPTG, no toxin): all binders ~uniform
        (growth isn't toxin-challenged, so binder identity shouldn't matter)
      * toxin_only (Ara, no IPTG so binder isn't expressed): also ~uniform,
        close to the no_inducer baseline (nothing rescues without induction)
      * toxin_plus_library (Ara + IPTG): true binders enriched relative to
        non-binders, more so at 1 mM than 0.1 mM IPTG
    plus multinomial sampling noise and a small amount of extra rep-to-rep
    variability, so the simulated data exercises the QC/consistency steps
    in analysis/DMS_analysis.Rmd too.

This is a teaching aid, not a generative model of real sequencing data -
reads have no flanking cloning sequence and errors are a simple per-base
substitution rate, which is enough to make alignment/counting non-trivial.

Usage: scripts/simulate_test_data.py [--n-binders 96] [--n-true-binders 6]
         [--reads-per-sample 60000] [--seed 0]
"""
import argparse
import gzip
import os
import random

import numpy as np

BASES = "ACGT"
SAMPLES_TSV = "config/samples.tsv"
DEFAULT_REF = "data/binders.fasta"
DEMUX_DIR = "data/demux"
READ_LEN = 150


def random_seq(length, rng):
    return "".join(rng.choice(BASES) for _ in range(length))


def revcomp(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def make_binders_fasta(path, n_binders, rng):
    with open(path, "w") as f:
        for i in range(1, n_binders + 1):
            length = rng.randint(300, 450)
            seq = random_seq(length, rng)
            f.write(f">binder_{i:03d}\n{seq}\n")
    print(f"Wrote {n_binders} synthetic binder sequences to {path}")


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


def read_samples_tsv(path):
    rows = []
    with open(path) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def relative_abundance(binder_ids, true_binder_ids, condition, iptg_mM, rng):
    """Return a probability vector over binder_ids for a given condition."""
    n = len(binder_ids)
    base = np.ones(n)

    is_toxin_plus_library = condition.startswith("toxin_plus_library")
    if is_toxin_plus_library:
        enrichment = 8.0 if float(iptg_mM) >= 1 else 3.5
        for idx, b in enumerate(binder_ids):
            if b in true_binder_ids:
                base[idx] *= enrichment

    # small multiplicative noise per binder to mimic biological variability
    noise = rng.lognormal(mean=0.0, sigma=0.15, size=n)
    weights = base * noise
    return weights / weights.sum()


def simulate_reads_for_binder(seq, n_reads, error_rate, rng, np_rng):
    """Yield (r1_seq, r2_seq) pairs sampled from the two ends of `seq`."""
    seq_len = len(seq)
    for _ in range(n_reads):
        r1 = seq[: min(READ_LEN, seq_len)].ljust(READ_LEN, "A")[:READ_LEN]
        r2_source = seq[max(0, seq_len - READ_LEN):]
        r2 = revcomp(r2_source).ljust(READ_LEN, "A")[:READ_LEN]
        r1 = introduce_errors(r1, error_rate, np_rng)
        r2 = introduce_errors(r2, error_rate, np_rng)
        yield r1, r2


def introduce_errors(seq, error_rate, np_rng):
    seq = list(seq)
    for i in range(len(seq)):
        if np_rng.random() < error_rate:
            seq[i] = np_rng.choice([b for b in BASES if b != seq[i]])
    return "".join(seq)


def write_fastq_pair(sample_dir, sample_id, pairs, rng):
    r1_path = os.path.join(sample_dir, f"{sample_id}_R1.fastq.gz")
    r2_path = os.path.join(sample_dir, f"{sample_id}_R2.fastq.gz")
    qual = "I" * READ_LEN
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for i, (r1, r2) in enumerate(pairs):
            read_name = f"@SIM:{sample_id}:{i}"
            f1.write(f"{read_name} 1\n{r1}\n+\n{qual}\n")
            f2.write(f"{read_name} 2\n{r2}\n+\n{qual}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ref", default=DEFAULT_REF)
    p.add_argument("--n-binders", type=int, default=96)
    p.add_argument("--n-true-binders", type=int, default=6)
    p.add_argument("--reads-per-sample", type=int, default=60000)
    p.add_argument("--error-rate", type=float, default=0.01)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    rng = random.Random(args.seed)
    np_rng = np.random.default_rng(args.seed)

    if not os.path.isfile(args.ref):
        make_binders_fasta(args.ref, args.n_binders, rng)
    else:
        print(f"Using existing reference at {args.ref} (delete it to regenerate)")

    binders = read_fasta(args.ref)
    binder_ids = sorted(binders.keys())
    true_binder_ids = set(rng.sample(binder_ids, min(args.n_true_binders, len(binder_ids))))
    print(f"Designated true binders for this simulation: {sorted(true_binder_ids)}")

    samples = read_samples_tsv(SAMPLES_TSV)

    for row in samples:
        sample_id = row["sample_id"]
        condition = row["condition"]
        iptg = row["iptg_mM"]

        probs = relative_abundance(binder_ids, true_binder_ids, condition, iptg, np_rng)
        counts = np_rng.multinomial(args.reads_per_sample, probs)

        os.makedirs(DEMUX_DIR, exist_ok=True)

        pairs = []
        for binder_id, n in zip(binder_ids, counts):
            if n == 0:
                continue
            pairs.extend(simulate_reads_for_binder(binders[binder_id], int(n), args.error_rate, rng, np_rng))
        np_rng.shuffle(np.arange(len(pairs)))  # cheap shuffle of order isn't essential; skip heavy shuffle for speed

        write_fastq_pair(DEMUX_DIR, sample_id, pairs, rng)
        print(f"  {sample_id} ({condition}, IPTG={iptg}mM): {len(pairs)} read pairs -> {DEMUX_DIR}/")

    print("\nSimulated data ready. You can now run:")
    print("  scripts/01_qc_trim.sh")
    print("  scripts/02_align_and_count.sh --ref", args.ref)
    print("then knit analysis/DMS_analysis.Rmd")


if __name__ == "__main__":
    main()
