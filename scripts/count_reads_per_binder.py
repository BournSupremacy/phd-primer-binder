#!/usr/bin/env python3
"""Tally per-binder read-pair counts from aligned BAM files.

For each sample BAM (aligned with bwa mem against data/binders.fasta), count
one "hit" per read pair whose primary alignment is properly paired, mapped
with MAPQ >= --min-mapq, and where R1 and R2 agree on the same reference
sequence (binder). Read pairs that don't meet these criteria are tallied
separately as QC counters (unmapped, low MAPQ, mate mismatch) so you can spot
samples with a lot of ambiguous/multi-mapping reads (e.g. if some binders
share long identical scaffold regions).

Output: a single long-format TSV with columns
    sample_id  binder_id  count
plus a companion *_qc.tsv summarising total pairs / usable pairs / dropped
reasons per sample, both written next to --out.
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

import pysam


def count_bam(bam_path, min_mapq):
    """Return (binder_id -> count, qc_dict) for one sample's BAM.

    Two passes: first collect each mate's assigned reference by read name,
    then reconcile pairs. This deliberately doesn't assume R1 appears before
    R2 in the coordinate-sorted BAM - that holds for long binders (R1 near
    the 5' end, R2 near the 3', so R1's POS < R2's POS), but isn't guaranteed
    for short binders where the mates can overlap and swap positional order.
    """
    binder_counts = defaultdict(int)
    qc = defaultdict(int)
    r1_ref = {}
    r2_ref = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            qc["total_primary_alignments"] += 1

            if read.is_unmapped:
                qc["unmapped"] += 1
                continue
            if not read.is_proper_pair:
                qc["not_properly_paired"] += 1
                continue
            if read.mapping_quality < min_mapq:
                qc["low_mapq"] += 1
                continue

            if read.is_read1:
                r1_ref[read.query_name] = read.reference_name
            else:
                r2_ref[read.query_name] = read.reference_name

    for qname in set(r1_ref) | set(r2_ref):
        ref1 = r1_ref.get(qname)
        ref2 = r2_ref.get(qname)
        if ref1 is not None and ref2 is not None:
            if ref1 == ref2:
                binder_counts[ref1] += 1
                qc["counted_pairs_agree"] += 1
            else:
                qc["mate_reference_mismatch"] += 1
        else:
            # Only one mate survived filtering (or was mapped) - still
            # credit the pair using whichever mate we have.
            binder_counts[ref1 or ref2] += 1
            qc["counted_single_mate"] += 1

    qc["total_binder_hits"] = sum(binder_counts.values())
    return binder_counts, qc


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True, help="config/samples.tsv")
    p.add_argument("--bam-dir", required=True, help="directory of <sample_id>.bam files")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--out", required=True, help="output long-format counts TSV")
    args = p.parse_args()

    sample_ids = []
    with open(args.samples) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample_ids.append(row["sample_id"])

    qc_rows = []
    with open(args.out, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["sample_id", "binder_id", "count"])

        for sample_id in sample_ids:
            bam_path = os.path.join(args.bam_dir, f"{sample_id}.bam")
            if not os.path.isfile(bam_path):
                print(f"WARNING: no BAM for sample '{sample_id}' at {bam_path} - skipping", file=sys.stderr)
                continue

            binder_counts, qc = count_bam(bam_path, args.min_mapq)
            for binder_id, count in sorted(binder_counts.items()):
                writer.writerow([sample_id, binder_id, count])

            qc["sample_id"] = sample_id
            qc_rows.append(qc)
            print(
                f"{sample_id}: {qc.get('total_binder_hits', 0)} usable read pairs "
                f"across {len(binder_counts)} binders "
                f"(unmapped={qc.get('unmapped', 0)}, low_mapq={qc.get('low_mapq', 0)}, "
                f"mate_mismatch={qc.get('mate_reference_mismatch', 0)})"
            )

    qc_out = os.path.splitext(args.out)[0] + "_qc.tsv"
    qc_fields = sorted({k for row in qc_rows for k in row if k != "sample_id"})
    with open(qc_out, "w", newline="") as qc_f:
        writer = csv.DictWriter(qc_f, delimiter="\t", fieldnames=["sample_id"] + qc_fields)
        writer.writeheader()
        for row in qc_rows:
            writer.writerow(row)

    print(f"\nWrote counts to {args.out}")
    print(f"Wrote per-sample QC summary to {qc_out}")


if __name__ == "__main__":
    main()
