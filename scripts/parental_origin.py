#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
from collections import defaultdict


REF_REGEX = re.compile(
    r"^(EHDV-[12])_s(\d+)$",
    re.IGNORECASE
)


# --------------------------------------------------
# Argument parsing
# --------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Infer per-segment parental origin (EHDV-1 vs EHDV-2) "
                    "from competitive alignment counts with confidence scoring."
    )
    p.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Input .per_reference.tsv files (one per sample)"
    )
    p.add_argument(
        "--min-total",
        type=int,
        default=100,
        help="Minimum total reads (EHDV-1 + EHDV-2) required per segment"
    )
    p.add_argument(
        "--min-diff",
        type=int,
        default=50,
        help="Minimum absolute read-count difference to call dominance"
    )
    p.add_argument(
        "--min-frac",
        type=float,
        default=0.65,
        help="Minimum fraction of reads required for dominant parent"
    )
    p.add_argument(
        "--out",
        required=True,
        help="Output TSV file"
    )
    return p.parse_args()

# --------------------------------------------------
# Helpers
# --------------------------------------------------

def infer_sample_name(path):
    return os.path.basename(path).replace(".per_reference.tsv", "")


def parse_reference(ref):
    """
    Parse references like:
      EHDV-1_s2
      EHDV-2_s10
    """
    m = REF_REGEX.match(ref.strip())
    if not m:
        return None, None

    virus = m.group(1).upper()
    segment = f"s{m.group(2)}"
    return virus, segment


def confidence_label(total, diff, frac):
    if total >= 500 and diff >= 200 and frac >= 0.80:
        return "HIGH"
    if total >= 200 and diff >= 100 and frac >= 0.70:
        return "MEDIUM"
    return "LOW"

# --------------------------------------------------
# Main
# --------------------------------------------------

def main():
    args = parse_args()

    # sample -> segment -> virus -> count
    counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for path in args.inputs:
        sample = infer_sample_name(path)
        with open(path) as fh:
            for line in fh:
                if not line.strip():
                    continue
                ref, val = line.rstrip().split("\t")[:2]
                try:
                    val = int(val)
                except ValueError:
                    continue

                virus, segment = parse_reference(ref)
                if virus is None:
                    continue

                counts[sample][segment][virus] += val

    with open(args.out, "w") as out:
        out.write(
            "\t".join([
                "sample",
                "segment",
                "ehdv1_reads",
                "ehdv2_reads",
                "total_reads",
                "dominant_parent",
                "dominant_fraction",
                "read_difference",
                "confidence",
                "note"
            ]) + "\n"
        )

        for sample in sorted(counts):
            for segment in sorted(counts[sample]):
                v1 = counts[sample][segment].get("EHDV-1", 0)
                v2 = counts[sample][segment].get("EHDV-2", 0)
                total = v1 + v2
                diff = abs(v1 - v2)

                dominant = "AMBIGUOUS"
                frac = 0.0
                conf = "LOW"
                note = ""

                if total < args.min_total:
                    note = "below_min_total_reads"
                else:
                    if v1 == v2:
                        note = "equal_counts"
                    else:
                        if v1 > v2:
                            dominant = "EHDV-1"
                            frac = v1 / total
                        else:
                            dominant = "EHDV-2"
                            frac = v2 / total

                        if diff < args.min_diff or frac < args.min_frac:
                            dominant = "AMBIGUOUS"
                            note = "dominance_below_thresholds"
                        else:
                            conf = confidence_label(total, diff, frac)

                out.write(
                    "\t".join(map(str, [
                        sample,
                        segment,
                        v1,
                        v2,
                        total,
                        dominant,
                        f"{frac:.3f}",
                        diff,
                        conf,
                        note
                    ])) + "\n"
                )

if __name__ == "__main__":
    main()