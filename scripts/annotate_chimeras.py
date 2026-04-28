#!/usr/bin/env python3
# Indentation: 4 spaces only (no tabs)

import argparse
import os
import re
import subprocess
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate parental origin table with chimera + coverage metrics (coverage computed from BAMs)."
    )
    p.add_argument("--parental", required=True, help="Input parental_origin.tsv")
    p.add_argument("--bam-dir", required=True, help="Directory containing {sample}.sorted.bam")
    p.add_argument("--split", nargs="+", required=True, help="Split chimera read-ID TSVs")
    p.add_argument("--primer", nargs="+", required=True, help="Internal primer read-ID TSVs")
    p.add_argument("--out", required=True, help="Output annotated TSV")
    return p.parse_args()


def safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return default


def safe_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


def extract_sample_name(path):
    """
    Extract sample name from files like:
      sample.split_reads.tsv
      sample.internal_primers.tsv
    """
    fname = os.path.basename(path)
    fname = re.sub(r"\.(split_reads|internal_primers)\.tsv$", "", fname)
    if fname.endswith(".tsv"):
        fname = fname[:-4]
    return fname


def load_read_ids(files):
    """
    Load read IDs per sample; first token per line is treated as read ID.
    """
    out = defaultdict(set)
    for f in files:
        sample = extract_sample_name(f)
        with open(f) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                out[sample].add(line.split()[0])
    return out


def map_reads_to_segments(bam_path, read_ids):
    """
    Map read IDs back to segments using samtools view.
    Returns counts[segment] = number of reads mapping to that segment.
    """
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"Missing BAM: {bam_path}")

    counts = defaultdict(int)

    proc = subprocess.Popen(
        ["samtools", "view", bam_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    stdout, _ = proc.communicate()
    stdout = stdout or ""

    for line in stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        qname = fields[0]
        rname = fields[2]

        if qname not in read_ids:
            continue
        if rname == "*":
            continue

        m = re.search(r"(s\d+)", rname)
        if not m:
            continue
        segment = m.group(1)
        counts[segment] += 1

    return counts


def classify_primer_fraction(frac):
    """
    Primer-driven chimera classification.
    """
    if frac < 0.02:
        return "low"
    if frac < 0.10:
        return "moderate"
    return "high"


def compute_coverage_from_bam(bam_path):
    """
    Run 'samtools coverage' on the BAM and parse:
      cov[virus][segment] = (pct_covered, mean_depth)

    virus: EHDV-1 or EHDV-2 inferred from rname (if present), else 'unknown'
    segment: s1..s10 inferred from rname using regex.
    """
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"Missing BAM: {bam_path}")

    cov = defaultdict(dict)

    proc = subprocess.Popen(
        ["samtools", "coverage", bam_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    stdout, stderr = proc.communicate()
    stdout = stdout or ""
    stderr = stderr or ""

    # Parse header-aware if present
    header = None
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue

        # Some versions prefix header with '#'
        if line.startswith("#"):
            header = line.lstrip("#").split("\t")
            continue

        parts = line.split("\t")
        if header and len(parts) == len(header):
            rec = dict(zip(header, parts))
            rname = rec.get("rname", rec.get("#rname", rec.get("ref", "")))
            pct = safe_float(rec.get("coverage", "0"), default=0.0)
            depth = safe_float(rec.get("meandepth", "0"), default=0.0)
        else:
            # Typical order:
            # rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
            rname = parts[0] if len(parts) >= 1 else ""
            pct = safe_float(parts[5], default=0.0) if len(parts) >= 6 else 0.0
            depth = safe_float(parts[6], default=0.0) if len(parts) >= 7 else 0.0

        if not rname or rname == "*":
            continue

        seg_m = re.search(r"(s\d+)", rname)
        if not seg_m:
            continue
        segment = seg_m.group(1)

        virus_m = re.search(r"(EHDV-1|EHDV-2)", rname)
        virus = virus_m.group(1) if virus_m else "unknown"

        cov[virus][segment] = (pct, depth)

    return cov


def main():
    args = parse_args()

    # Load chimera read IDs
    split_reads = load_read_ids(args.split)
    primer_reads = load_read_ids(args.primer)

    # Read parental file once to get list of samples we must compute coverage for
    parental_rows = []
    with open(args.parental) as inp:
        base_header = inp.readline().rstrip("\n").split("\t")
        for line in inp:
            line = line.rstrip("\n")
            if not line:
                continue
            vals = line.split("\t")
            if len(vals) < len(base_header):
                vals = vals + [""] * (len(base_header) - len(vals))
            row = dict(zip(base_header, vals))
            parental_rows.append((line, row))

    parental_samples = set(r[1].get("sample", "") for r in parental_rows if r[1].get("sample", ""))
    chimera_samples = set(split_reads) | set(primer_reads)
    all_samples = parental_samples | chimera_samples

    # Precompute chimera counts (segment-resolved)
    split_counts = {}
    primer_counts = {}
    coverage = {}

    for sample in all_samples:
        bam = os.path.join(args.bam_dir, f"{sample}.sorted.bam")

        if sample in split_reads:
            split_counts[sample] = map_reads_to_segments(bam, split_reads.get(sample, set()))
        else:
            split_counts[sample] = {}

        if sample in primer_reads:
            primer_counts[sample] = map_reads_to_segments(bam, primer_reads.get(sample, set()))
        else:
            primer_counts[sample] = {}

        coverage[sample] = compute_coverage_from_bam(bam)

    # Write output
    extra_cols = [
        "chimera_split_reads",
        "chimera_primer_reads",
        "primer_fraction",
        "primer_per_1k_assigned",
        "chimera_flag",
        "pct_ref_covered_dom",
        "mean_depth_dom",
        "pct_ref_covered_ehdv1",
        "mean_depth_ehdv1",
        "pct_ref_covered_ehdv2",
        "mean_depth_ehdv2",
    ]

    out_header = base_header + extra_cols

    with open(args.out, "w") as out:
        out.write("\t".join(out_header) + "\n")

        for raw_line, row in parental_rows:
            sample = row.get("sample", "")
            segment = row.get("segment", "")
            total = safe_int(row.get("total_reads", "0"), default=0)

            split_n = split_counts.get(sample, {}).get(segment, 0)
            primer_n = primer_counts.get(sample, {}).get(segment, 0)

            primer_frac = (primer_n / total) if total > 0 else 0.0
            primer_per_1k = primer_frac * 1000.0
            chimera_flag = classify_primer_fraction(primer_frac)

            dom = row.get("dominant_parent", "")
            dom_virus = dom if dom in ("EHDV-1", "EHDV-2") else None

            # Pull coverage for both viruses
            c1 = coverage.get(sample, {}).get("EHDV-1", {}).get(segment, None)
            c2 = coverage.get(sample, {}).get("EHDV-2", {}).get(segment, None)

            if c1:
                pct1, depth1 = c1
                pct1_s = f"{pct1:.2f}"
                depth1_s = f"{depth1:.2f}"
            else:
                pct1_s, depth1_s = "NA", "NA"

            if c2:
                pct2, depth2 = c2
                pct2_s = f"{pct2:.2f}"
                depth2_s = f"{depth2:.2f}"
            else:
                pct2_s, depth2_s = "NA", "NA"

            # Coverage for dominant parent
            if dom_virus == "EHDV-1" and c1:
                pctd_s, depthd_s = pct1_s, depth1_s
            elif dom_virus == "EHDV-2" and c2:
                pctd_s, depthd_s = pct2_s, depth2_s
            else:
                pctd_s, depthd_s = "NA", "NA"

            out.write(
                raw_line + "\t" +
                f"{split_n}\t{primer_n}\t"
                f"{primer_frac:.4f}\t"
                f"{primer_per_1k:.1f}\t"
                f"{chimera_flag}\t"
                f"{pctd_s}\t{depthd_s}\t"
                f"{pct1_s}\t{depth1_s}\t"
                f"{pct2_s}\t{depth2_s}\n"
            )


if __name__ == "__main__":
    main()