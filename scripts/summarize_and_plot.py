#!/usr/bin/env python3

import re
import argparse
import os
import gzip
import pandas as pd

UNMAPPED_LABEL = "UNMAPPED"


def natural_key(s):
    return tuple(
        int(text) if text.isdigit() else text.lower()
        for text in re.split(r"(\d+)", str(s))
    )


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize alignment counts and generate percentage + raw-count heatmaps (with optional unmapped accounting)."
    )
    p.add_argument(
        "--mode",
        required=True,
        choices=["competitive", "per_reference"],
        help="Input mode: competitive (TSV per sample) or per_reference (TXT per sample+ref).",
    )
    p.add_argument("--top-n", type=int, default=20, help="Top N references by mean percent for heatmap columns.")
    p.add_argument(
        "--reads-dir",
        default=None,
        help="Directory containing input reads named {sample}.fastq.gz (or suffix below). If provided, total reads and unmapped will be computed.",
    )
    p.add_argument("--reads-suffix", default=".fastq.gz", help="Suffix for read files in reads-dir (default: .fastq.gz).")
    p.add_argument(
        "--include-unmapped",
        action="store_true",
        help="If set and totals are available, add an UNMAPPED pseudo-reference row per sample and include it in heatmaps.",
    )
    p.add_argument("--out-summary", required=True, help="Output TSV summary path.")
    p.add_argument("--out-heatmap-pct", required=True, help="Output percentage heatmap path (PNG).")
    p.add_argument("--out-heatmap-counts", required=True, help="Output raw-count heatmap path (PNG).")
    p.add_argument("inputs", nargs="+", help="Input count files.")
    return p.parse_args()


def ensure_outdir(path: str):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)


def count_reads_in_fastq_gz(path: str) -> int:
    """Count reads in a FASTQ(.gz) by counting lines / 4 (streaming, memory-safe)."""
    line_count = 0
    opener = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"

    with opener(path, mode) as fh:
        for _ in fh:
            line_count += 1

    return line_count // 4


def infer_samples_from_inputs(mode: str, files):
    samples = []
    if mode == "competitive":
        # reports/counts/{sample}.per_reference.tsv
        for f in files:
            base = os.path.basename(f)
            samples.append(base.replace(".per_reference.tsv", ""))
    else:
        # reports/counts/{sample}.{ref}.txt (sample may include dots; split at first dot)
        for f in files:
            base = os.path.basename(f).replace(".txt", "")
            if "." not in base:
                raise ValueError(f"Unexpected count filename (expected sample.ref.txt): {f}")
            sample, _ref = base.split(".", 1)
            samples.append(sample)

    return sorted(set(samples))


def read_counts_competitive(files):
    # Each file: {sample}.per_reference.tsv; format: reference<TAB>count
    rows = []
    for f in files:
        sample = os.path.basename(f).replace(".per_reference.tsv", "")
        if os.path.getsize(f) == 0:
            continue
        df = pd.read_csv(f, sep="\t", header=None, names=["reference", "aligned_reads"])
        df["sample"] = sample
        rows.append(df[["sample", "reference", "aligned_reads"]])

    if rows:
        out = pd.concat(rows, ignore_index=True)
    else:
        out = pd.DataFrame(columns=["sample", "reference", "aligned_reads"])

    out["aligned_reads"] = pd.to_numeric(out["aligned_reads"], errors="coerce").fillna(0).astype(int)
    return out


def read_counts_per_reference(files):
    # Each file: {sample}.{ref}.txt containing integer count
    rows = []
    for f in files:
        base = os.path.basename(f).replace(".txt", "")
        sample, ref = base.split(".", 1)
        with open(f) as fh:
            val = fh.read().strip()
        aligned = int(val) if val else 0
        rows.append({"sample": sample, "reference": ref, "aligned_reads": aligned})
    return pd.DataFrame(rows)


def compute_total_reads(samples, reads_dir, reads_suffix):
    totals = {}
    for s in samples:
        fp = os.path.join(reads_dir, f"{s}{reads_suffix}")
        if not os.path.exists(fp):
            raise FileNotFoundError(
                f"Could not find reads for sample '{s}' at: {fp}. Check --reads-dir and --reads-suffix."
            )
        totals[s] = count_reads_in_fastq_gz(fp)
    return totals


def add_unmapped_and_percents(df, samples, total_reads, include_unmapped, mode):
    """
    Adds:
      - mapped_reads (sum of aligned_reads per sample)
      - total_reads / unmapped_reads if totals provided
      - percent_of_mapped
      - percent_of_total (if totals provided)
    Optionally adds UNMAPPED pseudo-reference rows if totals provided and include_unmapped is True.

    NOTE: In per_reference mode, mapped_reads may double-count reads across references.
    """
    if df.empty:
        return pd.DataFrame(columns=[
            "sample", "reference", "aligned_reads",
            "mapped_reads", "total_reads", "unmapped_reads",
            "percent_of_mapped", "percent_of_total", "note"
        ])

    out = df.copy()
    mapped_by_sample = out.groupby("sample")["aligned_reads"].sum().to_dict()
    out["mapped_reads"] = out["sample"].map(mapped_by_sample).fillna(0).astype(int)

    if total_reads is not None:
        out["total_reads"] = out["sample"].map(total_reads).astype(int)
        out["unmapped_reads"] = (out["total_reads"] - out["mapped_reads"]).clip(lower=0).astype(int)
    else:
        out["total_reads"] = pd.NA
        out["unmapped_reads"] = pd.NA

    out["percent_of_mapped"] = out.apply(
        lambda r: (r["aligned_reads"] / r["mapped_reads"] * 100) if r["mapped_reads"] > 0 else 0.0,
        axis=1
    )

    def pct_total(r):
        if pd.isna(r["total_reads"]) or r["total_reads"] == 0:
            return 0.0
        return (r["aligned_reads"] / r["total_reads"] * 100)

    out["percent_of_total"] = out.apply(pct_total, axis=1)

    out["percent_of_mapped"] = out["percent_of_mapped"].round(2)
    out["percent_of_total"] = out["percent_of_total"].round(2)

    if include_unmapped and total_reads is not None:
        extra = []
        for s in samples:
            mapped = mapped_by_sample.get(s, 0)
            tot = total_reads.get(s, 0)
            unmap = max(tot - mapped, 0)
            extra.append({
                "sample": s,
                "reference": UNMAPPED_LABEL,
                "aligned_reads": unmap,
                "mapped_reads": mapped,
                "total_reads": tot,
                "unmapped_reads": unmap,
                "percent_of_mapped": 0.0,
                "percent_of_total": (unmap / tot * 100) if tot > 0 else 0.0,
                "note": ""
            })
        out = pd.concat([out, pd.DataFrame(extra)], ignore_index=True)

    if mode == "competitive":
        out["note"] = "mapped_reads are primary mapped reads (competitive, best-hit)"
    else:
        out["note"] = "mapped_reads may double-count across references (non-competitive)"

    return out


def make_heatmap(summary_df, value_col, out_png, top_n, title, cbar_label):
    import matplotlib
    matplotlib.use("Agg")  # headless-safe
    import matplotlib.pyplot as plt
    import seaborn as sns

    if summary_df.empty:
        plt.figure(figsize=(6, 2))
        plt.text(0.5, 0.5, "No mapped reads to plot", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(out_png, dpi=300)
        plt.close()
        return

    mat = summary_df.pivot(index="sample", columns="reference", values=value_col).fillna(0)

    has_unmapped = UNMAPPED_LABEL in mat.columns
    mat_no_unmapped = mat.drop(columns=[UNMAPPED_LABEL], errors="ignore")

    # select top refs by mean of chosen value_col; keep UNMAPPED at end if present
    top_refs = mat_no_unmapped.mean(axis=0).sort_values(ascending=False).head(top_n).index.tolist()
    
    # Alphanumeric sorting
    top_refs = sorted(top_refs, key=natural_key)
    
    # If UNMAPPED exists, keep it at the end
    if has_unmapped:
        top_refs.append(UNMAPPED_LABEL)

    top_refs = [c for c in top_refs if c in mat.columns]
    mat = mat[top_refs]

    fig_w = max(6, 1 + 0.5 * len(top_refs))
    fig_h = max(4, 0.35 * len(mat.index) + 2)

    plt.figure(figsize=(fig_w, fig_h))
    sns.heatmap(
        mat,
        cmap="viridis",
        cbar_kws={"label": cbar_label},
        linewidths=0.3,
        linecolor="white"
    )
    plt.title(title)
    plt.xlabel("Reference")
    plt.ylabel("Sample")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    args = parse_args()

    samples = infer_samples_from_inputs(args.mode, args.inputs)

    if args.mode == "competitive":
        counts_df = read_counts_competitive(args.inputs)
    else:
        counts_df = read_counts_per_reference(args.inputs)

    total_reads = None
    if args.reads_dir is not None:
        total_reads = compute_total_reads(samples, args.reads_dir, args.reads_suffix)

    summary = add_unmapped_and_percents(counts_df, samples, total_reads, args.include_unmapped, args.mode)

    if not summary.empty:
        summary = summary.sort_values(
          ["sample", "reference"],
          key=lambda col: col.map(natural_key) if col.name == "reference" else col
        ).reset_index(drop=True)

    # Write summary TSV
    ensure_outdir(args.out_summary)
    summary.to_csv(args.out_summary, sep="\t", index=False)

    # Percent heatmap column selection:
    # (Reverted behavior) use percent_of_total only when totals exist AND unmapped rows are included;
    # otherwise use percent_of_mapped.
    pct_col = "percent_of_total" if (total_reads is not None and args.include_unmapped) else "percent_of_mapped"

    ensure_outdir(args.out_heatmap_pct)
    ensure_outdir(args.out_heatmap_counts)

    make_heatmap(
        summary_df=summary,
        value_col=pct_col,
        out_png=args.out_heatmap_pct,
        top_n=args.top_n,
        title=f"Alignment heatmap (percent; {pct_col}) - mode: {args.mode}",
        cbar_label="Percent of reads"
    )

    make_heatmap(
        summary_df=summary,
        value_col="aligned_reads",
        out_png=args.out_heatmap_counts,
        top_n=args.top_n,
        title=f"Alignment heatmap (raw counts) - mode: {args.mode}",
        cbar_label="Aligned reads"
    )


if __name__ == "__main__":
    main()
