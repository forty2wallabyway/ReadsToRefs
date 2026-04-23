#!/usr/bin/env python3

"""
Final publication-quality visualization for EHDV competition analysis.

Inputs
------
- alignment_summary_percentages.tsv
  Output from summarize_and_plot.py

Outputs
-------
- Figure1_EHDV_competition.png
  Three-panel figure:
    A) Genome-level dominance (fractions)
    B) Within-sample competition (log2 ratios)
    C) Absolute viral burden (log scale)

Notes
-----
- Uses a headless-safe matplotlib backend (Agg).
- Designed to run on HPC / CI / containerized environments.
"""

# -----------------------------
# Headless-safe plotting
# -----------------------------

import matplotlib
matplotlib.use("Agg")  # REQUIRED for non-interactive environments

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# -----------------------------
# Configuration
# -----------------------------

INPUT_TSV = Path("alignment_summary_percentages.tsv")
OUTPUT_PNG = Path("Figure1_EHDV_competition.png")

# Color-blind friendly palette
COLOR_EHDV1 = "#4C72B0"  # blue
COLOR_EHDV2 = "#DD8452"  # orange

# -----------------------------
# Helper functions
# -----------------------------

def classify_condition(sample_name: str) -> str:
    """
    Assign high-level experimental condition from sample name.
    Adjust mapping if sample naming conventions change.
    """
    if "and-2" in sample_name:
        return "Mixed"
    if "-to-" in sample_name or "Stag" in sample_name:
        return "Serial passage"
    if sample_name.startswith("EHDV-1"):
        return "Single EHDV-1"
    if sample_name.startswith("EHDV-2"):
        return "Single EHDV-2"
    return "Other"


def extract_virus(reference: str):
    """Infer virus identity from reference name."""
    if reference.startswith("EHDV-1"):
        return "EHDV-1"
    if reference.startswith("EHDV-2"):
        return "EHDV-2"
    return None


# -----------------------------
# Load and preprocess data
# -----------------------------

if not INPUT_TSV.exists():
    raise FileNotFoundError(f"Required input not found: {INPUT_TSV}")

df = pd.read_csv(INPUT_TSV, sep="\t")

# Drop unmapped rows
df = df[df["reference"] != "UNMAPPED"].copy()

# Assign virus and condition
df["virus"] = df["reference"].apply(extract_virus)
df = df.dropna(subset=["virus"])
df["condition"] = df["sample"].apply(classify_condition)

# Aggregate aligned reads per sample x virus
summary = (
    df.groupby(["sample", "condition", "virus"], as_index=False)
      .agg(aligned_reads=("aligned_reads", "sum"))
)

# Pivot for plotting convenience
pivot = (
    summary.pivot(index=["sample", "condition"],
                  columns="virus",
                  values="aligned_reads")
    .fillna(0)
    .reset_index()
)

# Ensure expected columns exist
for col in ["EHDV-1", "EHDV-2"]:
    if col not in pivot.columns:
        pivot[col] = 0

# Totals and derived metrics
pivot["total_mapped"] = pivot["EHDV-1"] + pivot["EHDV-2"]
pivot["frac_EHDV1"] = pivot["EHDV-1"] / pivot["total_mapped"].replace(0, np.nan)
pivot["frac_EHDV2"] = pivot["EHDV-2"] / pivot["total_mapped"].replace(0, np.nan)

pivot["log2_ratio"] = np.log2(
    (pivot["EHDV-1"] + 1) / (pivot["EHDV-2"] + 1)
)

# -----------------------------
# Figure layout
# -----------------------------

fig = plt.figure(figsize=(14, 5))
gs = fig.add_gridspec(1, 3, width_ratios=[2.5, 1.5, 2.0])

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[0, 2])

# -----------------------------
# Panel A: Stacked fractions
# -----------------------------

order = (
    pivot.groupby("sample")["condition"]
    .first()
    .map({
        "Single EHDV-1": 0,
        "Single EHDV-2": 1,
        "Mixed": 2,
        "Serial passage": 3,
        "Other": 4
    })
    .sort_values()
    .index
)

dataA = pivot.set_index("sample").loc[order]

x = np.arange(len(dataA))
axA.bar(x, dataA["frac_EHDV1"], color=COLOR_EHDV1, label="EHDV-1")
axA.bar(
    x,
    dataA["frac_EHDV2"],
    bottom=dataA["frac_EHDV1"],
    color=COLOR_EHDV2,
    label="EHDV-2",
)

axA.set_ylabel("Fraction of mapped reads")
axA.set_title("A. Genome-level dominance")
axA.set_xticks(x)
axA.set_xticklabels(dataA.index, rotation=90, fontsize=8)
axA.legend(frameon=False)

# -----------------------------
# Panel B: Log2 ratios (mixed)
# -----------------------------

mixed = pivot[pivot["condition"] == "Mixed"]

axB.axhline(0, color="black", lw=1, linestyle="--")

if not mixed.empty:
    jitter = np.random.normal(0, 0.03, size=len(mixed))
    axB.scatter(jitter, mixed["log2_ratio"],
                s=60, color=COLOR_EHDV1, alpha=0.7)

axB.set_xlim(-0.5, 0.5)
axB.set_xticks([0])
axB.set_xticklabels(["Mixed"])
axB.set_ylabel("log2(EHDV-1 / EHDV-2)")
axB.set_title("B. Within-sample competition")

# -----------------------------
# Panel C: Absolute burden
# -----------------------------

melted = pivot.melt(
    id_vars=["sample", "condition"],
    value_vars=["EHDV-1", "EHDV-2"],
    var_name="virus",
    value_name="aligned_reads",
)

for virus, color, xoff in [
    ("EHDV-1", COLOR_EHDV1, -0.15),
    ("EHDV-2", COLOR_EHDV2, 0.15),
]:
    subset = melted[melted["virus"] == virus]
    axC.scatter(
        np.full(len(subset), xoff),
        subset["aligned_reads"],
        s=50,
        color=color,
        alpha=0.7,
        label=virus,
    )

axC.set_yscale("log")
axC.set_xticks([0])
axC.set_xticklabels(["All conditions"])
axC.set_ylabel("Total aligned reads (log10)")
axC.set_title("C. Absolute viral burden")
axC.legend(frameon=False)

# -----------------------------
# Save figure
# -----------------------------

plt.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=300)
plt.close(fig)

print(f"[final_processing] Saved figure to: {OUTPUT_PNG.resolve()}")
