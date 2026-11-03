"""
argnet_short_reads_summary.py

Aggregate per-sample ARGNet short-read predictions into a category-level summary.
These values are intentionally read-level metrics rather than cpg-normalized
abundance estimates, so they can sit alongside the AMRFinder short-read outputs
without implying direct equivalence.
"""

import os

import pandas as pd

samples = snakemake.params.samples
output_dir = snakemake.params.output_dir
out_file = snakemake.output.tsv


def safe_read(path: str) -> pd.DataFrame:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame()


rows: list[dict[str, object]] = []

for sample in samples:
    path = os.path.join(output_dir, "data", "argnet_short_reads", f"{sample}_argnet.tsv")
    df = safe_read(path)
    if df.empty:
        rows.append(
            {
                "sample": sample,
                "resistance_category": "none_detected",
                "arg_read_count": 0,
                "total_reads_scored": 0,
                "total_arg_reads": 0,
                "non_arg_reads": 0,
                "short_fragment_reads": 0,
                "arg_fraction": 0.0,
                "arg_reads_per_million": 0.0,
                "mean_probability": 0.0,
                "max_probability": 0.0,
            }
        )
        continue

    prediction = df.get("ARG_prediction", pd.Series(dtype=str)).fillna("")
    total_reads_scored = int(len(df))
    total_arg_reads = int((prediction == "ARG").sum())
    non_arg_reads = int((prediction == "non-ARG").sum())
    short_fragment_reads = int((prediction == "<25aa").sum())

    arg_df = df[prediction == "ARG"].copy()
    if arg_df.empty:
        rows.append(
            {
                "sample": sample,
                "resistance_category": "none_detected",
                "arg_read_count": 0,
                "total_reads_scored": total_reads_scored,
                "total_arg_reads": total_arg_reads,
                "non_arg_reads": non_arg_reads,
                "short_fragment_reads": short_fragment_reads,
                "arg_fraction": 0.0,
                "arg_reads_per_million": 0.0,
                "mean_probability": 0.0,
                "max_probability": 0.0,
            }
        )
        continue

    arg_df["resistance_category"] = arg_df["resistance_category"].fillna("").replace("", "unclassified")
    arg_df["probability"] = pd.to_numeric(arg_df["probability"], errors="coerce")

    for category, group in arg_df.groupby("resistance_category", sort=True):
        arg_read_count = int(len(group))
        rows.append(
            {
                "sample": sample,
                "resistance_category": category,
                "arg_read_count": arg_read_count,
                "total_reads_scored": total_reads_scored,
                "total_arg_reads": total_arg_reads,
                "non_arg_reads": non_arg_reads,
                "short_fragment_reads": short_fragment_reads,
                "arg_fraction": round(arg_read_count / total_reads_scored, 6) if total_reads_scored else 0.0,
                "arg_reads_per_million": round(arg_read_count * 1_000_000 / total_reads_scored, 2) if total_reads_scored else 0.0,
                "mean_probability": round(float(group["probability"].mean()), 6) if group["probability"].notna().any() else 0.0,
                "max_probability": round(float(group["probability"].max()), 6) if group["probability"].notna().any() else 0.0,
            }
        )

columns = [
    "sample",
    "resistance_category",
    "arg_read_count",
    "total_reads_scored",
    "total_arg_reads",
    "non_arg_reads",
    "short_fragment_reads",
    "arg_fraction",
    "arg_reads_per_million",
    "mean_probability",
    "max_probability",
]

out_df = pd.DataFrame(rows, columns=columns)
os.makedirs(os.path.dirname(out_file), exist_ok=True)
out_df.to_csv(out_file, sep="\t", index=False)
