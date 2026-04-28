"""
argnet_mag_summary.py

Aggregate per-sample ARGNet MAG-protein predictions into a single summary table.
"""

import os
import pandas as pd

samples = snakemake.params.samples
output_dir = snakemake.params.output_dir
skip_gtdbtk = snakemake.params.skip_gtdbtk
out_file = snakemake.output.tsv


def safe_read_tsv(path: str, **kwargs) -> pd.DataFrame:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str, **kwargs)
    except Exception:
        return pd.DataFrame()


def parse_abundance(path: str) -> dict[str, str]:
    df = safe_read_tsv(path)
    if df.empty:
        return {}
    cov_col = next((c for c in df.columns[1:] if "Trimmed Mean" in c), None)
    if cov_col is None:
        return {}
    return dict(zip(df.iloc[:, 0].astype(str), df[cov_col].astype(str)))


def parse_gtdbtk(path: str) -> dict[str, str]:
    df = safe_read_tsv(path)
    if df.empty:
        return {}
    cols = {c.lower(): c for c in df.columns}
    user_col = cols.get("user_genome")
    class_col = cols.get("classification")
    if not user_col or not class_col:
        return {}
    return dict(zip(df[user_col].astype(str), df[class_col].astype(str)))


def parse_argnet(path: str) -> list[dict[str, str]]:
    df = safe_read_tsv(path)
    if df.empty:
        return []
    hits = []
    for _, row in df.iterrows():
        if str(row.get("ARG_prediction", "")).strip().upper() != "ARG":
            continue
        test_id = str(row.get("test_id", "")).strip()
        if "|" not in test_id:
            continue
        bin_id, protein_id = test_id.split("|", 1)
        contig_id = protein_id.rsplit("_", 1)[0] if "_" in protein_id else protein_id
        hits.append(
            {
                "bin_id": bin_id,
                "protein_id": protein_id,
                "contig_id": contig_id,
                "resistance_category": str(row.get("resistance_category", "")).strip(),
                "probability": str(row.get("probability", "")).strip(),
            }
        )
    return hits


rows: list[dict[str, str]] = []

for sample in samples:
    bins_dir = os.path.join(output_dir, "data", "bins", sample)
    abund_map = parse_abundance(os.path.join(bins_dir, "mag_abundance.tsv"))
    gtdbtk_map = {} if skip_gtdbtk else parse_gtdbtk(os.path.join(bins_dir, "gtdbtk.summary.tsv"))
    for hit in parse_argnet(os.path.join(bins_dir, "mag_argnet.tsv")):
        row = {
            "sample": sample,
            "bin_id": hit["bin_id"],
            "protein_id": hit["protein_id"],
            "contig_id": hit["contig_id"],
            "abundance_trimmed_mean": abund_map.get(hit["bin_id"], ""),
            "resistance_category": hit["resistance_category"],
            "probability": hit["probability"],
        }
        if not skip_gtdbtk:
            row["gtdbtk_classification"] = gtdbtk_map.get(hit["bin_id"], "")
        rows.append(row)

columns = [
    "sample",
    "bin_id",
    "protein_id",
    "contig_id",
    "abundance_trimmed_mean",
    "resistance_category",
    "probability",
]
if not skip_gtdbtk:
    columns.append("gtdbtk_classification")

out_df = pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)
os.makedirs(os.path.dirname(out_file), exist_ok=True)
out_df.to_csv(out_file, sep="\t", index=False)
