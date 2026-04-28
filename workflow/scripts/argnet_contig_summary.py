"""
argnet_contig_summary.py

Aggregate per-sample ARGNet contig-protein predictions into a single summary
table with contig abundance and context metadata.
"""

import os
import pandas as pd

samples = snakemake.params.samples
output_dir = snakemake.params.output_dir
out_file = snakemake.output.tsv


def safe_read(path: str, **kwargs) -> pd.DataFrame:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str, **kwargs)
    except Exception:
        return pd.DataFrame()


def parse_genomad_molecule_type(plas_sum_path: str, vir_sum_path: str) -> dict[str, str]:
    plas = safe_read(plas_sum_path)
    vir = safe_read(vir_sum_path)
    mol_map: dict[str, str] = {}
    if not plas.empty and "seq_name" in plas.columns:
        for cid in plas["seq_name"]:
            mol_map[str(cid)] = "plasmid"
    if not vir.empty and "seq_name" in vir.columns:
        for cid in vir["seq_name"]:
            mol_map.setdefault(str(cid), "phage")
    return mol_map


def parse_lca(lca_path: str) -> dict[str, str]:
    df = safe_read(lca_path, header=None)
    if df.empty or len(df.columns) < 4:
        return {}
    return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 3].astype(str)))


def parse_abundance(abund_path: str) -> dict[str, str]:
    df = safe_read(abund_path)
    if df.empty or "contig_id" not in df.columns:
        return {}
    return dict(zip(df["contig_id"].astype(str), df["trimmed_mean"].astype(str)))


def protein_to_contig(protein_id: str) -> str:
    token = str(protein_id).split()[0]
    return token.rsplit("_", 1)[0] if "_" in token else token


def parse_argnet(argnet_path: str) -> list[dict[str, str]]:
    df = safe_read(argnet_path)
    if df.empty:
        return []
    result = []
    for _, row in df.iterrows():
        if str(row.get("ARG_prediction", "")).strip().upper() != "ARG":
            continue
        protein_id = str(row.get("test_id", "")).strip()
        if not protein_id:
            continue
        result.append(
            {
                "protein_id": protein_id,
                "contig_id": protein_to_contig(protein_id),
                "resistance_category": str(row.get("resistance_category", "")).strip(),
                "probability": str(row.get("probability", "")).strip(),
            }
        )
    return result


rows: list[dict[str, str]] = []

for sample in samples:
    plas_sum = os.path.join(
        output_dir,
        "data",
        "genomad",
        sample,
        f"{sample}.contigs_summary",
        f"{sample}.contigs_plasmid_summary.tsv",
    )
    vir_sum = os.path.join(
        output_dir,
        "data",
        "genomad",
        sample,
        f"{sample}.contigs_summary",
        f"{sample}.contigs_virus_summary.tsv",
    )
    lca_path = os.path.join(output_dir, "data", "mmseqs", f"{sample}_lca.tsv")
    abund_path = os.path.join(output_dir, "data", "contig_abundance", f"{sample}_contig_abundance.tsv")
    argnet_path = os.path.join(output_dir, "data", "argnet_contigs", f"{sample}_argnet.tsv")

    mol_map = parse_genomad_molecule_type(plas_sum, vir_sum)
    tax_map = parse_lca(lca_path)
    abund_map = parse_abundance(abund_path)

    for hit in parse_argnet(argnet_path):
        cid = hit["contig_id"]
        rows.append(
            {
                "sample": sample,
                "protein_id": hit["protein_id"],
                "contig_id": cid,
                "trimmed_mean": abund_map.get(cid, ""),
                "molecule_type": mol_map.get(cid, "chromosome"),
                "taxonomy": tax_map.get(cid, ""),
                "resistance_category": hit["resistance_category"],
                "probability": hit["probability"],
            }
        )

columns = [
    "sample",
    "protein_id",
    "contig_id",
    "trimmed_mean",
    "molecule_type",
    "taxonomy",
    "resistance_category",
    "probability",
]

out_df = pd.DataFrame(rows, columns=columns) if rows else pd.DataFrame(columns=columns)
os.makedirs(os.path.dirname(out_file), exist_ok=True)
out_df.to_csv(out_file, sep="\t", index=False)
