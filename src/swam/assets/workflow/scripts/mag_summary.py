"""
mag_summary.py

Snakemake script (called via script: directive).

Aggregates per-sample MAG-level annotations into a single summary table with
one row per (sample, bin):

    sample | bin_id | abundance_trimmed_mean |
    n_amr_genes | amr_genes | n_mge | mge_elements |
    completeness | contamination | quality_score |   # CheckM2 (optional)
    gtdbtk_classification                            # GTDB-Tk  (optional)

Optional columns are included only when the corresponding tool was not skipped.
"""

import os
import pandas as pd

# ---------------------------------------------------------------------------
# Snakemake-injected variables
# ---------------------------------------------------------------------------
samples        = snakemake.params.samples
output_dir     = snakemake.params.output_dir
skip_checkm2   = snakemake.params.skip_checkm2
skip_gtdbtk    = snakemake.params.skip_gtdbtk
out_file       = snakemake.output.tsv


def safe_read_tsv(path, **kwargs):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str, **kwargs)
    except Exception:
        return pd.DataFrame()


def safe_read_csv(path, **kwargs):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, dtype=str, **kwargs)
    except Exception:
        return pd.DataFrame()


def parse_abundance(path):
    """Return dict {bin_id -> trimmed_mean} from CoverM genome output."""
    df = safe_read_tsv(path)
    if df.empty:
        return {}
    # CoverM outputs: Genome | {sample}_contigs Trimmed Mean
    # The bin_id column is named "Genome"
    bin_col = df.columns[0]
    cov_col = next((c for c in df.columns[1:] if "Trimmed Mean" in c), None)
    if cov_col is None:
        return {}
    return dict(zip(df[bin_col].astype(str), pd.to_numeric(df[cov_col], errors="coerce")))


def parse_amr(path):
    """Return dict {bin_id -> list of gene symbols} from mag_amr.tsv."""
    df = safe_read_tsv(path)
    if df.empty:
        return {}
    cols = {c.lower().replace(" ", "_"): c for c in df.columns}

    def gcol(*candidates):
        for c in candidates:
            if c in cols:
                return cols[c]
        return None

    bin_col  = gcol("contig_id", "contig", "sequence_name", "sequence")
    gene_col = gcol("gene_symbol", "element_symbol", "gene", "name")
    result: dict[str, list] = {}
    for _, row in df.iterrows():
        # bin_id is the contig name; strip the trailing gene index to get bin
        # AMRFinderPlus on bins uses the bin FASTA name as sequence identifier
        raw = str(row[bin_col]) if bin_col else ""
        # bin name is typically "bin.1", "bin.2", etc. (from MetaBAT2)
        # the contig_id in mag context is just the bin fasta name
        bid = raw
        gene = str(row[gene_col]) if gene_col else ""
        result.setdefault(bid, []).append(gene)
    return result


def parse_mge(path):
    """Return dict {bin_id -> list of MGE names} from mag_mge.tsv (CSV with comments)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    try:
        df = pd.read_csv(path, sep=",", comment="#", dtype=str)
    except Exception:
        return {}
    if df.empty:
        return {}
    cols_lower = {c.lower(): c for c in df.columns}

    def gcol(*candidates):
        for c in candidates:
            if c in cols_lower:
                return cols_lower[c]
        return None

    seq_col  = gcol("contig", "sequence", "seqid")
    name_col = gcol("name", "element_symbol", "type", "element")
    result: dict[str, list] = {}
    for _, row in df.iterrows():
        raw = str(row[seq_col]).split()[0] if seq_col else ""
        name = str(row[name_col]) if name_col else ""
        result.setdefault(raw, []).append(name)
    return result


def parse_checkm2(path):
    """Return dict {bin_id -> {completeness, contamination, quality_score}}."""
    df = safe_read_tsv(path)
    if df.empty:
        return {}
    cols = {c.lower(): c for c in df.columns}
    name_col  = cols.get("name", cols.get("bin", None))
    comp_col  = cols.get("completeness", None)
    cont_col  = cols.get("contamination", None)
    score_col = cols.get("completeness_model_used", None)
    result = {}
    for _, row in df.iterrows():
        if not name_col:
            break
        bid = str(row[name_col])
        result[bid] = {
            "completeness":   str(row[comp_col])  if comp_col  else "",
            "contamination":  str(row[cont_col])  if cont_col  else "",
        }
    return result


def parse_gtdbtk(path):
    """Return dict {bin_id -> classification string} from gtdbtk.summary.tsv."""
    df = safe_read_tsv(path)
    if df.empty:
        return {}
    cols = {c.lower(): c for c in df.columns}
    user_col  = cols.get("user_genome", None)
    class_col = cols.get("classification", None)
    if not user_col or not class_col:
        return {}
    return dict(zip(df[user_col].astype(str), df[class_col].astype(str)))


# ---------------------------------------------------------------------------
# Main aggregation
# ---------------------------------------------------------------------------
all_rows = []

for sample in samples:
    bins_dir      = os.path.join(output_dir, "data", "bins", sample)
    abund_path    = os.path.join(bins_dir, "mag_abundance.tsv")
    amr_path      = os.path.join(bins_dir, "mag_amr.tsv")
    mge_path      = os.path.join(bins_dir, "mag_mge.tsv")
    checkm2_path  = os.path.join(bins_dir, "checkm2_quality.tsv")
    gtdbtk_path   = os.path.join(bins_dir, "gtdbtk.summary.tsv")

    abund_map   = parse_abundance(abund_path)
    amr_map     = parse_amr(amr_path)
    mge_map     = parse_mge(mge_path)
    checkm2_map = {} if skip_checkm2 else parse_checkm2(checkm2_path)
    gtdbtk_map  = {} if skip_gtdbtk  else parse_gtdbtk(gtdbtk_path)

    all_bins = set(abund_map.keys()) | set(amr_map.keys()) | set(mge_map.keys())

    for bid in sorted(all_bins):
        amr_genes  = sorted(set(amr_map.get(bid, [])))
        mge_elems  = sorted(set(mge_map.get(bid, [])))
        qc         = checkm2_map.get(bid, {})

        row = {
            "sample":                sample,
            "bin_id":                bid,
            "abundance_trimmed_mean": abund_map.get(bid, ""),
            "n_amr_genes":           len(amr_genes),
            "amr_genes":             ",".join(amr_genes),
            "n_mge":                 len(mge_elems),
            "mge_elements":          ",".join(mge_elems),
        }
        if not skip_checkm2:
            row["completeness"]  = qc.get("completeness",  "")
            row["contamination"] = qc.get("contamination", "")
        if not skip_gtdbtk:
            row["gtdbtk_classification"] = gtdbtk_map.get(bid, "")
        all_rows.append(row)

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
BASE_COLS = [
    "sample", "bin_id", "abundance_trimmed_mean",
    "n_amr_genes", "amr_genes", "n_mge", "mge_elements",
]
opt_cols = []
if not skip_checkm2:
    opt_cols += ["completeness", "contamination"]
if not skip_gtdbtk:
    opt_cols += ["gtdbtk_classification"]

COLS = BASE_COLS + opt_cols

out_df = pd.DataFrame(all_rows, columns=COLS) if all_rows else pd.DataFrame(columns=COLS)
out_df = out_df.fillna("")
os.makedirs(os.path.dirname(out_file), exist_ok=True)
out_df.to_csv(out_file, sep="\t", index=False)
