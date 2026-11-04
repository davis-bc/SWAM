"""
contig_summary.py

Snakemake script (called via script: directive).

Aggregates per-sample contig-level annotations into one long-format TSV:

    sample | contig_id | trimmed_mean |
    molecule_type | taxonomy | feature_type | gene | start | stop | strand

One row per annotated feature (AMR gene or MGE). Contig abundance comes from
CoverM contig mode (trimmed_mean only).
"""

import os
import pandas as pd

# ---------------------------------------------------------------------------
# Snakemake-injected variables
# ---------------------------------------------------------------------------
samples    = snakemake.params.samples
output_dir = snakemake.params.output_dir
out_file   = snakemake.output.tsv

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def safe_read(path, **kwargs):
    """Read TSV, return empty DataFrame on missing / empty file."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", **kwargs)
    except Exception:
        return pd.DataFrame()


def parse_genomad_molecule_type(plas_sum_path, vir_sum_path):
    """Return dict {contig_id -> molecule_type} from geNomad summary TSVs."""
    plas = safe_read(plas_sum_path)
    vir  = safe_read(vir_sum_path)
    mol_map: dict[str, str] = {}
    if not plas.empty and "seq_name" in plas.columns:
        for cid in plas["seq_name"]:
            mol_map[str(cid)] = "plasmid"
    if not vir.empty and "seq_name" in vir.columns:
        for cid in vir["seq_name"]:
            mol_map.setdefault(str(cid), "phage")
    return mol_map  # missing == chromosome


def parse_lca(lca_path):
    """Return dict {contig_id -> taxonomy_string} from MMseqs2 easy-taxonomy LCA TSV."""
    df = safe_read(lca_path, header=None)
    if df.empty or len(df.columns) < 4:
        return {}
    # columns: query, taxid, rank, taxname, n_frags, n_direct, n_classified, fraction, lineage
    df.columns = list(df.columns)  # keep int index
    return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 3].astype(str)))


def parse_abundance(abund_path):
    """Return dict {contig_id -> trimmed_mean} from CoverM contig output."""
    df = safe_read(abund_path)
    if df.empty or "contig_id" not in df.columns:
        return {}
    result = {}
    for _, row in df.iterrows():
        cid = str(row["contig_id"])
        result[cid] = row.get("trimmed_mean", "")
    return result


def parse_amr(amr_path):
    """
    Parse AMRFinderPlus output (nucleotide+protein+GFF mode).

    Relevant columns when run with -n/-p/-g --plus:
      'Contig id', 'Gene symbol', 'Type', 'Start', 'Stop', 'Strand'
    The 'Type' column carries AMR / STRESS / VIRULENCE / ANTIFUNGAL.
    Column names vary slightly across versions; fall back gracefully.
    """
    df = safe_read(amr_path)
    if df.empty:
        return []

    cols = {c.lower().replace(" ", "_"): c for c in df.columns}

    def gcol(*candidates):
        for c in candidates:
            if c in cols:
                return cols[c]
        return None

    contig_col = gcol("contig_id", "contig", "sequence_name", "sequence")
    gene_col   = gcol("gene_symbol", "gene", "element_symbol", "name")
    type_col   = gcol("type")
    start_col  = gcol("start")
    stop_col   = gcol("stop")
    strand_col = gcol("strand")

    features = []
    for _, row in df.iterrows():
        # Contig ID: AMRFinderPlus nucleotide mode puts the contig name directly
        # in 'Contig id'. In protein-only mode it's in 'Protein identifier'
        # and we strip the trailing _<gene_num>.
        cid = ""
        if contig_col:
            cid = str(row[contig_col])
        else:
            prot_col = gcol("protein_identifier")
            if prot_col:
                prot_id = str(row[prot_col])
                cid = prot_id.rsplit("_", 1)[0]

        feature_type = str(row[type_col]).strip() if type_col else "AMR"

        features.append({
            "contig_id":    cid,
            "feature_type": feature_type,
            "gene":         str(row[gene_col])   if gene_col   else "",
            "start":        str(row[start_col])  if start_col  else "",
            "stop":         str(row[stop_col])   if stop_col   else "",
            "strand":       str(row[strand_col]) if strand_col else "",
        })
    return features


def parse_mge(mge_path):
    """
    Parse MobileElementFinder output.

    Observed format (mge_finder ≥ 1.1):
      CSV file with leading comment lines starting with '#'.
      Columns: mge_no, name, synonyms, prediction, type, allele_len, depth,
               e_value, identity, coverage, gaps, substitution, contig, start, end, cigar

    The 'contig' column often contains extra metadata after the first space
    (e.g. "mock1-k127_19 flag=1 multi=7.0000 len=5958"); only the first
    space-delimited token is the actual contig identifier.
    """
    if not os.path.exists(mge_path) or os.path.getsize(mge_path) == 0:
        return []
    try:
        df = pd.read_csv(mge_path, sep=",", comment="#", dtype=str)
    except Exception:
        return []
    if df.empty:
        return []

    cols_lower = {c.lower(): c for c in df.columns}

    def gcol(*candidates):
        for c in candidates:
            if c in cols_lower:
                return cols_lower[c]
        return None

    seq_col    = gcol("contig", "sequence", "seqid", "#sequence")
    name_col   = gcol("name", "element_symbol", "type", "element")
    start_col  = gcol("start", "pos_from", "pos_start")
    stop_col   = gcol("end", "pos_to", "stop", "pos_end")
    strand_col = gcol("strand")

    features = []
    for _, row in df.iterrows():
        # Strip extra content appended after first space in contig column
        raw_cid = str(row[seq_col]) if seq_col else ""
        cid = raw_cid.split()[0] if raw_cid else ""
        features.append({
            "contig_id":    cid,
            "feature_type": "MGE",
            "gene":         str(row[name_col])   if name_col   else "",
            "start":        str(row[start_col])  if start_col  else "",
            "stop":         str(row[stop_col])   if stop_col   else "",
            "strand":       str(row[strand_col]) if strand_col else "",
        })
    return features


# ---------------------------------------------------------------------------
# Main aggregation loop
# ---------------------------------------------------------------------------

all_rows = []

for sample in samples:
    plas_sum = os.path.join(output_dir, "data", "genomad", sample,
                            f"{sample}.contigs_summary",
                            f"{sample}.contigs_plasmid_summary.tsv")
    vir_sum  = os.path.join(output_dir, "data", "genomad", sample,
                            f"{sample}.contigs_summary",
                            f"{sample}.contigs_virus_summary.tsv")
    lca_path   = os.path.join(output_dir, "data", "mmseqs",           f"{sample}_lca.tsv")
    abund_path = os.path.join(output_dir, "data", "contig_abundance",  f"{sample}_contig_abundance.tsv")
    amr_path   = os.path.join(output_dir, "data", "amr_contigs",       f"{sample}_contig_amr.tsv")
    mge_path   = os.path.join(output_dir, "data", "mge_contigs",       f"{sample}_mge.tsv")

    mol_map   = parse_genomad_molecule_type(plas_sum, vir_sum)
    tax_map   = parse_lca(lca_path)
    abund_map = parse_abundance(abund_path)
    features  = parse_amr(amr_path) + parse_mge(mge_path)

    # Collect all known contig IDs from abundance + features
    all_contigs = set(abund_map.keys()) | {f["contig_id"] for f in features}

    if features:
        for feat in features:
            cid = feat["contig_id"]
            all_contigs.discard(cid)  # will be covered below
            abund = abund_map.get(cid, "")
            all_rows.append({
                "sample":        sample,
                "contig_id":     cid,
                "trimmed_mean":  abund,
                "molecule_type": mol_map.get(cid, "chromosome"),
                "taxonomy":      tax_map.get(cid, ""),
                "feature_type":  feat["feature_type"],
                "gene":          feat["gene"],
                "start":         feat["start"],
                "stop":          feat["stop"],
                "strand":        feat["strand"],
            })

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
COLS = ["sample", "contig_id", "trimmed_mean",
        "molecule_type", "taxonomy", "feature_type", "gene", "start", "stop", "strand"]

out_df = pd.DataFrame(all_rows, columns=COLS) if all_rows else pd.DataFrame(columns=COLS)
out_df.to_csv(out_file, sep="\t", index=False)
