"""
amr_risk_score.py

Snakemake script (called via script: directive).

Computes per-sample AMR risk scores using two complementary strategies from
unified_amr_risk_score_strategies.txt:

  Strategy 1 — Weighted Additive Composite Score
    FeatureScore_i = AbundanceWeight_i * (0.30*R + 0.30*M + 0.25*H + 0.15*E)
    SampleRisk_Additive = sum(FeatureScore_i), normalized to 0-100 study-wide

  Strategy 2 — Hazard-Pathway Multiplicative Score
    FeatureScore_i = AbundanceWeight_i * (R^0.30) * (M^0.30) * (H^0.20) * (E^0.20)
    SampleRisk_Multiplicative = sum(FeatureScore_i), normalized to 0-100 study-wide

Four component scores (0-1 each):
  R — Resistance Hazard    (from AMRFinderPlus subclass/class, WHO criticality tiers)
  M — Mobility             (from geNomad molecule_type + MGE co-location on same contig)
  H — Host/Pathogenicity   (from MMseqs2 LCA taxonomy string)
  E — Exposure/Colonization (from pBI143 / crAss001 cpg, sample-level)

AbundanceWeight = log10(1 + cpg), min-max normalized 0-1 study-wide.

Produces:
  AMR_abundance_summary.csv  — per-sample table with abundance + risk columns:
    sample, AMR_total_cpg, pBI143_cpg, crAss001_cpg,
    E_exposure, R_mean, M_mean, H_mean,
    amr_risk_additive_raw, amr_risk_multiplicative_raw,
    amr_risk_additive, amr_risk_multiplicative
"""

import math
import os

import pandas as pd

# ---------------------------------------------------------------------------
# Snakemake-injected variables
# ---------------------------------------------------------------------------

short_reads_csv  = snakemake.input.get("short_reads",   None)
_ctg             = snakemake.input.get("contig_summary", [])
contig_tsv       = str(_ctg[0]) if _ctg else None
markers_csv      = snakemake.input.get("markers_cpg",   None)
out_file         = snakemake.output[0]
samples          = snakemake.params.samples

# ---------------------------------------------------------------------------
# Component R — Resistance Hazard
# WHO Critical / Highly Important / Important / Lower-priority tiers
# ---------------------------------------------------------------------------

# Subclass takes priority over class when available.
# Values sourced from WHO CIA 2019 + AMRFinderPlus catalog field naming.

SUBCLASS_R: dict[str, float] = {
    # ---- WHO Critical: last-resort / no alternatives ----
    "CARBAPENEM":       1.00,
    "COLISTIN":         1.00,
    "OXAZOLIDINONE":    1.00,   # linezolid
    "GLYCOPEPTIDE":     1.00,   # vancomycin, teicoplanin
    "DAPTOMYCIN":       1.00,
    "FOSFOMYCIN":       1.00,
    # ---- WHO High: highly important, limited alternatives ----
    "FLUOROQUINOLONE":  0.80,
    "AMINOGLYCOSIDE":   0.80,
    "CEPHALOSPORIN":    0.80,   # all generations lumped in AFP; subtype would refine
    "MACROLIDE":        0.80,
    "TETRACYCLINE":     0.80,
    "PENICILLIN":       0.80,
    "METHICILLIN":      0.80,
    "STREPTOMYCIN":     0.80,
    "AZITHROMYCIN":     0.80,
    # ---- WHO Important: still significant in human medicine ----
    "SULFONAMIDE":      0.60,
    "TRIMETHOPRIM":     0.60,
    "CHLORAMPHENICOL":  0.60,
    "LINCOSAMIDE":      0.60,
    "RIFAMYCIN":        0.60,
    "FUSIDIC ACID":     0.60,
    "MUPIROCIN":        0.60,
}

CLASS_R: dict[str, float] = {
    "FLUOROQUINOLONE":  1.00,   # included in WHO Critical for enteric pathogens
    "GLYCOPEPTIDE":     1.00,
    "OXAZOLIDINONE":    1.00,
    "AMINOGLYCOSIDE":   0.80,
    "BETA-LACTAM":      0.80,
    "TETRACYCLINE":     0.80,
    "MACROLIDE":        0.80,
    "LINCOSAMIDE":      0.80,
    "SULFONAMIDE":      0.60,
    "CHLORAMPHENICOL":  0.60,
    "TRIMETHOPRIM":     0.60,
    "RIFAMYCIN":        0.60,
    "MDR":              0.80,   # multi-drug efflux pumps
    "MULTIDRUG":        0.80,
}

R_DEFAULT = 0.40  # unclassified or lower-priority


def resistance_score(subclass: str, cls: str) -> float:
    """Return R ∈ [0.40, 1.00] based on subclass then class."""
    for val in (subclass, cls):
        if val:
            v = str(val).strip().upper()
            if v in SUBCLASS_R:
                return SUBCLASS_R[v]
            if v in CLASS_R:
                return CLASS_R[v]
    return R_DEFAULT


# ---------------------------------------------------------------------------
# Component M — Mobility
# Derived from contig_summary.tsv: molecule_type column + MGE co-location
# ---------------------------------------------------------------------------

M_PLASMID_MGE    = 1.00
M_PLASMID        = 0.85
M_PHAGE          = 0.75
M_CHROMOSOME_MGE = 0.70
M_CHROMOSOME     = 0.30
M_UNKNOWN        = 0.10  # short_reads_only: no contig context


def mobility_score(molecule_type: str, has_nearby_mge: bool) -> float:
    mt = str(molecule_type).strip().lower()
    if mt == "plasmid":
        return M_PLASMID_MGE if has_nearby_mge else M_PLASMID
    if mt == "phage":
        return M_PHAGE
    if mt == "chromosome":
        return M_CHROMOSOME_MGE if has_nearby_mge else M_CHROMOSOME
    return M_UNKNOWN


# ---------------------------------------------------------------------------
# Component H — Host / Pathogenicity
# Derived from taxonomy string in contig_summary.tsv
# ---------------------------------------------------------------------------

# Keyword groups (case-insensitive substring match)
ESKAPE_KW = {
    "enterococcus faecium", "staphylococcus aureus",
    "klebsiella", "acinetobacter", "pseudomonas aeruginosa",
    "enterobacter",
}
ENTEROBACTERIACEAE_KW = {
    "escherichia", "salmonella", "shigella", "citrobacter",
    "proteus", "morganella", "serratia", "cronobacter",
}
HUMAN_ASSOC_KW = {
    "campylobacter", "clostridioides", "clostridium difficile",
    "streptococcus", "staphylococcus", "enterococcus",
    "haemophilus", "neisseria", "helicobacter",
    "bacteroides", "listeria", "vibrio", "yersinia",
}

H_ESKAPE         = 1.00
H_ENTEROBACTERIA = 0.80
H_HUMAN_ASSOC    = 0.60
H_ENVIRONMENTAL  = 0.35
H_UNKNOWN        = 0.15  # empty / unclassified / no rank


def host_score(taxonomy: str) -> float:
    """Return H ∈ [0.15, 1.00] based on taxonomy string."""
    if not taxonomy or str(taxonomy).strip() in ("", "nan", "unknown", "unclassified", "no rank"):
        return H_UNKNOWN
    tax_lower = str(taxonomy).strip().lower()
    for kw in ESKAPE_KW:
        if kw in tax_lower:
            return H_ESKAPE
    for kw in ENTEROBACTERIACEAE_KW:
        if kw in tax_lower:
            return H_ENTEROBACTERIA
    for kw in HUMAN_ASSOC_KW:
        if kw in tax_lower:
            return H_HUMAN_ASSOC
    # Non-empty but not in known human-relevant groups
    return H_ENVIRONMENTAL


# ---------------------------------------------------------------------------
# Component E — Exposure / Colonization (sample-level)
# ---------------------------------------------------------------------------

E_STRONG   = 1.00   # max marker cpg ≥ 1.0
E_MODERATE = 0.75   # max marker cpg ≥ 0.1
E_WEAK     = 0.40   # max marker cpg ≥ 0.01
E_MINIMAL  = 0.10   # absent / below detection


def exposure_score(pbi143_cpg: float, crass001_cpg: float) -> float:
    max_cpg = max(float(pbi143_cpg or 0), float(crass001_cpg or 0))
    if max_cpg >= 1.0:
        return E_STRONG
    if max_cpg >= 0.1:
        return E_MODERATE
    if max_cpg >= 0.01:
        return E_WEAK
    return E_MINIMAL


# ---------------------------------------------------------------------------
# Safe readers
# ---------------------------------------------------------------------------

def safe_read_csv(path, sep=",", **kwargs):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep=sep, dtype=str, **kwargs)
    except Exception:
        return pd.DataFrame()


# ---------------------------------------------------------------------------
# Build per-gene contig context lookup from contig_summary.tsv
#
# For each (sample, gene) pair, we want:
#   - best molecule_type (most mobile)
#   - best taxonomy (most specific non-empty value)
#   - whether that contig also carries an MGE feature
#
# Algorithm:
#   1. For each AMR row in contig_summary, record contig_id + molecule_type + taxonomy.
#   2. For each MGE row in contig_summary, record the set of contig_ids that carry an MGE.
#   3. For each (sample, gene), join and compute M and H.  Take max M across all contigs.
# ---------------------------------------------------------------------------

def build_contig_context(contig_tsv_path: str):
    """
    Returns dict: { (sample, gene) -> {"M": float, "H": float} }
    Falls back to (M_UNKNOWN, H_UNKNOWN) if gene not found in contigs.
    """
    df = safe_read_csv(contig_tsv_path, sep="\t")
    if df.empty or "feature_type" not in df.columns:
        return {}

    # Normalize column names
    df["sample"]   = df.get("sample", "").fillna("")
    df["gene"]     = df.get("gene",   "").fillna("")
    df["contig_id"] = df.get("contig_id", "").fillna("")
    df["molecule_type"] = df.get("molecule_type", "chromosome").fillna("chromosome")
    df["taxonomy"]  = df.get("taxonomy", "").fillna("")
    df["feature_type"] = df.get("feature_type", "").fillna("")

    # Contigs that carry at least one MGE
    mge_contigs: set[tuple[str, str]] = set()
    mge_df = df[df["feature_type"].str.upper() == "MGE"]
    for _, row in mge_df.iterrows():
        mge_contigs.add((str(row["sample"]), str(row["contig_id"])))

    # AMR rows
    amr_df = df[df["feature_type"].str.upper() == "AMR"]

    context: dict[tuple[str, str], dict[str, float]] = {}

    for (sample, gene), grp in amr_df.groupby(["sample", "gene"]):
        best_M = M_UNKNOWN
        best_H = H_UNKNOWN

        for _, row in grp.iterrows():
            cid = str(row["contig_id"])
            has_mge = (str(row["sample"]), cid) in mge_contigs
            m = mobility_score(row["molecule_type"], has_mge)
            h = host_score(row["taxonomy"])
            best_M = max(best_M, m)
            # Prefer the taxonomy that gives a higher (more specific) H score
            if h > best_H:
                best_H = h

        context[(str(sample), str(gene))] = {"M": best_M, "H": best_H}

    return context


# ---------------------------------------------------------------------------
# Main scoring
# ---------------------------------------------------------------------------

def main():
    # ---- load inputs ----
    sr_df      = safe_read_csv(short_reads_csv)
    markers_df = safe_read_csv(markers_csv)
    context    = build_contig_context(contig_tsv) if contig_tsv else {}

    # ---- build per-sample marker cpg ----
    marker_map: dict[str, tuple[float, float]] = {}  # sample -> (pBI143_cpg, crAss001_cpg)
    if not markers_df.empty and "sample" in markers_df.columns:
        for _, row in markers_df.iterrows():
            s = str(row["sample"])
            pbi = float(row.get("pBI143_cpg",  0) or 0)
            cra = float(row.get("crAss001_cpg", 0) or 0)
            marker_map[s] = (pbi, cra)

    # ---- pre-compute exposure per sample ----
    E_map: dict[str, float] = {}
    for s in samples:
        pbi, cra = marker_map.get(s, (0.0, 0.0))
        E_map[s] = exposure_score(pbi, cra)

    # ---- build feature table (one row per gene × sample) ----
    feature_rows = []

    if not sr_df.empty and "sample" in sr_df.columns:
        for _, row in sr_df.iterrows():
            s    = str(row.get("sample", ""))
            gene = str(row.get("allele",  row.get("gene_family", "")))
            cpg  = float(row.get("cpg", 0) or 0)
            subclass = str(row.get("subclass", "") or "")
            cls      = str(row.get("class",    "") or "")
            ftype    = str(row.get("type",     "") or "").upper()

            # Only score AMR-type genes (skip STRESS, VIRULENCE, ANTIFUNGAL)
            if ftype and ftype not in ("AMR", ""):
                continue

            R = resistance_score(subclass, cls)
            ctx = context.get((s, gene), {})
            M = ctx.get("M", M_UNKNOWN)
            H = ctx.get("H", H_UNKNOWN)
            E = E_map.get(s, E_MINIMAL)

            feature_rows.append({
                "sample": s,
                "gene":   gene,
                "cpg":    cpg,
                "R": R, "M": M, "H": H, "E": E,
            })

    feat_df = pd.DataFrame(feature_rows) if feature_rows else pd.DataFrame(
        columns=["sample", "gene", "cpg", "R", "M", "H", "E"]
    )

    # ---- Abundance weight: log10(1 + cpg), normalized 0-1 study-wide ----
    if not feat_df.empty:
        feat_df["cpg"] = pd.to_numeric(feat_df["cpg"], errors="coerce").fillna(0)
        feat_df["log_cpg"] = feat_df["cpg"].apply(lambda x: math.log10(1 + max(x, 0)))
        log_min = feat_df["log_cpg"].min()
        log_max = feat_df["log_cpg"].max()
        if log_max > log_min:
            feat_df["AbundanceWeight"] = (feat_df["log_cpg"] - log_min) / (log_max - log_min)
        else:
            feat_df["AbundanceWeight"] = 0.5  # all equal
    else:
        feat_df["log_cpg"] = 0.0
        feat_df["AbundanceWeight"] = 0.0

    # ---- Per-feature scores ----
    if not feat_df.empty:
        R = feat_df["R"]; M = feat_df["M"]
        H = feat_df["H"]; E = feat_df["E"]
        W = feat_df["AbundanceWeight"]

        # Strategy 1 — additive
        feat_df["fs_additive"] = W * (0.30 * R + 0.30 * M + 0.25 * H + 0.15 * E)

        # Strategy 2 — multiplicative (guard against zero base: use 1e-6 floor)
        R_s = R.clip(lower=1e-6); M_s = M.clip(lower=1e-6)
        H_s = H.clip(lower=1e-6); E_s = E.clip(lower=1e-6)
        feat_df["fs_multi"] = W * (R_s ** 0.30) * (M_s ** 0.30) * (H_s ** 0.20) * (E_s ** 0.20)

    # ---- Aggregate to per-sample ----
    sample_rows = []
    for s in samples:
        pbi, cra = marker_map.get(s, (0.0, 0.0))
        E_val = E_map.get(s, E_MINIMAL)

        if not feat_df.empty:
            sg = feat_df[feat_df["sample"] == s]
        else:
            sg = pd.DataFrame()

        amr_total = float(sg["cpg"].sum()) if not sg.empty else 0.0
        R_mean = float(sg["R"].mean()) if not sg.empty else 0.0
        M_mean = float(sg["M"].mean()) if not sg.empty else 0.0
        H_mean = float(sg["H"].mean()) if not sg.empty else 0.0

        raw_add  = float(sg["fs_additive"].sum()) if not sg.empty else 0.0
        raw_mult = float(sg["fs_multi"].sum())     if not sg.empty else 0.0

        sample_rows.append({
            "sample":               s,
            "AMR_total_cpg":        round(amr_total, 4),
            "pBI143_cpg":           round(pbi, 4),
            "crAss001_cpg":         round(cra, 4),
            "E_exposure":           round(E_val,  2),
            "R_mean":               round(R_mean, 3),
            "M_mean":               round(M_mean, 3),
            "H_mean":               round(H_mean, 3),
            "amr_risk_additive_raw":      round(raw_add,  6),
            "amr_risk_multiplicative_raw": round(raw_mult, 6),
        })

    out_df = pd.DataFrame(sample_rows)

    # ---- Study-level 0-100 normalization ----
    for raw_col, norm_col in [
        ("amr_risk_additive_raw",       "amr_risk_additive"),
        ("amr_risk_multiplicative_raw", "amr_risk_multiplicative"),
    ]:
        vals = out_df[raw_col]
        lo, hi = vals.min(), vals.max()
        if hi > lo:
            out_df[norm_col] = ((vals - lo) / (hi - lo) * 100).round(2)
        else:
            # All samples identical → assign 50 (relative middle)
            out_df[norm_col] = 50.0

    # ---- Write output ----
    os.makedirs(os.path.dirname(os.path.abspath(out_file)), exist_ok=True)
    out_df.to_csv(out_file, index=False)


main()
