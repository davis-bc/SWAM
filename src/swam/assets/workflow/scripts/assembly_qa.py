"""
assembly_qa.py

Snakemake script (called via script: directive).

Aggregates per-sample assembly quality metrics into one wide-format TSV
(one row per sample), analogous to fastp_summary.csv for the read-QC stage.

Columns:
    sample | n_contigs | total_length_bp | n50_bp | longest_contig_bp |
    shortest_contig_bp | mean_contig_length_bp |
    reads_total | reads_mapped | pct_reads_mapped

Sources:
  - MEGAHIT log: "N contigs, total X bp, min Y bp, max Z bp, avg A bp, N50 B bp"
  - samtools flagstat on the contig BAM from contig_abundance
"""

import re
import subprocess
import pandas as pd

samples    = snakemake.params.samples
output_dir = snakemake.params.output_dir
out_file   = snakemake.output.tsv

bam_map = {
    s: p for s, p in zip(samples, snakemake.input.bam_files)
}
log_map = {
    s: p for s, p in zip(samples, snakemake.input.megahit_logs)
}

# Regex for the MEGAHIT summary line:
#   "2026-03-24 17:40:20 - 331 contigs, total 8278319 bp, min 1014 bp,
#    max 313433 bp, avg 25010 bp, N50 67653 bp"
_MEGAHIT_RE = re.compile(
    r"(\d+)\s+contigs,\s+total\s+(\d+)\s+bp,\s+min\s+(\d+)\s+bp,\s+"
    r"max\s+(\d+)\s+bp,\s+avg\s+(\d+)\s+bp,\s+N50\s+(\d+)\s+bp"
)

# Regex for samtools flagstat lines we care about:
#   "569956 + 0 in total ..."
#   "565871 + 0 primary mapped (99.32% : N/A)"
_TOTAL_RE   = re.compile(r"^(\d+) \+ \d+ in total")
_MAPPED_RE  = re.compile(r"^(\d+) \+ \d+ primary mapped \(([0-9.]+)%")


def parse_megahit_log(log_path):
    """Return dict of assembly stats from the MEGAHIT summary line."""
    try:
        with open(log_path) as fh:
            for line in fh:
                m = _MEGAHIT_RE.search(line)
                if m:
                    return {
                        "n_contigs":            int(m.group(1)),
                        "total_length_bp":      int(m.group(2)),
                        "shortest_contig_bp":   int(m.group(3)),
                        "longest_contig_bp":    int(m.group(4)),
                        "mean_contig_length_bp": int(m.group(5)),
                        "n50_bp":               int(m.group(6)),
                    }
    except Exception:
        pass
    return {k: "" for k in [
        "n_contigs", "total_length_bp", "shortest_contig_bp",
        "longest_contig_bp", "mean_contig_length_bp", "n50_bp",
    ]}


def parse_flagstat(bam_path):
    """Return dict of read-mapping stats from samtools flagstat."""
    result = {"reads_total": "", "reads_mapped": "", "pct_reads_mapped": ""}
    try:
        proc = subprocess.run(
            ["samtools", "flagstat", bam_path],
            capture_output=True, text=True, check=True,
        )
        for line in proc.stdout.splitlines():
            m_total = _TOTAL_RE.match(line)
            if m_total:
                result["reads_total"] = int(m_total.group(1))
            m_mapped = _MAPPED_RE.match(line)
            if m_mapped:
                result["reads_mapped"]     = int(m_mapped.group(1))
                result["pct_reads_mapped"] = float(m_mapped.group(2))
    except Exception:
        pass
    return result


rows = []
for sample in samples:
    row = {"sample": sample}
    row.update(parse_megahit_log(log_map[sample]))
    row.update(parse_flagstat(bam_map[sample]))
    rows.append(row)

COLS = [
    "sample",
    "n_contigs", "total_length_bp", "n50_bp",
    "longest_contig_bp", "shortest_contig_bp", "mean_contig_length_bp",
    "reads_total", "reads_mapped", "pct_reads_mapped",
]
out_df = pd.DataFrame(rows, columns=COLS)
out_df.to_csv(out_file, sep="\t", index=False)
