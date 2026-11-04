#!/usr/bin/env python3
"""
generate_mock_data.py  –  SWAM mock dataset (v2)

Uses real, full-length reference genomes from a local database directory to
create two mock metagenomic samples that:

  • Produce >= 1 MAG per sample with MetaBAT2 (10x bacterial coverage)
  • Have overlapping composition (E. coli + Carjivirus in both samples) so
    CoverM multi-sample abundance mapping is exercised
  • Contain a phage (Carjivirus) for geNomad virus classification
  • Contain an RNA virus (Influenza A) in mock1 for geNomad
  • Contain bacteria with plasmid-carrying contigs for MobMess

Sample compositions
-------------------
  mock1: E. coli (10x) + E. faecalis (10x) + Carjivirus (3x) + Influenza A (3x)
  mock2: E. coli (10x) + Salmonella (10x)  + Campylobacter (10x) + Carjivirus (3x)

SCG / AMR database compatibility
---------------------------------
The existing test/dbs/scg/ and test/dbs/amrfinder/ databases embed 40 synthetic
CDSes (seed 52) and blaTEM-1 into the E. coli reference.  This script re-embeds
the EXACT same sequences into the first contig of the new E. coli assembly so
DIAMOND and KMA still find hits without rebuilding any database.

Run from the repository root:
    python test/scripts/generate_mock_data.py [--ref-dir /path/to/test-genomes]

Default ref-dir: /home/benda/Projects/databases/test-genomes

Prerequisites (available in .snakemake/conda/*/bin/ after snakemake --use-conda):
    wgsim
"""

import argparse
import os
import sys
import gzip
import glob
import random
import shutil
import textwrap
import tempfile
import subprocess

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ROOT      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEST_DATA = os.path.join(ROOT, "test", "data")
TEST_DBS  = os.path.join(ROOT, "test", "dbs")
REF_DIR_DEFAULT = "/home/benda/Projects/databases/test-genomes"

# ---------------------------------------------------------------------------
# Reference genome filenames (within ref_dir)
# ---------------------------------------------------------------------------

GENOME_FILES = {
    "ecoli":         "GCA_042156675.1_PDT002390328.1_genomic.fna.gz",   # E. coli 5.18 Mb
    "efaecalis":     "GCA_015304995.1_PDT000704954.1_genomic.fna.gz",   # E. faecalis 3.31 Mb
    "salmonella":    "GCA_052022995.1_PDT002880621.1_genomic.fna.gz",   # Salmonella 4.71 Mb
    "campylobacter": "GCF_013072625.1_ASM1307262v1_genomic.fna.gz",     # Campylobacter 1.94 Mb
    "carjivirus":    "GCF_009734245.1_ASM973424v1_genomic.fna.gz",      # Carjivirus 97 kb
    "influenza":     "GCA_038998445.1_ASM3899844v1_genomic.fna.gz",     # Influenza A 13.6 kb
}

# Coverage (fold) per genome component
COVERAGE = {
    "ecoli":         10,
    "efaecalis":     10,
    "salmonella":    10,
    "campylobacter": 10,
    "carjivirus":     3,
    "influenza":      3,
    # Small synthetic plasmid carrying blaTEM-1 at high coverage to ensure
    # reliable MEGAHIT assembly and contig-level AMRFinderPlus annotation.
    # At ~3 kb total and 100x, this adds only ~1,000 read pairs per sample.
    "arg_spike":    100,
}

# Sample compositions
SAMPLES = {
    "mock1": ["ecoli", "efaecalis", "carjivirus", "influenza", "arg_spike"],
    "mock2": ["ecoli", "salmonella", "campylobacter", "carjivirus", "arg_spike"],
}

# Per-sample base wgsim seeds (offset applied per genome within sample)
SAMPLE_SEEDS = {"mock1": 1000, "mock2": 2000}

READ_LEN    = 150
INSERT_SIZE = 400    # outer distance (fragment length)
INSERT_SD   = 40
ERROR_RATE  = 0.005

# ---------------------------------------------------------------------------
# Synthetic SCG CDSes (identical to those embedded by the original script)
# SEED+10 = 52 ensures the same 40 protein sequences as test/dbs/scg/SCGs_40.fasta
# ---------------------------------------------------------------------------

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
NON_STOP_CODONS = [c for c, aa in CODON_TABLE.items() if aa != "*"]

# blaTEM-1 (861 bp; same as original script)
BLA_TEM1 = (
    "ATGAGTATTCAACATTTTCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCT"
    "CACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAA"
    "CTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACT"
    "TTTAAAGTTCTGCTATGTGGTGCGGTATTATCCCGTGTTGACGCCGGGCAAGAGCAACTCGGTCGCCGC"
    "ATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATG"
    "ACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCTGCCAACTTACTTCTGACA"
    "ACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGAT"
    "CGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGCAGCAATG"
    "GCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGAC"
    "TGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCT"
    "GATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCC"
    "TCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCT"
    "GAGATAGGTGCCTCACTGATTAAGCATTGGTAA"
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd, check=True):
    print(f"  $ {cmd[:120]}")
    r = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, executable="/bin/bash"
    )
    if check and r.returncode != 0:
        print(f"  STDOUT: {r.stdout[-800:]}")
        print(f"  STDERR: {r.stderr[-800:]}")
        sys.exit(f"Command failed (rc={r.returncode}): {cmd[:120]}")
    return r


def find_tool(name):
    """Find a tool in snakemake conda envs or system PATH."""
    pattern = os.path.join(ROOT, ".snakemake", "conda", "*", "bin", name)
    hits = sorted(glob.glob(pattern))
    if hits:
        return hits[0]
    r = subprocess.run(f"which {name}", shell=True, capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


def parse_fasta(path):
    """Read gzipped or plain FASTA, return list of (header_line, seq)."""
    opener = gzip.open if path.endswith(".gz") else open
    records, cur_hdr, cur_seq = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_hdr is not None:
                    records.append((cur_hdr, "".join(cur_seq).upper()))
                cur_hdr, cur_seq = line[1:], []
            else:
                cur_seq.append(line)
    if cur_hdr is not None:
        records.append((cur_hdr, "".join(cur_seq).upper()))
    return records


def write_fasta(path, records):
    """Write list of (header, seq) to plain FASTA (uncompressed)."""
    with open(path, "w") as fh:
        for hdr, seq in records:
            fh.write(f">{hdr}\n")
            for chunk in textwrap.wrap(seq, 60):
                fh.write(chunk + "\n")


def make_scg_cdses(n=40, seed=52):
    """Generate 40 synthetic CDSes using the same seed as the original script."""
    rng = random.Random(seed)
    cdses = []
    for _ in range(n):
        body = [rng.choice(NON_STOP_CODONS) for _ in range(299)]
        cds = "ATG" + "".join(body)
        cdses.append(cds)
    return cdses


def embed(seq_list, insert, pos):
    """Overwrite seq_list[pos:pos+len(insert)] in-place."""
    for i, b in enumerate(insert):
        seq_list[pos + i] = b


# ---------------------------------------------------------------------------
# Build modified E. coli reference with embedded SCGs + blaTEM-1
# ---------------------------------------------------------------------------

def make_modified_ecoli(ref_dir, scg_cdses=None):
    """
    Load the E. coli reference, embed 40 synthetic SCG CDSes and blaTEM-1
    into the first contig, and return modified records.

    SCG positions: 5000, 15000, ..., 5000 + 39*10000 = 395000
    blaTEM-1 position: 200000  (between SCG19 at 195000 and SCG20 at 205000)
    All positions fit within the first contig (467182 bp).
    """
    gz_path = os.path.join(ref_dir, GENOME_FILES["ecoli"])
    records  = parse_fasta(gz_path)
    if scg_cdses is None:
        scg_cdses = make_scg_cdses()

    # Modify the first contig in-place
    hdr, seq = records[0]
    seg = list(seq)

    for i, cds in enumerate(scg_cdses):
        pos = 5000 + i * 10000
        embed(seg, cds, pos)

    embed(seg, BLA_TEM1, 200000)

    records[0] = (hdr, "".join(seg))
    print(f"  Embedded 40 SCG CDSes (positions 5000-395900) and blaTEM-1 "
          f"(position 200000) into {hdr[:60]}")
    return records


# ---------------------------------------------------------------------------
# Save reference FASTAs to test/dbs/genomes/
# ---------------------------------------------------------------------------

def save_reference_fastas(ref_dir, ecoli_records):
    genome_dir = os.path.join(TEST_DBS, "genomes")
    os.makedirs(genome_dir, exist_ok=True)

    # E. coli — save the modified version
    write_fasta(os.path.join(genome_dir, "ecoli.fa"), ecoli_records)
    print(f"  Saved ecoli.fa  ({len(ecoli_records)} contigs)")

    # All other genomes — save as-is
    others = {k: v for k, v in GENOME_FILES.items() if k != "ecoli"}
    for name, fname in others.items():
        gz_path = os.path.join(ref_dir, fname)
        records  = parse_fasta(gz_path)
        out_path = os.path.join(genome_dir, f"{name}.fa")
        write_fasta(out_path, records)
        total_bp = sum(len(s) for _, s in records)
        print(f"  Saved {name}.fa  ({len(records)} contigs, {total_bp:,} bp)")


# ---------------------------------------------------------------------------
# Read simulation
# ---------------------------------------------------------------------------

# Synthetic flanking sequence flanking blaTEM-1 in the ARG spike plasmid.
# Provides assembly context without relying on any external reference.
# 1000 bp each side (repeating ATCG pattern with variation to avoid tandem-
# repeat collapse in assembly).
_FLANK_SEED = 42
_FLANK_LEN  = 1000

def _make_flank(length: int, seed: int) -> str:
    rng = __import__("random").Random(seed)
    bases = "ATCG"
    return "".join(rng.choice(bases) for _ in range(length))

ARG_SPIKE_LEFT  = _make_flank(_FLANK_LEN, _FLANK_SEED)
ARG_SPIKE_RIGHT = _make_flank(_FLANK_LEN, _FLANK_SEED + 1)
ARG_SPIKE_SEQ   = ARG_SPIKE_LEFT + BLA_TEM1 + ARG_SPIKE_RIGHT


def make_arg_spike_fasta(path: str) -> int:
    """Write the ARG spike plasmid FASTA and return its length."""
    with open(path, "w") as fh:
        fh.write(">arg_spike_blaTEM1\n")
        for i in range(0, len(ARG_SPIKE_SEQ), 70):
            fh.write(ARG_SPIKE_SEQ[i:i + 70] + "\n")
    return len(ARG_SPIKE_SEQ)


def simulate_wgsim(wgsim, fa_path, n_pairs, seed, r1_out, r2_out):
    cmd = (
        f"{wgsim} -N {n_pairs} -1 {READ_LEN} -2 {READ_LEN} "
        f"-d {INSERT_SIZE} -s {INSERT_SD} "
        f"-S {seed} -e {ERROR_RATE} -r 0 -R 0 "
        f"{fa_path} {r1_out} {r2_out} > /dev/null"
    )
    run(cmd)


def make_sample(wgsim, sample, components, base_seed, ref_dir, ecoli_records, tmpdir):
    """
    Simulate reads for all genome components of a sample, combine, and gzip.
    Returns dict with per-component read counts and total.
    """
    all_r1, all_r2 = [], []
    stats = {}

    for offset, name in enumerate(components):
        # Get the reference FASTA for this component
        if name == "ecoli":
            fa_path = os.path.join(tmpdir, "ecoli_modified.fa")
            if not os.path.exists(fa_path):
                write_fasta(fa_path, ecoli_records)
            genome_len = sum(len(s) for _, s in ecoli_records)
        elif name == "arg_spike":
            fa_path = os.path.join(tmpdir, "arg_spike.fa")
            if not os.path.exists(fa_path):
                genome_len = make_arg_spike_fasta(fa_path)
            else:
                genome_len = len(ARG_SPIKE_SEQ)
        else:
            gz_path = os.path.join(ref_dir, GENOME_FILES[name])
            fa_path  = os.path.join(tmpdir, f"{name}.fa")
            if not os.path.exists(fa_path):
                records = parse_fasta(gz_path)
                write_fasta(fa_path, records)
                genome_len = sum(len(s) for _, s in records)
            else:
                # Re-use cached; compute length
                records = parse_fasta(fa_path)
                genome_len = sum(len(s) for _, s in records)

        cov      = COVERAGE[name]
        n_pairs  = max(200, int(genome_len * cov / (2 * READ_LEN)))
        seed     = base_seed + offset * 100

        r1_tmp = os.path.join(tmpdir, f"{sample}_{name}_R1.fastq")
        r2_tmp = os.path.join(tmpdir, f"{sample}_{name}_R2.fastq")
        simulate_wgsim(wgsim, fa_path, n_pairs, seed, r1_tmp, r2_tmp)

        all_r1.append(r1_tmp)
        all_r2.append(r2_tmp)
        stats[name] = {"coverage": cov, "genome_bp": genome_len, "pairs": n_pairs}
        print(f"    {sample}/{name}: {n_pairs:,} pairs  ({cov}x of {genome_len:,} bp)")

    # Combine and gzip
    r1_out = os.path.join(TEST_DATA, f"{sample}_R1.fastq.gz")
    r2_out = os.path.join(TEST_DATA, f"{sample}_R2.fastq.gz")

    with gzip.open(r1_out, "wt") as out:
        for f in all_r1:
            with open(f) as inp:
                shutil.copyfileobj(inp, out)
    with gzip.open(r2_out, "wt") as out:
        for f in all_r2:
            with open(f) as inp:
                shutil.copyfileobj(inp, out)

    total_pairs = sum(v["pairs"] for v in stats.values())
    print(f"  -> {r1_out}  ({total_pairs:,} total pairs)")
    stats["_total_pairs"] = total_pairs
    return stats


# ---------------------------------------------------------------------------
# Summary file
# ---------------------------------------------------------------------------

def write_summary(all_stats):
    out_path = os.path.join(TEST_DATA, "mock_data_summary.txt")
    lines = [
        "SWAM mock dataset v2 — generation summary",
        "=" * 60,
        "",
        "Reference genome source",
        "-----------------------",
        f"  {REF_DIR_DEFAULT}",
        "",
        "Genome sizes",
        "------------",
    ]
    genome_info = {
        "ecoli":         ("Escherichia coli 2024GO-0438",               "GCA_042156675.1", 5175654, 198, "chromosome + plasmid contigs (plasmid replication genes in GFF)"),
        "efaecalis":     ("Enterococcus faecalis ST6",                   "GCA_015304995.1", 3313510,  89, "chromosome + plasmid contigs (MobC mobilization proteins in GFF)"),
        "salmonella":    ("Salmonella enterica ENWS000537-005",          "GCA_052022995.1", 4705311,  43, "chromosome + plasmid contigs (RepB replication genes in GFF)"),
        "campylobacter": ("Campylobacter coli 18QD2YX31C",              "GCF_013072625.1", 1939122,   5, "chromosome + plasmid contigs (RepL replication genes in GFF)"),
        "carjivirus":    ("Carjivirus communis (phage)",                 "GCF_009734245.1",   97065,   1, "complete phage genome"),
        "influenza":     ("Influenza A H3N2 A/USA/PV00499/2017",        "GCA_038998445.1",   13628,   8, "8 segments (RNA virus)"),
    }
    for name, (species, accession, total_bp, ncontigs, notes) in genome_info.items():
        lines.append(f"  {name:<14}  {species}")
        lines.append(f"               Accession: {accession}")
        lines.append(f"               {total_bp:,} bp  |  {ncontigs} contigs  |  {notes}")
        lines.append("")

    lines += [
        "Sample compositions and read simulation",
        "----------------------------------------",
        "  Tool: wgsim",
        f"  Read length: {READ_LEN} bp (paired-end)",
        f"  Fragment size: {INSERT_SIZE} bp outer distance (SD {INSERT_SD} bp)",
        f"  Error rate: {ERROR_RATE}",
        "",
    ]
    for sample, components in SAMPLES.items():
        lines.append(f"  {sample}:")
        sample_stats = all_stats[sample]
        for name in components:
            s = sample_stats[name]
            lines.append(f"    {name:<14}  {s['coverage']}x  →  {s['pairs']:,} read pairs")
        lines.append(f"    {'TOTAL':<14}  {sample_stats['_total_pairs']:,} read pairs")
        lines.append("")

    lines += [
        "Overlap between samples",
        "-----------------------",
        "  E. coli and Carjivirus appear in both mock1 and mock2.",
        "  This tests CoverM genome-mode abundance mapping across samples.",
        "",
        "Expected workflow outcomes",
        "--------------------------",
        "  geNomad:        Carjivirus -> virus; E. faecalis/Salmonella/Campylobacter plasmid contigs -> plasmid",
        "  MetaBAT2 bins:  >= 1 MAG per sample (10x coverage of multi-Mb genomes)",
        "  CheckM2:        Completeness + contamination estimates per MAG",
        "  GTDB-tk:        Taxonomy for bacterial MAGs (requires large DB; skipped in test mode)",
        "  KMA/AMRFinder:  blaTEM-1 detected in short reads and contig AMR (embedded in E. coli first contig)",
        "  DIAMOND SCGs:   40 synthetic CDSes embedded in E. coli first contig -> n_genomes > 0",
        "",
        "SCG / AMR database embedding",
        "------------------------------",
        "  The 40 synthetic SCG CDSes (seed=52, 900 bp each, 300 aa translated) are embedded",
        "  in the E. coli first contig at positions 5000, 15000, ..., 395000.",
        "  blaTEM-1 (861 bp) is embedded at position 200000.",
        "  These are IDENTICAL to the sequences in test/dbs/scg/SCGs_40.fasta and",
        "  test/dbs/amrfinder/AMR_CDS.fa, so those databases do NOT need to be rebuilt.",
        "",
        "Files generated",
        "---------------",
        "  test/data/mock1_R1.fastq.gz",
        "  test/data/mock1_R2.fastq.gz",
        "  test/data/mock2_R1.fastq.gz",
        "  test/data/mock2_R2.fastq.gz",
        "  test/dbs/genomes/ecoli.fa           (modified: SCGs + blaTEM-1 embedded)",
        "  test/dbs/genomes/efaecalis.fa",
        "  test/dbs/genomes/salmonella.fa",
        "  test/dbs/genomes/campylobacter.fa",
        "  test/dbs/genomes/carjivirus.fa",
        "  test/dbs/genomes/influenza.fa",
        "  test/data/mock_data_summary.txt     (this file)",
        "  test/dbs/scg/SCGs_40.fasta          (40 synthetic 300-aa proteins, seed=52)",
        "  test/dbs/uniref50/uniref50_mini.fasta",
        "  test/dbs/uniref50/tax_mapping.tsv",
        "  test/dbs/uniref50/nodes.dmp",
        "  test/dbs/uniref50/names.dmp",
        "  test/dbs/uniref50/merged.dmp",
    ]
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  Written: {out_path}")


# ---------------------------------------------------------------------------
# Regenerate protein database source files (text FASTAs; binary DBs rebuilt by workflow)
# ---------------------------------------------------------------------------

def translate(dna):
    aa = []
    for i in range(0, len(dna) - 2, 3):
        a = CODON_TABLE.get(dna[i:i+3], "X")
        if a == "*":
            break
        aa.append(a)
    return "".join(aa)


def update_protein_databases(scg_cdses):
    """
    Write test/dbs/scg/SCGs_40.fasta and test/dbs/uniref50/ text files.
    These are the source files from which the workflow builds binary databases
    (DIAMOND SCG db via initiate_dbs; MMseqs2 UniRef50 db via init_mmseqs_db).
    """
    import textwrap as _tw

    scg_dir    = os.path.join(TEST_DBS, "scg")
    uniref_dir = os.path.join(TEST_DBS, "uniref50")
    os.makedirs(scg_dir, exist_ok=True)
    os.makedirs(uniref_dir, exist_ok=True)

    # Translate CDS -> protein
    proteins = [(f"SCG_{i:02d}", translate(cds)) for i, cds in enumerate(scg_cdses)]

    # SCGs_40.fasta
    scg_fasta = os.path.join(scg_dir, "SCGs_40.fasta")
    with open(scg_fasta, "w") as fh:
        for name, prot in proteins:
            fh.write(f">{name}\n")
            for chunk in _tw.wrap(prot, 60):
                fh.write(chunk + "\n")
    print(f"  Written {scg_fasta}  ({len(proteins)} proteins)")

    # UniRef50 mini FASTA: SCG proteins + blaTEM-1 protein
    bla_prot = translate(BLA_TEM1)
    acc_taxid = []
    uniref_records = []
    for i, (_, prot) in enumerate(proteins):
        acc = f"UniRef50_SCG{i:02d}"
        uniref_records.append((
            f"{acc} n=1 Tax=Escherichia coli TaxID=562 RepID=SCG{i:02d}_ECOLI",
            prot,
        ))
        acc_taxid.append((acc, 562))
    bla_acc = "UniRef50_BLATEM1"
    uniref_records.append((
        f"{bla_acc} n=1 Tax=Escherichia coli TaxID=562 RepID=BLATEM1_ECOLI",
        bla_prot,
    ))
    acc_taxid.append((bla_acc, 562))

    uniref_fasta = os.path.join(uniref_dir, "uniref50_mini.fasta")
    with open(uniref_fasta, "w") as fh:
        for hdr, seq in uniref_records:
            fh.write(f">{hdr}\n")
            for chunk in _tw.wrap(seq, 60):
                fh.write(chunk + "\n")

    # tax_mapping.tsv
    with open(os.path.join(uniref_dir, "tax_mapping.tsv"), "w") as fh:
        for acc, taxid in acc_taxid:
            fh.write(f"{acc}\t{taxid}\n")

    # Minimal NCBI taxonomy dump covering E. coli lineage
    fake_taxa = {
        1:      ("root",               1,      "no rank"),
        131567: ("cellular organisms", 1,      "no rank"),
        2:      ("Bacteria",           131567, "superkingdom"),
        1224:   ("Proteobacteria",     2,      "phylum"),
        1236:   ("Gammaproteobacteria",1224,   "class"),
        91347:  ("Enterobacterales",   1236,   "order"),
        543:    ("Enterobacteriaceae", 91347,  "family"),
        561:    ("Escherichia",        543,    "genus"),
        562:    ("Escherichia coli",   561,    "species"),
    }
    with open(os.path.join(uniref_dir, "nodes.dmp"), "w") as fh:
        for taxid, (_, parent, rank) in fake_taxa.items():
            fh.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n")
    with open(os.path.join(uniref_dir, "names.dmp"), "w") as fh:
        for taxid, (name, _, _) in fake_taxa.items():
            fh.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    open(os.path.join(uniref_dir, "merged.dmp"), "w").close()

    print(f"  Written {uniref_fasta}  ({len(uniref_records)} sequences)")
    print(f"  Written UniRef50 taxonomy files (nodes.dmp, names.dmp, merged.dmp, tax_mapping.tsv)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate SWAM mock test data v2")
    parser.add_argument("--ref-dir", default=REF_DIR_DEFAULT,
                        help=f"Directory containing reference genome .fna.gz files (default: {REF_DIR_DEFAULT})")
    args = parser.parse_args()
    ref_dir = args.ref_dir

    # Verify reference files exist
    missing = [f for f in GENOME_FILES.values()
               if not os.path.exists(os.path.join(ref_dir, f))]
    if missing:
        sys.exit(f"ERROR: missing reference files in {ref_dir}:\n  " + "\n  ".join(missing))

    wgsim = find_tool("wgsim")
    if not wgsim:
        sys.exit("ERROR: wgsim not found.  Run: snakemake --use-conda first.")

    print("=" * 60)
    print("SWAM: generating mock test data v2")
    print(f"  ROOT:    {ROOT}")
    print(f"  ref_dir: {ref_dir}")
    print(f"  wgsim:   {wgsim}")
    print("=" * 60)

    os.makedirs(TEST_DATA, exist_ok=True)

    print("\n[1] Building modified E. coli reference (embedding SCGs + blaTEM-1)...")
    scg_cdses     = make_scg_cdses()
    ecoli_records = make_modified_ecoli(ref_dir, scg_cdses=scg_cdses)

    print("\n[2] Updating protein database source FASTAs...")
    update_protein_databases(scg_cdses)

    print("\n[3] Saving reference FASTAs to test/dbs/genomes/...")
    save_reference_fastas(ref_dir, ecoli_records)

    print("\n[4] Simulating reads...")
    all_stats = {}
    with tempfile.TemporaryDirectory() as tmpdir:
        for sample, components in SAMPLES.items():
            print(f"\n  --- {sample} ---")
            stats = make_sample(
                wgsim, sample, components,
                SAMPLE_SEEDS[sample], ref_dir, ecoli_records, tmpdir
            )
            all_stats[sample] = stats

    print("\n[5] Writing summary...")
    write_summary(all_stats)

    print("\n" + "=" * 60)
    print("Done. Mock test data ready.")
    print("  FASTQs:      test/data/")
    print("  Ref FASTAs:  test/dbs/genomes/")
    print("  Summary:     test/data/mock_data_summary.txt")
    print()
    print("NOTE: After running this script, also run:")
    print("  python test/scripts/inject_marker_reads.py")
    print("to append pBI143 and crAss001 reads to the test FASTQs.")
    print("This is required for the AMR risk scoring module (Exposure component).")
    print("=" * 60)


if __name__ == "__main__":
    main()
