#!/usr/bin/env python3
"""
inject_marker_reads.py  —  SWAM test data helper

Generates synthetic paired-end reads from the anthropogenic marker sequences
(pBI143 and crAss001) and appends them to the mock test FASTQs.

Without this step, KMA finds no reads matching the markers in the mock data
(which is composed of bacteria/phage that are not pBI143 or crAss001), so
markers_cpg.csv shows 0 for both markers and the Exposure component (E) of
the AMR risk scoring is always at its minimum.

Target coverages (chosen to give realistic wastewater-like cpg values):
  pBI143  (2.7 kb) :  200 read pairs  →  ~10x coverage
  crAss001 (103 kb):  500 read pairs  →  ~1x coverage

Run from the repository root (no dependencies beyond Python stdlib):
    python test/scripts/inject_marker_reads.py [--dry-run]
"""

import argparse
import gzip
import os
import random
import sys

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ROOT      = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TEST_DATA = os.path.join(ROOT, "test", "data")
TEST_DBS  = os.path.join(ROOT, "test", "dbs")

MARKER_FASTAS = {
    "pBI143":   os.path.join(TEST_DBS, "markers", "pBI143.fasta"),
    "crAss001": os.path.join(TEST_DBS, "markers", "crAss001.fasta"),
}

MOCK_SAMPLES = ["mock1", "mock2"]

# Target number of read pairs per marker per sample
N_READS = {
    "pBI143":   200,
    "crAss001": 500,
}

READ_LEN     = 150
INSERT_MEAN  = 350   # mean insert size (outer distance)
INSERT_SD    = 50    # stdev of insert size
BASE_QUAL    = 37    # uniform Phred score (FASTQ char 'F' = 38 — high quality)
SEED         = 99    # reproducibility

# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def read_fasta(path: str) -> str:
    """Return concatenated sequence from a FASTA file (ignores header lines)."""
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith(">"):
                seq_parts.append(line.upper())
    return "".join(seq_parts)

# ---------------------------------------------------------------------------
# Read generation
# ---------------------------------------------------------------------------

def random_insert_size(rng: random.Random) -> int:
    """Draw an insert size from a truncated normal, minimum = 2 * READ_LEN."""
    size = int(rng.gauss(INSERT_MEAN, INSERT_SD))
    return max(size, 2 * READ_LEN)


def simulate_pairs(seq: str, n_pairs: int, read_len: int, rng: random.Random):
    """
    Yield (r1_bases, r2_bases) pairs sampled uniformly from seq.
    r2 is the reverse complement of the downstream fragment end.
    """
    seq_len = len(seq)
    RC = str.maketrans("ACGTN", "TGCAN")

    generated = 0
    attempts = 0
    max_attempts = n_pairs * 100

    while generated < n_pairs and attempts < max_attempts:
        attempts += 1
        insert = random_insert_size(rng)
        if insert > seq_len:
            insert = seq_len
        if insert < 2 * read_len:
            continue

        start = rng.randint(0, seq_len - insert)
        fragment = seq[start : start + insert]

        r1 = fragment[:read_len]
        r2_raw = fragment[-read_len:]
        r2 = r2_raw.translate(RC)[::-1]

        # Replace any non-ACGT with N
        r1 = "".join(b if b in "ACGT" else "N" for b in r1)
        r2 = "".join(b if b in "ACGT" else "N" for b in r2)

        yield r1, r2
        generated += 1

    if generated < n_pairs:
        # Sequence is shorter than ideal; tile the entire sequence
        fragment = seq[:min(read_len * 2, seq_len)]
        r1 = fragment[:read_len].ljust(read_len, "N")
        r2 = fragment[-read_len:].ljust(read_len, "N")
        for _ in range(n_pairs - generated):
            yield r1, r2

# ---------------------------------------------------------------------------
# FASTQ writer
# ---------------------------------------------------------------------------

def append_reads_gz(
    r1_path: str,
    r2_path: str,
    pairs,
    marker_name: str,
    read_offset: int,
):
    """Append FASTQ read pairs to existing gzipped files."""
    qual = chr(BASE_QUAL + 33) * READ_LEN  # Phred + 33 offset for FASTQ

    with gzip.open(r1_path, "at") as fh1, gzip.open(r2_path, "at") as fh2:
        for i, (r1, r2) in enumerate(pairs):
            read_id = f"@{marker_name}_inj_{read_offset + i}/1"
            fh1.write(f"{read_id}\n{r1}\n+\n{qual[:len(r1)]}\n")
            read_id2 = f"@{marker_name}_inj_{read_offset + i}/2"
            fh2.write(f"{read_id2}\n{r2}\n+\n{qual[:len(r2)]}\n")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true", help="Report what would be done without writing files")
    parser.add_argument("--seed", type=int, default=SEED, help=f"Random seed (default {SEED})")
    args = parser.parse_args()

    rng = random.Random(args.seed)

    for marker_name, fasta_path in MARKER_FASTAS.items():
        if not os.path.exists(fasta_path):
            print(f"ERROR: {fasta_path} not found", file=sys.stderr)
            sys.exit(1)

        seq = read_fasta(fasta_path)
        n_pairs = N_READS[marker_name]
        print(f"  {marker_name}: {len(seq):,} bp, generating {n_pairs} read pairs per sample")

        if len(seq) < READ_LEN:
            print(f"  WARNING: {marker_name} sequence ({len(seq)} bp) shorter than read length ({READ_LEN} bp) — reads will be padded with N")

        # Pre-generate pairs once per marker; re-use same list for each sample
        # (seeded, so reproducible regardless of sample order)
        rng_marker = random.Random(args.seed ^ hash(marker_name) & 0xFFFF)
        pairs = list(simulate_pairs(seq, n_pairs, READ_LEN, rng_marker))

        for sample in MOCK_SAMPLES:
            r1_path = os.path.join(TEST_DATA, f"{sample}_R1.fastq.gz")
            r2_path = os.path.join(TEST_DATA, f"{sample}_R2.fastq.gz")

            if not os.path.exists(r1_path) or not os.path.exists(r2_path):
                print(f"  ERROR: {r1_path} or {r2_path} not found", file=sys.stderr)
                sys.exit(1)

            print(f"    → appending {len(pairs)} pairs to {sample} FASTQs", end="")
            if args.dry_run:
                print(" [dry-run, skipped]")
            else:
                append_reads_gz(r1_path, r2_path, pairs, marker_name, read_offset=0)
                print()

    if not args.dry_run:
        print("\nDone. Re-run the pipeline (forcerun short_reads) to see updated marker cpg values.")
        print("Expected: pBI143_cpg ~5-15, crAss001_cpg ~0.3-1.5 in both mock samples.")


if __name__ == "__main__":
    main()
