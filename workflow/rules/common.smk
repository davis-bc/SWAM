import glob
import os
import re

# ---------------------------------------------------------------------------
#   Test-mode path overrides
#   Activated with:  snakemake --config test=True
# ---------------------------------------------------------------------------

_TEST   = config.get("test", False)
_DB_TAG = ".test" if _TEST else ""
_REPO   = os.path.dirname(os.path.dirname(os.path.dirname(workflow.snakefile)))

if _TEST:
    input_dir  = os.path.join(_REPO, "test", "data")
    output_dir = os.path.join(_REPO, "test", "output")
    # GTDB-tk, METABOLIC, and CheckM2 are skipped in test mode
    config["skip_gtdbtk"]    = True
    config["skip_metabolic"] = True
    config["skip_checkm2"]   = True
    # Increase geNomad splits to reduce peak memory on low-RAM machines
    config.setdefault("genomad_splits", 4)
else:
    input_dir  = config.get("in_dir", "")
    output_dir = config.get("out_dir", "")

# ---------------------------------------------------------------------------
#   Database paths
#   All production databases live in SWAM-meta/dbs/ under one subdirectory per
#   tool or reference bundle. Legacy paths are kept here so init rules can
#   migrate older layouts in-place without forcing unnecessary re-downloads.
# ---------------------------------------------------------------------------

_DBS_DIR = os.path.join(_REPO, "dbs")

_DB_AMR_KMA_DIR  = os.path.join(_DBS_DIR, "amrfinderplus_kma")
_DB_AMR_DIR      = os.path.join(_DBS_DIR, "amrfinderplus")
_DB_HUMAN_DIR    = os.path.join(_DBS_DIR, "human")
_DB_MARKERS_DIR  = os.path.join(_DBS_DIR, "markers")
_DB_SCG_DIR      = os.path.join(_DBS_DIR, "scg")
_DB_GENOMAD_DIR  = os.path.join(_DBS_DIR, "genomad")
_DB_UNIREF50_DIR = os.path.join(_DBS_DIR, "uniref50")
_DB_CHECKM2_DIR  = os.path.join(_DBS_DIR, "checkm2")
_DB_METABOLIC_DIR = os.path.join(_DBS_DIR, "metabolic")
_DB_ARGNET_DIR   = os.path.join(_DBS_DIR, "argnet")

_LEGACY_SR_DBS_DIR     = os.path.join(_DBS_DIR, "short_reads")
_LEGACY_GENOMAD_DB_DIR = os.path.join(_DBS_DIR, "genomad_db")
_LEGACY_AFP_DB_DIR     = os.path.join(_DBS_DIR, "amrfinderplus_db")

_SCG_SOURCE_FASTA = (
    os.path.join(_REPO, "test", "dbs", "scg", "SCGs_40.fasta")
    if _TEST else os.path.join(_REPO, "workflow", "resources", "SCGs_40_All.fasta")
)

_SCG_DB = os.path.join(_DB_SCG_DIR, "scg_db")
_UNIREF50_DB = os.path.join(_DB_UNIREF50_DIR, "uniref50_mmseqs")
_AFP_DB_DIR = _DB_AMR_DIR
_CHECKM2_DB = os.path.join(_DB_CHECKM2_DIR, "CheckM2_database", "uniref100.KO.1.dmnd")
_METABOLIC_DIR = os.path.join(_DB_METABOLIC_DIR, "METABOLIC")
_ARGNET_REPO_DIR = os.path.join(_DB_ARGNET_DIR, "ARGNet")

# ---------------------------------------------------------------------------
#   Output path helpers
# ---------------------------------------------------------------------------

FASTP_SUMMARY = os.path.join(output_dir, "fastp_summary.csv")
SHORT_READS_OUTPUT = os.path.join(output_dir, "short_reads_output.csv")
MARKERS_CPG_OUTPUT = os.path.join(output_dir, "data", "QAQC", "markers_cpg.csv")
ASSEMBLY_QA_OUTPUT = os.path.join(output_dir, "assembly_qa.tsv")
CONTIG_SUMMARY_OUTPUT = os.path.join(output_dir, "contig_summary.tsv")
AMR_ABUNDANCE_SUMMARY_OUTPUT = os.path.join(output_dir, "AMR_abundance_summary.csv")
MAG_SUMMARY_OUTPUT = os.path.join(output_dir, "mag_summary.tsv")
ARGNET_CONTIG_SUMMARY_OUTPUT = os.path.join(output_dir, "ARGnet_contig_summary.tsv")
ARGNET_MAG_SUMMARY_OUTPUT = os.path.join(output_dir, "ARGnet_mag_summary.tsv")


def sample_log(rule_name):
    return os.path.join(output_dir, "logs", rule_name, "{sample}.log")


def rule_log(rule_name):
    return os.path.join(output_dir, "logs", f"{rule_name}.log")

# ---------------------------------------------------------------------------
#   Resource helper: returns test value when running in test mode
# ---------------------------------------------------------------------------

def res(production, test=4000):
    """Return test-mode resource value when config test=True, else production."""
    return test if _TEST else production


def attempt_value(values, attempt):
    return values[min(attempt, len(values)) - 1]

# ---------------------------------------------------------------------------
#   Module flags — control which pipeline stages are active
# ---------------------------------------------------------------------------

_RUN_SR   = config.get("run_short_reads", True)
_RUN_CTG  = config.get("run_contigs",     True)
_RUN_MAGS = config.get("run_mags",        True)

# ---------------------------------------------------------------------------
#   Config validation for external-input scenarios
# ---------------------------------------------------------------------------

from snakemake.exceptions import WorkflowError

_AMR_MODE = str(config.get("amr_mode", "operational")).strip().lower()
if _AMR_MODE not in {"operational", "discovery"}:
    raise WorkflowError(
        "amr_mode must be either 'operational' or 'discovery' in config.yaml."
    )
_RUN_ARGNET = _AMR_MODE == "discovery"

if not _RUN_SR and (_RUN_CTG or _RUN_MAGS) and not config.get("clean_reads_dir"):
    raise WorkflowError(
        "run_short_reads: False requires 'clean_reads_dir' in config.yaml "
        "(directory with pre-filtered {sample}*R1*.fastq* read pairs)."
    )

if not _RUN_CTG and _RUN_MAGS and not config.get("contigs_dir"):
    raise WorkflowError(
        "run_contigs: False with run_mags: True requires 'contigs_dir' in config.yaml "
        "(directory with pre-assembled {sample}.contigs.fa files)."
    )

# ---------------------------------------------------------------------------
#   Input helper functions for dynamic stage routing
# ---------------------------------------------------------------------------

def get_clean_r1(wildcards):
    """Return path to host-filtered R1 reads."""
    if _RUN_SR:
        return os.path.join(output_dir, "data", "clean_reads", f"{wildcards.sample}_R1.clean.fastq.gz")
    for pattern in [
        os.path.join(config["clean_reads_dir"], f"{wildcards.sample}*R1*.fastq*"),
        os.path.join(config["clean_reads_dir"], f"{wildcards.sample}*_1.fastq*"),
    ]:
        hits = glob.glob(pattern)
        if hits:
            return sorted(hits)[0]
    raise WorkflowError(
        f"No R1 reads found for sample '{wildcards.sample}' in clean_reads_dir: {config['clean_reads_dir']}"
    )

def get_clean_r2(wildcards):
    """Return path to host-filtered R2 reads."""
    if _RUN_SR:
        return os.path.join(output_dir, "data", "clean_reads", f"{wildcards.sample}_R2.clean.fastq.gz")
    for pattern in [
        os.path.join(config["clean_reads_dir"], f"{wildcards.sample}*R2*.fastq*"),
        os.path.join(config["clean_reads_dir"], f"{wildcards.sample}*_2.fastq*"),
    ]:
        hits = glob.glob(pattern)
        if hits:
            return sorted(hits)[0]
    raise WorkflowError(
        f"No R2 reads found for sample '{wildcards.sample}' in clean_reads_dir: {config['clean_reads_dir']}"
    )

def get_contigs(wildcards):
    """Return path to assembled contigs (MEGAHIT output or user-provided)."""
    if _RUN_CTG:
        return os.path.join(output_dir, "data", "megahit", f"{wildcards.sample}.contigs.fa")
    return os.path.join(config["contigs_dir"], f"{wildcards.sample}.contigs.fa")

# ----------------------------------------------
#   Build a mapping from sample R1 and R2 files
# ----------------------------------------------

# Helper: extract sample name from fastq filename
def extract_sample_name(filename):
    base = os.path.basename(filename)
    if base.endswith('.gz'):
        base = base[:-3]
    m = re.match(r"(.+?)(?:_R?1(?:_001)?|_1)(?:\.fastq)?$", base)
    if m:
        return m.group(1)
    else:
        return re.sub(r"(_R?1(?:_001)?|_1)(\.fastq)?(\.gz)?$", "", base)

# Find R1/R2 files
r1_files = glob.glob(os.path.join(input_dir, "*R1*.fastq*")) + glob.glob(os.path.join(input_dir, "*_1.fastq*"))
r2_files = glob.glob(os.path.join(input_dir, "*R2*.fastq*")) + glob.glob(os.path.join(input_dir, "*_2.fastq*"))

# Map sample names to paired files
sample_to_files = {}
for r1 in r1_files:
    sample = extract_sample_name(r1)
    possible_r2s = [r1.replace('R1', 'R2'), r1.replace('_1', '_2')]
    found_r2 = None
    for r2_candidate in possible_r2s:
        if os.path.exists(r2_candidate):
            found_r2 = r2_candidate
            break
    if not found_r2:
        for r2 in r2_files:
            if extract_sample_name(r2) == sample:
                found_r2 = r2
                break
    if found_r2:
        sample_to_files[sample] = {'r1': r1, 'r2': found_r2}

samples = sorted(sample_to_files.keys())

def get_r1(wildcards):
    return sample_to_files[wildcards.sample]['r1']

def get_r2(wildcards):
    return sample_to_files[wildcards.sample]['r2']
