#!/usr/bin/env python3

import argparse
import os
import re
from collections import Counter
from pathlib import Path


SAMPLE_REGEX = re.compile(r"(.+?)(?:_R?[12](?:_001)?|_[12])(?:\.fastq)?$")

RULES = [
    {
        "name": "short_reads",
        "outputs": [
            "data/QAQC/fastp_reports/{sample}.json",
            "data/clean_reads/{sample}_R1.clean.fastq.gz",
            "data/clean_reads/{sample}_R2.clean.fastq.gz",
            "data/QAQC/nonpareil/{sample}.npo",
            "data/alignments/{sample}.afp.res",
            "data/alignments/{sample}.scgs",
            "data/alignments/{sample}.markers.res",
        ],
    },
    {"name": "contigs", "outputs": ["data/megahit/{sample}.contigs.fa"]},
    {
        "name": "genomad",
        "outputs": [
            "data/genomad/{sample}/{sample}.contigs_summary/{sample}.contigs_plasmid.fna",
            "data/genomad/{sample}/{sample}.contigs_summary/{sample}.contigs_plasmid_summary.tsv",
            "data/genomad/{sample}/{sample}.contigs_summary/{sample}.contigs_virus_summary.tsv",
        ],
    },
    {
        "name": "mobmess",
        "outputs": [
            "data/mobmess/{sample}.contigs_circular.txt",
            "data/mobmess/{sample}.contigs_plasmid.bam",
            "data/mobmess/{sample}-mobmess_contigs.txt",
        ],
    },
    {
        "name": "prodigal",
        "outputs": [
            "data/prodigal/{sample}.faa",
            "data/prodigal/{sample}.gff",
        ],
    },
    {"name": "contig_amr", "outputs": ["data/amr_contigs/{sample}_contig_amr.tsv"]},
    {"name": "mge_annotation", "outputs": ["data/mge_contigs/{sample}_mge.tsv"]},
    {"name": "mmseqs_taxonomy", "outputs": ["data/mmseqs/{sample}_lca.tsv"]},
    {
        "name": "contig_abundance",
        "outputs": [
            "data/contig_abundance/{sample}_contigs.bam",
            "data/contig_abundance/{sample}_contigs.bam.bai",
            "data/contig_abundance/{sample}_contig_abundance.tsv",
        ],
    },
    {"name": "bin", "outputs": ["data/bins/{sample}/.binning.done"]},
    {"name": "mag_prodigal", "outputs": ["data/bins/{sample}/.prodigal.done"]},
    {"name": "mag_amr", "outputs": ["data/bins/{sample}/mag_amr.tsv"]},
    {"name": "mag_mge", "outputs": ["data/bins/{sample}/mag_mge.tsv"]},
    {"name": "mag_taxonomy", "outputs": ["data/bins/{sample}/gtdbtk.summary.tsv"]},
    {"name": "mag_metabolism", "outputs": ["data/bins/{sample}/.metabolic.done"]},
    {"name": "mag_qc", "outputs": ["data/bins/{sample}/checkm2_quality.tsv"]},
    {"name": "mag_abundance", "outputs": ["data/bins/{sample}/mag_abundance.tsv"]},
]

RETRYABLE_CAUSES = {"SLURM", "WALLTIME", "TRANSIENT_IO", "FILESYSTEM", "OOM"}


def extract_sample_name(filename: str) -> str | None:
    base = os.path.basename(filename)
    if base.endswith(".gz"):
        base = base[:-3]
    match = SAMPLE_REGEX.match(base)
    return match.group(1) if match else None


def discover_sample_pairs(input_dir: Path) -> dict[str, dict[str, str]]:
    r1_files = sorted(list(input_dir.glob("*R1*.fastq*")) + list(input_dir.glob("*_1.fastq*")))
    r2_files = sorted(list(input_dir.glob("*R2*.fastq*")) + list(input_dir.glob("*_2.fastq*")))

    pairs: dict[str, dict[str, str]] = {}
    for r1 in r1_files:
        sample = extract_sample_name(r1.name)
        if sample is None:
            continue

        r2_candidates = [
            str(r1).replace("_R1", "_R2"),
            str(r1).replace("_1", "_2"),
        ]
        found_r2 = next((candidate for candidate in r2_candidates if Path(candidate) in r2_files), None)
        if found_r2 is None:
            for r2 in r2_files:
                if extract_sample_name(r2.name) == sample:
                    found_r2 = str(r2)
                    break

        if found_r2 is not None:
            pairs[sample] = {"r1": str(r1), "r2": found_r2}

    return dict(sorted(pairs.items()))


def classify_failure(log_tail: str) -> tuple[str, str]:
    text = log_tail.lower()
    if not text.strip():
        return "UNKNOWN", "Inspect the rule log and Snakemake driver log."
    if "out of memory" in text or "oom" in text or "cannot allocate memory" in text or "killed" in text:
        return "OOM", "Increase mem_mb for this rule or reduce group-components before rerunning."
    if "time limit" in text or "due to time limit" in text or "timelimit" in text:
        return "WALLTIME", "Increase runtime for this rule and rerun the affected samples."
    if "slurmstepd" in text or "preempt" in text or "cancelled" in text:
        return "SLURM", "Retry the affected samples; this looks like a scheduler-side interruption."
    if "stale file handle" in text or "input/output error" in text or "temporarily unavailable" in text:
        return "FILESYSTEM", "Retry the affected samples; this looks like shared-filesystem instability."
    if "connection reset" in text or "connection timed out" in text or "network is unreachable" in text:
        return "TRANSIENT_IO", "Retry the affected samples; this looks transient."
    if "no such file or directory" in text or "missing input" in text:
        return "MISSING_INPUT", "Check upstream outputs, database staging, or sample naming before rerunning."
    if "command not found" in text or "environmentlocationnotfound" in text or "conda" in text:
        return "ENVIRONMENT", "Check the conda environment and tool installation before rerunning."
    if "dbs/" in text and ("missing" in text or "not found" in text):
        return "DATABASE", "Re-run the database init rules on the head node and rerun affected samples."
    return "TOOL_ERROR", "Inspect the rule log for the underlying tool error before rerunning."


def log_path_for(output_dir: Path, rule_name: str, sample: str) -> Path:
    return output_dir / "logs" / rule_name / f"{sample}.log"


def outputs_for(output_dir: Path, rule: dict, sample: str) -> list[Path]:
    return [output_dir / template.format(sample=sample) for template in rule["outputs"]]


def has_rule_evidence(output_dir: Path, samples: list[str], rule: dict) -> bool:
    for sample in samples:
        if log_path_for(output_dir, rule["name"], sample).exists():
            return True
        if any(path.exists() for path in outputs_for(output_dir, rule, sample)):
            return True
    return False


def determine_rule_status(output_dir: Path, sample: str, rule: dict, blocked: bool) -> dict:
    log_path = log_path_for(output_dir, rule["name"], sample)
    outputs = outputs_for(output_dir, rule, sample)
    outputs_exist = all(path.exists() for path in outputs)
    log_tail = ""
    if log_path.exists():
        try:
            log_tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-40:])
        except OSError:
            log_tail = ""

    if outputs_exist:
        status = "SUCCESS"
        cause = ""
        action = ""
    elif log_path.exists():
        status = "FAILED"
        cause, action = classify_failure(log_tail)
    elif blocked:
        status = "BLOCKED"
        cause = ""
        action = ""
    else:
        status = "NOT_RUN"
        cause = ""
        action = ""

    return {
        "sample": sample,
        "rule": rule["name"],
        "status": status,
        "likely_cause": cause,
        "suggested_action": action,
        "log_path": str(log_path),
    }


def summarize_sample(rule_rows: list[dict]) -> dict:
    first_failure = next((row for row in rule_rows if row["status"] == "FAILED"), None)
    if all(row["status"] == "SUCCESS" for row in rule_rows):
        overall_status = "COMPLETED"
        failed_rule = ""
        likely_cause = ""
        suggested_action = "No action needed."
    elif first_failure:
        failed_rule = first_failure["rule"]
        likely_cause = first_failure["likely_cause"]
        suggested_action = first_failure["suggested_action"]
        overall_status = "NEEDS_RETRY" if likely_cause in RETRYABLE_CAUSES else "NEEDS_REVIEW"
    elif any(row["status"] == "SUCCESS" for row in rule_rows):
        overall_status = "PARTIAL"
        failed_rule = ""
        likely_cause = ""
        suggested_action = "Resume the run after reviewing blocked or not-run stages."
    else:
        overall_status = "NOT_STARTED"
        failed_rule = ""
        likely_cause = ""
        suggested_action = "Submit the run for this sample."

    completed_rules = sum(row["status"] == "SUCCESS" for row in rule_rows)
    last_completed_rule = ""
    for row in rule_rows:
        if row["status"] == "SUCCESS":
            last_completed_rule = row["rule"]
        else:
            break

    return {
        "sample": rule_rows[0]["sample"] if rule_rows else "",
        "overall_status": overall_status,
        "completed_rules": str(completed_rules),
        "total_rules": str(len(rule_rows)),
        "last_completed_rule": last_completed_rule,
        "failed_rule": failed_rule,
        "likely_cause": likely_cause,
        "suggested_action": suggested_action,
    }


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(fieldnames) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(field, "")) for field in fieldnames) + "\n")


def build_text_report(report_path: Path, sample_rows: list[dict], failed_rows: list[dict], driver_log: str) -> None:
    counts = Counter(row["overall_status"] for row in sample_rows)
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("SWAM run report\n")
        handle.write("===================\n\n")
        if driver_log:
            handle.write(f"Driver log: {driver_log}\n\n")

        handle.write("Sample summary\n")
        handle.write("--------------\n")
        for status in ("COMPLETED", "NEEDS_RETRY", "NEEDS_REVIEW", "PARTIAL", "NOT_STARTED"):
            handle.write(f"{status}: {counts.get(status, 0)}\n")
        handle.write("\n")

        handle.write("Top failed sample-step pairs\n")
        handle.write("----------------------------\n")
        if not failed_rows:
            handle.write("No failed sample-step pairs detected.\n")
        else:
            for row in failed_rows[:15]:
                handle.write(
                    f"{row['sample']}: {row['rule']} -> {row['likely_cause'] or 'UNKNOWN'}\n"
                    f"  log: {row['log_path']}\n"
                    f"  action: {row['suggested_action']}\n"
                )
        handle.write("\n")

        handle.write("Completed samples\n")
        handle.write("-----------------\n")
        completed = [row["sample"] for row in sample_rows if row["overall_status"] == "COMPLETED"]
        handle.write(", ".join(completed) + "\n" if completed else "None\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a post-run status report for SWAM.")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--driver-log", default="")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    run_status_dir = output_dir / "logs" / "run_status"
    samples = list(discover_sample_pairs(input_dir).keys())
    active_rules = [rule for rule in RULES if has_rule_evidence(output_dir, samples, rule)]
    if not active_rules:
        active_rules = RULES

    rule_rows: list[dict] = []
    sample_rows: list[dict] = []

    for sample in samples:
        sample_rule_rows = []
        blocked = False
        for rule in active_rules:
            row = determine_rule_status(output_dir, sample, rule, blocked)
            sample_rule_rows.append(row)
            if row["status"] == "FAILED":
                blocked = True
        rule_rows.extend(sample_rule_rows)
        sample_rows.append(summarize_sample(sample_rule_rows))

    failed_rows = [row for row in rule_rows if row["status"] == "FAILED"]
    write_tsv(
        run_status_dir / "sample_rule_status.tsv",
        rule_rows,
        ["sample", "rule", "status", "likely_cause", "suggested_action", "log_path"],
    )
    build_text_report(run_status_dir / "run_report.txt", sample_rows, failed_rows, args.driver_log)


if __name__ == "__main__":
    main()
