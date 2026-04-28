"""
run_argnet.py

Thin wrapper around an ARGNet checkout. It imports the appropriate upstream
module from the cloned repository, runs inference, then copies the results file
to the requested output path.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import sys
import tempfile
import uuid

from Bio import SeqIO


MODULE_MAP = {
    ("aa", "argnet-l"): ("argnet_lsaa_speed_sgpu.py", "argnet_lsaa"),
    ("nt", "argnet-l"): ("argnet_lsnt.py", "argnet_lsnt"),
    ("aa", "argnet-s"): ("argnet_ssaa_chunk.py", "argnet_ssaa"),
    ("nt", "argnet-s"): ("argnet_ssnt_new_chunk.py", "argnet_ssnt"),
}

AA_ALLOWED = set("ACDEFGHIKLMNPQRSTVWXYBJZX")


def load_module(module_path: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load ARGNet module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True, help="Path to the cloned ARGNet repository")
    parser.add_argument("--input", required=True, help="Input FASTA path")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--seq-type", required=True, choices=["aa", "nt"])
    parser.add_argument("--model", required=True, choices=["argnet-s", "argnet-l"])
    return parser.parse_args()


def prepare_input_fasta(input_path: str, seq_type: str) -> str:
    if seq_type != "aa":
        return input_path

    fd, temp_path = tempfile.mkstemp(suffix=".faa", prefix="argnet_sanitized_")
    os.close(fd)

    cleaned_records = []
    for record in SeqIO.parse(input_path, "fasta"):
        seq = str(record.seq).upper().replace("*", "")
        seq = "".join(base if base in AA_ALLOWED else "X" for base in seq)
        if not seq:
            continue
        record.seq = record.seq.__class__(seq)
        cleaned_records.append(record)

    SeqIO.write(cleaned_records, temp_path, "fasta")
    return temp_path


def main() -> None:
    args = parse_args()
    key = (args.seq_type, args.model)
    if key not in MODULE_MAP:
        raise RuntimeError(f"Unsupported ARGNet mode: {key}")

    os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

    script_name, function_name = MODULE_MAP[key]
    scripts_dir = os.path.join(args.repo, "scripts")
    module_path = os.path.join(scripts_dir, script_name)
    if not os.path.exists(module_path):
        raise FileNotFoundError(f"ARGNet script not found: {module_path}")

    os.makedirs(os.path.join(args.repo, "results"), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    sys.path.insert(0, scripts_dir)
    module = load_module(module_path, f"copilot_{os.path.splitext(script_name)[0]}")
    run_fn = getattr(module, function_name, None)
    if run_fn is None:
        raise AttributeError(f"{script_name} does not define {function_name}")

    out_name = f"copilot-{uuid.uuid4().hex}.tsv"
    prepared_input = prepare_input_fasta(args.input, args.seq_type)
    try:
        run_fn(prepared_input, out_name)
    finally:
        if prepared_input != args.input and os.path.exists(prepared_input):
            os.remove(prepared_input)

    source = os.path.join(args.repo, "results", out_name)
    if not os.path.exists(source):
        raise FileNotFoundError(f"ARGNet did not create the expected results file: {source}")

    shutil.move(source, args.output)


if __name__ == "__main__":
    main()
