#!/usr/bin/env bash
# Backwards-compatible shim: forwards to the packaged CLI on the Slurm profile.
#   bash run_swam.sh [input_dir] [output_dir] [extra snakemake args...]
# Equivalent to:
#   swam run --profile slurm --input <input_dir> --output <output_dir> [extra args...]
# Configure the profile first with:  swam config --profile slurm

set -euo pipefail

if ! command -v swam >/dev/null 2>&1; then
    echo "swam CLI not found. Install it from the repo with:  python -m pip install -e ." >&2
    exit 1
fi

INPUT="${1:-$(pwd)/input}"
OUTPUT="${2:-$(pwd)/output}"
shift $(( $# >= 1 ? 1 : 0 ))
shift $(( $# >= 1 ? 1 : 0 ))

exec swam run --profile slurm --input "$INPUT" --output "$OUTPUT" "$@"
