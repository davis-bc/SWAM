from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import SUPPORTED_PROFILES, app_config_dir, default_dbs_dir, ensure_config_path, runtime_root

CONFIG_THEMES = ("osrs", "classic")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="swam",
        description=(
            "SWAM command-line interface.\n\n"
            "Use `swam config` to edit a user-owned profile in the terminal UI,\n"
            "then run the workflow with `swam run`."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="command")

    config_parser = subparsers.add_parser(
        "config",
        help="Open the SWAM configuration TUI.",
        description=(
            "Create or edit a SWAM profile (pipeline stages, paths, Slurm, Snakemake).\n\n"
            "Examples:\n"
            "  swam config\n"
            "  swam config --profile local\n"
            "  swam config --print-path --profile slurm"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    config_parser.add_argument("--profile", choices=SUPPORTED_PROFILES, default="slurm", help="Profile to edit.")
    config_parser.add_argument("--path", type=Path, help="Edit a specific profile file instead of the user profile.")
    config_parser.add_argument("--theme", choices=CONFIG_THEMES, default="osrs", help="TUI theme.")
    config_parser.add_argument("--print-path", action="store_true", help="Print the resolved profile path and exit.")
    config_parser.set_defaults(handler=handle_config)

    run_parser = subparsers.add_parser(
        "run",
        help="Run the packaged SWAM workflow.",
        description=(
            "Run SWAM from the packaged workflow staged in a writable runtime cache.\n"
            "Unrecognised arguments are passed through to Snakemake.\n\n"
            "Examples:\n"
            "  swam run --profile local --input /data/fastq --output /data/out -n\n"
            "  swam run --profile slurm\n"
            "  swam run --profile local --test\n"
            "  swam run --profile local --test --config run_mags=False -n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    run_parser.add_argument("--profile", choices=SUPPORTED_PROFILES, default="slurm", help="Profile to run with.")
    run_parser.add_argument("--input", dest="input_dir", type=Path, help="FASTQ directory (default: profile in_dir, else ./input).")
    run_parser.add_argument("--output", dest="output_dir", type=Path, help="Output directory (default: profile out_dir, else ./output).")
    run_parser.add_argument("--test", action="store_true", help="Run on the bundled mock samples (outputs to <test_dir>/output).")
    run_parser.add_argument("--test-dir", type=Path, help="Path to a checkout's test/ directory (auto-detected for editable installs).")
    run_parser.add_argument("--print-runtime-root", action="store_true", help="Print the staged runtime path and exit.")
    run_parser.set_defaults(handler=handle_run)
    return parser


def handle_config(args: argparse.Namespace) -> int:
    config_path = ensure_config_path(args.profile, args.path)
    if args.print_path:
        print(config_path)
        return 0
    from .tui import launch_config_tui

    return launch_config_tui(config_path=config_path, profile=args.profile, theme_name=args.theme)


def handle_run(args: argparse.Namespace) -> int:
    if args.print_runtime_root:
        print(runtime_root())
        return 0
    from .runtime import run_workflow

    return run_workflow(
        profile=args.profile,
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        extra_args=args.snakemake_args,
        test=args.test,
        test_dir=args.test_dir,
    )


def main(argv: list[str] | None = None) -> int:
    raw_argv = list(argv) if argv is not None else sys.argv[1:]
    parser = build_parser()
    args, unknown_args = parser.parse_known_args(raw_argv)
    if not hasattr(args, "handler"):
        parser.print_help()
        print(f"\nUser config root:       {app_config_dir()}")
        print(f"Runtime workspace root: {runtime_root()}")
        print(f"Default database dir:   {default_dbs_dir()}")
        return 0
    if args.command == "run":
        args.snakemake_args = unknown_args
    elif unknown_args:
        parser.error(f"unrecognized arguments: {' '.join(unknown_args)}")
    return args.handler(args)
