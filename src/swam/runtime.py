from __future__ import annotations

import os
import shutil
import subprocess
import sys
from datetime import datetime
from hashlib import sha256
from importlib.resources import files
from pathlib import Path

from .config import default_dbs_dir, ensure_config_path, load_profile, runtime_root

RUNTIME_SENTINEL = ".swam-runtime-version"
CONFIG_FLAGS = ("--config", "-C")


def _staged_resources():
    """(resource, relative destination) pairs copied into the runtime root."""
    assets = files("swam.assets")
    return [
        (assets.joinpath("workflow"), Path("workflow")),
        (assets.joinpath("config", "config.yaml"), Path("config") / "config.yaml"),
    ]


def _copy_resource(resource, destination: Path) -> None:
    if resource.is_dir():
        destination.mkdir(parents=True, exist_ok=True)
        for child in resource.iterdir():
            if child.name == "__pycache__":
                continue
            _copy_resource(child, destination / child.name)
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(resource.read_bytes())


def _fingerprint(resource, digest, prefix: str) -> None:
    if resource.is_dir():
        digest.update(f"dir:{prefix}/".encode("utf-8"))
        for child in sorted(resource.iterdir(), key=lambda item: item.name):
            if child.name == "__pycache__":
                continue
            _fingerprint(child, digest, f"{prefix}/{child.name}")
    else:
        digest.update(f"file:{prefix}".encode("utf-8"))
        digest.update(resource.read_bytes())


def ensure_runtime_workspace() -> Path:
    """Stage the packaged workflow into a writable, versioned runtime cache.

    Re-staged whenever the packaged assets change (fingerprint mismatch). Reference
    databases are NOT kept here (see dbs_dir) because this directory is wiped.
    """
    root = runtime_root()
    sentinel = root / RUNTIME_SENTINEL
    digest = sha256()
    for resource, rel in _staged_resources():
        _fingerprint(resource, digest, str(rel))
    expected = digest.hexdigest()
    staged = sentinel.read_text(encoding="utf-8").strip() if sentinel.exists() else ""
    if not (root / "workflow" / "Snakefile").exists() or staged != expected:
        # Only the staged assets are replaced; .snakemake/ (conda envs, metadata) is kept.
        for _resource, rel in _staged_resources():
            target = root / rel
            if target.is_dir():
                shutil.rmtree(target)
            elif target.exists():
                target.unlink()
        root.mkdir(parents=True, exist_ok=True)
        for resource, rel in _staged_resources():
            _copy_resource(resource, root / rel)
        sentinel.write_text(expected, encoding="utf-8")
    return root


def source_test_dir() -> Path | None:
    """test/ of a source checkout (editable install); None for a packaged install."""
    candidate = Path(__file__).resolve().parents[2] / "test"
    return candidate if (candidate / "data").is_dir() else None


def profile_pipeline_config(profile: str) -> dict:
    data = load_profile(ensure_config_path(profile)) or {}
    block = data.get("config") or {}
    return dict(block)


def split_config_args(extra_args: list[str]) -> tuple[list[str], list[str]]:
    """Pull `--config k=v ...` out of passthrough args.

    Snakemake keeps only the last --config flag, so user overrides are merged into
    the single --config we emit instead of silently dropping in_dir/out_dir.
    """
    config_items: list[str] = []
    remaining: list[str] = []
    index = 0
    while index < len(extra_args):
        token = extra_args[index]
        if token in CONFIG_FLAGS:
            index += 1
            while index < len(extra_args) and not extra_args[index].startswith("-"):
                config_items.append(extra_args[index])
                index += 1
            continue
        if token.startswith("--config="):
            config_items.append(token.split("=", 1)[1])
        else:
            remaining.append(token)
        index += 1
    return config_items, remaining


def resolve_run_paths(
    profile_config: dict,
    input_dir: str | os.PathLike[str] | None,
    output_dir: str | os.PathLike[str] | None,
) -> tuple[Path, Path]:
    cwd = Path.cwd()
    chosen_input = input_dir or profile_config.get("in_dir") or cwd / "input"
    chosen_output = output_dir or profile_config.get("out_dir") or cwd / "output"
    return Path(chosen_input).expanduser().resolve(), Path(chosen_output).expanduser().resolve()


def snakemake_command(
    runtime_dir: Path,
    profile: str,
    config_items: list[str],
    extra_args: list[str],
) -> list[str]:
    command = [
        "snakemake",
        "--snakefile",
        str(runtime_dir / "workflow" / "Snakefile"),
        "--directory",
        str(runtime_dir),
        "--profile",
        str(ensure_config_path(profile).parent),
        "--scheduler",
        "greedy",
        "--config",
        *config_items,
    ]
    if profile == "slurm":
        command.extend(["--local-cores", "1", "--quiet", "rules"])
    command.extend(extra_args)
    return command


def build_report(runtime_dir: Path, input_dir: Path, output_dir: Path, driver_log: Path) -> int:
    script_path = runtime_dir / "workflow" / "scripts" / "build_run_report.py"
    command = [
        sys.executable,
        str(script_path),
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--driver-log",
        str(driver_log),
    ]
    return subprocess.run(command, cwd=runtime_dir).returncode


def run_workflow(
    profile: str,
    input_dir: str | os.PathLike[str] | None,
    output_dir: str | os.PathLike[str] | None,
    extra_args: list[str],
    test: bool = False,
    test_dir: str | os.PathLike[str] | None = None,
) -> int:
    if shutil.which("snakemake") is None:
        print("snakemake was not found on PATH. Install Snakemake in the active environment before running SWAM.", file=sys.stderr)
        return 1

    profile_config = profile_pipeline_config(profile)
    user_config_items, passthrough = split_config_args(extra_args)

    if test:
        resolved_test = Path(test_dir).expanduser().resolve() if test_dir else source_test_dir()
        if resolved_test is None or not (resolved_test / "data").is_dir():
            print(
                "Test mode needs the mock data in a source checkout's test/ directory. "
                "Pass --test-dir /path/to/SWAM/test.",
                file=sys.stderr,
            )
            return 1
        # common.smk hard-wires test input/output under test_dir.
        resolved_input, resolved_output = resolved_test / "data", resolved_test / "output"
    else:
        resolved_input, resolved_output = resolve_run_paths(profile_config, input_dir, output_dir)

    dbs_dir = Path(profile_config.get("dbs_dir") or default_dbs_dir()).expanduser().resolve()
    config_items = [f"in_dir={resolved_input}", f"out_dir={resolved_output}", f"dbs_dir={dbs_dir}"]
    if test:
        config_items += ["test=True", f"test_dir={resolved_test}"]
    config_items += user_config_items

    runtime_dir = ensure_runtime_workspace()
    dbs_dir.mkdir(parents=True, exist_ok=True)
    run_status_dir = resolved_output / "logs" / "run_status"
    run_status_dir.mkdir(parents=True, exist_ok=True)
    driver_log = run_status_dir / f"snakemake.{datetime.now().strftime('%Y%m%d-%H%M%S')}.log"

    command = snakemake_command(runtime_dir, profile, config_items, passthrough)
    print(f"Input:      {resolved_input}")
    print(f"Output:     {resolved_output}")
    print(f"Databases:  {dbs_dir}")
    print(f"Profile:    {ensure_config_path(profile).parent}/")
    print(f"Runtime:    {runtime_dir}")
    print(f"Driver log: {driver_log}")
    print("")

    with driver_log.open("w", encoding="utf-8") as handle:
        process = subprocess.Popen(
            command,
            cwd=runtime_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            handle.write(line)
        process.wait()

    report_rc = build_report(runtime_dir, resolved_input, resolved_output, driver_log)
    print("")
    print(f"Run report:        {run_status_dir / 'run_report.txt'}")
    print(f"Rule status table: {run_status_dir / 'sample_rule_status.tsv'}")
    if process.returncode != 0:
        return process.returncode
    return report_rc
