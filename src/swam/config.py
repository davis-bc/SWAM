from __future__ import annotations

import os
import shutil
from collections.abc import Mapping
from datetime import datetime
from importlib.resources import files
from pathlib import Path

from ruamel.yaml import YAML

from . import __version__

APP_NAME = "swam"
SUPPORTED_PROFILES = ("slurm", "local")
OPTIONAL_MAPPING_PROFILE_KEYS = (
    "resources",
    "default-resources",
    "set-threads",
    "set-resources",
    "groups",
    "group-components",
)


def _validate_profile(profile: str) -> str:
    if profile not in SUPPORTED_PROFILES:
        raise ValueError(f"Unsupported profile '{profile}'. Expected one of: {', '.join(SUPPORTED_PROFILES)}.")
    return profile


def _xdg(var: str, fallback: str) -> Path:
    configured = os.environ.get(var)
    if configured:
        return Path(configured).expanduser()
    return Path.home() / fallback


def app_config_dir() -> Path:
    return _xdg("XDG_CONFIG_HOME", ".config") / APP_NAME


def app_cache_dir() -> Path:
    return _xdg("XDG_CACHE_HOME", ".cache") / APP_NAME


def app_data_dir() -> Path:
    return _xdg("XDG_DATA_HOME", ".local/share") / APP_NAME


def default_dbs_dir() -> Path:
    """Reference DBs live outside the versioned runtime cache (which is wiped on upgrade)."""
    return app_data_dir() / "dbs"


def profiles_dir() -> Path:
    return app_config_dir() / "profiles"


def backups_dir() -> Path:
    return app_config_dir() / "backups"


def runtime_root() -> Path:
    return app_cache_dir() / "runtime" / __version__


def default_profile_resource(profile: str):
    profile = _validate_profile(profile)
    return files("swam.assets").joinpath("config", profile, "config.yaml")


def user_profile_path(profile: str) -> Path:
    profile = _validate_profile(profile)
    return profiles_dir() / profile / "config.yaml"


def resolve_config_path(profile: str, explicit_path: str | os.PathLike[str] | None = None) -> Path:
    if explicit_path is None:
        return user_profile_path(profile)
    return Path(explicit_path).expanduser().resolve()


def seed_config(profile: str, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(default_profile_resource(profile).read_text(encoding="utf-8"), encoding="utf-8")
    return destination


def make_yaml() -> YAML:
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)
    return yaml


def load_profile(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        return make_yaml().load(handle)


def prune_empty_profile_sections(data: object) -> bool:
    if not isinstance(data, Mapping):
        return False
    changed = False
    for key in OPTIONAL_MAPPING_PROFILE_KEYS:
        value = data.get(key)
        if isinstance(value, Mapping) and not value:
            del data[key]
            changed = True
    return changed


def sanitize_profile_config(config_path: Path) -> bool:
    if not config_path.exists():
        return False
    data = load_profile(config_path)
    if not prune_empty_profile_sections(data):
        return False
    with config_path.open("w", encoding="utf-8") as handle:
        make_yaml().dump(data, handle)
    return True


def ensure_config_path(profile: str, explicit_path: str | os.PathLike[str] | None = None) -> Path:
    path = resolve_config_path(profile, explicit_path)
    if path.exists():
        sanitize_profile_config(path)
        return path
    return seed_config(profile, path)


def create_backup_for(config_path: Path, profile: str | None = None) -> Path:
    profile_dir = backups_dir() / (profile or "custom")
    profile_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    backup_path = profile_dir / f"{config_path.stem}.{stamp}{config_path.suffix}.bak"
    shutil.copy2(config_path, backup_path)
    return backup_path
