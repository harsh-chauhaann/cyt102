#!/usr/bin/env python3
"""
Shared configuration helpers for the reproducible delta-analysis pipeline.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MP_KEY_ENV = "MP_API_KEY"
DEFAULT_MP_KEY_FILE_ENV = "MP_API_KEY_FILE"
DEFAULT_DFTD3_CMD_ENV = "S_DFTD3_CMD"


def load_dotenv_if_present(env_path: Path | None = None) -> None:
    path = env_path or (REPO_ROOT / ".env")
    if not path.exists():
        return

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if not key or key in os.environ:
            continue
        value = value.strip().strip("'").strip('"')
        os.environ[key] = value


def env_str(name: str, default: str) -> str:
    value = os.environ.get(name)
    if value is None:
        return default
    value = value.strip()
    return value if value else default


def env_int(name: str, default: int) -> int:
    value = os.environ.get(name)
    if value is None or not value.strip():
        return default
    return int(value.strip())


def env_float(name: str, default: float) -> float:
    value = os.environ.get(name)
    if value is None or not value.strip():
        return default
    return float(value.strip())


def env_bool(name: str, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None or not value.strip():
        return default
    return value.strip().lower() in {"1", "true", "yes", "y", "on"}


def env_csv_list(name: str, default: Iterable[str]) -> list[str]:
    value = os.environ.get(name)
    if value is None or not value.strip():
        return list(default)
    return [item.strip() for item in value.split(",") if item.strip()]


def resolve_mp_key_file() -> Path:
    raw = os.environ.get(DEFAULT_MP_KEY_FILE_ENV)
    if raw and raw.strip():
        return Path(raw.strip()).expanduser()
    return Path.home() / ".mp_api_key"


def get_api_key() -> str | None:
    key = os.environ.get(DEFAULT_MP_KEY_ENV)
    if key and key.strip():
        return key.strip()

    key_file = resolve_mp_key_file()
    if key_file.exists():
        value = key_file.read_text(encoding="utf-8").strip()
        if value:
            return value
    return None


def get_dftd3_command(default: str = "s-dftd3") -> str:
    return env_str(DEFAULT_DFTD3_CMD_ENV, default)


load_dotenv_if_present()

