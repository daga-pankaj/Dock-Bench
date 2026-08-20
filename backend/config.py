"""Paths and runtime settings shared by the rest of the backend."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

FRONTEND_DIR = ROOT / "frontend"
RUNS_DIR = Path(os.environ.get("DOCKING_RUNS_DIR", ROOT / "runs")).resolve()
RECEPTOR_DIR = RUNS_DIR / "receptors"
JOB_DIR = RUNS_DIR / "jobs"

# Tool environment created by setup.sh (micromamba prefix).
TOOL_ENV = ROOT / ".env-tools"

# How many ligands may dock at the same time. Each smina call gets one core.
MAX_WORKERS = int(os.environ.get("DOCKING_WORKERS", max(1, (os.cpu_count() or 4) - 1)))

# Hard stop for a single ligand, seconds.
SMINA_TIMEOUT = int(os.environ.get("SMINA_TIMEOUT", 1800))

# Poses kept per ligand.
DEFAULT_NUM_MODES = 5


def find_smina() -> str | None:
    """Locate the smina executable: explicit override, local env, then PATH."""
    override = os.environ.get("SMINA_BIN")
    if override and Path(override).exists():
        return override
    local = TOOL_ENV / "bin" / "smina"
    if local.exists():
        return str(local)
    return shutil.which("smina")


def ensure_dirs() -> None:
    for d in (RUNS_DIR, RECEPTOR_DIR, JOB_DIR):
        d.mkdir(parents=True, exist_ok=True)
