"""Paths and runtime settings shared by the rest of the backend."""

from __future__ import annotations

import os
import secrets
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

# moebatch runs allowed at once across the whole server - every run checks
# out MOE license tokens shared with other users (6 tokens, 2 per run, at
# the time of writing). Ligands beyond this wait in a queue rather than
# failing on a license error, and tokens stay free for other users.
MOE_MAX_CONCURRENT = int(os.environ.get("MOE_MAX_CONCURRENT", 1))

# How long one moebatch run keeps retrying when the license server is out
# of tokens (someone else is using them), seconds, and how often it retries.
MOE_LICENSE_WAIT = int(os.environ.get("MOE_LICENSE_WAIT", 3600))
MOE_LICENSE_RETRY = 30

# Hard stop for a single ligand, seconds.
SMINA_TIMEOUT = int(os.environ.get("SMINA_TIMEOUT", 1800))

# Poses kept per ligand.
DEFAULT_NUM_MODES = 5

# Ligands allowed in a single job. Kept low by default since this is meant
# for a local/trusted-link setup, not a public multi-tenant service.
MAX_LIGANDS_PER_JOB = int(os.environ.get("DOCKING_MAX_LIGANDS", 50))

# HTTP Basic auth for sharing the app over a tunnel/link. Unset
# DOCKING_AUTH_PASS to leave the app open (fine for pure localhost use). Set
# it explicitly for a stable password across restarts, or leave it to
# auto-generate one printed to the console each start.
AUTH_USER = os.environ.get("DOCKING_AUTH_USER", "colleague")
AUTH_PASS = os.environ.get("DOCKING_AUTH_PASS")
AUTH_GENERATED = False
if AUTH_PASS is None and os.environ.get("DOCKING_REQUIRE_AUTH"):
    AUTH_PASS = secrets.token_urlsafe(9)
    AUTH_GENERATED = True


def _tool_candidates(tool: str) -> list[Path]:
    """Where setup.sh/setup.ps1 could have put a compiled tool in TOOL_ENV.

    Unix micromamba envs put compiled binaries at <prefix>/bin/<tool>; Windows
    conda-forge packages put them under <prefix>\\Library\\bin\\<tool>.exe.
    """
    if os.name == "nt":
        return [
            TOOL_ENV / "Library" / "bin" / f"{tool}.exe",
            TOOL_ENV / "Scripts" / f"{tool}.exe",
            TOOL_ENV / f"{tool}.exe",
        ]
    return [TOOL_ENV / "bin" / tool]


def find_smina() -> str | None:
    """Locate the smina executable: explicit override, local env, then PATH."""
    override = os.environ.get("SMINA_BIN")
    if override and Path(override).exists():
        return override
    for candidate in _tool_candidates("smina"):
        if candidate.exists():
            return str(candidate)
    return shutil.which("smina")


def find_gnina() -> str | None:
    """Locate the gnina executable: explicit override, local env, then PATH.

    gnina is a smina/Vina fork that can run its CNN scoring on a CUDA GPU
    (falling back to CPU automatically when none is present). Upstream ships
    no Windows build, so this only ever resolves there via SMINA_BIN-style
    override or PATH (e.g. a WSL passthrough shim).
    """
    override = os.environ.get("GNINA_BIN")
    if override and Path(override).exists():
        return override
    for candidate in _tool_candidates("gnina"):
        if candidate.exists():
            return str(candidate)
    return shutil.which("gnina")


def find_moebatch() -> str | None:
    """Locate moebatch executable for MOE docking.

    Checks common macOS install paths, versioned installs, then PATH.
    """
    override = os.environ.get("MOEBATCH_BIN")
    if override and Path(override).exists():
        return override

    # Common MOE install paths on macOS
    candidates = [
        Path("/Applications/moe2024.0604/bin/moebatch"),
        Path("/Applications/moe/bin/moebatch"),
    ]

    # Check versioned installs
    apps = Path("/Applications")
    if apps.is_dir():
        for moe_dir in sorted(apps.glob("moe*"), reverse=True):
            candidates.append(moe_dir / "bin" / "moebatch")

    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    return shutil.which("moebatch")


def find_moe_reactions() -> Path | None:
    """Locate MOE covalent reaction template directory."""
    moebatch = find_moebatch()
    if not moebatch:
        return None

    # moebatch at: /Applications/moe2024.0604/bin/moebatch
    # reactions at: /Applications/moe2024.0604/lib/reactions/covalent
    moe_root = Path(moebatch).parent.parent
    rxn_dir = moe_root / "lib" / "reactions" / "covalent"

    return rxn_dir if rxn_dir.is_dir() else None


def list_covalent_reactions() -> list[str]:
    """List available covalent reaction template names."""
    rxn_dir = find_moe_reactions()
    if not rxn_dir:
        return []

    return sorted([f.stem for f in rxn_dir.glob("*.rxn")])


def subprocess_env(exe: str) -> dict[str, str]:
    """Environment for running a bundled tool (smina/gnina) as a subprocess.

    Extends the current environment; when `exe` resolved from TOOL_ENV (our
    micromamba prefix, not a system install) this also points BABEL_LIBDIR at
    the bundled OpenBabel's plugin dir. smina/gnina link against OpenBabel for
    structure I/O (PDB/SDF/mol2/...) and it resolves its format plugins at
    runtime via that variable - without it every read/write silently fails
    ("could not open ... for reading") even though the file exists, since the
    PDB/SDF format handlers themselves never registered.
    """
    env = dict(os.environ)
    try:
        Path(exe).resolve().relative_to(TOOL_ENV.resolve())
    except ValueError:
        return env  # not our bundled install - leave BABEL_LIBDIR (if any) alone
    lib_root = TOOL_ENV / ("Library/lib/openbabel" if os.name == "nt" else "lib/openbabel")
    versions = sorted(lib_root.glob("*")) if lib_root.is_dir() else []
    if versions:
        env["BABEL_LIBDIR"] = str(versions[-1])
    return env


def has_nvidia_gpu() -> bool:
    """Best-effort local hint for the UI; gnina itself decides at run time."""
    return shutil.which("nvidia-smi") is not None


def ensure_dirs() -> None:
    for d in (RUNS_DIR, RECEPTOR_DIR, JOB_DIR):
        d.mkdir(parents=True, exist_ok=True)
