"""Thin wrapper around the smina command line."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem, RDLogger

from . import config

RDLogger.DisableLog("rdApp.*")

SCORING_FUNCTIONS = ("vina", "vinardo", "dkoes_scoring", "ad4_scoring")


class SminaNotFound(RuntimeError):
    pass


class SminaFailed(RuntimeError):
    pass


@dataclass
class Pose:
    rank: int
    affinity: float | None
    rmsd_lb: float | None
    rmsd_ub: float | None


def smina_version() -> str:
    exe = config.find_smina()
    if not exe:
        raise SminaNotFound("smina executable not found; run ./setup.sh")
    out = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=60)
    return (out.stdout or out.stderr).strip().splitlines()[0] if (out.stdout or out.stderr) else "smina"


def build_command(
    receptor: Path,
    ligand: Path,
    out_sdf: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    exhaustiveness: int = 8,
    num_modes: int = config.DEFAULT_NUM_MODES,
    energy_range: float = 4.0,
    seed: int = 42,
    scoring: str = "vina",
    cpu: int = 1,
) -> list[str]:
    exe = config.find_smina()
    if not exe:
        raise SminaNotFound("smina executable not found; run ./setup.sh")
    if scoring not in SCORING_FUNCTIONS:
        raise ValueError(f"unknown scoring function: {scoring}")
    return [
        exe,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--out", str(out_sdf),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", str(energy_range),
        "--seed", str(seed),
        "--scoring", scoring,
        "--cpu", str(cpu),
        "--quiet",
    ]


def parse_poses(sdf_path: Path) -> list[Pose]:
    """Read the scores smina writes as SDF properties on each pose."""
    poses: list[Pose] = []
    if not sdf_path.exists():
        return poses
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            continue

        def prop(name: str) -> float | None:
            if not mol.HasProp(name):
                return None
            try:
                return float(mol.GetProp(name))
            except ValueError:
                return None

        poses.append(
            Pose(
                rank=i,
                affinity=prop("minimizedAffinity"),
                rmsd_lb=prop("minimizedRMSD"),
                rmsd_ub=prop("RMSD"),
            )
        )
    return poses


def dock(
    receptor: Path,
    ligand_sdf: Path,
    out_sdf: Path,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    log_path: Path | None = None,
    **kwargs,
) -> list[Pose]:
    """Dock one ligand and return its poses, best first."""
    cmd = build_command(receptor, ligand_sdf, out_sdf, center, size, **kwargs)
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=config.SMINA_TIMEOUT
    )
    if log_path is not None:
        log_path.write_text(
            "$ " + " ".join(cmd) + "\n\n" + (proc.stdout or "") + "\n" + (proc.stderr or "")
        )
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()
        raise SminaFailed(tail[-1] if tail else f"smina exited {proc.returncode}")
    poses = parse_poses(out_sdf)
    if not poses:
        raise SminaFailed("smina produced no poses (ligand may not fit the box)")
    return poses
