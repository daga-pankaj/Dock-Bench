"""RMSD of each docked pose against the receptor's own bound ligand - the
standard "self-docking" sanity check: strip the co-crystallized ligand,
redock it, and see how far the top pose drifted from where it really sat.

Computed with obrms (David Koes's OpenBabel-based heavy-atom RMSD tool,
bundled next to smina/gnina in TOOL_ENV - see backend.config.find_obrms).
obrms works out its own graph-based, symmetry-aware atom correspondence
between the reference and each pose (so e.g. a flipped phenyl ring doesn't
inflate the RMSD), and does it straight from the reference's raw PDB
coordinates - no bond perception or template SMILES needed on our side. A
pose obrms can't match to the reference's heavy-atom graph at all (wrong
bound ligand picked as the reference, or a genuinely different molecule)
comes back reported as "inf", which this module turns into None rather than
failing the whole batch.
"""

from __future__ import annotations

import math
import subprocess
import tempfile
from pathlib import Path

from . import config

OBRMS_TIMEOUT = 60


def rmsds_to_reference(reference_pdb_block: str, poses_sdf_path: Path | str, n_poses: int) -> list[float | None]:
    """RMSD (Å) of each pose in `poses_sdf_path` (an SDF on disk, in the same
    order the job reports poses in) to `reference_pdb_block` (raw HETATM
    lines for the bound ligand - see prep.het_ligand_pdb_block). Always
    `n_poses` long: every entry is None if obrms isn't installed or the call
    fails outright, and an individual entry is None if obrms couldn't match
    that one pose to the reference at all.
    """
    exe = config.find_obrms()
    if not exe:
        return [None] * n_poses

    with tempfile.TemporaryDirectory() as tmp:
        ref_path = Path(tmp) / "reference.pdb"
        ref_path.write_text(reference_pdb_block)
        try:
            result = subprocess.run(
                [exe, "--firstonly", str(ref_path), str(poses_sdf_path)],
                capture_output=True,
                text=True,
                timeout=OBRMS_TIMEOUT,
                env=config.subprocess_env(exe),
            )
        except (subprocess.TimeoutExpired, OSError):
            return [None] * n_poses

    if result.returncode != 0:
        return [None] * n_poses

    values: list[float | None] = []
    for line in result.stdout.splitlines():
        if not line.startswith("RMSD "):
            continue
        try:
            v = float(line.rsplit(maxsplit=1)[-1])
        except ValueError:
            continue
        values.append(v if math.isfinite(v) else None)
    values += [None] * (n_poses - len(values))
    return values[:n_poses]
