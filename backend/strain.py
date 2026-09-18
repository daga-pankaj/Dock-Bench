"""Ligand-only strain energy for a docked pose.

The docking score alone can't tell you a pose forces the ligand into an
unrealistic geometry to fit the site - this can. Defined as the MMFF94
single-point energy of the pose's own conformation minus the lowest MMFF94
energy found across a small random-conformer search of the same molecule (an
approximation of its unconstrained global minimum, not a guarantee - see
N_CONFORMERS below). A negative result (the pose already beat the reference
search) is clipped to 0.

All poses of one ligand share the same molecular graph, so the reference
search only needs to run once per ligand (see strain_energies) rather than
once per pose.
"""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem

N_CONFORMERS = 20
_SEED = 0xF00D


def _mmff_energy(mol: Chem.Mol, conf_id: int = -1) -> float | None:
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        return None  # e.g. an element MMFF has no parameters for
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    if ff is None:
        return None
    return ff.CalcEnergy()


def _with_explicit_hs(mol: Chem.Mol) -> Chem.Mol | None:
    """smina's pose SDF output doesn't carry explicit Hs - MMFF needs them
    for a meaningful energy, so they're added back onto the existing pose
    conformer's heavy-atom positions."""
    out = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None
    return Chem.AddHs(out, addCoords=True)


def reference_minimum_energy(mol: Chem.Mol, n_conformers: int = N_CONFORMERS, seed: int = _SEED) -> float | None:
    """Lowest MMFF94 energy found across `n_conformers` random embeddings of
    `mol`'s own graph - the zero-strain baseline for pose_energy()."""
    prepped = _with_explicit_hs(mol)
    if prepped is None:
        return None
    prepped.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    conf_ids = list(AllChem.EmbedMultipleConfs(prepped, numConfs=n_conformers, params=params))
    if not conf_ids:
        return None
    try:
        results = AllChem.MMFFOptimizeMoleculeConfs(prepped, maxIters=1000)
    except Exception:
        return None
    energies = [energy for converged, energy in results if converged == 0]
    return min(energies) if energies else None


def pose_energy(mol: Chem.Mol) -> float | None:
    """Single-point MMFF94 energy of `mol`'s current (docked) conformation."""
    prepped = _with_explicit_hs(mol)
    if prepped is None:
        return None
    return _mmff_energy(prepped)


def strain_energies(pose_mols: list[Chem.Mol]) -> list[float | None]:
    """One strain value per mol in `pose_mols`, sharing a single reference
    search across all of them (same ligand, same graph, one docking job).
    None for a pose MMFF couldn't score (uncommon elements/valences)."""
    if not pose_mols:
        return []
    reference = reference_minimum_energy(pose_mols[0])
    if reference is None:
        return [None] * len(pose_mols)
    out: list[float | None] = []
    for mol in pose_mols:
        energy = pose_energy(mol)
        out.append(max(0.0, energy - reference) if energy is not None else None)
    return out
