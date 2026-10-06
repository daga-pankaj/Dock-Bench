"""Tests for backend.rmsd: RMSD of each docked pose against a receptor's own
bound ligand, computed with obrms straight off the reference's HETATM lines
and the on-disk poses SDF (no template, no bonds needed on our side)."""

from __future__ import annotations

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from backend import config, rmsd


def _mol(smiles: str, seed: int = 1) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.RemoveHs(mol)


def _pdb_line(serial: int, name: str, resname: str, chain: str, resseq: str, x: float, y: float, z: float, element: str) -> str:
    """One HETATM record with fields at the exact fixed-width columns a real
    receptor PDB uses (see prep._xyz/_resname/_element)."""
    line = [" "] * 80
    line[0:6] = list("HETATM")
    line[6:11] = list(f"{serial:>5}")
    line[12:16] = list(f" {name:<3}"[:4])
    line[17:20] = list(f"{resname:>3}")
    line[21] = chain
    line[22:26] = list(f"{resseq:>4}")
    line[30:38] = list(f"{x:8.3f}")
    line[38:46] = list(f"{y:8.3f}")
    line[46:54] = list(f"{z:8.3f}")
    line[76:78] = list(f"{element:>2}")
    return "".join(line)


def _het_pdb_block(mol: Chem.Mol, resname: str = "LIG", chain: str = "A", resseq: str = "1") -> str:
    conf = mol.GetConformer()
    lines = []
    for i, atom in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        lines.append(_pdb_line(i + 1, atom.GetSymbol(), resname, chain, resseq, p.x, p.y, p.z, atom.GetSymbol()))
    return "\n".join(lines) + "\nEND\n"


def _poses_sdf(tmp_path, mols: list[Chem.Mol]):
    path = tmp_path / "poses.sdf"
    with Chem.SDWriter(str(path)) as writer:
        for mol in mols:
            writer.write(mol)
    return path


@pytest.fixture(autouse=True)
def _require_obrms():
    if not config.find_obrms():
        pytest.skip("obrms not installed")


def test_rmsds_to_reference_near_zero_for_the_same_pose(tmp_path):
    mol = _mol("c1ccccc1C(=O)O")  # benzoic acid
    ref = _het_pdb_block(mol)
    poses = _poses_sdf(tmp_path, [mol])
    values = rmsd.rmsds_to_reference(ref, poses, 1)
    assert len(values) == 1
    assert values[0] is not None
    assert values[0] < 1e-2


def test_rmsds_to_reference_reflects_a_translation(tmp_path):
    mol = _mol("c1ccccc1C(=O)O")
    ref = _het_pdb_block(mol)
    translated = Chem.Mol(mol)
    conf = translated.GetConformer()
    for i in range(translated.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x + 2.0, p.y, p.z))
    poses = _poses_sdf(tmp_path, [translated])
    values = rmsd.rmsds_to_reference(ref, poses, 1)
    assert values[0] is not None
    assert abs(values[0] - 2.0) < 1e-2


def test_rmsds_to_reference_tolerates_reordered_atoms(tmp_path):
    mol = _mol("c1ccccc1C(=O)O")
    ref = _het_pdb_block(mol)
    reordered = Chem.RenumberAtoms(mol, list(range(mol.GetNumAtoms()))[::-1])
    poses = _poses_sdf(tmp_path, [reordered])
    values = rmsd.rmsds_to_reference(ref, poses, 1)
    assert values[0] is not None
    assert values[0] < 1e-2


def test_rmsds_to_reference_none_for_a_different_molecule(tmp_path):
    mol = _mol("c1ccccc1C(=O)O")
    ref = _het_pdb_block(mol)
    other = _mol("CCO")  # ethanol - not even the same heavy-atom composition
    poses = _poses_sdf(tmp_path, [other])
    values = rmsd.rmsds_to_reference(ref, poses, 1)
    assert values == [None]


def test_rmsds_to_reference_one_per_pose_same_order(tmp_path):
    mol = _mol("c1ccccc1C(=O)O")
    ref = _het_pdb_block(mol)
    translated = Chem.Mol(mol)
    conf = translated.GetConformer()
    for i in range(translated.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x + 5.0, p.y, p.z))
    poses = _poses_sdf(tmp_path, [mol, translated])
    values = rmsd.rmsds_to_reference(ref, poses, 2)
    assert len(values) == 2
    assert values[0] < values[1]


def test_rmsds_to_reference_missing_obrms_returns_all_none(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "find_obrms", lambda: None)
    mol = _mol("c1ccccc1C(=O)O")
    ref = _het_pdb_block(mol)
    poses = _poses_sdf(tmp_path, [mol, mol])
    assert rmsd.rmsds_to_reference(ref, poses, 2) == [None, None]
