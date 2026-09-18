"""Tests for backend.pharmacophore: geometric constraint checking."""

from __future__ import annotations

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from backend import pharmacophore as pharm


def _phenol() -> Chem.Mol:
    """A small molecule with one clear H-bond donor (OH) and one acceptor-only ring N."""
    mol = Chem.MolFromSmiles("Oc1ccncc1")  # 3-hydroxypyridine-ish: phenol O-H + pyridine N
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=1)
    Chem.SanitizeMol(mol)
    return mol


def test_constraint_rejects_unknown_kind():
    with pytest.raises(ValueError):
        pharm.Constraint(kind="nope", point=(0, 0, 0))


def test_from_dict_roundtrip():
    c = pharm.Constraint.from_dict(
        {"kind": "positional", "point": [1.0, 2.0, 3.0], "radius": 2.0, "required": False}
    )
    assert c.kind == "positional"
    assert c.point == (1.0, 2.0, 3.0)
    assert c.required is False


def test_donor_atom_found_near_its_own_position():
    mol = _phenol()
    conf = mol.GetConformer()
    donor_idx = next(iter(pharm._candidate_atoms(mol, "hbond_donor")))
    p = conf.GetAtomPosition(donor_idx)
    c = pharm.Constraint(kind="hbond_donor", point=(p.x, p.y, p.z), radius=0.5)
    results = pharm.evaluate(mol, [c])
    assert results[0].satisfied
    assert results[0].distance == pytest.approx(0.0, abs=1e-6)


def test_constraint_far_away_is_unsatisfied_but_scored():
    mol = _phenol()
    c = pharm.Constraint(kind="hbond_donor", point=(100.0, 100.0, 100.0), radius=3.5, weight=1.0)
    results = pharm.evaluate(mol, [c])
    assert results[0].satisfied is False
    assert results[0].distance > 3.5
    assert pharm.penalty(results, [c]) > 0
    assert pharm.all_required_met(results) is False


def test_positional_ignores_hydrogens():
    mol = _phenol()
    idxs = pharm._candidate_atoms(mol, "positional")
    assert all(mol.GetAtomWithIdx(i).GetAtomicNum() > 1 for i in idxs)


def test_no_matching_atom_gives_none_distance():
    # A bare methane has no H-bond donor or acceptor atom at all.
    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    c = pharm.Constraint(kind="hbond_acceptor", point=(0, 0, 0))
    results = pharm.evaluate(mol, [c])
    assert results[0].distance is None
    assert results[0].satisfied is False
    assert pharm.penalty(results, [c]) == pytest.approx(c.weight * c.radius * 2)


def _benzene_ring() -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    return mol


def test_aromatic_matches_ring_centroid():
    mol = _benzene_ring()
    conf = mol.GetConformer()
    centroid = pharm._ring_centroids(mol, conf)[0]
    c = pharm.Constraint(kind="aromatic", point=centroid, radius=0.5)
    results = pharm.evaluate(mol, [c])
    assert results[0].satisfied
    assert results[0].distance == pytest.approx(0.0, abs=1e-6)


def test_aromatic_finds_no_ring_on_acyclic_molecule():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    c = pharm.Constraint(kind="aromatic", point=(0, 0, 0))
    results = pharm.evaluate(mol, [c])
    assert results[0].distance is None


def test_metal_coordination_reuses_acceptor_atoms():
    mol = _phenol()
    conf = mol.GetConformer()
    idx = next(iter(pharm._candidate_atoms(mol, "metal_coordination")))
    p = conf.GetAtomPosition(idx)
    c = pharm.Constraint(kind="metal_coordination", point=(p.x, p.y, p.z), radius=0.5)
    results = pharm.evaluate(mol, [c])
    assert results[0].satisfied


def test_smarts_constraint_matches_pattern():
    mol = _phenol()
    conf = mol.GetConformer()
    # Match the phenolic oxygen directly by pattern rather than the built-in kinds.
    c = pharm.Constraint(kind="smarts", smarts="[OX2H]", point=(0, 0, 0), radius=100)
    results = pharm.evaluate(mol, [c])
    assert results[0].distance is not None


def test_smarts_constraint_requires_pattern():
    with pytest.raises(ValueError):
        pharm.Constraint(kind="smarts", point=(0, 0, 0))


def test_smarts_constraint_rejects_invalid_pattern():
    with pytest.raises(ValueError):
        pharm.Constraint(kind="smarts", point=(0, 0, 0), smarts="not a smarts(")
