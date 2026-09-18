"""Tests for backend.strain: ligand-only MMFF strain energy for a docked pose."""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem

from backend import strain


def _unoptimized_mol(smiles: str, seed: int = 7) -> Chem.Mol:
    """A deliberately un-relaxed embedding - not a local minimum, so it has
    somewhere to strain-relax down from."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    AllChem.EmbedMolecule(mol, params)
    return mol


def test_pose_energy_returns_a_float():
    mol = _unoptimized_mol("CCCCC(=O)Oc1ccccc1C(=O)O")  # aspirin
    energy = strain.pose_energy(mol)
    assert isinstance(energy, float)


def test_reference_minimum_energy_returns_a_float():
    mol = _unoptimized_mol("CCCCC(=O)Oc1ccccc1C(=O)O")
    energy = strain.reference_minimum_energy(mol, n_conformers=5)
    assert isinstance(energy, float)


def test_strain_energies_is_positive_for_an_unrelaxed_pose():
    mol = _unoptimized_mol("CCCCC(=O)Oc1ccccc1C(=O)O")
    values = strain.strain_energies([mol])
    assert len(values) == 1
    assert values[0] is not None
    assert values[0] > 0


def test_strain_energies_shares_one_reference_across_poses_of_the_same_ligand():
    mols = [_unoptimized_mol("CCCCC(=O)Oc1ccccc1C(=O)O", seed=s) for s in (1, 2, 3)]
    values = strain.strain_energies(mols)
    assert len(values) == 3
    assert all(v is not None and v >= 0 for v in values)


def test_strain_energies_empty_input_returns_empty_list():
    assert strain.strain_energies([]) == []
