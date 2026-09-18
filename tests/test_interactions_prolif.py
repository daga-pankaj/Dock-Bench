"""Tests for backend.interactions_prolif: ProLIF-based interaction detection.

The receptor/ligand fixture below is a small, deterministic excerpt (7
residues within contact range of one ligand) taken from a real PDB entry
(1EVE, an acetylcholinesterase/E20 complex) - real enough for ProLIF's own
chemistry perception to behave normally, but compact enough to run fast and
stay self-contained (no network fetch at test time).
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from backend import interactions_prolif as ip

RECEPTOR_PDB = """\
ATOM    636  N   TRP A  84       0.030  73.938  64.893  1.00 18.27           N
ATOM    637  CA  TRP A  84      -0.120  72.702  64.142  1.00 18.35           C
ATOM    638  C   TRP A  84      -1.546  72.166  64.252  1.00 19.37           C
ATOM    639  O   TRP A  84      -1.942  71.281  63.494  1.00 23.73           O
ATOM    640  CB  TRP A  84       0.881  71.657  64.643  1.00 20.76           C
ATOM    641  CG  TRP A  84       2.313  72.117  64.597  1.00 21.50           C
ATOM    642  CD1 TRP A  84       3.105  72.450  65.665  1.00 22.33           C
ATOM    643  CD2 TRP A  84       3.131  72.271  63.431  1.00 20.27           C
ATOM    644  NE1 TRP A  84       4.364  72.797  65.231  1.00 18.82           N
ATOM    645  CE2 TRP A  84       4.410  72.696  63.866  1.00 21.26           C
ATOM    646  CE3 TRP A  84       2.910  72.091  62.057  1.00 15.77           C
ATOM    647  CZ2 TRP A  84       5.461  72.942  62.976  1.00 16.84           C
ATOM    648  CZ3 TRP A  84       3.958  72.339  61.173  1.00 15.19           C
ATOM    649  CH2 TRP A  84       5.216  72.759  61.638  1.00 19.00           C
ATOM    941  N   TYR A 121      -2.241  60.928  61.351  1.00 18.85           N
ATOM    942  CA  TYR A 121      -3.362  60.985  62.293  1.00 18.91           C
ATOM    943  C   TYR A 121      -3.864  62.423  62.443  1.00 19.80           C
ATOM    944  O   TYR A 121      -4.841  62.682  63.150  1.00 18.02           O
ATOM    945  CB  TYR A 121      -2.948  60.436  63.671  1.00 16.47           C
ATOM    946  CG  TYR A 121      -2.050  61.359  64.470  1.00 16.74           C
ATOM    947  CD1 TYR A 121      -2.589  62.341  65.307  1.00 18.93           C
ATOM    948  CD2 TYR A 121      -0.664  61.277  64.363  1.00 18.12           C
ATOM    949  CE1 TYR A 121      -1.764  63.226  66.013  1.00 12.92           C
ATOM    950  CE2 TYR A 121       0.169  62.150  65.061  1.00 11.04           C
ATOM    951  CZ  TYR A 121      -0.384  63.123  65.883  1.00 19.23           C
ATOM    952  OH  TYR A 121       0.451  63.985  66.564  1.00 16.41           O
ATOM   2154  N   TRP A 279      -5.471  56.010  70.172  1.00 25.18           N
ATOM   2155  CA  TRP A 279      -4.158  55.604  70.665  1.00 20.31           C
ATOM   2156  C   TRP A 279      -4.113  54.264  71.397  1.00 23.46           C
ATOM   2157  O   TRP A 279      -3.050  53.654  71.520  1.00 21.20           O
ATOM   2158  CB  TRP A 279      -3.592  56.707  71.555  1.00 15.41           C
ATOM   2159  CG  TRP A 279      -3.619  58.045  70.890  1.00 19.11           C
ATOM   2160  CD1 TRP A 279      -4.591  58.997  71.008  1.00 16.99           C
ATOM   2161  CD2 TRP A 279      -2.648  58.572  69.976  1.00 14.91           C
ATOM   2162  NE1 TRP A 279      -4.287  60.084  70.224  1.00 19.81           N
ATOM   2163  CE2 TRP A 279      -3.101  59.852  69.578  1.00 16.43           C
ATOM   2164  CE3 TRP A 279      -1.440  58.089  69.452  1.00 18.23           C
ATOM   2165  CZ2 TRP A 279      -2.389  60.657  68.679  1.00 19.10           C
ATOM   2166  CZ3 TRP A 279      -0.729  58.890  68.556  1.00 19.14           C
ATOM   2167  CH2 TRP A 279      -1.209  60.160  68.180  1.00 21.45           C
ATOM   2231  N   PHE A 288       6.813  58.252  68.754  1.00 17.67           N
ATOM   2232  CA  PHE A 288       6.953  58.359  67.302  1.00 17.28           C
ATOM   2233  C   PHE A 288       6.297  57.171  66.592  1.00 17.22           C
ATOM   2234  O   PHE A 288       6.844  56.607  65.640  1.00 18.34           O
ATOM   2235  CB  PHE A 288       8.427  58.495  66.897  1.00 16.03           C
ATOM   2236  CG  PHE A 288       8.655  59.399  65.706  1.00 19.36           C
ATOM   2237  CD1 PHE A 288       7.732  59.465  64.665  1.00 17.46           C
ATOM   2238  CD2 PHE A 288       9.797  60.187  65.632  1.00 14.85           C
ATOM   2239  CE1 PHE A 288       7.942  60.300  63.572  1.00 13.44           C
ATOM   2240  CE2 PHE A 288      10.015  61.024  64.540  1.00 21.65           C
ATOM   2241  CZ  PHE A 288       9.079  61.078  63.506  1.00 14.68           C
ATOM   2562  N   PHE A 330      10.840  66.722  69.637  1.00 23.81           N
ATOM   2563  CA  PHE A 330       9.385  66.885  69.715  1.00 20.58           C
ATOM   2564  C   PHE A 330       8.681  65.689  70.350  1.00 18.35           C
ATOM   2565  O   PHE A 330       7.762  65.847  71.153  1.00 21.68           O
ATOM   2566  CB  PHE A 330       8.809  67.097  68.307  1.00 23.95           C
ATOM   2567  CG  PHE A 330       7.876  68.275  68.193  1.00 25.29           C
ATOM   2568  CD1 PHE A 330       6.999  68.600  69.227  1.00 27.09           C
ATOM   2569  CD2 PHE A 330       7.884  69.069  67.050  1.00 20.97           C
ATOM   2570  CE1 PHE A 330       6.147  69.702  69.122  1.00 25.76           C
ATOM   2571  CE2 PHE A 330       7.038  70.174  66.935  1.00 20.51           C
ATOM   2572  CZ  PHE A 330       6.169  70.491  67.971  1.00 23.59           C
ATOM   2573  N   PHE A 331       9.149  64.496  70.002  1.00 18.28           N
ATOM   2574  CA  PHE A 331       8.565  63.248  70.475  1.00 17.88           C
ATOM   2575  C   PHE A 331       8.984  62.832  71.876  1.00 20.54           C
ATOM   2576  O   PHE A 331       8.240  62.136  72.568  1.00 29.60           O
ATOM   2577  CB  PHE A 331       8.844  62.153  69.446  1.00 13.39           C
ATOM   2578  CG  PHE A 331       8.438  62.546  68.058  1.00 16.82           C
ATOM   2579  CD1 PHE A 331       9.286  63.323  67.268  1.00  7.68           C
ATOM   2580  CD2 PHE A 331       7.171  62.236  67.579  1.00  9.84           C
ATOM   2581  CE1 PHE A 331       8.875  63.792  66.028  1.00 11.86           C
ATOM   2582  CE2 PHE A 331       6.750  62.699  66.341  1.00 12.22           C
ATOM   2583  CZ  PHE A 331       7.603  63.482  65.563  1.00 17.31           C
ATOM   2600  N   TYR A 334       6.864  64.813  74.364  1.00 21.52           N
ATOM   2601  CA  TYR A 334       5.496  64.337  74.512  1.00 20.21           C
ATOM   2602  C   TYR A 334       5.287  62.959  75.132  1.00 24.56           C
ATOM   2603  O   TYR A 334       4.178  62.654  75.569  1.00 35.64           O
ATOM   2604  CB  TYR A 334       4.764  64.404  73.169  1.00 18.50           C
ATOM   2605  CG  TYR A 334       4.228  65.774  72.835  1.00 17.65           C
ATOM   2606  CD1 TYR A 334       5.090  66.827  72.529  1.00 19.70           C
ATOM   2607  CD2 TYR A 334       2.858  66.023  72.837  1.00 13.12           C
ATOM   2608  CE1 TYR A 334       4.600  68.094  72.235  1.00 10.68           C
ATOM   2609  CE2 TYR A 334       2.358  67.286  72.544  1.00 16.55           C
ATOM   2610  CZ  TYR A 334       3.231  68.318  72.244  1.00 15.14           C
ATOM   2611  OH  TYR A 334       2.729  69.571  71.958  1.00 12.97           O
"""

# The E20 ligand (donepezil) HETATM block from the same structure, at its
# original crystallographic pose relative to RECEPTOR_PDB above.
LIGAND_PDB = """\
HETATM 4326  C1  E20 A2001       1.057  60.056  71.240  1.00 26.85           C
HETATM 4327  C2  E20 A2001       0.188  61.079  71.628  1.00 22.90           C
HETATM 4328  C3  E20 A2001       0.347  62.415  71.171  1.00 21.66           C
HETATM 4329  C4  E20 A2001       1.447  62.720  70.270  1.00 15.86           C
HETATM 4330  C5  E20 A2001       2.306  61.665  69.881  1.00 18.21           C
HETATM 4331  C6  E20 A2001       2.140  60.371  70.362  1.00 24.60           C
HETATM 4332  C7  E20 A2001       3.365  62.094  69.016  1.00 25.58           C
HETATM 4333  C8  E20 A2001       3.086  63.590  68.740  1.00 21.53           C
HETATM 4334  C9  E20 A2001       1.839  63.992  69.558  1.00 17.21           C
HETATM 4335  C10 E20 A2001       4.105  64.576  69.245  1.00 18.02           C
HETATM 4336  C11 E20 A2001       4.307  65.848  68.533  1.00 18.84           C
HETATM 4337  C12 E20 A2001       5.006  65.637  67.195  1.00 15.13           C
HETATM 4338  C13 E20 A2001       5.161  66.919  66.425  1.00 19.67           C
HETATM 4339  N14 E20 A2001       3.854  67.635  66.237  1.00 18.93           N
HETATM 4340  C15 E20 A2001       3.177  67.997  67.535  1.00 13.92           C
HETATM 4341  C16 E20 A2001       3.004  66.708  68.390  1.00 13.21           C
HETATM 4342  C17 E20 A2001       3.743  68.943  65.483  1.00 20.75           C
HETATM 4343  C18 E20 A2001       3.992  68.847  63.988  1.00 17.17           C
HETATM 4344  C19 E20 A2001       2.864  68.675  63.126  1.00 18.68           C
HETATM 4345  C20 E20 A2001       3.052  68.432  61.723  1.00 21.92           C
HETATM 4346  C21 E20 A2001       4.345  68.390  61.195  1.00 14.62           C
HETATM 4347  C22 E20 A2001       5.479  68.598  62.046  1.00 19.28           C
HETATM 4348  C23 E20 A2001       5.310  68.869  63.440  1.00 13.50           C
HETATM 4349  O24 E20 A2001       4.325  61.375  68.615  1.00 28.32           O
HETATM 4350  O25 E20 A2001       0.890  58.762  71.644  1.00 29.40           O
HETATM 4351  C26 E20 A2001       1.842  57.647  71.265  1.00 21.15           C
HETATM 4352  O27 E20 A2001      -0.925  61.025  72.321  1.00 22.73           O
HETATM 4353  C28 E20 A2001      -1.434  59.862  72.928  1.00 31.57           C
"""


def _ligand_mol() -> Chem.Mol:
    mol = Chem.MolFromPDBBlock(LIGAND_PDB + "END\n", sanitize=False, removeHs=False)
    Chem.SanitizeMol(mol, catchErrors=True)
    return mol


def test_detect_interactions_finds_known_contacts():
    mol = _ligand_mol()
    hits = ip.detect_interactions(RECEPTOR_PDB, mol)

    kinds = {h.kind for h in hits}
    assert "hydrophobic" in kinds

    for hit in hits:
        assert hit.distance > 0
        assert hit.reschain == "A"
        if hit.ligand_atom_idx is not None:
            assert 0 <= hit.ligand_atom_idx < mol.GetNumAtoms()


def test_detect_interactions_returns_empty_when_nothing_is_nearby():
    mol = _ligand_mol()
    far_residue = (
        "ATOM      1  N   ALA A   1     500.000 500.000 500.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1     501.500 500.000 500.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1     502.000 501.500 500.000  1.00  0.00           C\n"
        "ATOM      4  O   ALA A   1     501.500 502.500 500.000  1.00  0.00           O\n"
        "ATOM      5  CB  ALA A   1     502.200 499.000 501.000  1.00  0.00           C\n"
    )
    hits = ip.detect_interactions(far_residue, mol)
    assert hits == []


def test_extract_pose_out_of_range_raises(tmp_path):
    sdf_path = tmp_path / "poses.sdf"
    block = Chem.MolToMolBlock(_ligand_mol())
    sdf_path.write_text(block + "$$$$\n")
    assert ip.extract_pose(sdf_path, 1) is not None
    with pytest.raises(ValueError):
        ip.extract_pose(sdf_path, 2)


def test_extract_pose_empty_file_raises(tmp_path):
    sdf_path = tmp_path / "empty.sdf"
    sdf_path.write_text("")
    with pytest.raises(ValueError):
        ip.extract_pose(sdf_path, 1)


def test_extract_pose_rejects_zero_index(tmp_path):
    sdf_path = tmp_path / "poses.sdf"
    sdf_path.write_text("")
    with pytest.raises(ValueError):
        ip.extract_pose(sdf_path, 0)
