"""Tests for backend.covalent_check and the MOE score plumbing around it."""

from pathlib import Path

from rdkit import Chem

from backend import covalent_check as cc
from backend import docking, moe_dock

# Ser + carboxylic acid -> ester, as drawn by MOESketch (note the off-by-one
# "M  ALS" columns, which covalent_check has to normalise).
ESTER_RXN = """$RXN
ester formation O
      MOESketch100520261500
acylation carboxylic acid serine ester
  2  1
$MOL

MOESketch           2D

  5  4  0  0  1  0            999 V2000
   22.9159    4.8820    0.0000 O   0  0  0  0  0  0  0  0  0  5  0  0
   21.6163    4.1329    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0
   21.6163    2.6232    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
   20.3172    4.8828    0.0000 L   0  0  0  0  0  0  0  0  0  3  0  0
   21.6163    5.6330    0.0000 H   0  0  0  0  0  0  0  0  0  4  0  0
  2  1  1  0  0  0  0
  3  2  1  0  0  0  0
  2  4  1  1  0  0  0
  2  5  1  6  0  0  0
M  SUB  2   1  -2   3   3
M  ALS   4  2 F    H   C
M  END
$MOL

MOESketch           2D

  4  3  0  0  1  0            999 V2000
   28.2721    4.9561    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
   29.2020    5.8858    0.0000 O   0  0  0  0  0  0  0  0  0  9  0  0
   29.5000    4.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   26.9572    4.9561    0.0000 C   0  0  0  0  0  0  0  0  0  7  0  0
  1  2  2  0  0  0  0
  1  3  1  0  0  0  4
  1  4  1  0  0  0  0
M  SUB  1   3  -2
M  END
$MOL

MOESketch           2D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"""


def _molblock(smiles: str) -> str:
    return Chem.MolToMolBlock(Chem.AddHs(Chem.MolFromSmiles(smiles)))


def test_reactant_count_and_warning(tmp_path: Path):
    rxn = tmp_path / "ester.rxn"
    rxn.write_text(ESTER_RXN)
    query = cc.ligand_reactant_query(rxn)
    assert query is not None

    def count(smi):
        return cc.count_reactive_sites(Chem.AddHs(Chem.MolFromSmiles(smi)), query)

    assert count("CCC(=O)O") == 1
    assert count("CCC(=O)[O-]") == 1
    assert count("CCC=O") == 0
    assert count("OC(=O)CCC(=O)O") == 2

    assert cc.reactant_warning(_molblock("CCC(=O)O"), rxn) is None
    assert "no group" in cc.reactant_warning(_molblock("CCC=O"), rxn)
    assert "2 groups" in cc.reactant_warning(_molblock("OC(=O)CCC(=O)O"), rxn)
    assert cc.reactant_warning(_molblock("CCC=O"), None) is None


def test_unreadable_rxn_skips_check(tmp_path: Path):
    rxn = tmp_path / "bad.rxn"
    rxn.write_text("$RXN\nnot really\n")
    assert cc.ligand_reactant_query(rxn) is None
    assert cc.reactant_warning(_molblock("CCC(=O)O"), rxn) is None


LINKED_PDB = """\
LINK         OG  SER A  70                 C1  ALD A 301     1555   1555  1.40
ATOM      1  OG  SER A  70       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  C1  ALD A 301       1.400   0.000   0.000  1.00 10.00           C
HETATM    3  O1  ALD A 301       2.000   1.050   0.000  1.00 10.00           O
HETATM    4  C2  ALD A 301       2.150  -1.300   0.000  1.00 10.00           C
HETATM    5  C3  ALD A 301       3.650  -1.300   0.000  1.00 10.00           C
END
"""


def test_acid_form_from_link_rebuilds_acyl_adduct():
    assert cc.acid_form_from_link(LINKED_PDB, "ALD:A:301", "CCC=O") == "CCC(=O)O"


def test_acid_form_from_link_leaves_unlinked_or_non_acyl_alone():
    no_link = "\n".join(l for l in LINKED_PDB.splitlines() if not l.startswith("LINK"))
    assert cc.acid_form_from_link(no_link, "ALD:A:301", "CCC=O") is None
    # linked through a carbon that's already an acid - nothing to rebuild
    assert cc.acid_form_from_link(LINKED_PDB, "ALD:A:301", "CCC(=O)O") is None
    assert cc.acid_form_from_link(LINKED_PDB, "bad-key", "CCC=O") is None


def _covalent_pose() -> Chem.Mol:
    # residue fragment (CB, OG) first, then ligand acyl C(=O)CC
    mol = Chem.MolFromSmiles("COC(=O)CC")
    mol.SetProp("useAtoms", "0 0 1 2 3 4")
    return mol


def test_pose_attachment():
    assert cc.pose_attachment(_covalent_pose()) == "*C(=O)CC"
    plain = Chem.MolFromSmiles("CCO")
    assert cc.pose_attachment(plain) is None
    plain.SetProp("useAtoms", "1 2 3")
    assert cc.pose_attachment(plain) is None


def test_parse_moe_sdf_keeps_every_moe_score(tmp_path: Path):
    mol = _covalent_pose()
    fields = {"S": "-7.9", "E_place": "-7.1", "E_score1": "-13.4", "E_refine": "inf",
              "E_score2": "-7.9", "E_conf": "-142.2", "rmsd_refine": "1.08"}
    for k, v in fields.items():
        mol.SetProp(k, v)
    sdf = tmp_path / "poses.sdf"
    with Chem.SDWriter(str(sdf)) as w:
        w.write(mol)

    (pose,) = docking._parse_moe_sdf(sdf)
    assert pose.affinity == -7.9
    # inf (a failed refinement) is dropped rather than breaking JSON
    assert pose.moe_scores == {"S": -7.9, "E_place": -7.1, "E_score1": -13.4,
                               "E_score2": -7.9, "E_conf": -142.2, "rmsd_refine": 1.08}
    assert pose.covalent_attachment == "*C(=O)CC"


def test_svl_uses_moe_dock_panel_defaults(tmp_path: Path):
    script = moe_dock.build_svl_script(
        tmp_path / "rec.moe", tmp_path / "lig.sdf", tmp_path / "o.mdb", tmp_path / "o.sdf",
        (0, 0, 0), (20, 20, 20), num_poses=5,
    )
    assert "scoring: 'London dG'" in script
    assert "rescoring: 'GBVI/WSA dG'" in script
    assert "maxpose: 30," in script
    assert "remaxpose: 5," in script
    assert docking.MOE_SCORING_FUNCTIONS[0] == "GBVI/WSA dG"
