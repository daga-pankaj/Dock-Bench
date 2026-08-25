"""Tests for backend.prep: receptor cleanup, HET-ligand detection, ligand parsing."""

from __future__ import annotations

from backend import prep

RECEPTOR_PDB = """\
ATOM      1 N    ALA A   1      11.104  13.207  10.000                       N
ATOM      2 CA   ALA A   1      12.560  13.207  10.000                       C
HETATM    3 O    HOH A 101      20.000  20.000  20.000                       O
HETATM    4 ZN    ZN A 102      15.000  15.000  15.000                      ZN
HETATM    5 C1   LIG A 200       9.000   9.000   9.000                       C
HETATM    6 C2   LIG A 200      11.000  11.000  11.000                       C
HETATM    7 C3   LIG A 200      10.000  12.000  10.000                       C
HETATM    8 C4   LIG A 200      10.000  10.000  12.000                       C
HETATM    9 NA    NA A 300       5.000   5.000   5.000                      NA
END
"""


def test_list_het_ligands_excludes_water_and_orders_by_size():
    hets = prep.list_het_ligands(RECEPTOR_PDB)
    keys = [h.key for h in hets]
    assert "HOH:A:101" not in keys
    # LIG has 4 atoms and isn't an ion -> ranked ahead of the ions.
    assert keys[0] == "LIG:A:200"
    lig = next(h for h in hets if h.key == "LIG:A:200")
    assert lig.n_atoms == 4
    assert lig.is_ion is False
    assert lig.center == [10.0, 10.5, 10.5]


def test_list_het_ligands_keeps_ions_even_below_min_atoms():
    hets = prep.list_het_ligands(RECEPTOR_PDB)
    keys = {h.key: h for h in hets}
    assert "ZN:A:102" in keys
    assert keys["ZN:A:102"].is_ion is True
    assert "NA:A:300" in keys


def test_clean_receptor_drops_waters_by_default():
    cleaned = prep.clean_receptor(RECEPTOR_PDB)
    assert "HOH" not in cleaned
    assert "ATOM" in cleaned
    assert cleaned.strip().endswith("END")


def test_clean_receptor_keep_waters_true_keeps_them():
    cleaned = prep.clean_receptor(RECEPTOR_PDB, keep_waters=True)
    assert "HOH" in cleaned


def test_clean_receptor_drops_het_not_in_keep_list():
    cleaned = prep.clean_receptor(RECEPTOR_PDB, keep_het=())
    assert "LIG" not in cleaned


def test_clean_receptor_keeps_requested_het():
    cleaned = prep.clean_receptor(RECEPTOR_PDB, keep_het=("LIG",))
    assert "LIG" in cleaned


def test_clean_receptor_keeps_metals_by_default():
    cleaned = prep.clean_receptor(RECEPTOR_PDB)
    assert " ZN " in cleaned or "ZN A 102" in cleaned


TWO_CHAIN_PDB = RECEPTOR_PDB.replace("END\n", "") + """\
ATOM     10 N    ALA B   1     111.000 113.000 110.000                       N
ATOM     11 CA   ALA B   1     112.000 113.000 110.000                       C
END
"""


def test_clean_receptor_keeps_all_chains_without_site_center():
    cleaned = prep.clean_receptor(TWO_CHAIN_PDB)
    assert "ALA B" in cleaned


def test_clean_receptor_drops_chains_far_from_site_center():
    # LIG:A:200 sits around (10, 10.5, 10.5); chain B is ~170 A away.
    cleaned = prep.clean_receptor(TWO_CHAIN_PDB, keep_het=("LIG",), site_center=(10.0, 10.5, 10.5))
    assert "ALA A" in cleaned
    assert "LIG" in cleaned
    assert "ALA B" not in cleaned


def test_chains_near_point_finds_only_close_chains():
    near = prep.chains_near_point(TWO_CHAIN_PDB, center=(10.0, 10.5, 10.5), cutoff=15.0)
    assert near == {"A"}


def test_parse_smiles_block_names_and_dedupes():
    text = "CCO ethanol\nCCO ethanol\nc1ccccc1\n"
    specs = prep.parse_smiles_block(text)
    assert [s.name for s in specs] == ["ethanol", "ethanol_2", "lig_3"]
    assert specs[0].smiles == "CCO"


def test_parse_smiles_block_skips_blank_and_comment_lines():
    text = "# a comment\n\nCCO first\n"
    specs = prep.parse_smiles_block(text)
    assert len(specs) == 1
    assert specs[0].name == "first"


def test_parse_ligand_file_dispatches_on_csv_extension():
    text = "name,smiles\naspirin,CC(=O)Oc1ccccc1C(=O)O\n"
    specs = prep.parse_ligand_file("ligands.csv", text)
    assert len(specs) == 1
    assert specs[0].name == "aspirin"
    assert specs[0].smiles == "CC(=O)Oc1ccccc1C(=O)O"


def test_parse_ligand_file_defaults_to_smiles_for_unknown_extension():
    specs = prep.parse_ligand_file("ligands.txt", "CCO x\n")
    assert specs[0].name == "x"


def test_build_3d_from_smiles_produces_molblock_with_conformer():
    spec = prep.LigandSpec(name="ethanol", smiles="CCO")
    block = prep.build_3d(spec, seed=1)
    assert "V2000" in block or "V3000" in block
    assert block.count("\n") > 3


def test_build_3d_rejects_bad_smiles():
    spec = prep.LigandSpec(name="bad", smiles="not-a-smiles(((")
    try:
        prep.build_3d(spec)
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_build_3d_rejects_empty_spec():
    spec = prep.LigandSpec(name="empty")
    try:
        prep.build_3d(spec)
        assert False, "expected ValueError"
    except ValueError:
        pass
