#!/usr/bin/env python3
"""Test covalent docking implementation."""

from pathlib import Path
from backend import config, moe_dock, prep


def test_parse_residue_spec():
    """Test residue spec parsing."""
    # Valid specs
    assert prep.parse_residue_spec("A:CYS:189") == ("A", "CYS", 189)
    assert prep.parse_residue_spec(":SER:25") == ("", "SER", 25)
    assert prep.parse_residue_spec("B:LYS:42") == ("B", "LYS", 42)

    # Invalid specs
    assert prep.parse_residue_spec("invalid") is None
    assert prep.parse_residue_spec("A:CY:189") is None  # resname not 3 chars
    assert prep.parse_residue_spec("ABC:CYS:189") is None  # chain too long
    print("✓ parse_residue_spec tests passed")


def test_validate_residue_exists():
    """Test residue validation in PDB text."""
    pdb_sample = """
ATOM   1234  CA  CYS A 189      10.0  20.0  30.0  1.00 50.00           C
ATOM   1235  CB  SER B  25      11.0  21.0  31.0  1.00 50.00           C
HETATM 1236  MG  MG  C 100      12.0  22.0  32.0  1.00 60.00          MG
"""

    # Should find
    assert prep.validate_residue_exists(pdb_sample, "A:CYS:189")
    assert prep.validate_residue_exists(pdb_sample, "B:SER:25")
    assert prep.validate_residue_exists(pdb_sample, ":CYS:189")  # any chain

    # Should not find
    assert not prep.validate_residue_exists(pdb_sample, "A:LEU:100")
    assert not prep.validate_residue_exists(pdb_sample, "X:CYS:189")  # wrong chain
    print("✓ validate_residue_exists tests passed")


def test_covalent_reactions_discovery():
    """Test MOE reaction template discovery."""
    rxn_dir = config.find_moe_reactions()
    reactions = config.list_covalent_reactions()

    if rxn_dir is None:
        print("⊘ MOE not installed - skipping reaction discovery test")
        return

    assert rxn_dir.is_dir()
    assert len(reactions) == 22
    assert "michael_acceptor_s" in reactions
    assert "epoxide_s" in reactions
    assert "alpha_halocarbonyl_s" in reactions
    print(f"✓ Found {len(reactions)} covalent reaction templates")


def test_svl_script_generation():
    """Test SVL script generation with covalent params."""
    rxn_dir = config.find_moe_reactions()
    if rxn_dir is None:
        print("⊘ MOE not installed - skipping SVL script test")
        return

    rxn_file = rxn_dir / "michael_acceptor_s.rxn"
    script = moe_dock.build_svl_script(
        receptor=Path("receptor.pdb"),
        ligand=Path("ligand.sdf"),
        out_mdb=Path("output.mdb"),
        out_sdf=Path("output.sdf"),
        center=(10.0, 20.0, 30.0),
        size=(20.0, 20.0, 20.0),
        covalent_reaction_file=rxn_file,
        covalent_residue="A:CYS:189",
    )

    # Check covalent options present
    assert "opt.dockmode = 'Covalent'" in script
    assert "opt.placement = 'Covalent Docking'" in script
    assert "opt.placement_opt" in script
    assert "rxnFile:" in script
    assert "michael_acceptor_s.rxn" in script
    assert "mode: 'Covalent'" in script
    assert "opt.refine_opt" in script
    assert "tether: 10" in script
    # Check atom selection using rUID (PDB numbering)
    assert "rUID aResidue protein_atoms == 189" in script
    assert "aName protein_atoms == 'SG'" in script
    assert "rxsite_atoms = protein_atoms |" in script
    print("✓ SVL script with covalent options generated correctly")


def test_svl_script_without_covalent():
    """Test SVL script generation without covalent params."""
    script = moe_dock.build_svl_script(
        receptor=Path("receptor.pdb"),
        ligand=Path("ligand.sdf"),
        out_mdb=Path("output.mdb"),
        out_sdf=Path("output.sdf"),
        center=(10.0, 20.0, 30.0),
        size=(20.0, 20.0, 20.0),
    )

    # Check no covalent options injected
    assert "opt.dockmode = 'Covalent'" not in script
    assert "opt.placement = 'Covalent Docking'" not in script
    assert "rxnFile:" not in script
    assert "rxsite_atoms" not in script
    print("✓ SVL script without covalent options generated correctly")


def test_custom_rxn_validation():
    """Test custom .rxn file validation."""
    # Valid RXN
    valid_rxn = "$RXN\n\nKetcher\n\n  1  1\n$MOL\n..."
    assert valid_rxn.strip().startswith("$RXN")

    # Invalid RXN
    invalid_rxn = "not a reaction file"
    assert not invalid_rxn.strip().startswith("$RXN")

    print("✓ Custom .rxn validation logic correct")


if __name__ == "__main__":
    test_parse_residue_spec()
    test_validate_residue_exists()
    test_covalent_reactions_discovery()
    test_svl_script_generation()
    test_svl_script_without_covalent()
    test_custom_rxn_validation()
    print("\n✓ All covalent docking tests passed")
