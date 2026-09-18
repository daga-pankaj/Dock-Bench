"""Tests for backend.interaction_render: SVG diagram rendering (no ProLIF import here)."""

from __future__ import annotations

import xml.etree.ElementTree as ET

from rdkit import Chem
from rdkit.Chem import AllChem

from backend import interaction_render as render
from backend.interactions_prolif import Interaction


def _mol() -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1O"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    return mol


def test_render_placeholder_svg_is_valid_xml():
    svg = render.render_placeholder_svg("no interactions detected")
    ET.fromstring(svg)
    assert "no interactions detected" in svg


def test_render_svg_with_no_interactions_returns_placeholder():
    svg = render.render_svg(_mol(), [])
    ET.fromstring(svg)
    assert "no interactions detected" in svg


def test_render_svg_produces_valid_svg_with_nodes_and_legend():
    mol = _mol()
    hits = [
        Interaction("hbond", "ASP", 25, "A", 2.9, ligand_atom_idx=0),
        Interaction("hydrophobic", "PHE", 88, "A", 3.7, ligand_atom_idx=1),
        Interaction("pistacking", "TRP", 90, "A", 4.0),  # no ligand_atom_idx -> centroid fallback
    ]
    svg = render.render_svg(mol, hits)
    root = ET.fromstring(svg)
    assert root.tag.endswith("svg")
    assert "ASP25A" in svg
    assert "PHE88A" in svg
    assert "TRP90A" in svg
    assert "H-bond" in svg
    assert "hydrophobic" in svg
    assert "pi-stacking" in svg


def test_render_svg_escapes_malicious_residue_name():
    mol = _mol()
    malicious = Interaction("hbond", '"><script>alert(1)</script>', 1, "A", 3.0)
    svg = render.render_svg(mol, [malicious])
    assert "<script>" not in svg
    assert "&lt;script&gt;" in svg
    ET.fromstring(svg)  # must still be well-formed XML
