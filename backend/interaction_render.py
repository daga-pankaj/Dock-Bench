"""Render a LigPlot-style flattened 2D interaction schematic as SVG.

Takes only an RDKit mol and the interaction list from interactions_prolif.py -
no ProLIF import here at all, so detection and rendering stay swappable
independently (e.g. pointing this at a real LigPlot+ render, or a different
detection engine, later only touches this file).

Layout: the ligand's own 2D depiction sits in the center; each residue that
contacts it gets one labeled node placed on a circle around it, connected by
a dashed line colored/styled by interaction type (this mirrors classic
LigPlot diagrams without reproducing their exact geometry).
"""

from __future__ import annotations

import html
import math
import re

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from .interactions_prolif import Interaction

WIDTH = 720
HEIGHT = 720
LIGAND_BOX = 360  # square the ligand depiction is drawn into, centered on the canvas
RESIDUE_RADIUS = 300  # circle radius the residue nodes sit on, around the canvas center
NODE_RADIUS = 20

INTERACTION_STYLE = {
    "hbond": {"color": "#37c98b", "dash": "6,3", "label": "H-bond"},
    "hydrophobic": {"color": "#8b97a6", "dash": "3,3", "label": "hydrophobic"},
    "pistacking": {"color": "#4da3ff", "dash": "1,0", "label": "pi-stacking"},
    "pication": {"color": "#4da3ff", "dash": "4,2", "label": "pi-cation"},
    "saltbridge": {"color": "#f0b429", "dash": "6,3", "label": "salt bridge"},
    "halogen": {"color": "#ff8bd6", "dash": "6,3", "label": "halogen bond"},
    "waterbridge": {"color": "#4da3ff", "dash": "2,4", "label": "water bridge"},
    "metal": {"color": "#ff6b6b", "dash": "1,0", "label": "metal contact"},
}
DEFAULT_STYLE = {"color": "#e6ebf1", "dash": "4,2", "label": "contact"}

_SVG_INNER_RE = re.compile(r"<svg[^>]*>(.*)</svg>", re.DOTALL)


def render_placeholder_svg(message: str) -> str:
    cx, cy = WIDTH / 2, HEIGHT / 2
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}">'
        f'<rect width="100%" height="100%" fill="#171c23"/>'
        f'<text x="{cx}" y="{cy}" text-anchor="middle" fill="#8b97a6" '
        f'font-family="sans-serif" font-size="16">{html.escape(message)}</text>'
        "</svg>"
    )


def _ligand_svg_inner(svg_text: str) -> str:
    match = _SVG_INNER_RE.search(svg_text)
    return match.group(1) if match else svg_text


def render_svg(mol: Chem.Mol, interactions: list[Interaction]) -> str:
    if not interactions:
        return render_placeholder_svg("no interactions detected for this pose")

    # Draw on a copy so computing 2D coords never disturbs the 3D pose conformer
    # the caller may still be using elsewhere. Atoms are NOT stripped of
    # explicit Hs here: `ligand_atom_idx` on each Interaction was computed
    # against this exact atom ordering (see interactions_prolif.py), and
    # removing/reordering atoms would silently break that correspondence.
    depict_mol = Chem.Mol(mol)
    AllChem.Compute2DCoords(depict_mol)

    off_x = (WIDTH - LIGAND_BOX) / 2
    off_y = (HEIGHT - LIGAND_BOX) / 2
    drawer = rdMolDraw2D.MolDraw2DSVG(LIGAND_BOX, LIGAND_BOX)
    drawer.DrawMolecule(depict_mol)
    drawer.FinishDrawing()
    ligand_svg = drawer.GetDrawingText()

    anchor: dict[int, tuple[float, float]] = {}
    for atom in depict_mol.GetAtoms():
        idx = atom.GetIdx()
        point = drawer.GetDrawCoords(idx)
        anchor[idx] = (point.x + off_x, point.y + off_y)
    ligand_centroid = (
        sum(x for x, _ in anchor.values()) / len(anchor),
        sum(y for _, y in anchor.values()) / len(anchor),
    )

    by_residue: dict[tuple[str, int, str], list[Interaction]] = {}
    for hit in interactions:
        by_residue.setdefault((hit.restype, hit.resnr, hit.reschain), []).append(hit)
    residues = sorted(by_residue.items(), key=lambda kv: (kv[0][2], kv[0][1]))

    cx, cy = WIDTH / 2, HEIGHT / 2
    n = len(residues)
    node_pos: dict[tuple[str, int, str], tuple[float, float]] = {}
    for i, (key, _) in enumerate(residues):
        angle = 2 * math.pi * i / n - math.pi / 2
        node_pos[key] = (cx + RESIDUE_RADIUS * math.cos(angle), cy + RESIDUE_RADIUS * math.sin(angle))

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="#171c23"/>',
    ]

    kinds_present: set[str] = set()
    for key, hits in residues:
        rx, ry = node_pos[key]
        for hit in hits:
            kinds_present.add(hit.kind)
            style = INTERACTION_STYLE.get(hit.kind, DEFAULT_STYLE)
            start = anchor.get(hit.ligand_atom_idx, ligand_centroid) if hit.ligand_atom_idx is not None else ligand_centroid
            parts.append(
                f'<line x1="{start[0]:.1f}" y1="{start[1]:.1f}" x2="{rx:.1f}" y2="{ry:.1f}" '
                f'stroke="{style["color"]}" stroke-width="1.6" stroke-dasharray="{style["dash"]}"/>'
            )

    parts.append(f'<g transform="translate({off_x:.1f},{off_y:.1f})">{_ligand_svg_inner(ligand_svg)}</g>')

    for key, hits in residues:
        restype, resnr, reschain = key
        rx, ry = node_pos[key]
        style = INTERACTION_STYLE.get(hits[0].kind, DEFAULT_STYLE)
        parts.append(
            f'<circle cx="{rx:.1f}" cy="{ry:.1f}" r="{NODE_RADIUS}" fill="#1e252e" '
            f'stroke="{style["color"]}" stroke-width="2"/>'
        )
        # restype/reschain ultimately trace back to ATOM/HETATM fields of a
        # user-uploaded receptor PDB (POST /api/receptor/upload) - escape
        # before landing in this SVG's <text> content.
        label = html.escape(f"{restype}{resnr}{reschain}")
        parts.append(
            f'<text x="{rx:.1f}" y="{ry + 4:.1f}" text-anchor="middle" fill="#e6ebf1" '
            f'font-family="sans-serif" font-size="11">{label}</text>'
        )

    legend_y = HEIGHT - 16 * len(kinds_present) - 10
    for kind in sorted(kinds_present):
        style = INTERACTION_STYLE.get(kind, DEFAULT_STYLE)
        parts.append(
            f'<line x1="16" y1="{legend_y}" x2="44" y2="{legend_y}" stroke="{style["color"]}" '
            f'stroke-width="2" stroke-dasharray="{style["dash"]}"/>'
        )
        parts.append(
            f'<text x="50" y="{legend_y + 4}" fill="#8b97a6" font-family="sans-serif" '
            f'font-size="11">{html.escape(style["label"])}</text>'
        )
        legend_y += 16

    parts.append("</svg>")
    return "".join(parts)
