"""Tests for the MOE Ligand Interactions diagram plumbing in backend.moe_dock."""

from __future__ import annotations

from pathlib import Path

from backend import moe_dock


def test_dock_script_draws_one_diagram_per_pose_next_to_the_poses():
    script = moe_dock.build_svl_script(
        Path("/r.moe"), Path("/l.sdf"), Path("/j/poses/sita_poses.mdb"),
        Path("/j/poses/sita_poses.sdf"), (0, 0, 0), (20, 20, 20),
    )
    assert "local function WriteLigandInteractions" in script
    assert "twrite ['/j/poses/sita_poses_lid_{}.svg', lid_k]" in script
    assert script.index("db_ExportSD") < script.index("WriteLigandInteractions [lid_lig")


def test_interactions_svg_path_matches_what_the_script_writes():
    assert moe_dock.interactions_svg_path(Path("/j/poses/sita_poses.sdf"), 3) == Path("/j/poses/sita_poses_lid_3.svg")


def test_browser_svg_swaps_font_and_adds_white_background():
    svg = '<?xml version="1.0"?><svg width="10" height="10"><text font-family="sf">N</text></svg>'
    out = moe_dock.browser_svg(svg)
    assert 'font-family="sf"' not in out
    assert '<svg width="10" height="10"><rect width="100%" height="100%" fill="#ffffff"/>' in out
