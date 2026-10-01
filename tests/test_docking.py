"""Tests for backend.docking: command building and pose parsing."""

from __future__ import annotations

from pathlib import Path

import pytest

from backend import docking
from backend import pharmacophore as pharm



def test_parse_poses_reads_affinity_property(tmp_path):
    sdf = tmp_path / "poses.sdf"
    sdf.write_text(
        "mol\n     RDKit          3D\n\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n"
        ">  <minimizedAffinity>\n"
        "-7.5\n\n"
        "$$$$\n"
    )
    poses = docking.parse_poses(sdf)
    assert len(poses) == 1
    assert poses[0].affinity == -7.5
    assert poses[0].rank == 1


def _co_pose(o_xyz: tuple[float, float, float], affinity: float) -> str:
    """One C-O 'pose' molblock; O gets an implicit H on sanitize, so it's a
    donor/acceptor candidate atom without needing an explicit H in the file."""
    ox, oy, oz = o_xyz
    return (
        "mol\n     RDKit          3D\n\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        f"{0.0:10.4f}{0.0:10.4f}{0.0:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        f"{ox:10.4f}{oy:10.4f}{oz:10.4f} O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0\n"
        "M  END\n"
        ">  <minimizedAffinity>\n"
        f"{affinity}\n\n"
        "$$$$\n"
    )


def test_parse_poses_reranks_constraint_satisfying_pose_first(tmp_path):
    sdf = tmp_path / "poses.sdf"
    # Pose 1: best affinity, but its O sits far from the constraint point.
    # Pose 2: worse affinity, but its O lands right on the constraint point.
    sdf.write_text(_co_pose((100.0, 100.0, 100.0), -7.5) + _co_pose((5.0, 0.0, 0.0), -6.0))

    constraint = pharm.Constraint(kind="hbond_donor", point=(5.0, 0.0, 0.0), radius=1.0)
    poses = docking.parse_poses(sdf, constraints=[constraint])

    assert [p.affinity for p in poses] == [-6.0, -7.5]
    assert [p.rank for p in poses] == [1, 2]
    assert poses[0].constraint_satisfied is True
    assert poses[1].constraint_satisfied is False
    assert poses[0].constraint_detail[0]["distance"] == pytest.approx(0.0, abs=1e-6)


def test_parse_poses_without_constraints_keeps_original_order(tmp_path):
    sdf = tmp_path / "poses.sdf"
    sdf.write_text(_co_pose((100.0, 100.0, 100.0), -7.5) + _co_pose((5.0, 0.0, 0.0), -6.0))
    poses = docking.parse_poses(sdf)
    assert [p.affinity for p in poses] == [-7.5, -6.0]
    assert all(p.constraint_satisfied is None for p in poses)


def test_dock_raises_sminafailed_on_nonzero_exit(monkeypatch, tmp_path):
    monkeypatch.setattr(docking.config, "find_smina", lambda: "/usr/bin/smina")

    class FakeProc:
        returncode = 1
        stdout = ""
        stderr = "smina: fatal error, box outside receptor\n"

    monkeypatch.setattr(docking.subprocess, "run", lambda *a, **k: FakeProc())
    with pytest.raises(docking.SminaFailed):
        docking.dock(
            receptor=tmp_path / "r.pdb",
            ligand_sdf=tmp_path / "l.sdf",
            out_sdf=tmp_path / "out.sdf",
            center=(0, 0, 0),
            size=(10, 10, 10),
        )
