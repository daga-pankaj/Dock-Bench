"""Tests for backend.docking: command building and pose parsing."""

from __future__ import annotations

from pathlib import Path

import pytest

from backend import docking
from backend import pharmacophore as pharm


def test_build_command_includes_all_flags(monkeypatch):
    monkeypatch.setattr(docking.config, "find_smina", lambda: "/usr/bin/smina")
    cmd = docking.build_command(
        receptor=Path("r.pdb"),
        ligand=Path("l.sdf"),
        out_sdf=Path("out.sdf"),
        center=(1.0, 2.0, 3.0),
        size=(20.0, 20.0, 20.0),
        exhaustiveness=16,
        num_modes=9,
        energy_range=3.0,
        seed=7,
        scoring="vinardo",
        cpu=2,
    )
    assert cmd[0] == "/usr/bin/smina"
    assert "--center_x" in cmd and cmd[cmd.index("--center_x") + 1] == "1.000"
    assert "--exhaustiveness" in cmd and cmd[cmd.index("--exhaustiveness") + 1] == "16"
    assert "--scoring" in cmd and cmd[cmd.index("--scoring") + 1] == "vinardo"
    assert "--cpu" in cmd and cmd[cmd.index("--cpu") + 1] == "2"
    assert "--quiet" in cmd


def test_build_command_raises_when_smina_missing(monkeypatch):
    monkeypatch.setattr(docking.config, "find_smina", lambda: None)
    with pytest.raises(docking.SminaNotFound):
        docking.build_command(Path("r.pdb"), Path("l.sdf"), Path("o.sdf"), (0, 0, 0), (10, 10, 10))


def test_build_command_rejects_unknown_scoring(monkeypatch):
    monkeypatch.setattr(docking.config, "find_smina", lambda: "/usr/bin/smina")
    with pytest.raises(ValueError):
        docking.build_command(
            Path("r.pdb"), Path("l.sdf"), Path("o.sdf"), (0, 0, 0), (10, 10, 10), scoring="nope"
        )


def test_build_command_uses_gnina_when_selected(monkeypatch):
    monkeypatch.setattr(docking.config, "find_gnina", lambda: "/usr/bin/gnina")
    cmd = docking.build_command(
        Path("r.pdb"), Path("l.sdf"), Path("o.sdf"), (0, 0, 0), (10, 10, 10), engine="gnina"
    )
    assert cmd[0] == "/usr/bin/gnina"
    assert "--no_gpu" not in cmd


def test_build_command_gnina_cpu_only_adds_no_gpu_flag(monkeypatch):
    monkeypatch.setattr(docking.config, "find_gnina", lambda: "/usr/bin/gnina")
    cmd = docking.build_command(
        Path("r.pdb"), Path("l.sdf"), Path("o.sdf"), (0, 0, 0), (10, 10, 10),
        engine="gnina", use_gpu=False,
    )
    assert "--no_gpu" in cmd


def test_build_command_smina_ignores_use_gpu(monkeypatch):
    monkeypatch.setattr(docking.config, "find_smina", lambda: "/usr/bin/smina")
    cmd = docking.build_command(
        Path("r.pdb"), Path("l.sdf"), Path("o.sdf"), (0, 0, 0), (10, 10, 10), use_gpu=False
    )
    assert "--no_gpu" not in cmd


def test_build_command_raises_when_gnina_missing(monkeypatch):
    monkeypatch.setattr(docking.config, "find_gnina", lambda: None)
    with pytest.raises(docking.GninaNotFound):
        docking.build_command(
            Path("r.pdb"), Path("l.sdf"), Path("o.sdf"), (0, 0, 0), (10, 10, 10), engine="gnina"
        )


def test_build_command_rejects_unknown_engine(monkeypatch):
    monkeypatch.setattr(docking.config, "find_smina", lambda: "/usr/bin/smina")
    with pytest.raises(ValueError):
        docking.build_command(
            Path("r.pdb"), Path("l.sdf"), Path("o.sdf"), (0, 0, 0), (10, 10, 10), engine="nope"
        )


def test_parse_poses_missing_file_returns_empty():
    assert docking.parse_poses(Path("/nonexistent/path.sdf")) == []


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
