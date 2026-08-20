"""Tests for backend.docking: command building and pose parsing."""

from __future__ import annotations

from pathlib import Path

import pytest

from backend import docking


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
