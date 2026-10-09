"""End-to-end run of the Streamlit app with streamlit.testing (no browser):
load a receptor, parse ligands, dock with smina, see the results."""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from backend import config, services

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = Path(__file__).resolve().parent.parent / "streamlit_app" / "app.py"

# A tiny receptor: one alanine-ish fragment plus a bound "LIG" group.
PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.988  -0.773  -1.199  1.00  0.00           C
HETATM    6  C1  LIG A 101       5.000   5.000   5.000  1.00  0.00           C
HETATM    7  C2  LIG A 101       6.500   5.000   5.000  1.00  0.00           C
HETATM    8  C3  LIG A 101       7.000   6.400   5.000  1.00  0.00           C
HETATM    9  O1  LIG A 101       8.400   6.400   5.000  1.00  0.00           O
END
"""


@pytest.fixture
def isolated_runs_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "RUNS_DIR", tmp_path)
    monkeypatch.setattr(config, "RECEPTOR_DIR", tmp_path / "receptors")
    monkeypatch.setattr(config, "JOB_DIR", tmp_path / "jobs")
    config.ensure_dirs()
    return tmp_path


def _app() -> AppTest:
    at = AppTest.from_file(str(APP), default_timeout=120)
    at.run()
    assert not at.exception, at.exception
    return at


def test_app_renders_without_receptor(isolated_runs_dir):
    at = _app()
    assert any(c.value == "no receptor loaded" for c in at.caption)
    labels = {b.label for b in at.button}
    assert {"Fetch", "Run", "Load from box", "make cubic"} <= labels


def test_run_without_receptor_toasts_reason(isolated_runs_dir):
    at = _app()
    next(b for b in at.button if b.label == "Run").click().run()
    assert any("load a receptor first" in t.value for t in at.toast)


@pytest.mark.skipif(not config.find_smina(), reason="smina not installed")
def test_dock_through_the_ui(isolated_runs_dir):
    at = _app()
    meta = services.store_receptor("tiny.pdb", PDB)
    # what the upload callback does, minus the browser file picker
    from streamlit_app import logic

    state = {}
    logic.init_state(state)
    logic.load_receptor(state, meta, PDB)
    for k in ("receptor", "pdb_text", "active_het", "cx", "cy", "cz", "sx", "sy", "sz"):
        at.session_state[k] = state[k]
    at.session_state["exh"] = 1
    at.session_state["modes"] = 2
    at.run()
    assert any("tiny.pdb - 5 protein atoms - 1 het groups" in c.value for c in at.caption)

    at.text_area(key="smiles_box").input("CCO ethanol\nnot_smiles bad")
    next(b for b in at.button if b.label == "Load from box").click().run()
    assert any(c.value == "1 ligand(s) ready, 1 rejected" for c in at.caption)

    next(b for b in at.button if b.label == "Run").click().run()
    job_id = at.session_state["job_id"]
    assert job_id
    deadline = time.time() + 120
    while services.get_job(job_id).status != "done" and time.time() < deadline:
        time.sleep(0.5)
    assert services.get_job(job_id).ligands[0].status == "done"

    at.run()
    assert not at.exception
    at.session_state["selected_ligand"] = "ethanol"
    at.run()
    assert not at.exception
    assert any("<img" in m.value and "interaction diagram" in m.value for m in at.markdown)
