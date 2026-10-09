"""Tests for streamlit_app.logic: the Streamlit UI's decisions, ported from frontend/src/app.js."""

from __future__ import annotations

import pytest

from streamlit_app import logic

HET = {"key": "STI:A:201", "resname": "STI", "chain": "A", "resseq": "201", "n_atoms": 37,
       "center": [10.0, 20.0, 30.0], "size": [8.0, 6.0, 4.0], "is_ion": False}
ION = {"key": "ZN:A:301", "resname": "ZN", "chain": "A", "resseq": "301", "n_atoms": 1,
       "center": [0.0, 0.0, 0.0], "size": [0.0, 0.0, 0.0], "is_ion": True}
RECEPTOR = {"id": "0123456789abcdef", "name": "x.pdb", "n_protein_atoms": 100,
            "het_ligands": [HET, ION], "chains": [{"chain": "A"}, {"chain": "B"}], "prepared_from": None}


@pytest.fixture
def state():
    s = {}
    logic.init_state(s)
    return s


def atom(**kw):
    base = {"resn": "LEU", "hetflag": False, "chain": "A", "resi": 10, "atom": "CA", "x": 1.0, "y": 2.0, "z": 3.0,
            "ring": None}
    return {**base, **kw}


def test_load_receptor_boxes_biggest_non_ion_ligand(state):
    logic.load_receptor(state, RECEPTOR, "PDB")
    assert state["active_het"] == HET
    # padding 5 on every side, min edge 8
    assert logic.box_center(state) == [10.0, 20.0, 30.0]
    assert logic.box_size(state) == [18.0, 16.0, 14.0]
    assert state["focus"]["kind"] == "het"


def test_set_box_enforces_minimum_edge(state):
    logic.set_box(state, [1.234, 2, 3], [2, 30, 7.96])
    assert logic.box_size(state) == [8.0, 30.0, 8.0]
    assert state["cx"] == 1.23


def test_load_receptor_resets_previous_setup(state):
    logic.load_receptor(state, RECEPTOR, "PDB")
    logic.add_constraint(state)
    state["picks"] = [[0, 0, 0]]
    logic.load_receptor(state, {**RECEPTOR, "het_ligands": []}, "PDB")
    assert state["constraints"] == [] and state["picks"] == [] and state["active_het"] is None


def test_click_on_bound_ligand_sets_site(state):
    logic.load_receptor(state, {**RECEPTOR, "het_ligands": [ION, HET]}, "PDB")
    logic.clear_picks(state)
    logic.handle_click(state, atom(resn="sti", hetflag=True, resi="201"))
    assert state["active_het"] == HET


def test_click_in_pick_mode_builds_box_from_points(state):
    logic.load_receptor(state, RECEPTOR, "PDB")
    state["pick_mode"] = "pick atoms"
    logic.handle_click(state, atom(x=0.0, y=0.0, z=0.0))
    assert state["active_het"] is None
    assert logic.box_size(state) == [20.0, 20.0, 20.0]  # one pick: 20 Å minimum
    logic.handle_click(state, atom(x=10.0, y=0.0, z=0.0))
    assert logic.box_center(state) == [5.0, 0.0, 0.0]
    assert logic.box_size(state) == [20.0, 14.0, 14.0]  # 10 + 2*5 padding; 14 Å minimum after that


def test_water_click_is_a_pick_not_a_site(state):
    logic.load_receptor(state, RECEPTOR, "PDB")
    logic.handle_click(state, atom(resn="HOH", hetflag=True))
    assert len(state["picks"]) == 1


def test_covalent_pick_fills_residue_and_warns_on_odd_atom(state):
    state["covalent_pick_active"] = True
    logic.handle_click(state, atom(resn="CYS", resi=189, atom="CB"))
    assert state["cov_residue"] == "A:CYS:189"
    assert state["covalent_pick_active"] is False
    assert any(bad for _, bad in state["toasts"])


def test_covalent_pick_rejects_non_nucleophile(state):
    state["covalent_pick_active"] = True
    logic.handle_click(state, atom(resn="ALA"))
    assert state["cov_residue"] == ""
    assert state["covalent_pick_active"] is False


def test_hbond_constraint_takes_point_then_direction(state):
    logic.add_constraint(state)
    cid = state["constraints"][0]["id"]
    logic.toggle_constraint_pick(state, cid)
    logic.handle_click(state, atom(x=0.0, y=0.0, z=0.0))
    assert state["pending_pick_id"] == cid  # still armed for the direction click
    assert logic.constraint_status(logic.constraint_view(state, state["constraints"][0])) == "needs direction"
    logic.handle_click(state, atom(x=0.0, y=0.0, z=2.0))
    c = state["constraints"][0]
    assert (c["dx"], c["dy"], c["dz"]) == (0.0, 0.0, 1.0)
    assert state["pending_pick_id"] is None


def test_aromatic_constraint_snaps_to_ring(state):
    logic.add_constraint(state)
    cid = state["constraints"][0]["id"]
    state[logic.ckey(cid, "kind")] = "aromatic"
    logic.toggle_constraint_pick(state, cid)
    logic.handle_click(state, atom(ring={"centroid": [4.0, 5.0, 6.0], "normal": [0.0, 1.0, 0.0]}))
    view = logic.constraint_view(state, state["constraints"][0])
    assert (view["x"], view["y"], view["z"]) == (4.0, 5.0, 6.0)
    assert logic.constraint_status(view) == "✓ planar"
    assert logic.constraint_markers(state)[0]["ny"] == 1.0


def test_events_apply_once_per_nonce(state):
    logic.load_receptor(state, RECEPTOR, "PDB")
    state["pick_mode"] = "pick atoms"
    event = {"nonce": 7, "type": "click", "atom": atom()}
    logic.handle_event(state, event)
    logic.handle_event(state, event)
    assert len(state["picks"]) == 1
    logic.handle_event(state, {"nonce": 8, "type": "pose", "index": 2})
    assert state["active_pose"] == 2


def test_constraints_payload_requires_points_and_smarts(state):
    logic.add_constraint(state)
    cid = state["constraints"][0]["id"]
    assert logic.constraints_payload(state)[0] is None
    state[logic.ckey(cid, "x")], state[logic.ckey(cid, "y")], state[logic.ckey(cid, "z")] = 1.0, 2.0, 3.0
    state[logic.ckey(cid, "kind")] = "smarts"
    assert "SMARTS" in logic.constraints_payload(state)[1]
    state[logic.ckey(cid, "smarts")] = " [OX1] "
    payload, err = logic.constraints_payload(state)
    assert err is None
    assert payload == [{"kind": "smarts", "point": [1.0, 2.0, 3.0], "radius": 3.5, "required": True,
                        "label": "", "smarts": "[OX1]"}]


def test_build_job_request(state):
    assert logic.build_job_request(state) == (None, "load a receptor first")
    logic.load_receptor(state, RECEPTOR, "PDB")
    assert logic.build_job_request(state)[1].startswith("no ligands")
    logic.show_ligands(state, {"ligands": [{"name": "a", "smiles": "CCO", "molblock": None, "valid": True, "note": ""},
                                           {"name": "b", "smiles": "C1", "molblock": None, "valid": False, "note": ""}],
                               "n_valid": 1})
    state[logic.chain_key(state, "B")] = True
    state["validate_rmsd"] = True
    state["device"] = "CPU only"
    payload, err = logic.build_job_request(state)
    assert err is None
    assert payload["ligands"] == [{"name": "a", "smiles": "CCO", "molblock": None}]
    assert payload["exclude_chains"] == ["B"]
    assert payload["reference_ligand_key"] == "STI:A:201"
    assert payload["use_gpu"] is False
    assert payload["covalent_reaction"] is None and payload["covalent_residue"] is None


def test_rmsd_validation_needs_bound_ligand_site(state):
    logic.load_receptor(state, {**RECEPTOR, "het_ligands": []}, "PDB")
    logic.set_box(state, [0, 0, 0], [20, 20, 20])
    logic.show_ligands(state, {"ligands": [{"name": "a", "smiles": "C", "molblock": None, "valid": True, "note": ""}],
                               "n_valid": 1})
    state["validate_rmsd"] = True
    assert "RMSD validation" in logic.build_job_request(state)[1]


def test_swap_receptor_keeps_site(state):
    logic.load_receptor(state, RECEPTOR, "PDB")
    prepared = {**RECEPTOR, "id": "fedcba9876543210", "prepared_from": RECEPTOR["id"]}
    logic.swap_receptor(state, prepared, "PDB2")
    assert state["active_het"]["key"] == HET["key"] and state["receptor"] is prepared
    # no matching ligand in the new version: the box itself is kept
    logic.set_box(state, [1, 2, 3], [30, 30, 30])
    logic.swap_receptor(state, {**prepared, "het_ligands": []}, "PDB3")
    assert logic.box_center(state) == [1.0, 2.0, 3.0] and logic.box_size(state) == [30.0, 30.0, 30.0]


def test_engine_options_and_change(state):
    status = {"smina": None, "gnina": None, "moe_available": True,
              "scoring_functions": ["vina", "vinardo"], "moe_scoring_functions": ["GBVI/WSA dG", "London dG"]}
    assert logic.engine_options(status) == {"gnina": "gnina (not installed)", "moe": "MOE"}
    state["engine"], state["scoring"] = "moe", "vina"
    logic.on_engine_change(state, status)
    assert state["scoring"] == "GBVI/WSA dG"
    state["covalent_enabled"] = True
    state["engine"] = "smina"
    logic.on_engine_change(state, status)
    assert state["scoring"] == "vina" and state["covalent_enabled"] is False


def test_peptide_mode_bumps_search(state):
    state["peptide_mode"] = True
    logic.apply_peptide_mode(state)
    assert (state["exh"], state["erange"]) == (64, 10.0)
    state["peptide_mode"] = False
    logic.apply_peptide_mode(state)
    assert (state["exh"], state["erange"]) == (8, 4.0)


def test_sorted_results_and_rmsd_mismatch():
    ligs = [{"name": "b", "best_score": None}, {"name": "a", "best_score": None},
            {"name": "c", "best_score": -5.0}, {"name": "d", "best_score": -7.0}]
    assert [l["name"] for l in logic.sorted_results(ligs)] == ["d", "c", "a", "b"]
    job = {"params": {"reference_ligand_pdb": "HETATM"},
           "ligands": [{"status": "done", "rmsd_to_reference": [None, None]}, {"status": "failed", "rmsd_to_reference": []}]}
    assert logic.rmsd_mismatch(job)
    job["ligands"][0]["rmsd_to_reference"] = [1.2, None]
    assert not logic.rmsd_mismatch(job)
