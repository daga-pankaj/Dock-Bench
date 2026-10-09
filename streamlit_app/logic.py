"""UI decisions for the Streamlit app, kept free of Streamlit itself.

Every function takes `state` - st.session_state in the app, a plain dict in
tests - and reads/writes the keys documented in DEFAULTS. Ported from the
original web page's frontend/src/app.js so both front ends behave the same.
"""

from __future__ import annotations

import math
from collections.abc import MutableMapping

State = MutableMapping

# Docking-parameter defaults. Peptides/peptidomimetics are floppier and need
# a much wider conformer search than a typical small-molecule ligand.
STANDARD_DOCKING_PARAMS = {"exhaustiveness": 8, "energy_range": 4.0}
PEPTIDE_DOCKING_PARAMS = {"exhaustiveness": 64, "energy_range": 10.0}

NUCLEOPHILE_ATOMS = {"CYS": "SG", "SER": "OG", "THR": "OG1", "TYR": "OH", "LYS": "NZ"}
WATER_RESNAMES = {"HOH", "WAT", "DOD"}
HBOND_KINDS = ("hbond_donor", "hbond_acceptor")

CONSTRAINT_KIND_LABELS = {
    "hbond_donor": "H-bond donor",
    "hbond_acceptor": "H-bond acceptor",
    "positional": "positional (any atom)",
    "aromatic": "aromatic ring (pi-stacking)",
    "metal_coordination": "metal coordination",
    "smarts": "SMARTS match",
}

# Typical contact distance per kind, shown as a hint next to the radius
# field - not enforced, just a starting point for the ones that aren't 3.5A.
CONSTRAINT_RADIUS_HINT = {
    "hbond_donor": "typical 2.5-3.5 Å",
    "hbond_acceptor": "typical 2.5-3.5 Å",
    "positional": "",
    "aromatic": "ring-centroid distance; typical 4-5.5 Å",
    "metal_coordination": "metal-ligand bond; typical 2.0-2.6 Å",
    "smarts": "",
}

MOE_PLACEMENTS = ["Triangle Matcher", "Alpha PMI", "Pharmacophore", "Proxy Triangle", "Alpha Triangle"]
MOE_REFINEMENTS = ["Rigid Receptor", "Induced Fit", "Rigid Body", "None"]

ENGINE_INFO = {
    "smina": "CPU-based Vina fork - good for general docking",
    "gnina": "CNN scoring + GPU support - better pose prediction",
    "moe": "Commercial engine - supports covalent docking, pharmacophore constraints",
}

# Session keys and their starting values. Widget-backed keys (box fields,
# docking parameters, checkboxes) live here too so code can set them before
# their widget is drawn.
DEFAULTS = {
    "receptor": None,           # receptor meta from services.store_receptor
    "original_receptor": None,  # the pre-QuickPrep receptor, for "use original"
    "pdb_text": "",
    "qp_report": None,
    "qp_info": "",
    "active_het": None,         # het_ligands entry defining the site
    "picks": [],                # [x, y, z] click points in "pick atoms" mode
    "pick_mode": "click ligand",
    "constraints": [],          # {id, dx/dy/dz, nx/ny/nz}; editable fields live in c<id>_* widget keys
    "constraint_seq": 0,
    "pending_pick_id": None,    # constraint row waiting for a viewer click
    "covalent_pick_active": False,
    "ligands": [],
    "lig_preview": None,
    "ph4_id": None,
    "ph4_name": "",
    "custom_rxn_id": None,
    "custom_rxn_name": "",
    "job_id": None,
    "job_done_seen": None,      # job id whose completion toast was shown
    "selected_ligand": None,
    "active_pose": 0,
    "hide_ligand": False,
    "hide_protein": False,
    "focus": None,              # {nonce, kind, ...} - viewer zooms once per nonce
    "focus_nonce": 0,
    "last_event_nonce": None,
    "toasts": [],               # [(message, bad)] shown on the next render
    "residue_hits": None,       # (resname, [specs]) from "Find in receptor"
    # widgets
    "cx": None, "cy": None, "cz": None,
    "sx": 22.0, "sy": 22.0, "sz": 22.0,
    "padding": 5.0,
    "engine": "smina",
    "remove_ligand": True,
    "keep_waters": False,
    "trim_chains": False,
    "validate_rmsd": False,
    "exh": 8,
    "modes": 5,
    "erange": 4.0,
    "seed": 42,
    "scoring": "vina",
    "peptide_mode": False,
    "minimize_after": False,
    "device": "GPU if available",
    "moe_placement": MOE_PLACEMENTS[0],
    "moe_refinement": MOE_REFINEMENTS[0],
    "covalent_enabled": False,
    "covalent_reaction": "",
    "cov_residue": "",
    "diagram_style": "MOE",
    "show_diagram": True,
    "hide_moe_scores": False,
}


def init_state(state: State) -> None:
    for key, value in DEFAULTS.items():
        if key not in state:
            state[key] = list(value) if isinstance(value, list) else value


def toast(state: State, message: str, bad: bool = False) -> None:
    state["toasts"] = [*state["toasts"], (message, bad)]


def _finite(*values) -> bool:
    return all(isinstance(v, (int, float)) and math.isfinite(v) for v in values)


# ---------------------------------------------------------------------------
# docking box
# ---------------------------------------------------------------------------


def box_center(state: State) -> list[float] | None:
    c = [state["cx"], state["cy"], state["cz"]]
    return [float(v) for v in c] if _finite(*c) else None


def box_size(state: State) -> list[float]:
    return [state["sx"], state["sy"], state["sz"]]


def set_box(state: State, center, size) -> None:
    for key, v in zip(("cx", "cy", "cz"), center):
        state[key] = round(float(v), 2)
    for key, v in zip(("sx", "sy", "sz"), size):
        state[key] = round(max(8.0, float(v)), 1)


def make_cubic(state: State) -> None:
    sizes = [v for v in box_size(state) if _finite(v)]
    if sizes:
        edge = round(max(sizes), 1)
        state["sx"] = state["sy"] = state["sz"] = edge


def box_from_points(points, pad: float):
    mins = [min(p[i] for p in points) for i in range(3)]
    maxs = [max(p[i] for p in points) for i in range(3)]
    center = [(mins[i] + maxs[i]) / 2 for i in range(3)]
    size = [maxs[i] - mins[i] + 2 * pad for i in range(3)]
    return center, size


def site_info(state: State) -> str:
    s = box_size(state)
    if box_center(state) is None or not _finite(*s) or min(s) <= 0:
        return ""
    return f"box {' x '.join(f'{v:.1f}' for v in s)} Å  ({s[0] * s[1] * s[2]:.0f} Å³)"


def request_focus(state: State, focus: dict) -> None:
    state["focus_nonce"] += 1
    state["focus"] = {"nonce": state["focus_nonce"], **focus}


def het_ref(het: dict) -> dict:
    return {"resname": het["resname"], "chain": het["chain"], "resseq": het["resseq"]}


def use_het(state: State, het: dict) -> None:
    """Box the site around a bound ligand (padding added on every side)."""
    state["active_het"] = het
    state["picks"] = []
    pad = state["padding"] or 0
    set_box(state, het["center"], [v + 2 * pad for v in het["size"]])
    request_focus(state, {"kind": "het", "het": het_ref(het)})
    toast(state, f"site set from {het['resname']} {het['chain']}{het['resseq']}")


def clear_picks(state: State) -> None:
    state["picks"] = []
    state["active_het"] = None


# ---------------------------------------------------------------------------
# receptor
# ---------------------------------------------------------------------------


def load_receptor(state: State, meta: dict, pdb_text: str) -> None:
    """A new receptor replaces the old one: site, picks, chains to delete,
    constraints and the pose display all start over; the biggest non-ion
    bound group (almost always the co-crystallised ligand) becomes the site."""
    state["receptor"] = meta
    state["pdb_text"] = pdb_text
    state["selected_ligand"] = None
    state["picks"] = []
    state["active_het"] = None
    state["constraints"] = []
    state["pending_pick_id"] = None
    best = next((h for h in meta.get("het_ligands") or [] if not h["is_ion"]), None)
    if best:
        use_het(state, best)
    if not meta.get("prepared_from"):
        state["original_receptor"] = None
        state["qp_report"] = None
        state["qp_info"] = ""


def swap_receptor(state: State, meta: dict, pdb_text: str) -> None:
    """Swap the loaded receptor for another version of itself (QuickPrep'd or
    original) without losing the docking site the user had set up."""
    site_key = (state["active_het"] or {}).get("key")
    center, size = box_center(state), box_size(state)
    load_receptor(state, meta, pdb_text)
    het = next((h for h in meta.get("het_ligands") or [] if h["key"] == site_key), None) if site_key else None
    if het:
        use_het(state, het)
    elif center:
        state["active_het"] = None
        set_box(state, center, size)


def chain_key(state: State, chain: str) -> str:
    """Widget key for a chain's "delete chain" box - per receptor, so a new
    receptor starts with every box unticked."""
    return f"chain_{state['receptor']['id']}_{chain}"


def excluded_chains(state: State) -> list[str]:
    if not state["receptor"]:
        return []
    return [c["chain"] for c in state["receptor"].get("chains") or [] if state.get(chain_key(state, c["chain"]))]


def quickprep_summary(report: dict) -> str:
    return (
        f"+{report['hydrogens_added']} H, {report['heavy_atoms_added']} heavy atoms built, "
        f"{report['waters_removed']} waters removed, {report['n_moved']} heavy atoms moved >0.5 Å "
        f"(RMSD {report['rmsd_heavy']} Å, max {report['max_shift']} Å)"
    )


def residue_focus(residue_label: str) -> dict:
    """Focus request for a QuickPrep protonation entry like "HIS A57"."""
    chain_res = residue_label.split(" ")[1]
    return {"kind": "residue", "chain": chain_res[0], "resseq": int(chain_res[1:])}


# ---------------------------------------------------------------------------
# constraints
# ---------------------------------------------------------------------------


def ckey(cid: int, field: str) -> str:
    return f"c{cid}_{field}"


def add_constraint(state: State) -> None:
    state["constraint_seq"] += 1
    cid = state["constraint_seq"]
    state["constraints"] = [*state["constraints"], {"id": cid}]
    for field, value in (("kind", "hbond_donor"), ("x", None), ("y", None), ("z", None),
                         ("radius", 3.5), ("required", True), ("label", ""), ("smarts", "")):
        state[ckey(cid, field)] = value


def remove_constraint(state: State, cid: int) -> None:
    state["constraints"] = [c for c in state["constraints"] if c["id"] != cid]
    if state["pending_pick_id"] == cid:
        state["pending_pick_id"] = None


def constraint_view(state: State, c: dict) -> dict:
    """One constraint with its widget-held fields merged in."""
    return {**c, **{f: state.get(ckey(c["id"], f)) for f in
                    ("kind", "x", "y", "z", "radius", "required", "label", "smarts")}}


def constraint_status(view: dict) -> str:
    has_point = _finite(view["x"], view["y"], view["z"])
    if view["kind"] in HBOND_KINDS:
        if has_point and _finite(view.get("dx"), view.get("dy"), view.get("dz")):
            return "✓ directional"
        return "needs direction" if has_point else ""
    if view["kind"] == "aromatic" and _finite(view.get("nx"), view.get("ny"), view.get("nz")):
        return "✓ planar"
    return ""


def constraint_markers(state: State) -> list[dict]:
    """Placed constraints, in the shape the viewer component draws."""
    out = []
    for c in state["constraints"]:
        v = constraint_view(state, c)
        if not _finite(v["x"], v["y"], v["z"]):
            continue
        out.append({k: v.get(k) for k in ("id", "kind", "x", "y", "z", "nx", "ny", "nz", "dx", "dy", "dz")})
    return out


def toggle_constraint_pick(state: State, cid: int) -> None:
    state["pending_pick_id"] = None if state["pending_pick_id"] == cid else cid
    if state["pending_pick_id"] is not None:
        toast(state, "click an atom in the viewer to set this constraint's point")


def constraints_payload(state: State) -> tuple[list[dict] | None, str | None]:
    """Job-submission constraints, or (None, reason) if a row is incomplete -
    such a constraint can't be checked against anything."""
    views = [constraint_view(state, c) for c in state["constraints"]]
    if any(not _finite(v["x"], v["y"], v["z"]) for v in views):
        return None, "a pharmacophore constraint is missing its point - pick one or remove the row"
    if any(v["kind"] == "smarts" and not (v["smarts"] or "").strip() for v in views):
        return None, "a SMARTS constraint is missing its pattern"
    payload = []
    for v in views:
        item = {"kind": v["kind"], "point": [v["x"], v["y"], v["z"]], "radius": v["radius"],
                "required": v["required"], "label": v["label"] or ""}
        if v["kind"] == "smarts":
            item["smarts"] = v["smarts"].strip()
        payload.append(item)
    return payload, None


# ---------------------------------------------------------------------------
# viewer clicks
# ---------------------------------------------------------------------------


def _covalent_click(state: State, atom: dict) -> None:
    state["covalent_pick_active"] = False
    resname = (atom.get("resn") or "").upper()
    expected = NUCLEOPHILE_ATOMS.get(resname)
    if not expected:
        toast(state, f"⚠ {resname} not a typical nucleophile. Supported: CYS, SER, THR, TYR, LYS", True)
        return
    if atom.get("atom") != expected:
        toast(state, f"⚠ Picked {atom.get('atom')}, but {resname} nucleophile is typically {expected}. "
                     "Proceeding anyway.", True)
    spec = f"{atom.get('chain') or ''}:{resname}:{atom.get('resi')}"
    state["cov_residue"] = spec
    toast(state, f"✓ Covalent residue set: {spec} ({atom.get('atom')})")


def _constraint_click(state: State, atom: dict) -> None:
    c = next((c for c in state["constraints"] if c["id"] == state["pending_pick_id"]), None)
    if c is None:
        state["pending_pick_id"] = None
        return
    cid = c["id"]
    kind = state[ckey(cid, "kind")]
    needs_direction = kind in HBOND_KINDS
    has_point = _finite(state[ckey(cid, "x")], state[ckey(cid, "y")], state[ckey(cid, "z")])

    # H-bond: two-click mode - first click = point, second click = direction
    if needs_direction and has_point:
        d = [atom["x"] - state[ckey(cid, "x")], atom["y"] - state[ckey(cid, "y")], atom["z"] - state[ckey(cid, "z")]]
        mag = math.sqrt(sum(v * v for v in d))
        if mag < 0.1:
            toast(state, "direction atom too close to constraint point - pick a different atom", True)
            return
        c["dx"], c["dy"], c["dz"] = (v / mag for v in d)
        state["pending_pick_id"] = None
        toast(state, f"direction set from {atom['resn']}{atom['resi']} {atom['atom']}")
        return

    # Aromatic: snap to the ring centroid and keep its plane normal
    ring = atom.get("ring")
    if kind == "aromatic" and ring:
        state[ckey(cid, "x")], state[ckey(cid, "y")], state[ckey(cid, "z")] = ring["centroid"]
        if ring.get("normal"):
            c["nx"], c["ny"], c["nz"] = ring["normal"]
        state["pending_pick_id"] = None
        toast(state, f"aromatic plane set from {atom['resn']}{atom['resi']} ring")
        return

    state[ckey(cid, "x")], state[ckey(cid, "y")], state[ckey(cid, "z")] = atom["x"], atom["y"], atom["z"]
    if needs_direction:
        # stays armed for the second (direction) click
        toast(state, f"point set - now click a second atom to define "
                     f"{'donor' if kind == 'hbond_donor' else 'acceptor'} direction")
    else:
        state["pending_pick_id"] = None
        toast(state, f"constraint point set from {atom['resn']}{atom['resi']} {atom['atom']}")


def handle_click(state: State, atom: dict) -> None:
    """What a click on an atom in the viewer means, in priority order:
    covalent-residue pick, then a constraint waiting for its point, then
    the docking site (a bound ligand in "click ligand" mode, else a pick
    point - three or four picks define a site well)."""
    if not atom:
        return
    if state["covalent_pick_active"]:
        _covalent_click(state, atom)
        return
    if state["pending_pick_id"] is not None:
        _constraint_click(state, atom)
        return

    resn = (atom.get("resn") or "").upper()
    if state["pick_mode"] == "click ligand" and atom.get("hetflag") and resn not in WATER_RESNAMES and state["receptor"]:
        key = f"{resn}:{atom.get('chain') or '_'}:{atom.get('resi')}"
        het = next((h for h in state["receptor"].get("het_ligands") or [] if h["key"] == key), None)
        if het:
            use_het(state, het)
            return

    state["picks"] = [*state["picks"], [atom["x"], atom["y"], atom["z"]]]
    state["active_het"] = None
    center, size = box_from_points(state["picks"], state["padding"] or 0)
    min_edge = 20 if len(state["picks"]) == 1 else 14
    set_box(state, center, [max(v, min_edge) for v in size])
    toast(state, f"{len(state['picks'])} pick(s) - {atom.get('resn')}{atom.get('resi')} {atom.get('atom')}")


def handle_event(state: State, event: dict | None) -> None:
    """Apply a viewer component event once (events carry a unique nonce)."""
    if not event or event.get("nonce") == state["last_event_nonce"]:
        return
    state["last_event_nonce"] = event.get("nonce")
    if event.get("type") == "click":
        handle_click(state, event.get("atom") or {})
    elif event.get("type") == "pose":
        state["active_pose"] = int(event.get("index", 0))


# ---------------------------------------------------------------------------
# engines, ligands, jobs
# ---------------------------------------------------------------------------


def engine_options(status: dict) -> dict[str, str]:
    """engine -> label. gnina is always listed (marked when missing), smina
    and MOE only when installed - same as the web page's engine picker."""
    options = {}
    if status.get("smina"):
        options["smina"] = "smina (CPU)"
    options["gnina"] = "gnina (CPU/GPU)" if status.get("gnina") else "gnina (not installed)"
    if status.get("moe_available"):
        options["moe"] = "MOE"
    return options


def scoring_options(status: dict, engine: str) -> list[str]:
    key = "moe_scoring_functions" if engine == "moe" else "scoring_functions"
    return list(status.get(key) or [])


def on_engine_change(state: State, status: dict) -> None:
    """Scoring choices differ per engine; covalent docking is MOE-only."""
    options = scoring_options(status, state["engine"])
    if state["engine"] != "moe":
        state["covalent_enabled"] = False
        state["scoring"] = "vina" if "vina" in options else (options[0] if options else "vina")
    elif state["scoring"] not in options and options:
        state["scoring"] = options[0]


def apply_peptide_mode(state: State) -> None:
    params = PEPTIDE_DOCKING_PARAMS if state["peptide_mode"] else STANDARD_DOCKING_PARAMS
    state["exh"] = params["exhaustiveness"]
    state["erange"] = params["energy_range"]


def show_ligands(state: State, data: dict) -> None:
    state["ligands"] = [l for l in data["ligands"] if l["valid"]]
    state["lig_preview"] = data


def ligand_info(data: dict | None) -> str:
    if not data:
        return ""
    bad = len(data["ligands"]) - data["n_valid"]
    return f"{data['n_valid']} ligand(s) ready" + (f", {bad} rejected" if bad else "")


def build_job_request(state: State) -> tuple[dict | None, str | None]:
    """The POST /api/jobs body for the current setup, or (None, reason)."""
    if not state["receptor"]:
        return None, "load a receptor first"
    if not state["ligands"]:
        return None, "no ligands loaded - click 'Load from box'"
    center = box_center(state)
    if not center:
        return None, "define the active site first"
    constraints, err = constraints_payload(state)
    if err:
        return None, err
    if state["validate_rmsd"] and not state["active_het"]:
        return None, "RMSD validation needs a bound ligand selected above as the site"
    return {
        "receptor_id": state["receptor"]["id"],
        "center": center,
        "size": box_size(state),
        "ligands": [{"name": l["name"], "smiles": l["smiles"], "molblock": l["molblock"]} for l in state["ligands"]],
        "exhaustiveness": int(state["exh"]),
        "num_modes": int(state["modes"]),
        "energy_range": float(state["erange"]),
        "seed": int(state["seed"]),
        "scoring": state["scoring"],
        "engine": state["engine"],
        "use_gpu": state["device"] != "CPU only",
        "minimize_after_docking": state["minimize_after"],
        "keep_waters": state["keep_waters"],
        "remove_ligand": state["remove_ligand"],
        "exclude_chains": excluded_chains(state),
        "trim_distant_chains": state["trim_chains"],
        "constraints": constraints,
        "reference_ligand_key": state["active_het"]["key"] if state["validate_rmsd"] else None,
        "moe_placement": state["moe_placement"],
        "moe_refinement": state["moe_refinement"],
        "ph4_id": state["ph4_id"],
        "covalent_enabled": state["covalent_enabled"],
        "covalent_reaction": state["covalent_reaction"] or None,
        "covalent_residue": (state["cov_residue"] or "").strip() or None,
        "custom_rxn_id": state["custom_rxn_id"],
    }, None


def sorted_results(ligands: list[dict]) -> list[dict]:
    """Best score first; unscored ligands last, by name."""
    scored = sorted((l for l in ligands if l["best_score"] is not None), key=lambda l: l["best_score"])
    unscored = sorted((l for l in ligands if l["best_score"] is None), key=lambda l: l["name"])
    return scored + unscored


def rmsd_mismatch(job: dict) -> bool:
    """RMSD validation was requested but no docked pose matched the bound
    ligand's heavy-atom graph - obrms reports that silently as "inf"."""
    done = [l for l in job["ligands"] if l["status"] == "done"]
    return bool(job["params"].get("reference_ligand_pdb")) and bool(done) and all(
        not l["rmsd_to_reference"] or all(r is None for r in l["rmsd_to_reference"]) for l in done
    )


def constraint_text(lig: dict) -> str:
    if lig.get("constraint_satisfied") is None:
        return ""
    return "✓ met" if lig["constraint_satisfied"] else "✗ unmet"


def constraint_detail_text(lig: dict) -> str:
    """Each constraint's distance on the top pose, one per line."""
    lines = []
    for d in lig.get("constraint_detail") or []:
        dist = "no matching atom" if d["distance"] is None else f"{d['distance']:.2f} Å"
        lines.append(f"{d['label']}: {dist}{'' if d['required'] else ' (optional)'}")
    return "\n".join(lines)
