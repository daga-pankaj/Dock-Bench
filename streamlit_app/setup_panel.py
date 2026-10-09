"""Left column: receptor, active site, constraints, covalent, ligands, parameters."""

from __future__ import annotations

import streamlit as st

from backend import services
from backend.schemas import JobRequest
from streamlit_app import logic


def _call(state, fn, *args):
    """Run a service call, turning a ServiceError into an error toast."""
    try:
        return fn(*args)
    except services.ServiceError as exc:
        logic.toast(state, exc.detail, True)
    except Exception as exc:  # e.g. pydantic validation of a job request
        logic.toast(state, str(exc), True)
    return None


def _uploader_key(state, name: str) -> str:
    """File uploaders are re-keyed after each upload so the same file can be
    picked again (the web page cleared its <input> the same way)."""
    return f"{name}_{state.get(name + '_n', 0)}"


def _consume_upload(state, name: str):
    file = state.get(_uploader_key(state, name))
    state[name + "_n"] = state.get(name + "_n", 0) + 1
    return file


# ---------------------------------------------------------------------------
# 1 · Receptor
# ---------------------------------------------------------------------------


def _load(state, meta: dict | None) -> None:
    if meta:
        logic.load_receptor(state, meta, services.receptor_pdb(meta["id"]))


def _fetch(state) -> None:
    pdb_id = state["pdb_id"].strip()
    if pdb_id:
        _load(state, _call(state, services.fetch_receptor, pdb_id))


def _upload_receptor(state) -> None:
    file = _consume_upload(state, "pdb_upload")
    if file is not None:
        _load(state, _call(state, services.upload_receptor, file.name, file.getvalue()))


def _quickprep(state) -> None:
    original = state["receptor"]
    het = state["active_het"]
    meta = _call(state, services.quickprep_receptor, original["id"], het["key"] if het else None)
    if not meta:
        return
    meta["name"] = f"{original['name']} (QuickPrep)"
    logic.swap_receptor(state, meta, services.receptor_pdb(meta["id"]))
    state["original_receptor"] = original
    state["qp_report"] = meta.get("quickprep_report")
    het = state["active_het"]
    state["qp_info"] = (f"site prepared around {het['resname']} {het['chain']}{het['resseq']}" if het
                        else "no site ligand selected - receptor protonated but not minimized")


def _revert_quickprep(state) -> None:
    original = state["original_receptor"]
    if original:
        logic.swap_receptor(state, original, services.receptor_pdb(original["id"]))


def _receptor_card(state, status) -> None:
    with st.expander("**1 · Receptor**", expanded=True):
        c1, c2 = st.columns([3, 1], vertical_alignment="bottom")
        c1.text_input("PDB ID", placeholder="PDB ID (e.g. 1STP)", max_chars=4, key="pdb_id",
                      on_change=_fetch, args=(state,))
        c2.button("Fetch", on_click=_fetch, args=(state,), use_container_width=True)
        st.file_uploader("Upload .pdb", type=["pdb", "ent", "txt"], key=_uploader_key(state, "pdb_upload"),
                         on_change=_upload_receptor, args=(state,))
        receptor = state["receptor"]
        st.caption(
            f"{receptor['name']} - {receptor['n_protein_atoms']} protein atoms - "
            f"{len(receptor.get('het_ligands') or [])} het groups" if receptor else "no receptor loaded"
        )

        options = logic.engine_options(status)
        st.selectbox("Docking engine", list(options), format_func=options.get, key="engine",
                     on_change=logic.on_engine_change, args=(state, status))
        st.caption(logic.ENGINE_INFO.get(state["engine"], ""))
        if state["engine"] == "gnina" and not status["gnina"]:
            st.warning("gnina is not installed on this machine - run ./setup.sh")

        st.checkbox("delete existing ligand(s) before docking", key="remove_ligand")
        st.checkbox("keep waters in docking", key="keep_waters")
        st.checkbox("drop extra monomers far from site", key="trim_chains")

        if state["engine"] == "moe":
            prepared = bool(receptor and receptor.get("prepared_from"))
            b1, b2 = st.columns(2)
            b1.button("MOE QuickPrep", disabled=not receptor or prepared, on_click=_quickprep, args=(state,),
                      help="Structure Preparation, Protonate3D, distant-water removal and tethered "
                           "minimization with MOE's QuickPrep settings. The selected bound ligand defines "
                           "the site that stays flexible. The prepared receptor replaces the loaded one "
                           "and is what gets docked.")
            if prepared and state["original_receptor"]:
                b2.button("use original", on_click=_revert_quickprep, args=(state,))
            if prepared:
                if state["qp_info"]:
                    st.caption(state["qp_info"])
                report = state["qp_report"]
                if report:
                    st.caption(logic.quickprep_summary(report))
                    cols = st.columns(3)
                    # Site residues come first from the backend and are highlighted.
                    for i, p in enumerate(report["protonation"]):
                        cols[i % 3].button(
                            f"{p['residue']} {p['state']}", key=f"qp_{i}",
                            type="primary" if p["near_site"] else "secondary",
                            help="within 6 Å of the site ligand - click to zoom" if p["near_site"] else "click to zoom",
                            on_click=logic.request_focus, args=(state, logic.residue_focus(p["residue"])),
                        )
        st.caption("Surface/waters/style etc. are all in the viewer's own panel (right side) - "
                   "this is the full Mol* app, not a stripped-down view.")


# ---------------------------------------------------------------------------
# 2 · Active site
# ---------------------------------------------------------------------------


def _load_self_dock(state, smiles: str, label: str | None = None) -> None:
    """Load `smiles` as the ligand to dock (named after the selected bound
    ligand) and turn on RMSD validation against it."""
    het = state["active_het"]
    data = _call(state, services.parse_ligand_text, f"{smiles} {het['resname']}_bound")
    if data is None:
        return
    if not data["n_valid"]:
        logic.toast(state, "that SMILES didn't parse", True)
        return
    logic.show_ligands(state, data)
    state["validate_rmsd"] = True
    logic.toast(state, f"loaded {het['resname']}{f' ({label})' if label else ''} for self-docking - "
                       "run docking to get RMSD")


def _self_dock_rcsb(state) -> None:
    het = state["active_het"]
    if not het:
        logic.toast(state, "pick a bound ligand above first", True)
        return
    found = _call(state, services.pdb_component_smiles, het["resname"], state["receptor"]["id"], het["key"])
    if found:
        _load_self_dock(state, found["smiles"], found["name"])
        if found.get("covalent_note"):
            logic.toast(state, found["covalent_note"])


def _self_dock_manual(state) -> None:
    if not state["active_het"]:
        logic.toast(state, "pick a bound ligand above first", True)
        return
    smiles = state["self_dock_smiles"].strip()
    if not smiles:
        logic.toast(state, "paste the bound ligand's SMILES first", True)
        return
    _load_self_dock(state, smiles)


def _active_site_card(state) -> None:
    with st.expander("**2 · Active site**", expanded=True):
        st.caption("Click a bound ligand in the viewer (or the list below) to box it, or switch to "
                   "**pick atoms** and click residues to build a box from the picks.")
        c1, c2 = st.columns([3, 1], vertical_alignment="bottom")
        c1.segmented_control("Viewer click mode", ["click ligand", "pick atoms"], key="pick_mode",
                             label_visibility="collapsed")
        c2.button("clear picks", on_click=logic.clear_picks, args=(state,))

        receptor = state["receptor"]
        hets = (receptor or {}).get("het_ligands") or []
        with st.popover("Bound ligands (click to set site)", use_container_width=True):
            if not receptor:
                st.caption("load a receptor first")
            elif not hets:
                st.caption('no bound ligands found - use "pick atoms"')
            active_key = (state["active_het"] or {}).get("key")
            for h in hets:
                st.button(f"{h['resname']} {h['chain']}{h['resseq']}  ·  {h['n_atoms']} atoms", key=f"het_{h['key']}",
                          type="primary" if h["key"] == active_key else "secondary",
                          on_click=logic.use_het, args=(state, h), use_container_width=True)

        st.checkbox("validate: after docking, compute RMSD of each pose vs the selected bound ligand's "
                    "real position", key="validate_rmsd")
        st.button("Self-dock selected bound ligand (test)", on_click=_self_dock_rcsb, args=(state,),
                  help="looks up its SMILES from RCSB and loads it as the ligand to dock")
        c1, c2 = st.columns([3, 2], vertical_alignment="bottom")
        c1.text_input("Self-dock SMILES", key="self_dock_smiles", label_visibility="collapsed",
                      placeholder="or paste its SMILES (e.g. a novel/in-house ligand not on RCSB)",
                      on_change=_self_dock_manual, args=(state,))
        c2.button("Self-dock with this SMILES", on_click=_self_dock_manual, args=(state,))

        st.caption("A co-crystallized peptide/peptidomimetic ligand is its own chain, not a bound-ligand "
                   "group above - tick it here to delete it from the receptor before docking (its own "
                   "pocket would otherwise stay occupied).")
        chains = (receptor or {}).get("chains") or []
        if receptor and len(chains) < 2:
            st.caption("single chain - nothing extra to exclude")
        for c in chains if len(chains) >= 2 else []:
            st.checkbox(f"delete chain {c['chain']}  ·  {c['n_atoms']} atom(s), {c['n_hetatm']} het",
                        key=logic.chain_key(state, c["chain"]))

        cols = st.columns(3)
        for i, axis in enumerate("xyz"):
            cols[i].number_input(f"center {axis}", step=0.1, format="%.2f", key=f"c{axis}")
        for i, axis in enumerate("xyz"):
            cols[i].number_input(f"size {axis}", step=0.5, format="%.1f", key=f"s{axis}")
        c1, c2, c3 = st.columns([2, 2, 3], vertical_alignment="bottom")
        c1.number_input("padding (Å)", min_value=0.0, max_value=15.0, step=0.5, key="padding")
        c2.button("make cubic", on_click=logic.make_cubic, args=(state,))
        c3.caption(logic.site_info(state))


# ---------------------------------------------------------------------------
# 3 · Pharmacophore constraints / 3b · Covalent docking
# ---------------------------------------------------------------------------


def _upload_ph4(state) -> None:
    file = _consume_upload(state, "ph4_upload")
    data = file and _call(state, services.upload_pharmacophore, file.name, file.getvalue())
    if data:
        state["ph4_id"], state["ph4_name"] = data["ph4_id"], data["filename"]
        logic.toast(state, f"uploaded {data['filename']}")


def _constraint_row(state, c: dict) -> None:
    cid = c["id"]
    view = logic.constraint_view(state, c)
    hint = logic.CONSTRAINT_RADIUS_HINT.get(view["kind"]) or None
    with st.container(border=True):
        r1 = st.columns([4, 2, 2, 1], vertical_alignment="bottom")
        r1[0].selectbox("kind", list(logic.CONSTRAINT_KIND_LABELS), key=logic.ckey(cid, "kind"),
                        format_func=logic.CONSTRAINT_KIND_LABELS.get, help=hint, label_visibility="collapsed")
        armed = state["pending_pick_id"] == cid
        r1[1].button("picking…" if armed else "pick point", key=f"pick_{cid}", type="primary" if armed else "secondary",
                     on_click=logic.toggle_constraint_pick, args=(state, cid))
        r1[2].caption(logic.constraint_status(view))
        r1[3].button("✕", key=f"rm_{cid}", help="remove", on_click=logic.remove_constraint, args=(state, cid))
        r2 = st.columns(5, vertical_alignment="bottom")
        for i, axis in enumerate("xyz"):
            r2[i].number_input(axis, step=0.1, format="%.2f", key=logic.ckey(cid, axis))
        r2[3].number_input("radius Å", min_value=0.5, max_value=15.0, step=0.1, key=logic.ckey(cid, "radius"), help=hint)
        r2[4].checkbox("required", key=logic.ckey(cid, "required"))
        st.text_input("label", key=logic.ckey(cid, "label"), placeholder="label (optional)", label_visibility="collapsed")
        if view["kind"] == "smarts":
            st.text_input("SMARTS", key=logic.ckey(cid, "smarts"), placeholder="SMARTS pattern, e.g. [OX1]=[C,S]")


def _constraints_card(state) -> None:
    with st.expander("**3 · Pharmacophore constraints** (optional)", expanded=False):
        if state["engine"] == "moe":
            st.caption("**MOE pharmacophore (.ph4 file)** - upload a MOE pharmacophore query file to constrain docking.")
            st.file_uploader("Upload .ph4 file", type=["ph4"], key=_uploader_key(state, "ph4_upload"),
                             on_change=_upload_ph4, args=(state,))
            if state["ph4_name"]:
                st.caption(state["ph4_name"])
            return
        st.caption('Checked against each pose after docking (smina/gnina have no way to bias the search itself '
                   'toward a point) - poses meeting every "required" constraint are ranked first. Click "pick '
                   'point", then click an atom in the viewer.')
        for c in state["constraints"]:
            _constraint_row(state, c)
        st.button("+ add constraint", on_click=logic.add_constraint, args=(state,))


def _upload_rxn(state) -> None:
    file = _consume_upload(state, "rxn_upload")
    data = file and _call(state, services.upload_reaction, file.name, file.getvalue())
    if data:
        state["custom_rxn_id"], state["custom_rxn_name"] = data["rxn_id"], data["filename"]
        state["covalent_reaction"] = ""  # the custom file is used instead
        logic.toast(state, f"uploaded {data['filename']}")


def _toggle_covalent_pick(state) -> None:
    state["covalent_pick_active"] = not state["covalent_pick_active"]
    logic.toast(state, "Click a nucleophilic atom (CYS SG, SER OG, LYS NZ, etc.) in the viewer"
                if state["covalent_pick_active"] else "Atom picking cancelled")


def _find_residues(state) -> None:
    if not state["receptor"]:
        logic.toast(state, "load receptor first", True)
        return
    resname = state["residue_search"].strip().upper()
    if len(resname) != 3 or not resname.isalpha():
        logic.toast(state, "enter 3-letter residue code (e.g., CYS, SER, LYS)", True)
        return
    data = _call(state, services.list_residues, state["receptor"]["id"], resname)
    if data:
        state["residue_hits"] = (resname, data["residues"])


def _set_covalent_residue(state, spec: str) -> None:
    state["cov_residue"] = spec
    logic.toast(state, f"selected {spec}")


def _covalent_card(state, status) -> None:
    if state["engine"] != "moe":
        return
    with st.expander("**3b · Covalent docking** (MOE only)", expanded=False):
        st.checkbox("**Enable covalent docking**", key="covalent_enabled")
        if not state["covalent_enabled"]:
            return
        reactions = [""] + list(status.get("covalent_reactions") or [])
        custom = state["custom_rxn_name"]

        def label(rxn: str) -> str:
            if rxn:
                return rxn.replace("_", " ")
            return f"Custom: {custom}" if custom else "Select reaction..."

        c1, c2, c3 = st.columns([3, 3, 1], vertical_alignment="bottom")
        c1.selectbox("warhead type", reactions, format_func=label, key="covalent_reaction")
        c2.text_input("target residue", placeholder="A:CYS:189", key="cov_residue")
        c3.button("👆", help="Pick: click an atom in the viewer", on_click=_toggle_covalent_pick, args=(state,),
                  type="primary" if state["covalent_pick_active"] else "secondary")
        c1, c2 = st.columns([1, 2], vertical_alignment="bottom")
        c1.text_input("residue type", placeholder="CYS, SER, LYS...", max_chars=3, key="residue_search",
                      on_change=_find_residues, args=(state,))
        c2.button("Find in receptor", on_click=_find_residues, args=(state,))
        if state["residue_hits"]:
            resname, specs = state["residue_hits"]
            if not specs:
                st.caption(f"no {resname} residues found")
            else:
                st.caption(f"{len(specs)} {resname} residue(s) - click to use")
                cols = st.columns(4)
                for i, spec in enumerate(specs):
                    cols[i % 4].button(spec, key=f"res_{spec}", on_click=_set_covalent_residue, args=(state, spec))
        st.file_uploader("Upload custom .rxn (optional)", type=["rxn"], key=_uploader_key(state, "rxn_upload"),
                         on_change=_upload_rxn, args=(state,))
        if custom:
            st.caption(custom)
        st.caption("Format: `chain:resname:resnum` (e.g., A:CYS:189, B:SER:25, :LYS:42)  \n"
                   "Draw custom reactions in [Ketcher](https://lifescience.opensource.epam.com/ketcher/index.html), "
                   "export as .rxn, then upload.")


# ---------------------------------------------------------------------------
# 4 · Ligands
# ---------------------------------------------------------------------------


def _parse_box(state) -> None:
    text = state["smiles_box"].strip()
    if not text:
        logic.toast(state, "paste SMILES first", True)
        return
    data = _call(state, services.parse_ligand_text, text)
    if data:
        logic.show_ligands(state, data)


def _upload_ligands(state) -> None:
    file = _consume_upload(state, "lig_upload")
    data = file and _call(state, services.parse_ligand_file, file.name, file.getvalue())
    if data:
        logic.show_ligands(state, data)


def _ligands_card(state) -> None:
    with st.expander("**4 · Ligands**", expanded=True):
        st.text_area("SMILES", key="smiles_box", height=130, label_visibility="collapsed",
                     placeholder="One SMILES per line, optional name after a space:\n"
                                 "CC(=O)Oc1ccccc1C(=O)O aspirin\nCN1C=NC2=C1C(=O)N(C)C(=O)N2C caffeine")
        c1, c2 = st.columns([1, 3], vertical_alignment="center")
        c1.button("Load from box", on_click=_parse_box, args=(state,))
        c2.caption(logic.ligand_info(state["lig_preview"]))
        st.file_uploader("Upload .smi/.sdf/.csv", key=_uploader_key(state, "lig_upload"),
                         type=["smi", "smiles", "txt", "sdf", "sd", "mol", "csv", "tsv"],
                         on_change=_upload_ligands, args=(state,))
        data = state["lig_preview"]
        if data:
            rows = [{"name": l["name"], "smiles": l["smiles"] or "(3D block)", "note": l["note"],
                     "ok": "✓" if l["valid"] else "✗"} for l in data["ligands"][:500]]
            st.dataframe(rows, hide_index=True, use_container_width=True, height=min(250, 38 + 35 * len(rows)))
            if len(data["ligands"]) > 500:
                st.caption(f"... {len(data['ligands']) - 500} more")


# ---------------------------------------------------------------------------
# 5 · Docking parameters
# ---------------------------------------------------------------------------


def _run(state) -> None:
    payload, err = logic.build_job_request(state)
    if err:
        logic.toast(state, err, True)
        return
    job = _call(state, lambda: services.create_job(JobRequest(**payload)))
    if job:
        state["job_id"] = job["id"]
        state["selected_ligand"] = None
        logic.toast(state, f"job {job['id']} started - {job['n_ligands']} ligand(s)")


def _cancel(state) -> None:
    try:
        services.cancel_job(state["job_id"])
        logic.toast(state, "cancelling - running ligands finish first")
    except services.ServiceError as exc:
        logic.toast(state, exc.detail, True)


def job_running(state) -> bool:
    job = services.MANAGER.get(state["job_id"]) if state["job_id"] else None
    return bool(job and job.status not in ("done", "cancelled"))


def _params_card(state, status) -> None:
    with st.container(border=True):
        st.markdown("**5 · Docking parameters**")
        moe = state["engine"] == "moe"
        cols = st.columns(3)
        cols[0].number_input("exhaustiveness", min_value=1, max_value=64, step=1, key="exh")
        cols[1].number_input("poses kept", min_value=1, max_value=20, step=1, key="modes")
        cols[2].number_input("energy range", min_value=0.5, max_value=20.0, step=0.5, key="erange")
        cols[0].number_input("seed", step=1, key="seed")
        cols[1].selectbox("rescoring" if moe else "scoring", logic.scoring_options(status, state["engine"]),
                          key="scoring",
                          help="MOE rescoring after refinement - placement always uses London dG" if moe else None)
        if not moe:
            st.checkbox("peptide / peptidomimetic ligands (bump exhaustiveness & energy range)", key="peptide_mode",
                        on_change=logic.apply_peptide_mode, args=(state,),
                        help="raises exhaustiveness/energy range to search the larger, floppier conformer space "
                             "peptides and peptidomimetics need")
            st.checkbox("minimize each pose after docking", key="minimize_after",
                        help="extra smina/gnina --minimize pass on every pose after its own docking search "
                             "finishes - relaxes clashes a coarse exhaustiveness/energy range left behind")
        if state["engine"] == "gnina":
            c1, c2 = st.columns([1, 2], vertical_alignment="bottom")
            c1.selectbox("device", ["GPU if available", "CPU only"], key="device")
            c2.caption("NVIDIA GPU detected on this machine" if status["gpu_available"]
                       else "no NVIDIA GPU detected - gnina will run on CPU")
        if moe:
            c1, c2 = st.columns(2)
            c1.selectbox("placement", logic.MOE_PLACEMENTS, key="moe_placement")
            c2.selectbox("refinement", logic.MOE_REFINEMENTS, key="moe_refinement")
        running = job_running(state)
        c1, c2 = st.columns(2)
        c1.button("Run", type="primary", disabled=running, on_click=_run, args=(state,), use_container_width=True)
        c2.button("Cancel", disabled=not running, on_click=_cancel, args=(state,), use_container_width=True)


def render(state, status) -> None:
    _receptor_card(state, status)
    _active_site_card(state)
    _constraints_card(state)
    _covalent_card(state, status)
    _ligands_card(state)
    _params_card(state, status)
