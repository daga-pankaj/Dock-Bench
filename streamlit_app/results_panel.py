"""Right column: 3D viewer, 2D interaction diagram, results."""

from __future__ import annotations

import base64

import streamlit as st

from backend import services
from backend.docking import MOE_SCORE_FIELDS
from streamlit_app import logic
from streamlit_app.molstar_component import molstar_viewer

POLL_SECONDS = 1.5


@st.cache_data(show_spinner=False, max_entries=64)
def _depiction(job_id: str, name: str) -> str | None:
    try:
        return services.ligand_depiction_svg(job_id, name)
    except services.ServiceError:
        return None


@st.cache_data(show_spinner=False, max_entries=256)
def _diagram(job_id: str, name: str, pose: int, style: str) -> tuple[str | None, str | None, str | None]:
    """(svg, source, error) - cached so a failed MOE attempt isn't re-run on every rerun."""
    try:
        svg, source = services.interactions_2d(job_id, name, pose, style)
        return svg, source, None
    except services.ServiceError as exc:
        return None, None, exc.detail


def _svg_img(svg: str, alt: str, width: str = "100%") -> None:
    """SVG shown through <img>, never inlined: an <img>-rendered SVG can't run
    scripts even if escaping of PDB-derived text ever had a gap."""
    data = base64.b64encode(svg.encode()).decode()
    st.markdown(f'<img src="data:image/svg+xml;base64,{data}" alt="{alt}" style="width:{width};border-radius:6px">',
                unsafe_allow_html=True)


def _current_job(state) -> dict | None:
    if not state["job_id"]:
        return None
    try:
        return services.get_job(state["job_id"]).as_dict()
    except services.ServiceError:
        return None


def _poses_arg(state, job: dict | None) -> dict | None:
    name = state["selected_ligand"]
    lig = next((l for l in (job or {}).get("ligands", []) if l["name"] == name), None)
    if lig is None:
        return None
    try:
        sdf = services.ligand_poses_path(job["id"], name).read_text()
    except services.ServiceError:
        return None
    blocks = [b.strip() + "\n$$$$\n" for b in sdf.split("$$$$") if len(b.strip()) > 10]
    return {"key": f"{job['id']}/{name}", "name": name, "blocks": blocks, "scores": lig["scores"]}


# ---------------------------------------------------------------------------
# viewer
# ---------------------------------------------------------------------------


def _toggle(state, key: str) -> None:
    state[key] = not state[key]


def _zoom_site(state) -> None:
    if state["active_het"]:
        logic.request_focus(state, {"kind": "het", "het": logic.het_ref(state["active_het"])})
    elif state["picks"]:
        logic.toast(state, "zoom not available for free-hand picks yet - drag/scroll to frame the box", True)
    else:
        logic.toast(state, "set a docking box first", True)


def _viewer_card(state, job: dict | None) -> None:
    with st.container(border=True):
        head, tools = st.columns([1, 3], vertical_alignment="center")
        head.markdown("**Viewer**")
        t = tools.columns(4)
        t[0].button("hide ligand", type="primary" if state["hide_ligand"] else "secondary",
                    help="hide co-crystallized ligand from the complex",
                    on_click=_toggle, args=(state, "hide_ligand"), use_container_width=True)
        t[1].button("delete ligand", type="primary" if state["remove_ligand"] else "secondary",
                    help="delete the existing ligand(s) from the receptor before docking "
                         "(same as the checkbox in 1 · Receptor)",
                    on_click=_toggle, args=(state, "remove_ligand"), use_container_width=True)
        t[2].button("hide protein", type="primary" if state["hide_protein"] else "secondary",
                    help="hide the receptor protein from the complex",
                    on_click=_toggle, args=(state, "hide_protein"), use_container_width=True)
        t[3].button("site", help="zoom to docking box", on_click=_zoom_site, args=(state,),
                    use_container_width=True)

        if state["covalent_pick_active"]:
            st.info("click a nucleophilic atom (CYS SG, SER OG, LYS NZ, ...) in the viewer")
        elif state["pending_pick_id"] is not None:
            st.info("click an atom in the viewer to set the constraint point")

        receptor = state["receptor"]
        center = logic.box_center(state)
        poses = _poses_arg(state, job)
        svg = _depiction(job["id"], poses["name"]) if poses else None
        try:
            molstar_viewer(
                receptor={"id": receptor["id"], "name": receptor["name"], "pdb": state["pdb_text"]} if receptor else None,
                het_ligands=(receptor or {}).get("het_ligands") or [],
                style={
                    "hide_ligand": state["hide_ligand"],
                    "hide_protein": state["hide_protein"],
                    "active_het": logic.het_ref(state["active_het"]) if state["active_het"] else None,
                    "prepared": bool(receptor and receptor.get("prepared_from")),
                },
                box={"center": center, "size": logic.box_size(state)} if center else None,
                picks=state["picks"],
                constraints=logic.constraint_markers(state),
                poses=poses,
                depiction={"name": poses["name"], "svg": svg} if svg else None,
                focus=state["focus"],
                height=620,
                key="viewer",
            )
        except RuntimeError as exc:  # component not built yet
            st.error(str(exc))


# ---------------------------------------------------------------------------
# 2D interaction diagram
# ---------------------------------------------------------------------------


def _diagram_card(state, status, job: dict | None) -> None:
    with st.container(border=True):
        head, opts = st.columns([2, 2], vertical_alignment="center")
        head.markdown("**2D Interaction Diagram**")
        with opts:
            c1, c2 = st.columns([2, 1], vertical_alignment="center")
            if status["moe_available"]:
                c1.segmented_control("diagram style", ["MOE", "ProLIF"], key="diagram_style",
                                     label_visibility="collapsed",
                                     help="MOE's Ligand Interactions diagram, or the ProLIF-based one")
            c2.toggle("show", key="show_diagram")
        if not state["show_diagram"]:
            return
        name, pose = state["selected_ligand"], state["active_pose"]
        if not (job and name) or pose < 0:
            st.caption("select a pose to see its interaction diagram")
            return
        style = "moe" if status["moe_available"] and state["diagram_style"] == "MOE" else "prolif"
        with st.spinner("loading diagram…"):
            svg, source, error = _diagram(job["id"], name, pose + 1, style)
        if error:
            st.caption("diagram unavailable")
            st.error(f"interaction diagram: {error}")
            return
        if style == "moe" and source == "prolif":
            st.caption("MOE couldn't draw this diagram - showing the ProLIF one")
        _svg_img(svg, "2D protein-ligand interaction diagram")


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------


def _fmt(v, digits: int = 2) -> str:
    return "" if v is None else f"{v:.{digits}f}"


def _pose_breakdown(lig: dict) -> None:
    """Every pose's scores: MOE's full per-pose table, or score/strain/RMSD."""
    strain, rmsd = lig["strain_energies"], lig["rmsd_to_reference"]
    has_strain = any(v is not None for v in strain)
    has_rmsd = any(v is not None for v in rmsd)
    rows = []
    if lig["pose_scores"]:
        has_attach = any(lig["covalent_attachments"])
        for i, ps in enumerate(lig["pose_scores"]):
            row = {"pose": f"#{i + 1}", **{k: _fmt(ps.get(k)) for k in MOE_SCORE_FIELDS}}
            if has_strain:
                row["strain"] = _fmt(strain[i] if i < len(strain) else None, 1)
            if has_rmsd:
                row["RMSD"] = _fmt(rmsd[i] if i < len(rmsd) else None)
            if has_attach:
                row["bonded via"] = lig["covalent_attachments"][i] or ""
            rows.append(row)
        st.dataframe(rows, hide_index=True, use_container_width=True,
                     column_config={k: st.column_config.TextColumn(k, help=v) for k, v in MOE_SCORE_FIELDS.items()})
    else:
        for i, s in enumerate(lig["scores"]):
            row = {"pose": f"#{i + 1}", "kcal/mol": _fmt(s)}
            if has_strain:
                row["strain"] = _fmt(strain[i] if i < len(strain) else None, 1)
            if has_rmsd:
                row["RMSD Å"] = _fmt(rmsd[i] if i < len(rmsd) else None)
            rows.append(row)
        st.dataframe(rows, hide_index=True, use_container_width=True)


def _ligand_details(state, job: dict, lig: dict) -> None:
    is_moe = bool(lig["pose_scores"])
    has_breakdown = len(lig["scores"]) > 1 or is_moe or lig["warnings"]
    if not (has_breakdown or lig["error"] or lig["constraint_detail"]):
        return
    # MOE's per-pose score table is the point of a MOE result, so it starts
    # open (unless hidden from the Results header); other engines' start closed.
    with st.expander(f"{lig['name']} - details", expanded=is_moe and not state["hide_moe_scores"]):
        for w in lig["warnings"]:
            st.warning(w)
        if lig["error"]:
            st.error(lig["error"])
        if lig["constraint_detail"]:
            st.caption("constraints (top pose):  \n" + logic.constraint_detail_text(lig).replace("\n", "  \n"))
        if lig["scores"]:
            _pose_breakdown(lig)
        if lig["status"] == "done":
            c1, c2 = st.columns(2)
            try:
                path = services.ligand_poses_path(job["id"], lig["name"])
                c1.download_button("poses .sdf", path.read_bytes(), file_name=path.name,
                                   mime="chemical/x-mdl-sdfile", key=f"sdf_{job['id']}_{lig['name']}")
            except services.ServiceError:
                pass
            if c2.toggle("show log", key=f"log_{job['id']}_{lig['name']}"):
                try:
                    st.code(services.ligand_log(job["id"], lig["name"]), language=None)
                except services.ServiceError as exc:
                    st.caption(exc.detail)


def _view_ligand(state, name: str | None) -> None:
    """Show `name`'s poses. Widgets inside the polling fragment only rerun the
    fragment, so this flags a whole-page rerun for the viewer to pick it up."""
    if name and name != state["selected_ligand"]:
        state["selected_ligand"] = name
        state["active_pose"] = 0
        state["rerun_page"] = True


def _on_table_select(state, rows: list[dict]) -> None:
    selected = state["results_table"]["selection"]["rows"]
    if not selected:
        return
    lig = rows[selected[0]]
    if lig["status"] != "done":
        logic.toast(state, f"{lig['name']} has no poses ({lig['status']})", True)
    else:
        _view_ligand(state, lig["name"])


def _open_job(state) -> None:
    state["job_id"] = state["open_job"] or None
    state["selected_ligand"] = None


def _results_body(state) -> None:
    job = _current_job(state)
    if job is None:
        st.caption("No job yet. Results appear here as each ligand finishes.")
        return

    finished = job["status"] in ("done", "cancelled")
    st.progress(job["n_finished"] / job["n_ligands"] if job["n_ligands"] else 0.0,
                text=f"{job['id']} - {job['status']} - {job['n_finished']}/{job['n_ligands']}")
    if finished and state["job_done_seen"] != job["id"]:
        state["job_done_seen"] = job["id"]
        ok = sum(1 for l in job["ligands"] if l["status"] == "done")
        st.toast(f"job {job['status']}: {ok}/{job['n_ligands']} docked")
        if logic.rmsd_mismatch(job):
            st.toast("RMSD validation was on, but no docked pose matched the bound ligand's heavy-atom structure - "
                     "double check the SMILES you docked is really the same molecule as the crystal ligand",
                     icon="⚠️")
        st.rerun()  # re-enable Run and stop polling

    rows = logic.sorted_results(job["ligands"])
    table = [{
        "ligand": l["name"],
        "status": l["status"] + (" ⚠" if l["warnings"] else ""),
        "best kcal/mol": _fmt(l["best_score"]),
        "poses": str(l["n_poses"] or ""),
        "constraints": logic.constraint_text(l),
        # top-ranked pose vs the bound ligand - only for a validation run
        "RMSD Å": _fmt(l["rmsd_to_reference"][0] if l["rmsd_to_reference"] else None),
        "smiles": l["smiles"] or "",
    } for l in rows]
    st.dataframe(
        table, hide_index=True, use_container_width=True, key="results_table",
        on_select=lambda: _on_table_select(state, rows), selection_mode="single-row",
        column_config={"RMSD Å": st.column_config.TextColumn(help="top pose vs the bound ligand (≤ 2 Å is good)")},
    )
    done = [l["name"] for l in rows if l["status"] == "done"]
    state["view_ligand"] = state["selected_ligand"] if state["selected_ligand"] in done else ""
    st.selectbox("View poses of", [""] + done, key="view_ligand",
                 format_func=lambda n: n or "select a row above, or pick a ligand here",
                 on_change=lambda: _view_ligand(state, state["view_ligand"]))
    if state.pop("rerun_page", False):
        st.rerun()  # whole page, so the viewer loads the new poses

    if any(l["pose_scores"] for l in job["ligands"]):
        st.toggle("hide MOE scores", key="hide_moe_scores")
    for lig in rows:
        _ligand_details(state, job, lig)

    if finished:
        c1, c2 = st.columns(2)
        c1.download_button("Download all (ZIP)", services.results_zip(job["id"]), file_name=f"docking_{job['id']}.zip",
                           mime="application/zip", use_container_width=True)
        c2.download_button("CSV", services.results_csv(job["id"]), file_name=f"docking_{job['id']}_results.csv",
                           mime="text/csv", use_container_width=True)


def _results_card(state) -> None:
    with st.container(border=True):
        head, pick = st.columns([1, 2], vertical_alignment="center")
        head.markdown("**Results**")
        jobs = services.list_jobs()
        ids = [""] + [j["id"] for j in jobs]
        labels = {j["id"]: f"{j['id']} · {j['n_finished']}/{j['n_ligands']} · {j['status']}" for j in jobs}
        if state["job_id"] in ids:
            state["open_job"] = state["job_id"]
        pick.selectbox("previous jobs", ids, key="open_job", format_func=lambda i: labels.get(i, "open a previous job…"),
                       on_change=_open_job, args=(state,), label_visibility="collapsed")
        job = services.MANAGER.get(state["job_id"]) if state["job_id"] else None
        running = bool(job and job.status not in ("done", "cancelled"))
        st.fragment(_results_body, run_every=POLL_SECONDS if running else None)(state)


def render(state, status) -> None:
    job = _current_job(state)
    _viewer_card(state, job)
    _diagram_card(state, status, job)
    _results_card(state)
