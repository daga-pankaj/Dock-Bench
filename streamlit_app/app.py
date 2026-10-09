"""Docking Bench - Streamlit front end.

Same workflow as the Mol* web page (frontend/), served as one Python
process so it is easy to run on a shared machine:

    streamlit run streamlit_app/app.py

All docking work goes through backend.services; UI decisions live in
streamlit_app/logic.py; this file lays out the page.
"""

from __future__ import annotations

import hmac
import sys
from pathlib import Path

import streamlit as st

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from backend import config, services  # noqa: E402
from streamlit_app import logic, results_panel, setup_panel  # noqa: E402

st.set_page_config(page_title="Docking Bench", page_icon="🧪", layout="wide")
config.ensure_dirs()


@st.cache_resource
def announce_generated_password() -> None:
    """DOCKING_REQUIRE_AUTH makes up a password at startup - print it once
    per server process (this script itself reruns on every interaction)."""
    if config.AUTH_GENERATED:
        print(f"\n>>> auth enabled - user: {config.AUTH_USER}  password: {config.AUTH_PASS}\n", flush=True)


def require_password() -> None:
    """Gate the app behind DOCKING_AUTH_PASS when it's set (same variable
    the web page's HTTP Basic auth uses) - turn on before sharing a link."""
    if config.AUTH_PASS is None or st.session_state.get("authed"):
        return
    st.title("Docking Bench")
    with st.form("login"):
        user = st.text_input("User", value=config.AUTH_USER)
        password = st.text_input("Password", type="password")
        if st.form_submit_button("Sign in"):
            if hmac.compare_digest(user, config.AUTH_USER) and hmac.compare_digest(password, config.AUTH_PASS):
                st.session_state["authed"] = True
                st.rerun()
            st.error("wrong user or password")
    st.stop()


@st.cache_data(ttl=300, show_spinner=False)
def app_status() -> dict:
    return services.status()


def status_line(status: dict) -> tuple[str, bool]:
    parts = [
        status["smina_version"] or "smina" if status["smina"] else "smina not found",
        status["gnina_version"] or "gnina" if status["gnina"] else "gnina not found",
    ]
    if status["moe_available"]:
        parts.append("MOE available")
    ok = bool(status["smina"] or status["gnina"] or status["moe_available"])
    text = f"{' - '.join(parts)} - {status['workers']} workers{' - GPU detected' if status['gpu_available'] else ''}"
    return text, ok


def flush_toasts() -> None:
    for message, bad in st.session_state["toasts"]:
        st.toast(message, icon="⚠️" if bad else None)
    st.session_state["toasts"] = []


def main() -> None:
    announce_generated_password()
    require_password()
    state = st.session_state
    logic.init_state(state)
    # segmented controls can be clicked off; one option is always active
    state["pick_mode"] = state["pick_mode"] or "click ligand"
    state["diagram_style"] = state["diagram_style"] or "MOE"
    status = app_status()
    if state["engine"] not in logic.engine_options(status):
        state["engine"] = next(iter(logic.engine_options(status)))
        logic.on_engine_change(state, status)

    # The viewer's last event (a click or pose pick) is applied before any
    # widget is drawn, so it can update the box fields, constraint points
    # and covalent residue those widgets show.
    logic.handle_event(state, state.get("viewer"))
    flush_toasts()

    head, pill = st.columns([3, 2], vertical_alignment="center")
    head.markdown("## 🧪 Docking Bench  \n<span style='opacity:.7'>Local molecular docking workspace</span>",
                  unsafe_allow_html=True)
    text, ok = status_line(status)
    (pill.success if ok else pill.error)(text)

    left, right = st.columns([5, 7], gap="medium")
    with left:
        setup_panel.render(state, status)
    with right:
        results_panel.render(state, status)
    flush_toasts()


main()
