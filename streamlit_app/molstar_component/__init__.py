"""Streamlit component wrapping the full Mol* viewer (see README.md).

The frontend keeps one Mol* plugin alive across reruns and only rebuilds the
pieces whose args changed. It reports atom clicks and pose-bar selections back
as event dicts carrying a strictly increasing ``nonce`` - the last event is
re-returned on every rerun, so callers must handle each nonce only once.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Iterable

import streamlit.components.v1 as components

_FRONTEND_DIR = Path(__file__).resolve().parent / "frontend"
_BUILD_DIR = _FRONTEND_DIR / "dist"

_component_func = None


def _component():
    """Declare the component lazily so importing this module never fails."""
    global _component_func
    if _component_func is None:
        dev_url = os.environ.get("DOCKING_COMPONENT_DEV_URL")
        if dev_url:
            _component_func = components.declare_component("molstar_viewer", url=dev_url)
        else:
            if not (_BUILD_DIR / "index.html").exists():
                raise RuntimeError(
                    f"Mol* viewer component is not built - run `npm install && npm run build` in {_FRONTEND_DIR}"
                )
            _component_func = components.declare_component("molstar_viewer", path=str(_BUILD_DIR))
    return _component_func


def molstar_viewer(
    *,
    receptor: dict | None = None,
    het_ligands: Iterable[dict] = (),
    style: dict | None = None,
    box: dict | None = None,
    picks: Iterable[Iterable[float]] = (),
    constraints: Iterable[dict] = (),
    poses: dict | None = None,
    depiction: dict | None = None,
    focus: dict | None = None,
    height: int = 620,
    key: str | None = None,
) -> dict[str, Any] | None:
    """Render the viewer; returns the latest event dict, or None before any.

    receptor:    {"id", "name", "pdb"} - reloaded only when "id" changes
    het_ligands: the receptor's bound groups (prep.HetLigand.as_dict())
    style:       {"hide_ligand", "hide_protein", "active_het", "prepared"}
    box:         {"center": [x, y, z], "size": [x, y, z]}
    picks:       [[x, y, z], ...] pick-point markers
    constraints: [{"id", "kind", "x", "y", "z", optional "nx".."nz", "dx".."dz"}]
    poses:       {"key", "name", "blocks": [sdf block, ...], "scores": [...]} -
                 reloaded only when "key" changes
    depiction:   {"name", "svg"} 2D ligand overlay
    focus:       {"nonce", "kind": "het", "het": {...}} or
                 {"nonce", "kind": "residue", "chain", "resseq"} - zooms once per nonce

    Events: {"nonce", "type": "click", "atom": {...}} and
            {"nonce", "type": "pose", "index": int} (-1 = all poses).
    """
    return _component()(
        receptor=receptor,
        het_ligands=list(het_ligands),
        style=style,
        box=box,
        picks=[list(p) for p in picks],
        constraints=list(constraints),
        poses=poses,
        depiction=depiction,
        focus=focus,
        height=height,
        key=key,
        default=None,
    )


__all__ = ["molstar_viewer"]
