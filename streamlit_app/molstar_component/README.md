# Mol* viewer component

A Streamlit custom component (v1 API) that wraps the full [Mol*](https://molstar.org/)
viewer, ported from the original `frontend/src/app.js` viewer code. One Mol*
plugin lives for the whole session; every rerun diffs the args against what
was last applied and only rebuilds what changed (a 1 MB receptor isn't
reloaded because a checkbox flipped).

```python
from molstar_component import molstar_viewer

event = molstar_viewer(receptor={"id": rid, "name": name, "pdb": pdb_text},
                       box={"center": c, "size": s}, key="viewer")
```

## Args

| arg | shape | notes |
|---|---|---|
| `receptor` | `{id, name, pdb}` or `None` | reloaded only when `id` changes; clears poses and re-applies every overlay |
| `het_ligands` | `[prep.HetLigand.as_dict(), ...]` | excluded from the docked-pose selection for contacts |
| `style` | `{hide_ligand, hide_protein, active_het, prepared}` | cartoon / green bound ligands / orange active het / QuickPrep pocket sticks |
| `box` | `{center: [x,y,z], size: [x,y,z]}` | blue wireframe |
| `picks` | `[[x,y,z], ...]` | yellow pick markers |
| `constraints` | `[{id, kind, x, y, z, nx?..nz?, dx?..dz?}]` | colored markers; aromatic disk / H-bond arrow when the normal/direction is set |
| `poses` | `{key, name, blocks: [sdf block], scores: [float or None]}` | reloaded only when `key` changes; drives the pose bar |
| `depiction` | `{name, svg}` | 2D overlay top-left of the canvas (shown via `<img>`, never inlined) |
| `focus` | `{nonce, kind: "het", het}` / `{nonce, kind: "residue", chain, resseq}` | zooms once per new `nonce` |
| `height` | px | viewer height; the pose bar adds to the frame height |

## Events (return value)

The component returns the **latest** event, and Streamlit re-returns it on
every rerun, so handle each `nonce` once (keep the last handled nonce in
`st.session_state`).

- `{"nonce", "type": "click", "atom": {resn, hetflag, chain, resi, atom, x, y, z, ring}}`.
  `ring` is `{centroid, normal}` for an aromatic side-chain ring atom, else `None`.
  Clicks on the synthetic box/pick/constraint atoms are ignored.
- `{"nonce", "type": "pose", "index"}`. `index` is the selected pose (0-based, `-1` = all).
  It's sent when a new poses key loads (index 0) and when the user picks a pose in the bar.
  It's not sent during playback steps, only once when playback stops.

## Build

```bash
cd streamlit_app/molstar_component/frontend
npm install
npm run build   # type-check + bundle into dist/ (what the Python wrapper serves)
npm test        # vitest unit tests for the pure geometry / pose-stepping code
```

## Dev mode

```bash
npm run dev                                         # vite on http://localhost:3001
DOCKING_COMPONENT_DEV_URL=http://localhost:3001 streamlit run ...
```

## Layout

```
frontend/src/
  main.ts            Streamlit glue: arg diffing, events, depiction overlay
  viewer.ts          ViewerController - owns the Mol* plugin and every state ref
  poseBar.ts         pose bar + playback (1200 ms, "light" steps)
  synthetic.ts       pure geometry: PDB lines, box/marker points, ring centroid
  synthetic.test.ts  unit tests
  types.ts           arg / event types
```
