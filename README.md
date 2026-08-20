# Docking Bench

Local web app for small-molecule docking against a protein receptor, using
[smina](https://sourceforge.net/projects/smina/) as the docking engine and
[RDKit](https://www.rdkit.org/) for ligand 3D structure generation. Runs
entirely on your machine — no external services beyond fetching PDB files
from RCSB and a CDN-hosted 3Dmol.js viewer script.

## Features

- Fetch a receptor by PDB ID (RCSB) or upload a `.pdb` file
- Auto-detect bound HETATM ligands/cofactors to seed the docking box, or
  click atoms in the 3D viewer to build a box manually
- Load ligands from pasted SMILES, or upload `.smi`/`.sdf`/`.mol`/`.csv`/`.tsv`
- Dock many ligands concurrently (one smina process per ligand, thread pool)
- Live job progress, per-ligand poses viewable in-browser, CSV/ZIP export
- Job history survives a server restart (state reloaded from `runs/jobs/`)

## Setup

```bash
./setup.sh   # installs smina (via micromamba) into .env-tools/, and a
             # Python venv into .venv/ with fastapi, rdkit, etc.
```

Safe to re-run; existing installs are reused. Supported platforms: macOS
(arm64/x86_64), Linux (x86_64/aarch64).

## Run

```bash
./run.sh                 # serves on http://127.0.0.1:8000
PORT=9000 ./run.sh        # pick a different port
```

Open the printed URL in a browser. The top-right pill shows whether smina
was found; if not, re-run `./setup.sh`.

## Testing

```bash
.venv/bin/python -m pytest tests/ -q
```

Tests cover receptor cleanup/HET detection, ligand parsing, smina command
construction, and job persistence — they mock the smina subprocess call, so
they run without smina installed.

## Project layout

```
backend/
  main.py     FastAPI app: receptor/ligand/job HTTP routes
  jobs.py     Job store: one job -> many ligands, each docked in a thread
  docking.py  Thin wrapper around the smina CLI
  prep.py     Receptor cleanup, HET-ligand grouping, ligand parsing/3D embed
  config.py   Paths and runtime settings (env-var overridable)
frontend/
  index.html, app.js, style.css   Static single-page UI (3Dmol.js viewer)
tests/        pytest suite for the backend
runs/         receptor uploads + job outputs (gitignored; created at runtime)
```

## Notable env vars

| Variable            | Default                  | Purpose                              |
|----------------------|---------------------------|---------------------------------------|
| `PORT`               | `8000`                    | `run.sh` listen port                  |
| `DOCKING_RUNS_DIR`   | `<repo>/runs`             | where receptors/job outputs are stored |
| `DOCKING_WORKERS`    | `cpu_count - 1`           | max concurrent smina processes        |
| `SMINA_TIMEOUT`      | `1800` (seconds)          | per-ligand smina hard timeout         |
| `SMINA_BIN`          | (unset)                   | explicit path to a smina executable   |

## Known limitations

- mmCIF receptors are rejected — upload a `.pdb` file
- No authentication; intended for local/single-user use
- `runs/` is not pruned automatically — old jobs accumulate on disk
