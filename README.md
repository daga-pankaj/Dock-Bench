# Docking Bench

Local web app for small-molecule docking against a protein receptor, using
[smina](https://sourceforge.net/projects/smina/) or
[gnina](https://github.com/gnina/gnina) (GPU-capable, Linux only) as the
docking engine and [RDKit](https://www.rdkit.org/) for ligand 3D structure
generation. Runs entirely on your machine — no external services beyond
fetching PDB files from RCSB. The 3D viewer is [Mol*](https://molstar.org/),
bundled locally at setup time (see "Frontend viewer build" below).

## Features

- Fetch a receptor by PDB ID (RCSB) or upload a `.pdb` file
- Auto-detect binding sites using [fpocket](https://github.com/Discngine/fpocket)
  (geometric pocket finder), bound HETATM ligands/cofactors, or click atoms
  in the 3D viewer to build a box manually
- Load ligands from pasted SMILES, or upload `.smi`/`.sdf`/`.mol`/`.csv`/`.tsv`
- Dock many ligands concurrently (one docking process per ligand, thread pool)
- Pick smina or gnina as the engine per job; gnina can use a CUDA GPU when
  one's present (Linux/WSL2 only — no Windows build upstream)
- Optional pharmacophoric constraints (H-bond donor/acceptor, or any-atom
  "positional") checked against each pose after docking and used to rerank
  poses toward the ones that satisfy them — smina/gnina have no way to bias
  their own search this way, so this is post-hoc filtering/reranking, not a
  constrained search (see `backend/pharmacophore.py`)
- Self-docking validation: redock a receptor's own bound ligand and get
  each pose's RMSD to its real crystallographic position, read straight off
  the PDB's HETATM records and computed with `obrms` (bundled alongside
  smina/gnina — see `backend/rmsd.py`)
- Live job progress, per-ligand poses viewable in-browser, CSV/ZIP export
- Job history survives a server restart (state reloaded from `runs/jobs/`)

## Setup

Needs [Node.js](https://nodejs.org/) (any recent LTS) on `PATH` in addition to
Python 3 — it's used once, at setup time, to bundle the Mol* viewer.

macOS / Linux:
```bash
./setup.sh   # installs smina, fpocket (+ gnina on Linux) via micromamba
             # into .env-tools/, a Python venv into .venv/ with fastapi,
             # rdkit, etc., and builds frontend/app.js (see below).
```

Windows (native, CPU-only — no gnina/GPU build exists for Windows; run
inside WSL2 with `setup.sh` instead if you need GPU docking on a Windows
box):
```powershell
.\setup.ps1   # same idea, using micromamba's win-64 build and a
              # Scripts\python.exe venv
```

Safe to re-run; existing installs are reused. Supported platforms: macOS
(arm64/x86_64), Linux (x86_64/aarch64, gnina needs x86_64), Windows (x86_64,
smina only).

## Run

```bash
./run.sh                 # serves on http://127.0.0.1:8000
PORT=9000 ./run.sh        # pick a different port
```

```powershell
.\run.ps1                 # Windows equivalent
$env:PORT=9000; .\run.ps1  # pick a different port
```

Open the printed URL in a browser. The top-right pill shows whether smina
and/or gnina were found; if not, re-run `./setup.sh` (or `setup.ps1`).

## Sharing with a remote colleague

By default the app only listens on `127.0.0.1` — nothing outside your
machine can reach it. To hand a colleague a link:

1. Turn on auth so the tunnel link isn't wide open:
   ```bash
   DOCKING_REQUIRE_AUTH=1 ./run.sh
   ```
   Prints a generated username/password to the console on start. (Or set
   `DOCKING_AUTH_PASS` yourself for a password that survives restarts.)
2. In another terminal, open a tunnel — no router/port-forward needed:
   ```bash
   brew install cloudflared
   cloudflared tunnel --url http://127.0.0.1:8000
   ```
   Prints a `https://<random>.trycloudflare.com` URL. Share that + the
   username/password with your colleague; the browser will prompt for
   credentials on first load.
3. `Ctrl+C` the tunnel when you're done — the link stops working immediately.

## Testing

```bash
.venv/bin/python -m pytest tests/ -q
```

Tests cover receptor cleanup/HET detection, ligand parsing, smina/gnina
command construction, and job persistence — they mock the subprocess call,
so they run without smina or gnina installed.

## Frontend viewer build

The 3D viewer is the full [Mol*](https://molstar.org/) app itself (its own
sequence panel, right-hand controls, log, selection/measurement tools —
same as molstar.org, not a stripped-down embed), imported as a real ES
module in `frontend/src/app.js` and bundled with esbuild into
`frontend/app.js`, the plain script `index.html` loads. App-specific bits
(fetch/upload receptor, click-a-HETATM-to-box, ligand list, run docking,
job results) are layered on top via `plugin.behaviors.interaction.click`
and Mol*'s structure builders; representation/style/surface/waters are
left to Mol*'s own panel rather than duplicated in our UI. `./setup.sh` / `setup.ps1` run this for you;
to rebuild by hand after editing `src/app.js`:
```bash
cd frontend
npm install     # first time only
npm run build   # one-off bundle
npm run watch    # rebuild on every save while developing
```

## Project layout

```
backend/
  main.py           FastAPI app: receptor/ligand/job HTTP routes
  jobs.py           Job store: one job -> many ligands, each docked in a thread
  docking.py        Thin wrapper around the smina/gnina CLI
  prep.py           Receptor cleanup, HET-ligand grouping, ligand parsing/3D embed
  binding_sites.py  fpocket wrapper for automated binding site detection
  config.py         Paths and runtime settings (env-var overridable)
frontend/
  index.html, style.css     Static single-page UI shell
  src/app.js                Real source (imports Mol* directly) - edit this
  app.js, vendor/            Built by `npm run build`; gitignored, what index.html
                             actually loads (app.js is the esbuild bundle of
                             src/app.js; vendor/molstar.css is copied from the
                             molstar package)
  package.json, build.mjs   esbuild bundling (see "Frontend viewer build")
tests/        pytest suite for the backend
runs/         receptor uploads + job outputs (gitignored; created at runtime)
```

## Notable env vars

| Variable            | Default                  | Purpose                              |
|----------------------|---------------------------|---------------------------------------|
| `PORT`               | `8000`                    | `run.sh` listen port                  |
| `DOCKING_RUNS_DIR`   | `<repo>/runs`             | where receptors/job outputs are stored |
| `DOCKING_WORKERS`    | `cpu_count - 1`           | max concurrent docking processes      |
| `SMINA_TIMEOUT`      | `1800` (seconds)          | per-ligand docking hard timeout       |
| `SMINA_BIN`          | (unset)                   | explicit path to a smina executable   |
| `GNINA_BIN`          | (unset)                   | explicit path to a gnina executable   |
| `OBRMS_BIN`          | (unset)                   | explicit path to an obrms executable (self-docking RMSD) |
| `DOCKING_MAX_LIGANDS`| `50`                      | max ligands accepted per job          |
| `DOCKING_REQUIRE_AUTH`| (unset)                  | set to anything to turn on HTTP Basic auth with an auto-generated password (printed to console on start) |
| `DOCKING_AUTH_USER`  | `colleague`               | Basic auth username                   |
| `DOCKING_AUTH_PASS`  | (unset)                   | Basic auth password; set this instead of `DOCKING_REQUIRE_AUTH` for a stable password across restarts |

## Known limitations

- mmCIF receptors are rejected — upload a `.pdb` file
- No authentication; intended for local/single-user use
- `runs/` is not pruned automatically — old jobs accumulate on disk
