# Docking Bench

Local web app for small-molecule docking against a protein receptor. It comes
with two interchangeable front ends over the same backend:

- **Streamlit app** (`./run_streamlit.sh`) — the easiest way to run it on a
  shared machine and hand colleagues a URL. See "Streamlit app" below.
- **Original web page** (`./run.sh`) — FastAPI + a static Mol* page.

Three docking engines, picked per job:

- [smina](https://sourceforge.net/projects/smina/) (CPU)
- [gnina](https://github.com/gnina/gnina) (CNN scoring, CUDA GPU on Linux)
- [MOE](https://www.chemcomp.com/) Dock, driven headless through `moebatch`
  (needs a licensed MOE install; optional)

[RDKit](https://www.rdkit.org/) builds ligand 3D structures. Runs entirely on
your machine — no external services beyond fetching PDB files and chemical
component SMILES from RCSB. The 3D viewer is [Mol*](https://molstar.org/),
bundled locally at setup time (see "Frontend viewer build" below).

## Features

### All engines

- Fetch a receptor by PDB ID (RCSB) or upload a `.pdb` file
- Set the docking site from a bound HETATM ligand/cofactor, or click atoms in
  the 3D viewer to build a box manually
- Receptor cleanup before docking: delete the bound ligand(s), keep or drop
  waters, delete whole chains (e.g. a peptide ligand chain), drop extra
  monomers far from the site
- Load ligands from pasted SMILES, or upload `.smi`/`.sdf`/`.mol`/`.csv`/`.tsv`;
  "self-dock" a bound ligand straight from its PDB chemical component SMILES
- Optional pharmacophoric constraints (H-bond donor/acceptor, aromatic,
  metal, SMARTS or any-atom "positional") checked against each pose after
  docking and used to rerank poses — post-hoc filtering/reranking, not a
  constrained search (see `backend/pharmacophore.py`)
- Self-docking validation: each pose's RMSD to the receptor's own bound
  ligand, computed with `obrms` (see `backend/rmsd.py`)
- Ligand strain energy per pose (MMFF94, see `backend/strain.py`)
- 2D interaction diagram per pose: MOE's Ligand Interactions diagram when MOE
  is installed, otherwise a ProLIF-based one (switchable in the panel); 3D
  non-covalent contacts drawn in the viewer
- Live job progress, per-ligand poses viewable and steppable in-browser,
  CSV/ZIP export
- Job history survives a server restart (state reloaded from `runs/jobs/`)

### smina / gnina

- Many ligands dock concurrently (one process per ligand, thread pool)
- gnina can use a CUDA GPU when one's present (Linux/WSL2 only — no Windows
  build upstream)
- Optional extra `--minimize` pass on every pose after docking

### MOE

- Placement (Triangle Matcher, Alpha PMI, Pharmacophore, ...), refinement
  (Rigid Receptor, Induced Fit, ...) and rescoring (GBVI/WSA dG, London dG,
  ...) chosen per job, with MOE's Dock panel defaults; every pose's MOE
  scores (S, E_place, E_refine, ...) shown in a collapsible table
- **QuickPrep** button in the receptor panel: runs MOE QuickPrep (Structure
  Preparation, Protonate3D, distant-water removal, tethered minimization
  around the selected site ligand) with your saved MOE QuickPrep settings.
  The prepared receptor replaces the loaded one in the viewer, with a report
  of what changed (hydrogens added, atoms moved, protonation states chosen —
  site residues first), and is what gets docked
- Covalent docking against a chosen Cys/Ser/Thr/Tyr/Lys, using MOE's built-in
  reaction templates or an uploaded `.rxn` file; warns when a ligand isn't
  in the reactant form the reaction needs
- Pharmacophore-constrained placement from an uploaded MOE `.ph4` file
- Ligand Interactions diagrams are drawn during the docking run itself, so
  they cost no extra license checkout
- License-aware queue: MOE runs are capped server-wide
  (`MOE_MAX_CONCURRENT`, default 1) so ligands wait as "queued" instead of
  failing on license errors, and runs retry while other users hold the
  license tokens

MOE internals (how runs are scripted, QuickPrep, covalent docking, SVL
pitfalls, debugging against the GUI) are documented in [docs/MOE.md](docs/MOE.md).

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

MOE is not installed by setup. If you have it, the app finds `moebatch` under
`/Applications/moe*/bin/` or on `PATH`, or set `MOEBATCH_BIN` to its path.
Without MOE the MOE engine option is simply unavailable.

## Run

```bash
./run.sh                 # serves on http://127.0.0.1:8000
PORT=9000 ./run.sh        # pick a different port
```

```powershell
.\run.ps1                 # Windows equivalent
$env:PORT=9000; .\run.ps1  # pick a different port
```

Open the printed URL in a browser. The top-right pill shows which engines
(smina, gnina, MOE) were found; if smina/gnina are missing, re-run
`./setup.sh` (or `setup.ps1`).

The server doesn't reload code changes — restart it after pulling updates.

## Streamlit app

Same workflow as the web page — receptor, active site (click in the 3D
viewer), constraints, covalent docking, ligands, parameters, live results,
2D interaction diagrams, downloads — served by Streamlit:

```bash
./run_streamlit.sh                        # http://127.0.0.1:8501, this machine only
DOCKING_AUTH_PASS=secret HOST=0.0.0.0 ./run_streamlit.sh
                                          # colleagues open http://<this-machine>:8501
```

With `DOCKING_AUTH_PASS` (or `DOCKING_REQUIRE_AUTH=1`) set, the app asks for
the user (`DOCKING_AUTH_USER`, default `colleague`) and password before
showing anything. Every colleague gets their own session; jobs are shared on
the server (the Results card's "previous jobs" list shows them all) and
MOE runs queue across everyone (`MOE_MAX_CONCURRENT`). Settings such as the
upload limit and theme live in `.streamlit/config.toml`.

The 3D viewer is a Streamlit custom component wrapping the same full Mol*
app (`streamlit_app/molstar_component/`, built by `./setup.sh`; see its
README). Clicks in the viewer come back to Python, so site picking,
constraint points and the covalent residue all work by clicking atoms, as
on the web page.

How it fits together: `backend/services.py` holds every operation (shared by
both front ends), `streamlit_app/logic.py` holds the UI decisions (what a
click means, building the job request — unit-tested), and
`streamlit_app/app.py`, `setup_panel.py`, `results_panel.py` lay out the page.

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

## Hosting on Streamlit Community Cloud

Gives a permanent `https://<name>.streamlit.app` link guarded by a password,
with nothing running on your machine. smina docking only (no MOE licence and
no GPU in the cloud).

1. Push this repo to GitHub (private repos work). The cloud reads:
   - `requirements.txt` - Python packages
   - `packages.txt` - apt packages: `smina`, `openbabel` (for `obrms`)
   - `streamlit_app/molstar_component/frontend/dist/` - the prebuilt viewer.
     Rebuild it (`npm run build` in that folder) and commit after changing
     the viewer source; the cloud does not run npm.
2. Go to https://share.streamlit.io, sign in with GitHub, **Create app**:
   repo `daga-pankaj/Dock-Bench`, branch `remodel`, main file
   `streamlit_app/app.py`. Under **Advanced settings** pick Python 3.12 and
   paste the secrets (top-level secrets become environment variables):
   ```toml
   DOCKING_AUTH_USER = "colleague"
   DOCKING_AUTH_PASS = "pick-a-long-password"
   DOCKING_MAX_LIGANDS = "20"
   ```
   Always set `DOCKING_AUTH_PASS` - without it the app is open to anyone
   with the link.
3. Share the URL plus user/password. Change the password any time under the
   app's **Settings - Secrets**; the app restarts with the new one.

Limits: the free tier has about 1 GB RAM and a couple of shared CPUs, so
keep jobs small. `runs/` is not persistent - results vanish when the app
restarts or sleeps (after a few days idle), so download what you need.

## Testing

```bash
.venv/bin/python -m pytest tests/ -q
```

Tests cover receptor cleanup/HET detection, ligand parsing, smina/gnina
command construction, MOE script generation and license queueing, job
persistence, the Streamlit UI logic, and an end-to-end Streamlit run
(`tests/test_streamlit_app.py`, skipped without smina). Most mock the
subprocess call, so they run without smina, gnina or MOE installed.

The viewer component has its own unit tests:
`cd streamlit_app/molstar_component/frontend && npm test`.

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
  services.py             Every operation (receptors, ligands, jobs, diagrams,
                          downloads) - shared by both front ends
  schemas.py              Request models (job request, constraints, ...)
  main.py                 FastAPI app: thin HTTP routes over services.py
  jobs.py                 Job store: one job -> many ligands, each docked in a thread
  docking.py              Engine dispatch; smina/gnina CLI wrapper, MOE pose parsing
  moe_dock.py             moebatch wrapper: Dock, QuickPrep, Ligand Interactions,
                          license queue
  moe_dock_template.svl   SVL script template for a MOE docking run
  prep.py                 Receptor cleanup, HET-ligand grouping, ligand parsing/3D
                          embed, QuickPrep report
  covalent_check.py       Ligand reactant-form checks for covalent docking
  pharmacophore.py        Post-docking pharmacophoric constraints
  rmsd.py                 Self-docking RMSD via obrms
  strain.py               Ligand strain energy (MMFF94)
  interactions_prolif.py  ProLIF interaction detection
  interaction_render.py   ProLIF-based 2D interaction diagram (SVG)
  config.py               Paths and runtime settings (env-var overridable)
frontend/
  index.html, style.css     Static single-page UI shell
  src/app.js                Real source (imports Mol* directly) - edit this
  app.js, vendor/            Built by `npm run build`; gitignored, what index.html
                             actually loads (app.js is the esbuild bundle of
                             src/app.js; vendor/molstar.css is copied from the
                             molstar package)
  package.json, build.mjs   esbuild bundling (see "Frontend viewer build")
streamlit_app/
  app.py                    Streamlit entry point (page, auth, viewer events)
  setup_panel.py            Left column: receptor ... docking parameters
  results_panel.py          Right column: viewer, 2D diagram, results
  logic.py                  UI decisions, Streamlit-free and unit-tested
  molstar_component/        Mol* viewer as a Streamlit custom component
.streamlit/config.toml      Streamlit settings (theme, upload limit)
docs/         MOE.md (MOE integration notes); presentation/ slides
tests/        pytest suite for the backend and the Streamlit app
runs/         receptor uploads + job outputs (gitignored; created at runtime)
```

## Notable env vars

| Variable            | Default                  | Purpose                              |
|----------------------|---------------------------|---------------------------------------|
| `PORT`               | `8000` / `8501`           | `run.sh` / `run_streamlit.sh` listen port |
| `HOST`               | `127.0.0.1`               | `run_streamlit.sh` bind address (`0.0.0.0` to share on the network) |
| `DOCKING_RUNS_DIR`   | `<repo>/runs`             | where receptors/job outputs are stored |
| `DOCKING_WORKERS`    | `cpu_count - 1`           | max concurrent docking processes      |
| `SMINA_TIMEOUT`      | `1800` (seconds)          | per-ligand docking hard timeout       |
| `SMINA_BIN`          | (unset)                   | explicit path to a smina executable   |
| `GNINA_BIN`          | (unset)                   | explicit path to a gnina executable   |
| `MOEBATCH_BIN`       | (unset)                   | explicit path to MOE's `moebatch`     |
| `MOE_MAX_CONCURRENT` | `1`                       | MOE runs allowed at once, server-wide (each checks out license tokens) |
| `MOE_LICENSE_WAIT`   | `3600` (seconds)          | how long a MOE run keeps retrying while no license tokens are free |
| `OBRMS_BIN`          | (unset)                   | explicit path to an obrms executable (self-docking RMSD) |
| `DOCKING_MAX_LIGANDS`| `50`                      | max ligands accepted per job          |
| `DOCKING_REQUIRE_AUTH`| (unset)                  | set to anything to turn on HTTP Basic auth with an auto-generated password (printed to console on start) |
| `DOCKING_AUTH_USER`  | `colleague`               | Basic auth username                   |
| `DOCKING_AUTH_PASS`  | (unset)                   | Basic auth password (web page) / login password (Streamlit); set this instead of `DOCKING_REQUIRE_AUTH` for a stable password across restarts |

## Known limitations

- mmCIF receptors are rejected — upload a `.pdb` file
- Intended for local/trusted use; authentication is optional HTTP Basic only
  (see "Sharing with a remote colleague")
- Ligands are docked in the protonation state their SMILES/SDF gives (no
  Wash/pH 7 step) — write charged groups explicitly, e.g. `[NH3+]`
- `runs/` is not pruned automatically — old jobs accumulate on disk
