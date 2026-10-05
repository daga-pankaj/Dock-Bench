"""FastAPI app: receptor handling, ligand parsing, docking jobs, downloads.

Run with:  uvicorn backend.main:app --reload --port 8000
"""

from __future__ import annotations

import base64
import csv
import hmac
import io
import json
import os
import re
import urllib.error
import urllib.request
import uuid
import zipfile

# python.org's macOS Framework build ships without a CA bundle wired up, so
# urllib.request.urlopen fails with CERTIFICATE_VERIFY_FAILED until either
# "Install Certificates.command" is run or SSL_CERT_FILE points somewhere
# real. certifi carries its own bundle - use it as a fallback so RCSB fetch
# works out of the box on a fresh venv.
if not os.environ.get("SSL_CERT_FILE"):
    try:
        import certifi

        os.environ["SSL_CERT_FILE"] = certifi.where()
    except ImportError:
        pass

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from . import config, docking, prep
from . import interaction_render
from . import interactions_prolif as iprolif
from . import pharmacophore as pharm
from .jobs import MANAGER, JobParams

config.ensure_dirs()

if config.AUTH_GENERATED:
    print(f"\n>>> auth enabled - user: {config.AUTH_USER}  password: {config.AUTH_PASS}\n")

app = FastAPI(title="Docking Bench", version="1.0")

PDB_ID_RE = re.compile(r"^[0-9A-Za-z]{4}$")
COMP_ID_RE = re.compile(r"^[0-9A-Za-z]{1,5}$")
MAX_BOX = 60.0


@app.middleware("http")
async def basic_auth(request, call_next):
    """Gate every route behind HTTP Basic when DOCKING_AUTH_PASS is set.

    Off by default (plain localhost use); turn on before sharing a tunnel
    link with anyone outside your machine.
    """
    if config.AUTH_PASS is None:
        return await call_next(request)
    header = request.headers.get("authorization", "")
    ok = False
    if header.startswith("Basic "):
        try:
            user, _, pw = base64.b64decode(header[6:]).decode().partition(":")
            ok = hmac.compare_digest(user, config.AUTH_USER) and hmac.compare_digest(pw, config.AUTH_PASS)
        except Exception:
            ok = False
    if not ok:
        return Response(
            status_code=401,
            headers={"WWW-Authenticate": 'Basic realm="Docking Bench"'},
        )
    return await call_next(request)


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


class FetchRequest(BaseModel):
    pdb_id: str


class LigandInput(BaseModel):
    name: str
    smiles: str | None = None
    molblock: str | None = None


class ParseRequest(BaseModel):
    text: str


class ConstraintIn(BaseModel):
    """One pharmacophoric constraint: an H-bond donor/acceptor, an aromatic
    ring (pi-stacking proxy), a metal-coordinating atom, any heavy atom
    ("positional"), or any atom matching a user SMARTS ("smarts", needs the
    `smarts` field) - checked for proximity to `radius` of `point` (a
    receptor-atom or hand-picked coordinate). Checked against each pose after
    docking, not during smina/gnina's own search - see
    backend/pharmacophore.py."""

    kind: str
    point: tuple[float, float, float]
    radius: float = Field(3.5, ge=0.5, le=15.0)
    required: bool = True
    weight: float = Field(1.0, ge=0.0, le=20.0)
    label: str = ""
    smarts: str | None = None  # only used (and required) when kind == "smarts"


class JobRequest(BaseModel):
    receptor_id: str
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    ligands: list[LigandInput]
    exhaustiveness: int = Field(8, ge=1, le=64)
    num_modes: int = Field(config.DEFAULT_NUM_MODES, ge=1, le=20)
    energy_range: float = Field(4.0, ge=0.5, le=20.0)
    seed: int = 42
    scoring: str = "vina"
    engine: str = "smina"
    use_gpu: bool = True
    minimize_after_docking: bool = False
    keep_waters: bool = False
    remove_ligand: bool = True
    exclude_chains: list[str] = Field(default_factory=list)
    trim_distant_chains: bool = False
    chain_cutoff: float = Field(15.0, ge=5.0, le=60.0)
    constraints: list[ConstraintIn] = Field(default_factory=list)
    # Self-docking validation: the key (from het_ligands, resname:chain:resseq)
    # of a bound ligand to compare each docked pose against by RMSD (see
    # backend/rmsd.py). Read from this job's receptor before it's cleaned, so
    # it still works with remove_ligand=True stripping that same ligand out
    # of the receptor actually used for docking.
    reference_ligand_key: str | None = None
    # MOE-specific parameters
    moe_placement: str = "Triangle Matcher"
    moe_refinement: str = "Rigid Receptor"
    ph4_id: str | None = None  # MOE pharmacophore file ID from /api/pharmacophore/upload
    # Covalent docking (MOE-only)
    covalent_enabled: bool = False
    covalent_reaction: str | None = None  # e.g., "michael_acceptor_s" (built-in template)
    covalent_residue: str | None = None   # e.g., "A:CYS:189"
    custom_rxn_id: str | None = None      # Custom .rxn file ID from /api/reactions/upload


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _receptor_path(receptor_id: str):
    if not re.fullmatch(r"[0-9a-f]{8,32}", receptor_id):
        raise HTTPException(400, "bad receptor id")
    path = config.RECEPTOR_DIR / f"{receptor_id}.pdb"
    if not path.exists():
        raise HTTPException(404, "receptor not found; upload it again")
    return path


def _store_receptor(name: str, text: str) -> dict:
    if "ATOM" not in text and "HETATM" not in text:
        raise HTTPException(400, "file contains no ATOM/HETATM records")
    receptor_id = uuid.uuid4().hex[:16]
    (config.RECEPTOR_DIR / f"{receptor_id}.pdb").write_text(text)
    hets = prep.list_het_ligands(text)
    n_protein = sum(1 for tag, _ in prep._atom_records(text) if tag == "ATOM")
    return {
        "id": receptor_id,
        "name": name,
        "n_protein_atoms": n_protein,
        "het_ligands": [h.as_dict() for h in hets],
        "chains": [c.as_dict() for c in prep.list_chains(text)],
    }


def _safe_name(name: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9._+-]{1,120}", name):
        raise HTTPException(400, "bad ligand name")
    return name


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------


@app.get("/api/status")
def status() -> dict:
    def version_of(engine: str, exe: str | None) -> str | None:
        if not exe:
            return None
        try:
            return docking.engine_version(engine)
        except Exception as exc:
            return f"error: {exc}"

    smina_exe = config.find_smina()
    gnina_exe = config.find_gnina()
    moebatch_exe = config.find_moebatch()
    return {
        "smina": smina_exe,
        "smina_version": version_of("smina", smina_exe),
        "gnina": gnina_exe,
        "gnina_version": version_of("gnina", gnina_exe),
        "moebatch": moebatch_exe,
        "moe_available": moebatch_exe is not None,
        "gpu_available": config.has_nvidia_gpu(),
        "workers": config.MAX_WORKERS,
        "scoring_functions": list(docking.SCORING_FUNCTIONS),
        "moe_scoring_functions": list(docking.MOE_SCORING_FUNCTIONS),
        "default_num_modes": config.DEFAULT_NUM_MODES,
        "constraint_kinds": list(pharm.KINDS),
        "engines": list(docking.ENGINES),
        "covalent_reactions": config.list_covalent_reactions(),
    }


# ---------------------------------------------------------------------------
# Receptor
# ---------------------------------------------------------------------------


@app.post("/api/receptor/upload")
async def upload_receptor(file: UploadFile = File(...)) -> dict:
    raw = await file.read()
    if len(raw) > 60_000_000:
        raise HTTPException(413, "file too large")
    try:
        text = raw.decode("utf-8", errors="replace")
    except Exception:
        raise HTTPException(400, "cannot decode file as text")
    if file.filename and file.filename.lower().endswith(".cif"):
        raise HTTPException(400, "mmCIF not supported; upload a .pdb file")
    return _store_receptor(file.filename or "receptor.pdb", text)


@app.post("/api/pharmacophore/upload")
async def upload_pharmacophore(file: UploadFile = File(...)) -> dict:
    """Upload a MOE .ph4 pharmacophore file for docking constraints."""
    raw = await file.read()
    if len(raw) > 1_000_000:
        raise HTTPException(413, "file too large")
    try:
        text = raw.decode("utf-8", errors="replace")
    except Exception:
        raise HTTPException(400, "cannot decode file as text")

    # Basic validation: check for MOE .ph4 header
    if "#moe:ph4que" not in text and "#pharmacophore" not in text:
        raise HTTPException(400, "not a valid MOE .ph4 file")

    # Store in runs/pharmacophores/
    import uuid
    ph4_id = uuid.uuid4().hex[:12]
    ph4_dir = config.JOB_DIR.parent / "pharmacophores"
    ph4_dir.mkdir(exist_ok=True)
    ph4_path = ph4_dir / f"{ph4_id}.ph4"
    ph4_path.write_text(text)

    return {"ph4_id": ph4_id, "filename": file.filename or "pharmacophore.ph4"}


@app.post("/api/reactions/upload")
async def upload_reaction(file: UploadFile = File(...)) -> dict:
    """Upload a custom .rxn reaction file for covalent docking."""
    raw = await file.read()
    if len(raw) > 500_000:
        raise HTTPException(413, "file too large")
    try:
        text = raw.decode("utf-8", errors="replace")
    except Exception:
        raise HTTPException(400, "cannot decode file as text")

    # Basic validation: check for $RXN header (MDL RXN format)
    if not text.strip().startswith("$RXN"):
        raise HTTPException(400, "not a valid RXN file (missing $RXN header)")

    # Store in runs/custom_reactions/
    rxn_id = uuid.uuid4().hex[:12]
    rxn_dir = config.JOB_DIR.parent / "custom_reactions"
    rxn_dir.mkdir(exist_ok=True)
    rxn_path = rxn_dir / f"{rxn_id}.rxn"
    rxn_path.write_text(text)

    return {"rxn_id": rxn_id, "filename": file.filename or "custom.rxn"}


@app.post("/api/receptor/fetch")
def fetch_receptor(req: FetchRequest) -> dict:
    pdb_id = req.pdb_id.strip()
    if not PDB_ID_RE.match(pdb_id):
        raise HTTPException(400, "PDB ID must be 4 characters")
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            text = resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        raise HTTPException(404, f"RCSB returned {exc.code} for {pdb_id.upper()}")
    except Exception as exc:
        raise HTTPException(502, f"could not reach RCSB: {exc}")
    return _store_receptor(f"{pdb_id.upper()}.pdb", text)


@app.get("/api/receptor/{receptor_id}/pdb")
def receptor_pdb(receptor_id: str) -> Response:
    return PlainTextResponse(_receptor_path(receptor_id).read_text())


@app.get("/api/receptor/{receptor_id}/residues/{resname}")
def list_residues(receptor_id: str, resname: str) -> dict:
    """List all residues of given type (e.g., CYS, SER) in receptor."""
    if not re.fullmatch(r"[A-Z]{3}", resname.upper()):
        raise HTTPException(400, "resname must be 3-letter code")
    text = _receptor_path(receptor_id).read_text()
    residues = prep.list_residues_by_name(text, resname.upper())
    return {"resname": resname.upper(), "residues": residues}


@app.get("/api/pdb_components/{comp_id}/smiles")
def pdb_component_smiles(comp_id: str) -> dict:
    """A bound ligand's own PDB Chemical Component Dictionary SMILES, for the
    "self-dock" button - lets a HETATM group's residue code (e.g. STI) become
    a ready-to-dock ligand without needing its bond orders/aromaticity
    guessed back from bare coordinates (unreliable - see backend/rmsd.py's
    module docstring for how RMSD validation sidesteps that problem)."""
    comp_id = comp_id.strip().upper()
    if not COMP_ID_RE.match(comp_id):
        raise HTTPException(400, "bad chemical component id")
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise HTTPException(404, f"RCSB has no chemical component {comp_id} ({exc.code})")
    except Exception as exc:
        raise HTTPException(502, f"could not reach RCSB: {exc}")

    descriptors = data.get("pdbx_chem_comp_descriptor", [])
    smiles = next((d["descriptor"] for d in descriptors if d.get("type") == "SMILES_CANONICAL"), None)
    if smiles is None:
        smiles = next((d["descriptor"] for d in descriptors if d.get("type") == "SMILES"), None)
    if smiles is None:
        raise HTTPException(404, f"RCSB has no SMILES on file for {comp_id}")
    return {"comp_id": comp_id, "smiles": smiles, "name": data.get("chem_comp", {}).get("name", comp_id)}


# ---------------------------------------------------------------------------
# Ligands
# ---------------------------------------------------------------------------


def _preview(specs: list[prep.LigandSpec]) -> dict:
    from rdkit import Chem

    items = []
    for spec in specs:
        valid, note = True, ""
        if spec.molblock:
            note = "3D from file"
        elif spec.smiles:
            mol = Chem.MolFromSmiles(spec.smiles)
            if mol is None:
                valid, note = False, "invalid SMILES"
            else:
                note = f"{mol.GetNumAtoms()} heavy atoms"
        else:
            valid, note = False, "empty"
        items.append(
            {
                "name": spec.name,
                "smiles": spec.smiles,
                "molblock": spec.molblock,
                "valid": valid,
                "note": note,
            }
        )
    return {"ligands": items, "n_valid": sum(1 for i in items if i["valid"])}


@app.post("/api/ligands/parse_text")
def parse_ligand_text(req: ParseRequest) -> dict:
    text = req.text
    # An SDF pasted into the box still starts with a mol block, not a SMILES.
    specs = prep.parse_sdf_block(text) if "$$$$" in text else prep.parse_smiles_block(text)
    return _preview(specs)


@app.post("/api/ligands/parse_file")
async def parse_ligand_file(file: UploadFile = File(...)) -> dict:
    raw = await file.read()
    if len(raw) > 200_000_000:
        raise HTTPException(413, "file too large")
    specs = prep.parse_ligand_file(file.filename or "ligands.smi", raw.decode("utf-8", errors="replace"))
    return _preview(specs)


# ---------------------------------------------------------------------------
# Jobs
# ---------------------------------------------------------------------------


@app.post("/api/jobs")
def create_job(req: JobRequest) -> dict:
    if req.engine not in docking.ENGINES:
        raise HTTPException(400, f"unknown docking engine {req.engine}")

    # Check engine availability
    if req.engine == "moe":
        engine_exe = config.find_moebatch()
        if not engine_exe:
            raise HTTPException(503, "MOE not installed")
    elif req.engine == "gnina":
        engine_exe = config.find_gnina()
        if not engine_exe:
            raise HTTPException(503, "gnina not installed - run ./setup.sh")
    else:
        engine_exe = config.find_smina()
        if not engine_exe:
            raise HTTPException(503, "smina not installed - run ./setup.sh")

    if not req.ligands:
        raise HTTPException(400, "no ligands submitted")
    if len(req.ligands) > config.MAX_LIGANDS_PER_JOB:
        raise HTTPException(400, f"too many ligands in one job (max {config.MAX_LIGANDS_PER_JOB})")

    # Validate scoring function based on engine
    if req.engine == "moe":
        if req.scoring not in docking.MOE_SCORING_FUNCTIONS:
            raise HTTPException(400, f"unknown MOE scoring function {req.scoring}")
    else:
        if req.scoring not in docking.SCORING_FUNCTIONS:
            raise HTTPException(400, f"unknown scoring function {req.scoring}")
    for value in req.size:
        if not (4.0 <= value <= MAX_BOX):
            raise HTTPException(400, f"box edge {value} A outside 4-{MAX_BOX} A")
    for c in req.constraints:
        if c.kind not in pharm.KINDS:
            raise HTTPException(400, f"unknown constraint kind {c.kind}")
        try:
            pharm.Constraint.from_dict(c.model_dump())
        except ValueError as exc:
            raise HTTPException(400, str(exc)) from exc

    # Covalent docking validation
    if req.covalent_enabled:
        if req.engine != "moe":
            raise HTTPException(400, "covalent docking requires MOE engine")
        if not req.covalent_residue:
            raise HTTPException(400, "covalent residue required (e.g., 'A:CYS:189')")
        # Either built-in reaction or custom .rxn required
        if not req.covalent_reaction and not req.custom_rxn_id:
            raise HTTPException(400, "covalent reaction template or custom .rxn file required")
        # Validate built-in reaction if specified
        if req.covalent_reaction:
            available_reactions = config.list_covalent_reactions()
            if req.covalent_reaction not in available_reactions:
                raise HTTPException(400, f"unknown reaction template {req.covalent_reaction}")
        # Validate custom reaction file exists
        if req.custom_rxn_id:
            custom_rxn_path = config.RUNS_DIR / "custom_reactions" / f"{req.custom_rxn_id}.rxn"
            if not custom_rxn_path.exists():
                raise HTTPException(400, f"custom reaction file not found - re-upload")

    receptor_text = _receptor_path(req.receptor_id).read_text()

    # Validate covalent residue exists in receptor
    if req.covalent_enabled and req.covalent_residue:
        if not prep.validate_residue_exists(receptor_text, req.covalent_residue):
            # Extract resname to suggest available alternatives
            parsed = prep.parse_residue_spec(req.covalent_residue)
            suggestions = ""
            if parsed:
                _, resname, _ = parsed
                available = prep.list_residues_by_name(receptor_text, resname)
                if available:
                    suggestions = f" - available {resname} residues: {', '.join(available[:10])}"
            raise HTTPException(400, f"residue {req.covalent_residue} not found in receptor{suggestions}")

    reference_ligand_pdb = None
    if req.reference_ligand_key:
        reference_ligand_pdb = prep.het_ligand_pdb_block(receptor_text, req.reference_ligand_key)
        if reference_ligand_pdb is None:
            raise HTTPException(400, f"bound ligand {req.reference_ligand_key} not found in receptor")

    cleaned = prep.clean_receptor(
        receptor_text,
        keep_waters=req.keep_waters,
        strip_ligand=req.remove_ligand,
        exclude_chains=req.exclude_chains,
        site_center=tuple(req.center) if req.trim_distant_chains else None,
        chain_cutoff=req.chain_cutoff,
    )
    if not cleaned.strip().splitlines()[:-1]:
        raise HTTPException(400, "receptor is empty after cleanup")

    specs = [
        prep.LigandSpec(name=lig.name, smiles=lig.smiles, molblock=lig.molblock)
        for lig in req.ligands
    ]
    params = JobParams(
        center=tuple(req.center),
        size=tuple(req.size),
        exhaustiveness=req.exhaustiveness,
        num_modes=req.num_modes,
        energy_range=req.energy_range,
        seed=req.seed,
        scoring=req.scoring,
        engine=req.engine,
        use_gpu=req.use_gpu,
        minimize_after_docking=req.minimize_after_docking,
        constraints=[c.model_dump() for c in req.constraints],
        reference_ligand_pdb=reference_ligand_pdb,
        moe_placement=req.moe_placement,
        moe_refinement=req.moe_refinement,
        ph4_id=req.ph4_id,
        covalent_enabled=req.covalent_enabled,
        covalent_reaction=req.covalent_reaction,
        covalent_residue=req.covalent_residue,
        custom_rxn_id=req.custom_rxn_id,
    )
    job = MANAGER.submit(cleaned, f"receptor_{req.receptor_id}", specs, params)
    return job.as_dict()


@app.get("/api/jobs")
def list_jobs() -> dict:
    return {"jobs": MANAGER.list_jobs()}


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str) -> dict:
    job = MANAGER.get(job_id)
    if not job:
        raise HTTPException(404, "job not found")
    return job.as_dict()


@app.post("/api/jobs/{job_id}/cancel")
def cancel_job(job_id: str) -> dict:
    if not MANAGER.cancel(job_id):
        raise HTTPException(400, "job cannot be cancelled")
    return {"ok": True}


@app.get("/api/jobs/{job_id}/receptor.pdb")
def job_receptor(job_id: str) -> Response:
    job = MANAGER.get(job_id)
    if not job or not job.receptor_path.exists():
        raise HTTPException(404, "job receptor not found")
    return PlainTextResponse(job.receptor_path.read_text())


@app.get("/api/jobs/{job_id}/ligands/{name}/poses.sdf")
def ligand_poses(job_id: str, name: str) -> Response:
    job = MANAGER.get(job_id)
    if not job:
        raise HTTPException(404, "job not found")
    path = job.poses_path(_safe_name(name))
    if not path.exists():
        raise HTTPException(404, "no poses for this ligand")
    return FileResponse(path, media_type="chemical/x-mdl-sdfile", filename=path.name)


@app.get("/api/jobs/{job_id}/ligands/{name}/log")
def ligand_log(job_id: str, name: str) -> Response:
    job = MANAGER.get(job_id)
    if not job:
        raise HTTPException(404, "job not found")
    path = job.log_path(_safe_name(name))
    if not path.exists():
        raise HTTPException(404, "no log for this ligand")
    return PlainTextResponse(path.read_text())


@app.get("/api/jobs/{job_id}/ligands/{name}/poses/{pose_index}/interactions2d.svg")
def ligand_interactions_2d(job_id: str, name: str, pose_index: int) -> Response:
    """A flattened 2D interaction diagram (ProLIF detection + our own LigPlot-
    style render, see backend/interactions_prolif.py and backend/interaction_render.py)
    for one specific already-docked pose. pose_index is 1-based, matching the
    pose bar's labels and docking.Pose.rank."""
    job = MANAGER.get(job_id)
    if not job:
        raise HTTPException(404, "job not found")
    name = _safe_name(name)
    poses_path = job.poses_path(name)
    if not poses_path.exists():
        raise HTTPException(404, "no poses for this ligand")
    if pose_index < 1:
        raise HTTPException(400, "pose_index is 1-based")

    cache_path = job.dir / "interactions" / f"{name}_{pose_index}.svg"
    if cache_path.exists():
        return Response(cache_path.read_text(), media_type="image/svg+xml")

    try:
        mol, hits = iprolif.analyze_pose(job.receptor_path.read_text(), poses_path, pose_index)
        svg = interaction_render.render_svg(mol, hits)
    except ValueError as exc:
        raise HTTPException(404, str(exc)) from exc
    except Exception as exc:  # ProLIF failure - a broken diagram isn't a server outage
        svg = interaction_render.render_placeholder_svg(f"interaction detection failed: {exc}"[:200])

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(svg)
    return Response(svg, media_type="image/svg+xml")


def _results_csv(job) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(
        ["name", "smiles", "status", "best_affinity_kcal_per_mol", "n_poses",
         "all_scores", "all_strain_energies_kcal_per_mol", "all_rmsd_to_reference_angstrom",
         "constraints_satisfied", "error"]
    )
    for lig in job.ligands:
        writer.writerow(
            [
                lig.name,
                lig.smiles or "",
                lig.status,
                "" if lig.best_score is None else f"{lig.best_score:.2f}",
                lig.n_poses,
                ";".join(f"{s:.2f}" for s in lig.scores),
                ";".join("" if e is None else f"{e:.2f}" for e in lig.strain_energies),
                ";".join("" if r is None else f"{r:.2f}" for r in lig.rmsd_to_reference),
                "" if lig.constraint_satisfied is None else str(lig.constraint_satisfied),
                lig.error or "",
            ]
        )
    return buf.getvalue()


@app.get("/api/jobs/{job_id}/results.csv")
def results_csv(job_id: str) -> Response:
    job = MANAGER.get(job_id)
    if not job:
        raise HTTPException(404, "job not found")
    return PlainTextResponse(
        _results_csv(job),
        media_type="text/csv",
        headers={"Content-Disposition": f'attachment; filename="docking_{job_id}_results.csv"'},
    )


@app.get("/api/jobs/{job_id}/download")
def download_all(job_id: str) -> Response:
    """One ZIP: every ligand's top poses, the receptor used, and a summary CSV."""
    job = MANAGER.get(job_id)
    if not job:
        raise HTTPException(404, "job not found")
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for lig in job.ligands:
            path = job.poses_path(lig.name)
            if path.exists():
                zf.write(path, arcname=f"poses/{path.name}")
        if job.receptor_path.exists():
            zf.write(job.receptor_path, arcname="receptor.pdb")
        zf.writestr("results.csv", _results_csv(job))
        zf.writestr("job.json", JSONResponse(job.as_dict()).body.decode())
    buf.seek(0)
    return Response(
        buf.read(),
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="docking_{job_id}.zip"'},
    )


# Frontend last: it claims "/" and must not shadow the API routes above.
app.mount("/", StaticFiles(directory=str(config.FRONTEND_DIR), html=True), name="frontend")
