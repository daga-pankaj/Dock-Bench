"""FastAPI app: receptor handling, ligand parsing, docking jobs, downloads.

Run with:  uvicorn backend.main:app --reload --port 8000
"""

from __future__ import annotations

import csv
import io
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
from .jobs import MANAGER, JobParams

config.ensure_dirs()

app = FastAPI(title="Docking Bench", version="1.0")

PDB_ID_RE = re.compile(r"^[0-9A-Za-z]{4}$")
MAX_BOX = 60.0


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
    keep_waters: bool = False
    exclude_het: list[str] = Field(default_factory=list)


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
    exe = config.find_smina()
    version = None
    if exe:
        try:
            version = docking.smina_version()
        except Exception as exc:
            version = f"error: {exc}"
    return {
        "smina": exe,
        "smina_version": version,
        "workers": config.MAX_WORKERS,
        "scoring_functions": list(docking.SCORING_FUNCTIONS),
        "default_num_modes": config.DEFAULT_NUM_MODES,
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
    if not config.find_smina():
        raise HTTPException(503, "smina not installed - run ./setup.sh")
    if not req.ligands:
        raise HTTPException(400, "no ligands submitted")
    if req.scoring not in docking.SCORING_FUNCTIONS:
        raise HTTPException(400, f"unknown scoring function {req.scoring}")
    for value in req.size:
        if not (4.0 <= value <= MAX_BOX):
            raise HTTPException(400, f"box edge {value} A outside 4-{MAX_BOX} A")

    receptor_text = _receptor_path(req.receptor_id).read_text()
    cleaned = prep.clean_receptor(receptor_text, keep_waters=req.keep_waters)
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


def _results_csv(job) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["name", "smiles", "status", "best_affinity_kcal_per_mol", "n_poses", "all_scores", "error"])
    for lig in job.ligands:
        writer.writerow(
            [
                lig.name,
                lig.smiles or "",
                lig.status,
                "" if lig.best_score is None else f"{lig.best_score:.2f}",
                lig.n_poses,
                ";".join(f"{s:.2f}" for s in lig.scores),
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
