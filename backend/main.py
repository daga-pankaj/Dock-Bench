"""FastAPI app: the HTTP API behind the Mol* web page in frontend/.

All the work happens in backend.services (shared with the Streamlit app in
streamlit_app/); the routes here only translate HTTP to service calls and
ServiceError back to HTTP status codes.

Run with:  uvicorn backend.main:app --reload --port 8000
"""

from __future__ import annotations

import base64
import hmac

from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse, Response
from fastapi.staticfiles import StaticFiles

from . import config, services
from .schemas import FetchRequest, JobRequest, ParseRequest, QuickPrepRequest

config.ensure_dirs()

if config.AUTH_GENERATED:
    print(f"\n>>> auth enabled - user: {config.AUTH_USER}  password: {config.AUTH_PASS}\n")

app = FastAPI(title="Docking Bench", version="1.0")


@app.exception_handler(services.ServiceError)
async def service_error(_request: Request, exc: services.ServiceError) -> JSONResponse:
    return JSONResponse(status_code=exc.status, content={"detail": exc.detail})


@app.middleware("http")
async def no_cache_frontend(request, call_next):
    """Make the browser revalidate the frontend on every load. Without a
    Cache-Control header it heuristically caches app.js/style.css, so a
    rebuilt bundle can go unseen until a hard refresh."""
    response = await call_next(request)
    if not request.url.path.startswith("/api/"):
        response.headers.setdefault("Cache-Control", "no-cache")
    return response


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
# Status, receptors, uploads
# ---------------------------------------------------------------------------


@app.get("/api/status")
def status() -> dict:
    return services.status()


@app.post("/api/receptor/upload")
async def upload_receptor(file: UploadFile = File(...)) -> dict:
    return services.upload_receptor(file.filename, await file.read())


@app.post("/api/pharmacophore/upload")
async def upload_pharmacophore(file: UploadFile = File(...)) -> dict:
    """Upload a MOE .ph4 pharmacophore file for docking constraints."""
    return services.upload_pharmacophore(file.filename, await file.read())


@app.post("/api/reactions/upload")
async def upload_reaction(file: UploadFile = File(...)) -> dict:
    """Upload a custom .rxn reaction file for covalent docking."""
    return services.upload_reaction(file.filename, await file.read())


@app.post("/api/receptor/fetch")
def fetch_receptor(req: FetchRequest) -> dict:
    return services.fetch_receptor(req.pdb_id)


@app.post("/api/receptor/{receptor_id}/quickprep")
def quickprep_receptor(receptor_id: str, req: QuickPrepRequest) -> dict:
    return services.quickprep_receptor(receptor_id, req.site_ligand_key)


@app.get("/api/receptor/{receptor_id}/pdb")
def receptor_pdb(receptor_id: str) -> Response:
    return PlainTextResponse(services.receptor_pdb(receptor_id))


@app.get("/api/receptor/{receptor_id}/residues/{resname}")
def list_residues(receptor_id: str, resname: str) -> dict:
    """List all residues of given type (e.g., CYS, SER) in receptor."""
    return services.list_residues(receptor_id, resname)


@app.get("/api/pdb_components/{comp_id}/smiles")
def pdb_component_smiles(comp_id: str, receptor_id: str | None = None, het_key: str | None = None) -> dict:
    return services.pdb_component_smiles(comp_id, receptor_id, het_key)


# ---------------------------------------------------------------------------
# Ligands
# ---------------------------------------------------------------------------


@app.post("/api/ligands/parse_text")
def parse_ligand_text(req: ParseRequest) -> dict:
    return services.parse_ligand_text(req.text)


@app.post("/api/ligands/parse_file")
async def parse_ligand_file(file: UploadFile = File(...)) -> dict:
    return services.parse_ligand_file(file.filename, await file.read())


# ---------------------------------------------------------------------------
# Jobs
# ---------------------------------------------------------------------------


@app.post("/api/jobs")
def create_job(req: JobRequest) -> dict:
    return services.create_job(req)


@app.get("/api/jobs")
def list_jobs() -> dict:
    return {"jobs": services.list_jobs()}


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str) -> dict:
    return services.get_job(job_id).as_dict()


@app.post("/api/jobs/{job_id}/cancel")
def cancel_job(job_id: str) -> dict:
    services.cancel_job(job_id)
    return {"ok": True}


@app.get("/api/jobs/{job_id}/receptor.pdb")
def job_receptor(job_id: str) -> Response:
    return PlainTextResponse(services.job_receptor_pdb(job_id))


@app.get("/api/jobs/{job_id}/ligands/{name}/poses.sdf")
def ligand_poses(job_id: str, name: str) -> Response:
    path = services.ligand_poses_path(job_id, name)
    return FileResponse(path, media_type="chemical/x-mdl-sdfile", filename=path.name)


@app.get("/api/jobs/{job_id}/ligands/{name}/log")
def ligand_log(job_id: str, name: str) -> Response:
    return PlainTextResponse(services.ligand_log(job_id, name))


@app.get("/api/jobs/{job_id}/ligands/{name}/depiction.svg")
def ligand_depiction(job_id: str, name: str) -> Response:
    return Response(services.ligand_depiction_svg(job_id, name), media_type="image/svg+xml")


@app.get("/api/jobs/{job_id}/ligands/{name}/poses/{pose_index}/interactions2d.svg")
def ligand_interactions_2d(job_id: str, name: str, pose_index: int, style: str = "auto") -> Response:
    """The X-Diagram-Source response header says which renderer (moe/prolif)
    the diagram came from - see services.interactions_2d."""
    svg, source = services.interactions_2d(job_id, name, pose_index, style)
    return Response(svg, media_type="image/svg+xml", headers={"X-Diagram-Source": source})


@app.get("/api/jobs/{job_id}/results.csv")
def results_csv(job_id: str) -> Response:
    return PlainTextResponse(
        services.results_csv(job_id),
        media_type="text/csv",
        headers={"Content-Disposition": f'attachment; filename="docking_{job_id}_results.csv"'},
    )


@app.get("/api/jobs/{job_id}/download")
def download_all(job_id: str) -> Response:
    """One ZIP: every ligand's top poses, the receptor used, and a summary CSV."""
    return Response(
        services.results_zip(job_id),
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="docking_{job_id}.zip"'},
    )


# Frontend last: it claims "/" and must not shadow the API routes above.
app.mount("/", StaticFiles(directory=str(config.FRONTEND_DIR), html=True), name="frontend")
