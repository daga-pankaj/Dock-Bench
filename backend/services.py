"""Everything the app does, independent of how it's presented.

Both front ends call into this module: the HTTP API (backend.main, which
the original Mol* web page talks to) and the Streamlit app
(streamlit_app/app.py). Failures are raised as ServiceError carrying the
HTTP status the API should answer with; the Streamlit app shows the message.
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import urllib.error
import urllib.request
import uuid
import zipfile
from pathlib import Path

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

from rdkit import Chem

from . import config, covalent_check, docking, interaction_render, moe_dock, prep
from . import interactions_prolif as iprolif
from . import pharmacophore as pharm
from .jobs import MANAGER, Job, JobParams
from .schemas import JobRequest

PDB_ID_RE = re.compile(r"^[0-9A-Za-z]{4}$")
COMP_ID_RE = re.compile(r"^[0-9A-Za-z]{1,5}$")
MAX_BOX = 60.0

MAX_RECEPTOR_BYTES = 60_000_000
MAX_PH4_BYTES = 1_000_000
MAX_RXN_BYTES = 500_000
MAX_LIGAND_FILE_BYTES = 200_000_000

DIAGRAM_STYLES = ("auto", "moe", "prolif")


class ServiceError(Exception):
    """A request that can't be served: `status` is the HTTP code, `detail` the reason."""

    def __init__(self, status: int, detail: str):
        super().__init__(detail)
        self.status = status
        self.detail = detail


def _decode_upload(raw: bytes, limit: int) -> str:
    if len(raw) > limit:
        raise ServiceError(413, "file too large")
    return raw.decode("utf-8", errors="replace")


def _fetch_url(url: str, timeout: int) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read()


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------


def status() -> dict:
    """Which engines are installed (with versions) and the app's option lists."""
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
# Receptors
# ---------------------------------------------------------------------------


def receptor_path(receptor_id: str) -> Path:
    if not re.fullmatch(r"[0-9a-f]{8,32}", receptor_id):
        raise ServiceError(400, "bad receptor id")
    path = config.RECEPTOR_DIR / f"{receptor_id}.pdb"
    if not path.exists():
        raise ServiceError(404, "receptor not found; upload it again")
    return path


def prepared_from(receptor_id: str) -> str | None:
    """Id of the receptor a QuickPrep'd receptor was made from, else None."""
    meta = receptor_path(receptor_id).with_suffix(".json")  # validates the id
    return json.loads(meta.read_text())["prepared_from"] if meta.exists() else None


def source_text(receptor_id: str) -> str:
    """The as-uploaded PDB behind a receptor. For a QuickPrep'd one that is
    its source: the crystal ligand pose (RMSD reference) and the LINK
    records (covalent adducts) are only faithful there."""
    return receptor_path(prepared_from(receptor_id) or receptor_id).read_text()


def receptor_pdb(receptor_id: str) -> str:
    return receptor_path(receptor_id).read_text()


def store_receptor(name: str, text: str, prepared_from: str | None = None) -> dict:
    """Save a receptor PDB under a fresh id and describe it for the UI."""
    if "ATOM" not in text and "HETATM" not in text:
        raise ServiceError(400, "file contains no ATOM/HETATM records")
    receptor_id = uuid.uuid4().hex[:16]
    (config.RECEPTOR_DIR / f"{receptor_id}.pdb").write_text(text)
    if prepared_from:
        (config.RECEPTOR_DIR / f"{receptor_id}.json").write_text(
            json.dumps({"prepared_from": prepared_from}))
    hets = prep.list_het_ligands(text)
    n_protein = sum(1 for tag, _ in prep._atom_records(text) if tag == "ATOM")
    return {
        "id": receptor_id,
        "name": name,
        "n_protein_atoms": n_protein,
        "het_ligands": [h.as_dict() for h in hets],
        "chains": [c.as_dict() for c in prep.list_chains(text)],
        "prepared_from": prepared_from,
    }


def upload_receptor(filename: str | None, raw: bytes) -> dict:
    text = _decode_upload(raw, MAX_RECEPTOR_BYTES)
    if filename and filename.lower().endswith(".cif"):
        raise ServiceError(400, "mmCIF not supported; upload a .pdb file")
    return store_receptor(filename or "receptor.pdb", text)


def fetch_receptor(pdb_id: str) -> dict:
    """Download a PDB entry from RCSB and store it as a receptor."""
    pdb_id = pdb_id.strip()
    if not PDB_ID_RE.match(pdb_id):
        raise ServiceError(400, "PDB ID must be 4 characters")
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        text = _fetch_url(url, timeout=30).decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        raise ServiceError(404, f"RCSB returned {exc.code} for {pdb_id.upper()}")
    except Exception as exc:
        raise ServiceError(502, f"could not reach RCSB: {exc}")
    return store_receptor(f"{pdb_id.upper()}.pdb", text)


def quickprep_receptor(receptor_id: str, site_ligand_key: str | None) -> dict:
    """Run MOE QuickPrep on an uploaded receptor and store the result as a new
    receptor (hydrogens, protonation states, minimized site, ligands kept),
    so it can be inspected before docking and then docked into like any
    other. The raw upload goes in, header included: Structure Preparation
    uses SEQRES to build short missing loops, as in the GUI."""
    if not config.find_moebatch():
        raise ServiceError(503, "MOE not installed")
    if prepared_from(receptor_id):
        raise ServiceError(400, "receptor is already QuickPrep'd")
    src = receptor_path(receptor_id)
    text = src.read_text()

    site, site_coords = None, []
    if site_ligand_key:
        block = prep.het_ligand_pdb_block(text, site_ligand_key)
        if block is None:
            raise ServiceError(400, f"site ligand {site_ligand_key} not found in receptor")
        resname, _, resseq = site_ligand_key.split(":")
        site = (resname, int(re.sub(r"\D", "", resseq) or 0))
        site_coords = [prep._xyz(l) for l in block.splitlines() if l.startswith("HETATM")]

    out = config.RECEPTOR_DIR / f"{receptor_id}_quickprep.pdb"
    try:
        moe_dock.quickprep(src, out, site_ligand=site,
                           log_path=config.RECEPTOR_DIR / f"{receptor_id}_quickprep.log")
    except moe_dock.MoeDockingFailed as exc:
        raise ServiceError(500, str(exc)) from exc
    prepared = out.read_text()
    out.unlink()

    # The upload's file name isn't stored server-side; the UI relabels this.
    meta = store_receptor("QuickPrep", prepared, prepared_from=receptor_id)
    meta["quickprep_report"] = prep.quickprep_report(text, prepared, site_coords)
    return meta


def list_residues(receptor_id: str, resname: str) -> dict:
    """All residues of one type (e.g. CYS, SER) in a receptor."""
    if not re.fullmatch(r"[A-Z]{3}", resname.upper()):
        raise ServiceError(400, "resname must be 3-letter code")
    text = receptor_path(receptor_id).read_text()
    return {"resname": resname.upper(), "residues": prep.list_residues_by_name(text, resname.upper())}


def pdb_component_smiles(comp_id: str, receptor_id: str | None = None, het_key: str | None = None) -> dict:
    """A bound ligand's own PDB Chemical Component Dictionary SMILES, for the
    "self-dock" button - lets a HETATM group's residue code (e.g. STI) become
    a ready-to-dock ligand without needing its bond orders/aromaticity
    guessed back from bare coordinates (unreliable - see backend/rmsd.py's
    module docstring for how RMSD validation sidesteps that problem).

    Given the receptor and the bound group's key, a group LINKed to the
    protein as an acyl adduct comes back in its pre-reaction (carboxylic
    acid) form instead, since the CCD stores the cut adduct - see
    covalent_check.acid_form_from_link. `covalent_note` says when that
    happened."""
    comp_id = comp_id.strip().upper()
    if not COMP_ID_RE.match(comp_id):
        raise ServiceError(400, "bad chemical component id")
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
    try:
        data = json.loads(_fetch_url(url, timeout=15).decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise ServiceError(404, f"RCSB has no chemical component {comp_id} ({exc.code})")
    except Exception as exc:
        raise ServiceError(502, f"could not reach RCSB: {exc}")

    descriptors = data.get("pdbx_chem_comp_descriptor", [])
    smiles = next((d["descriptor"] for d in descriptors if d.get("type") == "SMILES_CANONICAL"), None)
    if smiles is None:
        smiles = next((d["descriptor"] for d in descriptors if d.get("type") == "SMILES"), None)
    if smiles is None:
        raise ServiceError(404, f"RCSB has no SMILES on file for {comp_id}")
    out = {"comp_id": comp_id, "smiles": smiles, "name": data.get("chem_comp", {}).get("name", comp_id)}
    if receptor_id and het_key:
        acid = covalent_check.acid_form_from_link(source_text(receptor_id), het_key, smiles)
        if acid:
            out.update(smiles=acid, ccd_smiles=smiles,
                       covalent_note=f"{comp_id} is LINKed to the protein as an acyl adduct - "
                                     "rebuilt as the carboxylic acid it reacts from")
    return out


# ---------------------------------------------------------------------------
# Uploads used by MOE jobs
# ---------------------------------------------------------------------------


def upload_pharmacophore(filename: str | None, raw: bytes) -> dict:
    """Store a MOE .ph4 pharmacophore file for docking constraints."""
    text = _decode_upload(raw, MAX_PH4_BYTES)
    if "#moe:ph4que" not in text and "#pharmacophore" not in text:
        raise ServiceError(400, "not a valid MOE .ph4 file")
    ph4_id = uuid.uuid4().hex[:12]
    path = config.pharmacophore_path(ph4_id)
    path.parent.mkdir(exist_ok=True)
    path.write_text(text)
    return {"ph4_id": ph4_id, "filename": filename or "pharmacophore.ph4"}


def upload_reaction(filename: str | None, raw: bytes) -> dict:
    """Store a custom MDL .rxn reaction file for covalent docking."""
    text = _decode_upload(raw, MAX_RXN_BYTES)
    if not text.strip().startswith("$RXN"):
        raise ServiceError(400, "not a valid RXN file (missing $RXN header)")
    rxn_id = uuid.uuid4().hex[:12]
    path = config.custom_reaction_path(rxn_id)
    path.parent.mkdir(exist_ok=True)
    path.write_text(text)
    return {"rxn_id": rxn_id, "filename": filename or "custom.rxn"}


# ---------------------------------------------------------------------------
# Ligands
# ---------------------------------------------------------------------------


def preview_ligands(specs: list[prep.LigandSpec]) -> dict:
    """Validity and a one-line note per parsed ligand."""
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
        items.append({"name": spec.name, "smiles": spec.smiles, "molblock": spec.molblock,
                      "valid": valid, "note": note})
    return {"ligands": items, "n_valid": sum(1 for i in items if i["valid"])}


def parse_ligand_text(text: str) -> dict:
    # An SDF pasted into the box still starts with a mol block, not a SMILES.
    specs = prep.parse_sdf_block(text) if "$$$$" in text else prep.parse_smiles_block(text)
    return preview_ligands(specs)


def parse_ligand_file(filename: str | None, raw: bytes) -> dict:
    text = _decode_upload(raw, MAX_LIGAND_FILE_BYTES)
    return preview_ligands(prep.parse_ligand_file(filename or "ligands.smi", text))


# ---------------------------------------------------------------------------
# Jobs
# ---------------------------------------------------------------------------


def _require_engine(engine: str) -> None:
    if engine == "moe":
        if not config.find_moebatch():
            raise ServiceError(503, "MOE not installed")
    elif engine == "gnina":
        if not config.find_gnina():
            raise ServiceError(503, "gnina not installed - run ./setup.sh")
    elif not config.find_smina():
        raise ServiceError(503, "smina not installed - run ./setup.sh")


def _validate_covalent(req: JobRequest) -> None:
    if req.engine != "moe":
        raise ServiceError(400, "covalent docking requires MOE engine")
    if not req.covalent_residue:
        raise ServiceError(400, "covalent residue required (e.g., 'A:CYS:189')")
    # Either built-in reaction or custom .rxn required
    if not req.covalent_reaction and not req.custom_rxn_id:
        raise ServiceError(400, "covalent reaction template or custom .rxn file required")
    if req.covalent_reaction and req.covalent_reaction not in config.list_covalent_reactions():
        raise ServiceError(400, f"unknown reaction template {req.covalent_reaction}")
    if req.custom_rxn_id and not config.custom_reaction_path(req.custom_rxn_id).exists():
        raise ServiceError(400, "custom reaction file not found - re-upload")


def _validate_covalent_residue(receptor_text: str, residue: str) -> None:
    if prep.validate_residue_exists(receptor_text, residue):
        return
    suggestions = ""
    parsed = prep.parse_residue_spec(residue)
    if parsed:
        _, resname, _ = parsed
        available = prep.list_residues_by_name(receptor_text, resname)
        if available:
            suggestions = f" - available {resname} residues: {', '.join(available[:10])}"
    raise ServiceError(400, f"residue {residue} not found in receptor{suggestions}")


def create_job(req: JobRequest) -> dict:
    """Validate a docking request, clean its receptor and queue every ligand."""
    if req.engine not in docking.ENGINES:
        raise ServiceError(400, f"unknown docking engine {req.engine}")
    _require_engine(req.engine)

    if not req.ligands:
        raise ServiceError(400, "no ligands submitted")
    if len(req.ligands) > config.MAX_LIGANDS_PER_JOB:
        raise ServiceError(400, f"too many ligands in one job (max {config.MAX_LIGANDS_PER_JOB})")

    if req.engine == "moe":
        if req.scoring not in docking.MOE_SCORING_FUNCTIONS:
            raise ServiceError(400, f"unknown MOE scoring function {req.scoring}")
    elif req.scoring not in docking.SCORING_FUNCTIONS:
        raise ServiceError(400, f"unknown scoring function {req.scoring}")
    for value in req.size:
        if not (4.0 <= value <= MAX_BOX):
            raise ServiceError(400, f"box edge {value} A outside 4-{MAX_BOX} A")
    for c in req.constraints:
        if c.kind not in pharm.KINDS:
            raise ServiceError(400, f"unknown constraint kind {c.kind}")
        try:
            pharm.Constraint.from_dict(c.model_dump())
        except ValueError as exc:
            raise ServiceError(400, str(exc)) from exc
    if req.covalent_enabled:
        _validate_covalent(req)

    receptor_text = receptor_path(req.receptor_id).read_text()
    prepared = prepared_from(req.receptor_id) is not None

    if req.covalent_enabled and req.covalent_residue:
        _validate_covalent_residue(receptor_text, req.covalent_residue)

    reference_ligand_pdb = None
    if req.reference_ligand_key:
        reference_ligand_pdb = prep.het_ligand_pdb_block(source_text(req.receptor_id), req.reference_ligand_key)
        if reference_ligand_pdb is None:
            raise ServiceError(400, f"bound ligand {req.reference_ligand_key} not found in receptor")

    cleaned = prep.clean_receptor(
        receptor_text,
        keep_waters=req.keep_waters,
        strip_ligand=req.remove_ligand,
        exclude_chains=req.exclude_chains,
        site_center=tuple(req.center) if req.trim_distant_chains else None,
        chain_cutoff=req.chain_cutoff,
        keep_hydrogens=prepared,
    )
    if not cleaned.strip().splitlines()[:-1]:
        raise ServiceError(400, "receptor is empty after cleanup")

    specs = [prep.LigandSpec(name=lig.name, smiles=lig.smiles, molblock=lig.molblock) for lig in req.ligands]
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
    return MANAGER.submit(cleaned, f"receptor_{req.receptor_id}", specs, params).as_dict()


def list_jobs() -> list[dict]:
    return MANAGER.list_jobs()


def get_job(job_id: str) -> Job:
    job = MANAGER.get(job_id)
    if not job:
        raise ServiceError(404, "job not found")
    return job


def cancel_job(job_id: str) -> None:
    if not MANAGER.cancel(job_id):
        raise ServiceError(400, "job cannot be cancelled")


def safe_ligand_name(name: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9._+-]{1,120}", name):
        raise ServiceError(400, "bad ligand name")
    return name


def job_receptor_pdb(job_id: str) -> str:
    job = MANAGER.get(job_id)
    if not job or not job.receptor_path.exists():
        raise ServiceError(404, "job receptor not found")
    return job.receptor_path.read_text()


def ligand_poses_path(job_id: str, name: str) -> Path:
    path = get_job(job_id).poses_path(safe_ligand_name(name))
    if not path.exists():
        raise ServiceError(404, "no poses for this ligand")
    return path


def ligand_log(job_id: str, name: str) -> str:
    path = get_job(job_id).log_path(safe_ligand_name(name))
    if not path.exists():
        raise ServiceError(404, "no log for this ligand")
    return path.read_text()


def ligand_depiction_svg(job_id: str, name: str) -> str:
    """Plain 2D structure of a docked ligand, for the viewer's corner overlay.
    Drawn from the submitted SMILES when there is one (cleaner bond orders
    than a docking engine's SDF), else from the top pose."""
    job = get_job(job_id)
    name = safe_ligand_name(name)
    lig = next((l for l in job.ligands if l.name == name), None)
    if lig is None:
        raise ServiceError(404, "ligand not found")

    mol = Chem.MolFromSmiles(lig.smiles) if lig.smiles else None
    if mol is None:
        poses_path = job.poses_path(name)
        if poses_path.exists():
            mol = next(iter(Chem.SDMolSupplier(str(poses_path), removeHs=False)), None)
    if mol is None:
        raise ServiceError(404, "no structure to depict")
    try:
        return interaction_render.render_depiction_svg(mol)
    except Exception as exc:
        raise ServiceError(500, f"depiction failed: {exc}") from exc


def _moe_interactions_svg(job: Job, poses_path: Path, pose_index: int) -> str | None:
    """MOE's own Ligand Interactions diagram for a pose. A MOE docking run
    draws these itself (moe_dock.LID_RENDER); for smina/gnina poses, or MOE
    jobs docked before that, it is drawn now (a short moebatch run in the
    license queue) and kept next to the poses. None if MOE can't draw it, so
    the caller falls back to the ProLIF diagram."""
    svg_path = moe_dock.interactions_svg_path(poses_path, pose_index)
    if not svg_path.exists():
        blocks = [b for b in poses_path.read_text().split("$$$$\n") if b.strip()]
        if pose_index > len(blocks) or not config.find_moebatch():
            return None
        pose_sdf = svg_path.with_suffix(".sdf")
        pose_sdf.write_text(blocks[pose_index - 1] + "$$$$\n")
        try:
            moe_dock.render_interactions(job.receptor_path, pose_sdf, svg_path)
        except Exception:
            return None
        finally:
            pose_sdf.unlink(missing_ok=True)
    return moe_dock.browser_svg(svg_path.read_text())


def interactions_2d(job_id: str, name: str, pose_index: int, style: str = "auto") -> tuple[str, str]:
    """(svg, source) - a flattened 2D interaction diagram for one already-docked
    pose, any engine. style "moe" (or "auto", when MOE is installed) is MOE's
    own Ligand Interactions diagram; "prolif" (or MOE failing) is ProLIF
    detection + our own LigPlot-style render, see backend/interactions_prolif.py
    and backend/interaction_render.py. `source` says which one came back.
    pose_index is 1-based, matching the pose bar's labels and docking.Pose.rank."""
    job = get_job(job_id)
    name = safe_ligand_name(name)
    poses_path = job.poses_path(name)
    if not poses_path.exists():
        raise ServiceError(404, "no poses for this ligand")
    if pose_index < 1:
        raise ServiceError(400, "pose_index is 1-based")
    if style not in DIAGRAM_STYLES:
        raise ServiceError(400, "style must be auto, moe or prolif")

    if style != "prolif":
        svg = _moe_interactions_svg(job, poses_path, pose_index)
        if svg is not None:
            return svg, "moe"

    cache_path = job.dir / "interactions" / f"{name}_{pose_index}.svg"
    if cache_path.exists():
        return cache_path.read_text(), "prolif"

    try:
        mol, hits = iprolif.analyze_pose(job.receptor_path.read_text(), poses_path, pose_index)
        svg = interaction_render.render_svg(mol, hits)
    except ValueError as exc:
        raise ServiceError(404, str(exc)) from exc
    except Exception as exc:  # ProLIF failure - a broken diagram isn't a server outage
        svg = interaction_render.render_placeholder_svg(f"interaction detection failed: {exc}"[:200])

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(svg)
    return svg, "prolif"


RESULTS_CSV_HEADER = [
    "name", "smiles", "status", "best_affinity_kcal_per_mol", "n_poses",
    "all_scores", "all_strain_energies_kcal_per_mol", "all_rmsd_to_reference_angstrom",
    "constraints_satisfied", "moe_pose_scores", "covalent_attachments", "warnings", "error",
]


def _fmt_list(values, sep: str = ";") -> str:
    return sep.join("" if v is None else f"{v:.2f}" for v in values)


def results_csv(job_id: str) -> str:
    """One row per ligand: scores, strain, RMSD, constraints, MOE fields."""
    job = get_job(job_id)
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(RESULTS_CSV_HEADER)
    for lig in job.ligands:
        writer.writerow([
            lig.name,
            lig.smiles or "",
            lig.status,
            "" if lig.best_score is None else f"{lig.best_score:.2f}",
            lig.n_poses,
            _fmt_list(lig.scores),
            _fmt_list(lig.strain_energies),
            _fmt_list(lig.rmsd_to_reference),
            "" if lig.constraint_satisfied is None else str(lig.constraint_satisfied),
            # one "field=value,..." group per pose, poses separated by ";"
            ";".join(",".join(f"{k}={v:.3f}" for k, v in ps.items()) for ps in lig.pose_scores),
            ";".join(a or "" for a in lig.covalent_attachments),
            " | ".join(lig.warnings),
            lig.error or "",
        ])
    return buf.getvalue()


def results_zip(job_id: str) -> bytes:
    """One ZIP: every ligand's top poses, the receptor used, and a summary CSV."""
    job = get_job(job_id)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for lig in job.ligands:
            path = job.poses_path(lig.name)
            if path.exists():
                zf.write(path, arcname=f"poses/{path.name}")
        if job.receptor_path.exists():
            zf.write(job.receptor_path, arcname="receptor.pdb")
        zf.writestr("results.csv", results_csv(job_id))
        zf.writestr("job.json", json.dumps(job.as_dict(), separators=(",", ":"), ensure_ascii=False))
    return buf.getvalue()
