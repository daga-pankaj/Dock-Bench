"""Job store: one docking job holds many ligands, each docked in its own thread.

Everything a job produces lands under runs/jobs/<job_id>/ so results survive a
server restart even though the live status lives in memory.
"""

from __future__ import annotations

import json
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field, fields
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem

from . import config, covalent_check, docking, prep, rmsd, strain
from . import pharmacophore as pharm


FINISHED = ("done", "failed", "cancelled")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _known_fields(cls, data: dict) -> dict:
    """The entries of `data` that are fields of dataclass `cls` - lets a
    results.json written by an older (fewer fields) or newer version load,
    missing fields falling back to their defaults."""
    names = {f.name for f in fields(cls)}
    return {k: v for k, v in data.items() if k in names}


@dataclass
class LigandRun:
    name: str
    smiles: str | None
    status: str = "queued"  # queued | running | done | failed | cancelled
    best_score: float | None = None
    scores: list[float] = field(default_factory=list)
    # Ligand-only MMFF strain per pose (see backend.strain), same order/index
    # as `scores`. None entries mean strain couldn't be scored for that pose
    # (e.g. an element MMFF has no parameters for) - never fails the job.
    strain_energies: list[float | None] = field(default_factory=list)
    # RMSD (Å) of each pose to the receptor's own bound-ligand reference (see
    # backend.rmsd), same order/index as scores. Empty unless the job was
    # submitted with a reference_ligand_key; None entries are poses obrms
    # couldn't match to the reference's heavy-atom graph at all (usually
    # because this ligand isn't the one that was bound).
    rmsd_to_reference: list[float | None] = field(default_factory=list)
    n_poses: int = 0
    error: str | None = None
    started: str | None = None
    finished: str | None = None
    # From the top-ranked pose only (parse_poses already sorted constraint-
    # satisfying poses first); None when the job had no constraints.
    constraint_satisfied: bool | None = None
    constraint_detail: list[dict] = field(default_factory=list)
    # MOE only, same order/index as scores: every score field MOE wrote per
    # pose (docking.MOE_SCORE_FIELDS), and for covalent jobs the ligand
    # environment each pose is bonded through ("*" = the residue atom).
    pose_scores: list[dict] = field(default_factory=list)
    covalent_attachments: list[str | None] = field(default_factory=list)
    # Non-fatal problems worth surfacing next to the result (e.g. the ligand
    # matching a covalent reaction zero or several times).
    warnings: list[str] = field(default_factory=list)

    @property
    def poses_file(self) -> str:
        return f"{self.name}_poses.sdf"


@dataclass
class JobParams:
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    exhaustiveness: int = 8
    num_modes: int = config.DEFAULT_NUM_MODES
    energy_range: float = 4.0
    seed: int = 42
    scoring: str = "vina"
    engine: str = "smina"
    use_gpu: bool = True
    # Extra `--minimize` pass (see docking.minimize) run on every pose after
    # its own docking search finishes - relaxes clashes a coarse
    # exhaustiveness/energy_range left behind, using the same scoring
    # function docking already used.
    minimize_after_docking: bool = False
    # Raw dicts (kind/point/radius/required/weight/label), not Constraint
    # objects, so this dataclass stays plain-JSON-serialisable via asdict().
    constraints: list[dict] = field(default_factory=list)
    # Raw HETATM lines for the receptor's own bound ligand, pulled from the
    # original (uncleaned) receptor - see backend.rmsd/prep.het_ligand_pdb_block.
    # None unless this job is a self-docking RMSD validation run.
    reference_ligand_pdb: str | None = None
    # MOE-specific parameters
    moe_placement: str = "Triangle Matcher"
    moe_refinement: str = "Rigid Receptor"
    ph4_id: str | None = None  # MOE pharmacophore file ID
    # Covalent docking (MOE-only)
    covalent_enabled: bool = False
    covalent_reaction: str | None = None  # e.g., "michael_acceptor_s" (built-in)
    covalent_residue: str | None = None   # e.g., "A:CYS:189"
    custom_rxn_id: str | None = None      # Custom .rxn file ID


class Job:
    def __init__(
        self,
        job_id: str,
        params: JobParams,
        ligands: list[LigandRun],
        receptor_name: str,
        created: str | None = None,
        status: str = "queued",
    ):
        self.id = job_id
        self.params = params
        self.ligands = ligands
        self.receptor_name = receptor_name
        self.created = created or _now()
        self.status = status  # queued | running | done | cancelled
        self.cancelled = threading.Event()
        self.lock = threading.Lock()
        self.dir = config.JOB_DIR / job_id

    # -- paths -------------------------------------------------------------
    @property
    def receptor_path(self) -> Path:
        return self.dir / "receptor.pdb"

    def ligand_input(self, name: str) -> Path:
        return self.dir / "ligands" / f"{name}.sdf"

    def poses_path(self, name: str) -> Path:
        return self.dir / "poses" / f"{name}_poses.sdf"

    def log_path(self, name: str) -> Path:
        return self.dir / "logs" / f"{name}.log"

    # -- serialisation -----------------------------------------------------
    @property
    def n_finished(self) -> int:
        return sum(1 for lig in self.ligands if lig.status in FINISHED)

    def summary(self) -> dict:
        """The short form GET /api/jobs lists."""
        return {
            "id": self.id,
            "status": self.status,
            "created": self.created,
            "receptor": self.receptor_name,
            "n_ligands": len(self.ligands),
            "n_finished": self.n_finished,
        }

    def as_dict(self) -> dict:
        params = asdict(self.params)
        params["center"] = list(self.params.center)
        params["size"] = list(self.params.size)
        return {
            "id": self.id,
            "status": self.status,
            "created": self.created,
            "receptor": self.receptor_name,
            "params": params,
            "n_ligands": len(self.ligands),
            "n_finished": self.n_finished,
            "ligands": [asdict(lig) for lig in self.ligands],
        }

    def save(self) -> None:
        self.dir.mkdir(parents=True, exist_ok=True)
        (self.dir / "results.json").write_text(json.dumps(self.as_dict(), indent=2))

    @classmethod
    def load(cls, job_dir: Path) -> "Job | None":
        """Rebuild a Job from the results.json a previous process wrote.

        Any ligand still "queued"/"running" belongs to a worker thread that no
        longer exists after a restart, so it's marked failed rather than left
        stuck forever.
        """
        results_path = job_dir / "results.json"
        if not results_path.exists():
            return None
        try:
            data = json.loads(results_path.read_text())
        except (json.JSONDecodeError, OSError):
            return None

        params = JobParams(**_known_fields(JobParams, data["params"]))
        params.center = tuple(params.center)
        params.size = tuple(params.size)
        ligands = []
        for lig in data["ligands"]:
            run = LigandRun(**_known_fields(LigandRun, lig))
            if run.status in ("queued", "running"):
                run.status, run.error = "failed", "server restarted before this ligand finished"
            ligands.append(run)
        job = cls(data["id"], params, ligands, data["receptor"], created=data["created"])
        job.status = "cancelled" if data["status"] == "cancelled" else "done"
        if params.engine == "moe":
            # MOE jobs saved before per-pose scores were recorded: recover
            # them from the poses SDF, which has always carried MOE's fields.
            for lig in job.ligands:
                sdf = job.poses_path(lig.name)
                if lig.status != "done" or lig.pose_scores or not sdf.exists():
                    continue
                try:
                    poses = docking._parse_moe_sdf(
                        sdf, [pharm.Constraint.from_dict(c) for c in params.constraints] or None
                    )
                except Exception:
                    continue
                lig.pose_scores = [p.moe_scores for p in poses]
                if params.covalent_enabled:
                    lig.covalent_attachments = [p.covalent_attachment for p in poses]
        job.save()
        return job


class JobManager:
    """Owns the thread pool and the in-memory job table."""

    def __init__(self) -> None:
        config.ensure_dirs()
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()
        self._pool = ThreadPoolExecutor(max_workers=config.MAX_WORKERS, thread_name_prefix="dock")
        # MOE ligands get their own pool, sized to the license budget: they
        # queue across all jobs (status "queued") instead of all starting
        # at once and failing for want of tokens, and they don't tie up the
        # smina/gnina workers while they wait.
        self._moe_pool = ThreadPoolExecutor(max_workers=config.MOE_MAX_CONCURRENT, thread_name_prefix="moe")
        self._load_from_disk()

    def _load_from_disk(self) -> None:
        """Repopulate the job table from runs/jobs/ so a restart doesn't lose history."""
        if not config.JOB_DIR.exists():
            return
        for job_dir in config.JOB_DIR.iterdir():
            if not job_dir.is_dir():
                continue
            job = Job.load(job_dir)
            if job:
                self._jobs[job.id] = job

    def get(self, job_id: str) -> Job | None:
        with self._lock:
            return self._jobs.get(job_id)

    def list_jobs(self) -> list[dict]:
        with self._lock:
            jobs = list(self._jobs.values())
        jobs.sort(key=lambda j: j.created, reverse=True)
        return [j.summary() for j in jobs]

    def cancel(self, job_id: str) -> bool:
        job = self.get(job_id)
        if not job or job.status in ("done", "cancelled"):
            return False
        job.cancelled.set()
        return True

    def submit(
        self,
        receptor_pdb: str,
        receptor_name: str,
        specs: list[prep.LigandSpec],
        params: JobParams,
    ) -> Job:
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
        job_id = f"{params.engine}_{timestamp}"
        ligands = [LigandRun(name=s.name, smiles=s.smiles) for s in specs]
        job = Job(job_id, params, ligands, receptor_name)

        for sub in ("ligands", "poses", "logs"):
            (job.dir / sub).mkdir(parents=True, exist_ok=True)
        job.receptor_path.write_text(receptor_pdb)
        job.save()

        with self._lock:
            self._jobs[job_id] = job

        # One pool task per ligand: a 500-ligand job streams through the pool
        # instead of blocking a worker for the whole batch.
        pool = self._moe_pool if params.engine == "moe" else self._pool
        for spec, run in zip(specs, ligands):
            pool.submit(self._run_one, job, spec, run)
        return job

    # -- worker ------------------------------------------------------------
    def _run_one(self, job: Job, spec: prep.LigandSpec, run: LigandRun) -> None:
        if job.cancelled.is_set():
            self._finish(job, run, status="cancelled", error="job cancelled")
            return

        with job.lock:
            if job.status == "queued":
                job.status = "running"
            run.status = "running"
            run.started = _now()

        try:
            molblock = prep.build_3d(spec, seed=job.params.seed)
            lig_path = job.ligand_input(run.name)
            lig_path.write_text(molblock + "$$$$\n")

            if job.params.engine == "moe" and job.params.covalent_enabled:
                warning = covalent_check.reactant_warning(
                    molblock,
                    docking.covalent_rxn_path(job.params.covalent_reaction, job.params.custom_rxn_id),
                )
                if warning:
                    with job.lock:
                        run.warnings.append(warning)

            # A job with fewer ligands than worker slots would otherwise leave
            # most cores idle - a 1-ligand job used to run smina pinned to a
            # single core even on a 9-worker box. Spread MAX_WORKERS cores
            # across however many ligands are actually running concurrently
            # (never more slots than ligands in the job), so a small job still
            # lets each smina process use several cores internally.
            n_concurrent = max(1, min(len(job.ligands), config.MAX_WORKERS))
            cpu_per_ligand = max(1, config.MAX_WORKERS // n_concurrent)

            constraint_list = [pharm.Constraint.from_dict(c) for c in job.params.constraints]
            constraints = constraint_list or None

            dock_params = self._dock_params(job, run, lig_path, constraints, cpu_per_ligand)
            poses = docking.dock(**dock_params)
            # MOE doesn't support separate minimize step via smina/gnina
            if job.params.minimize_after_docking and job.params.engine != "moe":
                poses_path = job.poses_path(run.name)
                # smina reads --ligand while it's writing --out - a distinct
                # temp path avoids reading and writing the same file at once,
                # then the temp file replaces the original once minimize
                # actually succeeds (a failed minimize leaves the un-
                # minimized poses in place rather than a half-written file).
                minimized_path = poses_path.with_name(poses_path.stem + ".minimized.sdf")
                poses = docking.minimize(
                    receptor=job.receptor_path,
                    poses_sdf=poses_path,
                    out_sdf=minimized_path,
                    center=job.params.center,
                    size=job.params.size,
                    log_path=job.log_path(run.name),
                    scoring=job.params.scoring,
                    engine=job.params.engine,
                    use_gpu=job.params.use_gpu,
                    cpu=cpu_per_ligand,
                    constraints=constraints,
                )
                minimized_path.replace(poses_path)
            # poses are already best-constrained-first when constraints were
            # given (see docking.parse_poses), so [0] is the one to report.
            scores = [p.affinity for p in poses if p.affinity is not None]
            run.scores = scores
            run.n_poses = len(poses)
            run.best_score = min(scores) if scores else None
            run.constraint_satisfied = poses[0].constraint_satisfied
            run.constraint_detail = poses[0].constraint_detail
            if job.params.engine == "moe":
                run.pose_scores = [p.moe_scores for p in poses]
                if job.params.covalent_enabled:
                    run.covalent_attachments = [p.covalent_attachment for p in poses]

            try:
                pose_mols = [
                    m for m in Chem.SDMolSupplier(str(job.poses_path(run.name)), removeHs=False, sanitize=False)
                    if m is not None
                ]
            except Exception:
                pose_mols = []

            try:  # strain is a nicety on top of the docking result, never worth failing the job over
                run.strain_energies = strain.strain_energies(pose_mols) if pose_mols else []
            except Exception:
                run.strain_energies = []

            try:  # RMSD to the receptor's own bound-ligand pose, see backend.rmsd
                if job.params.reference_ligand_pdb and poses:
                    run.rmsd_to_reference = rmsd.rmsds_to_reference(
                        job.params.reference_ligand_pdb, job.poses_path(run.name), len(poses)
                    )
            except Exception:
                run.rmsd_to_reference = []

            self._finish(job, run, status="done")
        except Exception as exc:  # one bad ligand must not sink the batch
            self._finish(job, run, status="failed", error=str(exc)[:400])

    @staticmethod
    def _dock_params(job: Job, run: LigandRun, lig_path: Path, constraints, cpu: int) -> dict:
        """docking.dock() keyword arguments - MOE takes different ones than smina/gnina."""
        p = job.params
        common = {
            "receptor": job.receptor_path,
            "ligand_sdf": lig_path,
            "out_sdf": job.poses_path(run.name),
            "center": p.center,
            "size": p.size,
            "log_path": job.log_path(run.name),
            "num_modes": p.num_modes,
            "scoring": p.scoring,
            "engine": p.engine,
            "constraints": constraints,
        }
        if p.engine == "moe":
            return {
                **common,
                "placement": p.moe_placement,
                "refinement": p.moe_refinement,
                "ph4_file": config.pharmacophore_path(p.ph4_id) if p.ph4_id else None,
                "covalent_enabled": p.covalent_enabled,
                "covalent_reaction": p.covalent_reaction,
                "covalent_residue": p.covalent_residue,
                "custom_rxn_id": p.custom_rxn_id,
            }
        return {
            **common,
            "exhaustiveness": p.exhaustiveness,
            "energy_range": p.energy_range,
            "seed": p.seed,
            "use_gpu": p.use_gpu,
            "cpu": cpu,
        }

    def _finish(self, job: Job, run: LigandRun, status: str, error: str | None = None) -> None:
        with job.lock:
            run.status = status
            run.error = error
            run.finished = _now()
            if all(lig.status in FINISHED for lig in job.ligands):
                job.status = "cancelled" if job.cancelled.is_set() else "done"
            job.save()


MANAGER = JobManager()
