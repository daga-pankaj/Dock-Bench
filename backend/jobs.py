"""Job store: one docking job holds many ligands, each docked in its own thread.

Everything a job produces lands under runs/jobs/<job_id>/ so results survive a
server restart even though the live status lives in memory.
"""

from __future__ import annotations

import json
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from . import config, docking, prep


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


@dataclass
class LigandRun:
    name: str
    smiles: str | None
    status: str = "queued"  # queued | running | done | failed | cancelled
    best_score: float | None = None
    scores: list[float] = field(default_factory=list)
    n_poses: int = 0
    error: str | None = None
    started: str | None = None
    finished: str | None = None

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
    def as_dict(self) -> dict:
        done = sum(1 for lig in self.ligands if lig.status in ("done", "failed", "cancelled"))
        return {
            "id": self.id,
            "status": self.status,
            "created": self.created,
            "receptor": self.receptor_name,
            "params": {
                "center": list(self.params.center),
                "size": list(self.params.size),
                "exhaustiveness": self.params.exhaustiveness,
                "num_modes": self.params.num_modes,
                "energy_range": self.params.energy_range,
                "seed": self.params.seed,
                "scoring": self.params.scoring,
                "engine": self.params.engine,
                "use_gpu": self.params.use_gpu,
            },
            "n_ligands": len(self.ligands),
            "n_finished": done,
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

        params = JobParams(
            center=tuple(data["params"]["center"]),
            size=tuple(data["params"]["size"]),
            exhaustiveness=data["params"]["exhaustiveness"],
            num_modes=data["params"]["num_modes"],
            energy_range=data["params"]["energy_range"],
            seed=data["params"]["seed"],
            scoring=data["params"]["scoring"],
            engine=data["params"].get("engine", "smina"),
            use_gpu=data["params"].get("use_gpu", True),
        )
        ligands = []
        for lig in data["ligands"]:
            status = lig["status"]
            error = lig["error"]
            if status in ("queued", "running"):
                status, error = "failed", "server restarted before this ligand finished"
            ligands.append(
                LigandRun(
                    name=lig["name"],
                    smiles=lig["smiles"],
                    status=status,
                    best_score=lig["best_score"],
                    scores=lig["scores"],
                    n_poses=lig["n_poses"],
                    error=error,
                    started=lig["started"],
                    finished=lig["finished"],
                )
            )
        job = cls(data["id"], params, ligands, data["receptor"], created=data["created"])
        job.status = "cancelled" if data["status"] == "cancelled" else "done"
        job.save()
        return job


class JobManager:
    """Owns the thread pool and the in-memory job table."""

    def __init__(self) -> None:
        config.ensure_dirs()
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()
        self._pool = ThreadPoolExecutor(max_workers=config.MAX_WORKERS, thread_name_prefix="dock")
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
        return [
            {
                "id": j.id,
                "status": j.status,
                "created": j.created,
                "receptor": j.receptor_name,
                "n_ligands": len(j.ligands),
                "n_finished": sum(
                    1 for lig in j.ligands if lig.status in ("done", "failed", "cancelled")
                ),
            }
            for j in jobs
        ]

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
        job_id = uuid.uuid4().hex[:12]
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
        for spec, run in zip(specs, ligands):
            self._pool.submit(self._run_one, job, spec, run)
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

            # A job with fewer ligands than worker slots would otherwise leave
            # most cores idle - a 1-ligand job used to run smina pinned to a
            # single core even on a 9-worker box. Spread MAX_WORKERS cores
            # across however many ligands are actually running concurrently
            # (never more slots than ligands in the job), so a small job still
            # lets each smina process use several cores internally.
            n_concurrent = max(1, min(len(job.ligands), config.MAX_WORKERS))
            cpu_per_ligand = max(1, config.MAX_WORKERS // n_concurrent)

            poses = docking.dock(
                receptor=job.receptor_path,
                ligand_sdf=lig_path,
                out_sdf=job.poses_path(run.name),
                center=job.params.center,
                size=job.params.size,
                log_path=job.log_path(run.name),
                exhaustiveness=job.params.exhaustiveness,
                num_modes=job.params.num_modes,
                energy_range=job.params.energy_range,
                seed=job.params.seed,
                scoring=job.params.scoring,
                engine=job.params.engine,
                use_gpu=job.params.use_gpu,
                cpu=cpu_per_ligand,
            )
            scores = [p.affinity for p in poses if p.affinity is not None]
            run.scores = scores
            run.n_poses = len(poses)
            run.best_score = min(scores) if scores else None
            self._finish(job, run, status="done")
        except Exception as exc:  # one bad ligand must not sink the batch
            self._finish(job, run, status="failed", error=str(exc)[:400])

    def _finish(self, job: Job, run: LigandRun, status: str, error: str | None = None) -> None:
        with job.lock:
            run.status = status
            run.error = error
            run.finished = _now()
            if all(lig.status in ("done", "failed", "cancelled") for lig in job.ligands):
                job.status = "cancelled" if job.cancelled.is_set() else "done"
            job.save()


MANAGER = JobManager()
