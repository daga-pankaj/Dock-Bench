"""Tests for backend.jobs: job persistence across a simulated restart."""

from __future__ import annotations

import importlib
import json

import pytest

from backend import config, jobs


@pytest.fixture
def isolated_runs_dir(tmp_path, monkeypatch):
    """Point config/jobs at a throwaway runs/ dir so tests never touch the real one."""
    monkeypatch.setattr(config, "RUNS_DIR", tmp_path)
    monkeypatch.setattr(config, "RECEPTOR_DIR", tmp_path / "receptors")
    monkeypatch.setattr(config, "JOB_DIR", tmp_path / "jobs")
    config.ensure_dirs()
    return tmp_path


def _write_job_dir(job_dir, *, job_status, ligand_status, ligand_error=None):
    job_dir.mkdir(parents=True)
    (job_dir / "results.json").write_text(
        json.dumps(
            {
                "id": job_dir.name,
                "status": job_status,
                "created": "2026-01-01T00:00:00+00:00",
                "receptor": "receptor_abc",
                "params": {
                    "center": [1.0, 2.0, 3.0],
                    "size": [20.0, 20.0, 20.0],
                    "exhaustiveness": 8,
                    "num_modes": 5,
                    "energy_range": 4.0,
                    "seed": 42,
                    "scoring": "vina",
                },
                "n_ligands": 1,
                "n_finished": 0 if ligand_status in ("queued", "running") else 1,
                "ligands": [
                    {
                        "name": "lig1",
                        "smiles": "CCO",
                        "status": ligand_status,
                        "best_score": None if ligand_status != "done" else -6.5,
                        "scores": [] if ligand_status != "done" else [-6.5],
                        "n_poses": 0 if ligand_status != "done" else 1,
                        "error": ligand_error,
                        "started": "2026-01-01T00:00:01+00:00",
                        "finished": None,
                    }
                ],
            }
        )
    )


def test_manager_reloads_finished_job_from_disk(isolated_runs_dir):
    _write_job_dir(config.JOB_DIR / "job1", job_status="done", ligand_status="done")
    manager = jobs.JobManager()
    job = manager.get("job1")
    assert job is not None
    assert job.status == "done"
    assert job.ligands[0].status == "done"
    assert job.ligands[0].best_score == -6.5


def test_manager_marks_interrupted_ligand_as_failed(isolated_runs_dir):
    _write_job_dir(config.JOB_DIR / "job2", job_status="running", ligand_status="running")
    manager = jobs.JobManager()
    job = manager.get("job2")
    assert job is not None
    # A ligand that was mid-flight when the process died can't be resumed.
    assert job.ligands[0].status == "failed"
    assert "restarted" in job.ligands[0].error
    assert job.status == "done"


def test_manager_preserves_cancelled_status(isolated_runs_dir):
    _write_job_dir(config.JOB_DIR / "job3", job_status="cancelled", ligand_status="cancelled")
    manager = jobs.JobManager()
    job = manager.get("job3")
    assert job.status == "cancelled"


def test_manager_ignores_dirs_without_results_json(isolated_runs_dir):
    (config.JOB_DIR / "empty_dir").mkdir(parents=True)
    manager = jobs.JobManager()
    assert manager.list_jobs() == []


def test_reload_round_trips_a_real_submitted_job(isolated_runs_dir, monkeypatch):
    """End-to-end: submit a job, let one ligand finish, restart the manager, check it survives."""
    from backend import prep

    manager = jobs.JobManager()
    captured = {}

    def fake_run_one(job, spec, run):
        run.status = "running"
        run.started = "2026-01-01T00:00:01+00:00"
        run.status, run.best_score, run.scores, run.n_poses = "done", -5.1, [-5.1], 1
        run.finished = "2026-01-01T00:00:02+00:00"
        job.status = "done"
        job.save()

    monkeypatch.setattr(manager, "_run_one", fake_run_one)
    specs = [prep.LigandSpec(name="lig1", smiles="CCO")]
    params = jobs.JobParams(center=(1, 2, 3), size=(20, 20, 20))
    job = manager.submit("ATOM\nEND\n", "receptor_x", specs, params)
    captured["id"] = job.id
    manager._pool.shutdown(wait=True)  # block until the (faked) docking task above has run

    # Simulate a process restart: fresh manager, same runs/ dir on disk.
    restarted = jobs.JobManager()
    reloaded = restarted.get(captured["id"])
    assert reloaded is not None
    assert reloaded.status == "done"
    assert reloaded.ligands[0].best_score == -5.1


def test_manager_submit_creates_job_and_lists_it(isolated_runs_dir, monkeypatch):
    manager = jobs.JobManager()
    from backend import prep

    monkeypatch.setattr(manager._pool, "submit", lambda fn, *a, **k: None)  # don't actually dock
    specs = [prep.LigandSpec(name="lig1", smiles="CCO")]
    params = jobs.JobParams(center=(0, 0, 0), size=(20, 20, 20))
    job = manager.submit("ATOM\nEND\n", "receptor_x", specs, params)
    assert manager.get(job.id) is job
    listed = manager.list_jobs()
    assert any(j["id"] == job.id for j in listed)


# -- covalent docking ---------------------------------------------------------


def _fake_dock(captured):
    def dock(**kwargs):
        captured.update(kwargs)
        from backend import docking
        return [docking.Pose(rank=1, affinity=-5.0, rmsd_lb=0.0, rmsd_ub=0.0)]
    return dock




