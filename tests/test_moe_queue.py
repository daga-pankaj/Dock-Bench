"""Tests for backend.moe_dock._run_moebatch: the server-wide moebatch cap and
the retry while the MOE license server is out of tokens."""

from __future__ import annotations

import subprocess
import threading
import time

from backend import config, moe_dock

LICENSE_FAIL = subprocess.CompletedProcess(
    [], 1, "", "moe: Can't check out license 'moe 2024.0604':\nNot enough tokens available (-4)")
OK = subprocess.CompletedProcess([], 0, "done", "")


def test_retries_until_license_is_free(monkeypatch):
    results = [LICENSE_FAIL, LICENSE_FAIL, OK]
    monkeypatch.setattr(moe_dock.subprocess, "run", lambda *a, **k: results.pop(0))
    monkeypatch.setattr(config, "MOE_LICENSE_RETRY", 0)
    assert moe_dock._run_moebatch(["moebatch"], timeout=10, cwd=".") is OK
    assert results == []


def test_gives_up_after_license_wait(monkeypatch):
    calls = []
    monkeypatch.setattr(moe_dock.subprocess, "run", lambda *a, **k: calls.append(1) or LICENSE_FAIL)
    monkeypatch.setattr(config, "MOE_LICENSE_WAIT", 0)
    assert moe_dock._run_moebatch(["moebatch"], timeout=10, cwd=".") is LICENSE_FAIL
    assert len(calls) == 1


def test_other_failures_are_not_retried(monkeypatch):
    fail = subprocess.CompletedProcess([], 1, "", "Error: atom SG not found")
    calls = []
    monkeypatch.setattr(moe_dock.subprocess, "run", lambda *a, **k: calls.append(1) or fail)
    assert moe_dock._run_moebatch(["moebatch"], timeout=10, cwd=".") is fail
    assert len(calls) == 1


def test_never_more_moebatch_runs_than_slots(monkeypatch):
    monkeypatch.setattr(moe_dock, "_LICENSE_SLOTS", threading.BoundedSemaphore(2))
    lock, running, peak = threading.Lock(), [0], [0]

    def fake_run(*a, **k):
        with lock:
            running[0] += 1
            peak[0] = max(peak[0], running[0])
        time.sleep(0.05)
        with lock:
            running[0] -= 1
        return OK

    monkeypatch.setattr(moe_dock.subprocess, "run", fake_run)
    threads = [threading.Thread(target=moe_dock._run_moebatch, args=(["moebatch"], 10, ".")) for _ in range(6)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert peak[0] == 2
