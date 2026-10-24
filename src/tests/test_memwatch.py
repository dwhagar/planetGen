# src/tests/test_memwatch.py

"""
The memory peaks of a job's steps and the nightly smoke run's budget (PERF.183): `util/memwatch.py` keeps the highest
figure of each group while a step runs, the job runner prints it after the step and keeps it in `memory.json`, and
`planetgen.cli.smoke` fails the run when a step peaked above its budget.
"""

import json
import os
import sys
import threading

import pytest

from planetgen import tuning
from planetgen.cli import smoke
from planetgen.util import memwatch
from planetgen.web import job_runner

HOLD_MB = 300
"""int: What the test step below holds, in megabytes."""

HOLD_STEP = (f"import time; held = b'x' * ({HOLD_MB} * 1024 * 1024); print('holding', flush=True); time.sleep(3)")


def _fake(figures, seen):
    """A sampler that gives the figures one after the other, then the last one again, and sets `seen` once all were given."""
    queue = list(figures)

    def sampler(_group):
        figure = queue.pop(0) if len(queue) > 1 else queue[0]
        if len(queue) == 1:
            seen.set()
        return figure
    return sampler


def _figures(mysql=0.0, workers=0.0, workers_n=0, step=0.0, host=0.0):
    return {"mysql": mysql, "workers": workers, "workers_n": workers_n, "step": step, "host": host}


def test_the_tracker_keeps_the_highest_figure_of_each_group():
    seen = threading.Event()
    sampler = _fake([_figures(mysql=100, step=50, host=900), _figures(mysql=400, workers=80, workers_n=2, step=20, host=700),
                     _figures(mysql=250, step=60, host=950)], seen)
    with memwatch.PeakTracker(interval=0.01, sampler=sampler) as watch:
        watch.set_group(1)
        assert seen.wait(10)
    peaks = watch.peaks
    assert (peaks["mysql"], peaks["workers"], peaks["workers_n"], peaks["step"], peaks["host"]) == (400, 80, 2, 60, 950)


def test_a_sampler_that_fails_never_stops_the_step():
    def broken(_group):
        raise RuntimeError("ps is gone")
    with memwatch.PeakTracker(interval=0.01, sampler=broken) as watch:
        pass
    assert watch.peaks["step"] == 0.0


def test_the_line_names_each_group_in_whole_megabytes():
    peaks = {"mysql": 884.4, "workers": 2310.2, "workers_n": 6, "step": 255.9, "host": 3100.0, "host_total": 15990.0}
    assert memwatch.format_line(peaks) == ("Memory peak: MySQL 884 MB, 6 workers 2,310 MB, the step itself 256 MB, "
                                           "the whole host 3,100 of 15,990 MB.")
    peaks.update(workers_n=1, host_total=None)
    assert "1 worker 2,310 MB" in memwatch.format_line(peaks)
    assert "the whole host 3,100 MB." in memwatch.format_line(peaks)


def test_the_sampler_sees_a_process_group_and_the_host():
    figures = memwatch.sample(os.getpgid(0))
    assert figures["step"] > 0
    assert figures["host"] > 0


@pytest.fixture
def held_job(tmp_path, monkeypatch):
    """A job whose one step holds HOLD_MB for three seconds, run by the job runner."""
    monkeypatch.setattr(memwatch, "INTERVAL_SECONDS", 0.2)
    job_dir = tmp_path / "jobs" / "job-1"
    job_dir.mkdir(parents=True)
    job = {"id": "job-1", "kind": "galaxy", "title": "Holds memory", "cwd": str(tmp_path), "env": {},
           "steps": [{"label": "Holds", "argv": [sys.executable, "-c", HOLD_STEP]}]}
    (job_dir / "job.json").write_text(json.dumps(job), encoding="utf-8")
    assert job_runner.run(str(job_dir)) == 0
    return job_dir


def test_the_job_runner_notes_each_steps_peak_in_the_log_and_in_memory_json(held_job):
    log = (held_job / "output.log").read_text(encoding="utf-8")
    assert "Memory peak: MySQL" in log
    entries = json.loads((held_job / job_runner.MEMORY_NAME).read_text(encoding="utf-8"))
    assert [entry["label"] for entry in entries] == ["Holds"]
    assert entries[0]["peaks"]["step"] >= HOLD_MB * 0.8
    assert entries[0]["seconds"] >= 3


def test_a_step_that_holds_far_more_memory_than_its_budget_fails_the_smoke_run(held_job):
    entries = json.loads((held_job / job_runner.MEMORY_NAME).read_text(encoding="utf-8"))
    budgets = smoke.parse_budgets([f"step={HOLD_MB // 2}", "host=1000000", "mysql=1000000", "workers=1000000"])   # the MySQL here is real
    over = smoke.over_budget(entries, budgets)
    assert [(label, group) for label, group, _peak, _budget in over] == [("Holds", "step")]
    assert "OVER BUDGET" in smoke.format_report(entries, budgets, over)


def test_a_step_within_its_budget_passes(held_job):
    entries = json.loads((held_job / job_runner.MEMORY_NAME).read_text(encoding="utf-8"))
    budgets = smoke.parse_budgets([f"step={HOLD_MB * 3}", "host=1000000", "mysql=1000000", "workers=1000000"])
    assert smoke.over_budget(entries, budgets) == []
    assert "stayed within its budget" in smoke.format_report(entries, budgets, [])


def test_budgets_start_from_the_tuning_defaults_and_refuse_nonsense():
    assert smoke.parse_budgets([]) == tuning.SMOKE_BUDGET_MB
    assert smoke.parse_budgets(["mysql=900"])["mysql"] == 900
    for bad in ("mysql", "mysql=lots", "mysql=0", "disk=100"):
        with pytest.raises(ValueError):
            smoke.parse_budgets([bad])


def test_a_step_with_a_budget_of_its_own_is_judged_by_it_unless_one_is_given_on_the_command_line():
    entry = {"label": "Scatter the phenomena", "seconds": 1,
             "peaks": {"mysql": 1000.0, "workers": 1000.0, "step": 3000.0, "host": 6000.0}}
    defaults = smoke.parse_budgets([])
    assert smoke.over_budget([entry], defaults) == []           # its own, larger, budget
    assert smoke.over_budget([dict(entry, label="Another step")], defaults) != []
    assert smoke.over_budget([dict(entry, peaks=dict(entry["peaks"], step=6400.0))], defaults) != []   # the old hold
    tight = smoke.parse_budgets(["step=2000"])
    assert [group for _label, group, _peak, _budget in smoke.over_budget([entry], tight)] == ["step"]
