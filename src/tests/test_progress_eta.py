# tests/test_progress_eta.py

"""
PERF.7: the decaying-average rate behind every generation progress bar's
ETA (`planetgen/queue/progress_rate.py`), the bar column that shows it, the
same rate and ETA in the web job's progress file and on the Generate
page, and the bright-star scatter drawn layer by layer on several
workers giving the same stars as one.

Tests that take the `mysql_config` fixture (see `conftest.py`) are
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import json

import pytest

import generate
from planetgen.queue import progress_file
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.queue.progress_rate import DecayingRate

from tests.test_bright_star_scatter import EDGE_PC, E_VALUE, EXTENTS, SHAPE, THRESHOLD, _plan_args, _seed_galaxy


class _Clock:
    def __init__(self):
        self.now = 0.0

    def __call__(self):
        return self.now


def _rate(tau=60.0):
    clock = _Clock()
    return DecayingRate(tau=tau, clock=clock), clock


def test_a_steady_pace_gives_that_rate_and_a_plain_eta():
    rate, clock = _rate()
    assert rate.eta(10) is None
    for _ in range(20):
        clock.now += 2.0
        rate.add(1)
    assert rate.rate == pytest.approx(0.5)
    assert rate.eta(10) == pytest.approx(20.0)
    assert rate.eta(0) == 0.0
    assert rate.eta(None) is None


def test_the_eta_counts_down_between_updates_and_holds_when_overdue():
    rate, clock = _rate()
    for _ in range(5):
        clock.now += 1.0
        rate.add(1)
    assert rate.eta(10) == pytest.approx(10.0)
    clock.now += 0.5
    assert rate.eta(10) == pytest.approx(9.5)
    clock.now += 30
    assert rate.eta(10) == pytest.approx(9.0)


def test_a_burst_moves_the_average_no_more_than_the_same_work_spread_out():
    burst, burst_clock = _rate()
    even, even_clock = _rate()
    for clock, rate in ((burst_clock, burst), (even_clock, even)):
        for _ in range(30):
            clock.now += 1.0
            rate.add(1)
    burst_clock.now += 4.0
    for _ in range(4):
        burst.add(1)
    for _ in range(4):
        even_clock.now += 1.0
        even.add(1)
    assert burst.rate == pytest.approx(even.rate, rel=0.02)


def test_the_average_follows_a_new_pace_at_the_time_constant():
    rate, clock = _rate(tau=10.0)
    for _ in range(100):
        clock.now += 1.0
        rate.add(1)
    for _ in range(10):
        clock.now += 1.0
        rate.add(3)
    # After one time constant at 3/s, about 63% of the way from 1 to 3.
    assert rate.rate == pytest.approx(1 + 2 * (1 - 2.718281828 ** -1), rel=0.01)


def test_the_remaining_column_shows_the_decaying_eta():
    progress = generate._generation_progress()
    task_id = progress.add_task("Sectors", total=10)
    task = progress.tasks[0]
    column = generate._DecayingRemainingColumn()
    assert column.render(task).plain == "-:--:--"
    clock = _Clock()
    rate = task.fields["rate"]
    rate.clock, rate.last = clock, 0.0
    for _ in range(4):
        clock.now += 30.0
        rate.add(1)
    progress.update(task_id, completed=4)
    assert column.render(task).plain.startswith("0:0")
    progress.update(task_id, completed=10)
    assert column.render(task).plain == "0:00:00"


def test_the_progress_file_carries_the_rate_and_eta(tmp_path, monkeypatch):
    path = tmp_path / "progress.json"
    monkeypatch.setenv(progress_file.ENV_VAR, str(path))
    progress_file.report(3, 10, "Sectors", force=True, rate=0.5, eta_s=14.0)
    body = json.loads(path.read_text())
    assert (body["completed"], body["total"], body["rate"], body["eta_s"]) == (3, 10, 0.5, 14.0)
    progress_file.report(0, 10, "Sectors", force=True)
    body = json.loads(path.read_text())
    assert body["rate"] is None and body["eta_s"] is None


def test_the_generate_page_counts_the_eta_down(monkeypatch):
    from planetgen.web import generate_page as generate_page

    monkeypatch.setattr(generate_page, "url_for", lambda *a, **k: "/job")
    monkeypatch.setattr(generate_page.time, "time", lambda: 1000.0)
    job = {"id": "x", "finished": False, "created_at": 900.0, "elapsed_s": 100,
           "progress": {"completed": 3, "total": 10, "eta_s": 200.0, "updated_at": 990.0}}
    assert generate_page._job_view(job)["remaining_text"] == "3 m 10 s"
    assert generate_page._job_view({**job, "finished": True})["remaining_text"] == ""
    assert generate_page._job_view({**job, "progress": {"completed": 0}})["remaining_text"] == ""


# ---------------------------------------------------------------------------
# The bright-star scatter, one layer per task
# ---------------------------------------------------------------------------

def test_layers_drawn_in_any_order_give_the_same_stars():
    together = brightStars.scatter(SHAPE, EXTENTS, EDGE_PC, E_VALUE, THRESHOLD, 21)
    apart = [row for layer, outer in reversed(EXTENTS)
             for row in brightStars.scatter_layer(SHAPE, layer, outer, EDGE_PC, E_VALUE, THRESHOLD, 21)]
    assert sorted(together) == sorted(apart)


def _stored_stars(config):
    conn = store.get_connection(config)
    try:
        columns = ", ".join(store.BRIGHT_STAR_COLUMNS)
        return sorted(tuple(row.values()) for row in conn.execute(f"SELECT {columns} FROM bright_stars").fetchall())
    finally:
        conn.close()


def test_a_parallel_scatter_places_the_same_stars_as_a_serial_one(mysql_config, monkeypatch):
    monkeypatch.setenv("PLANETGEN_CONTROL_DATABASE", mysql_config.database)
    # The galaxy's seed fixes the scatter's (GEN.39), so both draw alike.
    _seed_galaxy(mysql_config)
    serial = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--workers", "1"))
    serial_rows = _stored_stars(mysql_config)
    parallel = generate.scatter_bright_stars(_plan_args(mysql_config, "--bright-stars-only", "--workers", "2"))
    assert serial["total"] == parallel["total"] == len(serial_rows) > 20
    assert serial["counts"] == parallel["counts"]
    assert _stored_stars(mysql_config) == serial_rows
