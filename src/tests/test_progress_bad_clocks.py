# tests/test_progress_bad_clocks.py

"""
TEST.27: the progress bars' ETA (`progressRate.DecayingRate`) and the web
job's progress file (`progress_file.report`) under bad input -- a clock
stepping backwards, NaN or infinite amounts and times, many units
finishing at one instant, a tiny rate -- and the work queue's worker
processes never writing the progress file (only the run that owns the
bars does).
"""

import json
import math
import os

import pytest

from planetgen.queue import progress_file, work as workQueue
from planetgen.queue.progress_rate import DecayingRate


class _Clock:
    def __init__(self, now=0.0):
        self.now = now

    def __call__(self):
        return self.now


def _steady(rate, clock, steps=10, every=2.0, amount=1):
    for _ in range(steps):
        clock.now += every
        rate.add(amount)


def _rate(start=0.0):
    clock = _Clock(start)
    return DecayingRate(tau=60.0, clock=clock), clock


def _sane(value):
    return value is None or (math.isfinite(value) and value >= 0)


# ---------------------------------------------------------------------------
# A clock that goes backwards
# ---------------------------------------------------------------------------

def test_a_clock_stepping_back_keeps_the_rate_finite_and_positive():
    rate, clock = _rate(100.0)
    _steady(rate, clock)
    before = rate.rate
    clock.now -= 50.0
    rate.add(1)
    rate.add(1)
    assert math.isfinite(rate.rate) and rate.rate > 0
    # The units join the latest interval, as if they finished at its end.
    assert rate.rate >= before
    clock.now += 60.0
    rate.add(1)
    assert math.isfinite(rate.rate) and rate.rate > 0


def test_a_clock_stepping_back_never_raises_the_eta():
    rate, clock = _rate(100.0)
    _steady(rate, clock)
    plain = 10 / rate.rate
    assert rate.eta(10) == pytest.approx(plain)
    clock.now -= 3600.0
    assert rate.eta(10) == pytest.approx(plain)
    assert rate.eta(10, now=-1e9) == pytest.approx(plain)


def test_a_clock_that_never_moves_waits_for_a_real_interval():
    rate, clock = _rate(5.0)
    for _ in range(10):
        rate.add(1)
    assert rate.rate is None and rate.eta(10) is None
    clock.now = 10.0
    rate.add(1)
    # All eleven units over the five seconds since the bar started.
    assert rate.rate == pytest.approx(11 / 5.0)


# ---------------------------------------------------------------------------
# NaN and infinite amounts and times
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("amount", [math.nan, math.inf, -math.inf, 0, -3, None, "2"])
def test_a_bad_amount_is_ignored(amount):
    rate, clock = _rate()
    _steady(rate, clock)
    before, last = rate.rate, rate.last
    clock.now += 2.0
    rate.add(amount)
    assert (rate.rate, rate.last) == (before, last)
    clock.now += 2.0
    rate.add(1)
    assert math.isfinite(rate.rate) and rate.rate > 0
    assert _sane(rate.eta(10))


@pytest.mark.parametrize("when", [math.nan, math.inf, -math.inf])
def test_a_bad_time_is_ignored(when):
    rate, clock = _rate()
    _steady(rate, clock)
    before, last = rate.rate, rate.last
    rate.add(1, now=when)
    assert (rate.rate, rate.last) == (before, last)
    assert _sane(rate.eta(10, now=when))
    clock.now += 2.0
    rate.add(1)
    assert math.isfinite(rate.rate)


def test_a_nan_clock_never_poisons_the_rate():
    rate, clock = _rate()
    _steady(rate, clock)
    clock.now = math.nan
    rate.add(1)
    assert math.isfinite(rate.rate)
    assert _sane(rate.eta(10))


@pytest.mark.parametrize("remaining", [math.nan, math.inf])
def test_a_bad_remaining_count_has_no_eta(remaining):
    rate, clock = _rate()
    _steady(rate, clock)
    assert rate.eta(remaining) is None


# ---------------------------------------------------------------------------
# Many units at one instant, a tiny rate
# ---------------------------------------------------------------------------

def test_many_adds_at_one_instant_equal_one_add_of_their_sum():
    many, clock_many = _rate()
    one, clock_one = _rate()
    _steady(many, clock_many)
    _steady(one, clock_one)
    clock_many.now += 2.0
    clock_one.now += 2.0
    for _ in range(1000):
        many.add(0.5)
    one.add(500)
    assert many.rate == pytest.approx(one.rate)
    assert many.eta(1000) == pytest.approx(one.eta(1000))


def test_a_tiny_rate_gives_a_huge_but_finite_eta():
    rate, clock = _rate()
    clock.now += 1e9
    rate.add(1e-6)
    assert rate.rate == pytest.approx(1e-15)
    eta = rate.eta(1e6)
    assert math.isfinite(eta) and eta > 1e20


def test_a_rate_that_underflows_to_zero_has_no_eta():
    rate, clock = _rate()
    clock.now += 1e300
    rate.add(1e-300)
    assert rate.rate == 0.0
    assert rate.eta(5) is None


# ---------------------------------------------------------------------------
# The progress file
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_the_progress_file_never_holds_nan_or_infinity(tmp_path, monkeypatch, bad):
    path = tmp_path / "progress.json"
    monkeypatch.setenv(progress_file.ENV_VAR, str(path))
    progress_file.report(3, 10, "Sectors", force=True, rate=bad, eta_s=bad)
    text = path.read_text()
    # Strict JSON, as the browser's JSON.parse reads it.
    body = json.loads(text, parse_constant=lambda name: pytest.fail(f"progress file holds {name}"))
    assert body["rate"] is None and body["eta_s"] is None
    assert body["completed"] == 3


def test_the_progress_file_skips_a_write_it_cannot_encode(tmp_path, monkeypatch):
    path = tmp_path / "progress.json"
    monkeypatch.setenv(progress_file.ENV_VAR, str(path))
    progress_file.report(1, 10, "Sectors", force=True)
    progress_file.report(2, 10, "Sectors", force=True, detail={"eta_s": math.nan})
    assert json.loads(path.read_text())["completed"] == 1
    assert sorted(os.listdir(tmp_path)) == ["progress.json"]


# ---------------------------------------------------------------------------
# Workers never write the progress file
# ---------------------------------------------------------------------------

def _report_from_a_worker(payload):
    progress_file.report(payload, 10, "From a worker", force=True)
    return os.environ.get(progress_file.ENV_VAR)


def test_workers_never_write_the_progress_file(tmp_path, monkeypatch):
    path = tmp_path / "progress.json"
    monkeypatch.setenv(progress_file.ENV_VAR, str(path))
    seen = []
    with workQueue.WorkQueue("progress", workers=2) as queue:
        for n in range(4):
            queue.submit("report", n, _report_from_a_worker, n,
                         on_done=lambda result, _s, _w: seen.append(result))
    assert seen == [None] * 4
    assert not path.exists()
    # The run itself still writes it.
    assert os.environ[progress_file.ENV_VAR] == str(path)
    progress_file.report(4, 4, "Sectors", force=True)
    assert json.loads(path.read_text())["completed"] == 4
