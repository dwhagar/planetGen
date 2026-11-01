# tests/test_generator_battery.py

"""OPS.57: `generator_battery_digest()` is a database-free digest over eight canned cases that is stable from run to
run and from hash seed to hash seed, and changing one canned draw changes the digest of exactly one case."""

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from planetgen.generation import generator_battery as battery

SRC = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def baseline():
    return battery.generator_battery_digest()


def test_the_battery_has_the_eight_cases_in_order(baseline):
    assert tuple(baseline.cases) == battery.CASES
    assert len(battery.CASES) == 8
    assert all(len(digest) == 64 for digest in [baseline.overall] + [c.full for c in baseline.cases.values()])


def test_two_runs_give_the_same_digests(baseline):
    again = battery.generator_battery_digest()
    assert again == baseline
    assert again.first_difference(baseline) is None


def test_the_cases_are_not_empty():
    assert len(battery.case_value("bright-star-scatter")) > 0
    assert len(battery.case_value("population-stream")) == battery.POPULATION_STARS
    nebula = battery.case_value("nebula-shape")
    assert len(nebula["inside"]) == battery.NEBULA_POINTS and any(nebula["inside"]) and not all(nebula["inside"])
    assert len(battery.case_value("sector-a")["systems"]) == battery.SECTOR_SYSTEMS


def test_the_unit_seed_case_is_the_known_answer():
    galaxy_seed, kind, address = battery.KNOWN_UNIT_SEED
    expected = int.from_bytes(hashlib.sha256(galaxy_seed + f"{kind}:{'/'.join(map(str, address))}".encode()).digest(),
                              "big")
    assert battery.case_value("unit-seed") == expected


@pytest.mark.parametrize("name, constant, value", [
    ("unit-seed", "KNOWN_UNIT_SEED", (bytes(range(16)), "sector", (3, 1, 1))),
    ("nebula-shape", "NEBULA_SEED", 0x0EB1),
    ("sector-b", "SECTOR_ADDRESSES", ((7, 1, 2), (12, -2, 1))),
    ("population-stream", "POPULATION_SEED", 0x9092),
])
def test_changing_one_canned_draw_changes_exactly_one_case(baseline, monkeypatch, name, constant, value):
    monkeypatch.setattr(battery, constant, value)
    changed = battery.generator_battery_digest()
    assert [case for case in battery.CASES if changed.cases[case].full != baseline.cases[case].full] == [name]
    assert changed.overall != baseline.overall
    assert changed.first_difference(baseline) == name


def test_a_float_alone_changing_leaves_the_structural_digest(baseline, monkeypatch):
    monkeypatch.setattr(battery, "SCATTER_THRESHOLD_SOL", battery.SCATTER_THRESHOLD_SOL * 1.0000001)
    changed = battery.generator_battery_digest()
    case = changed.cases["population-stream"]
    assert case.full != baseline.cases["population-stream"].full
    assert case.structural == baseline.cases["population-stream"].structural
    assert case.floats != baseline.cases["population-stream"].floats


@pytest.mark.slow
def test_the_digest_does_not_depend_on_the_hash_seed(baseline):
    code = "from planetgen.generation import generator_battery as b; print(b.generator_battery_digest().overall)"
    outputs = set()
    for hash_seed in ("0", "1"):
        env = dict(os.environ, PYTHONHASHSEED=hash_seed,
                   PYTHONPATH=os.pathsep.join([str(SRC), os.environ.get("PYTHONPATH", "")]))
        done = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, timeout=600)
        assert done.returncode == 0, done.stderr[-2000:]
        outputs.add(done.stdout.strip().splitlines()[-1])
    assert outputs == {baseline.overall}


def test_tree_turns_objects_into_fields_and_sorts_sets():
    class Thing:
        def __init__(self):
            self.name = "x"
            self.values = {3, 1, 2}
            self.epoch_unix = 123.0
            self.array = np.array([1.5, 2.5])
            self.count = np.int64(4)
            self.self = self

    found = battery.tree(Thing())
    assert found["class"] == "Thing" and found["values"] == [1, 2, 3]
    assert "epoch_unix" not in found
    assert found["array"] == [1.5, 2.5] and found["count"] == 4 and type(found["count"]) is int
    assert found["self"] == "<Thing>"


def test_tree_leaves_the_clocks_out_of_dicts():
    assert battery.tree({"epoch_unix": 1.0, "next_update_due": 2.0, "kept": 3}) == {"kept": 3}


def test_the_digest_started_in_a_process_of_its_own_is_the_one_made_here(baseline, monkeypatch):
    """PERF.176: a plan works the digest out while it does other work; the digest it then stores is the same."""
    monkeypatch.setattr(battery, "_running", None)
    finish = battery.start_running_digest()
    assert battery._running is None, "nothing is kept until the caller waits for it"
    finish()
    assert battery._running == baseline
    assert battery.running_digest() is battery._running
    assert battery.start_running_digest()() is None, "a digest already kept starts nothing"


def test_a_digest_that_cannot_be_started_elsewhere_is_worked_out_where_it_is_asked_for(baseline, monkeypatch):
    def refuse(*args, **kwargs):
        raise OSError("no process")

    monkeypatch.setattr(battery, "_running", None)
    monkeypatch.setattr(battery.concurrent.futures, "ProcessPoolExecutor", refuse)
    battery.start_running_digest()()
    assert battery._running is None
    assert battery.running_digest() == baseline


def test_the_progress_callback_reports_every_case_and_changes_nothing(monkeypatch):
    """The update shows a bar: `progress(done, total)` before the first case and after each, with the digest unchanged."""
    monkeypatch.setattr(battery, "case_value", lambda name: {"case": name})
    seen = []
    with_bar = battery.generator_battery_digest(lambda done, total: seen.append((done, total)))
    total = len(battery.CASES)
    assert seen == [(done, total) for done in range(total + 1)]
    assert with_bar == battery.generator_battery_digest()


def test_the_running_digest_passes_the_progress_callback_on(monkeypatch):
    monkeypatch.setattr(battery, "_running", None)
    monkeypatch.setattr(battery, "case_value", lambda name: {"case": name})
    monkeypatch.delenv("PLANETGEN_BATTERY_DIGEST_CACHE", raising=False)
    seen = []
    battery.running_digest(lambda done, total: seen.append(done))
    assert seen[-1] == len(battery.CASES)
