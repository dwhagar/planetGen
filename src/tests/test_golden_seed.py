# tests/test_golden_seed.py

"""
TEST.77: a fixed seed makes what the pinned golden files say (`src/tests/golden/`, written by `tests/golden_seed.py`).

- Tier A, no database: the battery's digest of each case (OPS.57) against `battery.json`; the Rosetta text against
  `rosetta.txt` (a failure there is read as a diff).
- Tier B, database: a small galaxy of a fixed 128-bit seed (the plan, four sectors of four systems, the bright-star
  backfill) made at one worker: its digest of each sector, each section of it and the region equal `galaxy.json`. It is
  built once per run, and a rewrite of the golden files (`golden_update`) is that one build.
- Tier B at four workers, the end-of-version run only (`scripts/test.py full` sets `PLANETGEN_VERSION_GATE`): the same
  galaxy at four workers equals the one-worker galaxy, whatever the golden says. It builds the galaxy a second time,
  about as long again as the first, so an ordinary merge's gate leaves it to the version gate.
- Tier C, the hash-seed and locale probe, is `test_reproducible_draws.test_a_sector_does_not_depend_on_hash_seed_or_locale`.

A change of what the generator makes fails these until the golden files are rewritten with
`python -m tests.golden_update --reason "..."`, which also bumps the generator epoch (OPS.56); `bump_version.py
--check-pr` then wants the note's `generation-output: changed` line.
"""

import difflib
import os

import pytest

from planetgen.db import store

from tests import golden_seed
from tests.test_galaxy_gen import _plan_wide_galaxy
from tests.test_parallel_generation import _run, make_database  # noqa: F401 -- a fixture


def _check_or_write(path, actual, written=None):
    """Compares `actual` with the golden file at `path`, or writes the file (`written`, else `actual`) when this run
    rewrites the goldens; a missing file is an error."""
    if golden_seed.writing():
        golden_seed.write(path, actual if written is None else written)
        return None
    assert path.exists(), f"{path.name} is missing: run `python -m tests.golden_update --reason \"...\"`"
    return golden_seed.read(path)


def test_the_battery_makes_the_pinned_digests():
    actual = golden_seed.battery_record()
    expected = _check_or_write(golden_seed.BATTERY_FILE, actual, golden_seed.written_battery())
    if expected is None:
        return
    expected.pop("made_on", None)
    differences = golden_seed.describe_difference(expected, actual)
    assert not differences, "\n".join(differences)
    assert expected["overall"] == actual["overall"]


def test_the_rosetta_text_is_the_pinned_text():
    actual = golden_seed.rosetta_text()
    expected = _check_or_write(golden_seed.ROSETTA_FILE, actual)
    if expected is None:
        return
    diff = "".join(difflib.unified_diff(expected.splitlines(True), actual.splitlines(True), "rosetta.txt (pinned)",
                                        "rosetta.txt (made now)", n=1))
    assert not diff, "\n" + diff[:4000]


VERSION_GATE = "PLANETGEN_VERSION_GATE"
"""str: The environment variable `scripts/test.py full` sets: the end-of-version run, which also builds the golden galaxy
at four workers."""

_ONE_WORKER = {}
"""dict: The one-worker galaxy record of this process (built by the first test that needs it, then reused)."""


def _galaxy_record(make_database, monkeypatch, workers):
    config = make_database()
    _plan_wide_galaxy(config, galaxy_seed=golden_seed.GALAXY_SEED)
    _run("galaxy", [*golden_seed.GALAXY_ARGUMENTS, "--workers", str(workers)], config, monkeypatch)
    conn = store.get_connection(config)
    try:
        return golden_seed.galaxy_record(conn)
    finally:
        conn.close()


def _one_worker_record(make_database, monkeypatch):
    if "record" not in _ONE_WORKER:
        _ONE_WORKER["record"] = _galaxy_record(make_database, monkeypatch, 1)
    return _ONE_WORKER["record"]


def _sector_differences(expected, actual):
    """A line for each sector (and its named sections) whose digest differs, by address."""
    lines = []
    for label in sorted(expected["sectors"].keys() | actual["sectors"].keys()):
        want, got = expected["sectors"].get(label), actual["sectors"].get(label)
        if want is None or got is None:
            lines.append(f"sector {label}: only in the {'golden' if got is None else 'galaxy made now'}")
        elif want["digest"] != got["digest"]:
            sections = [name for name in sorted(want["sections"].keys() | got["sections"].keys())
                        if want["sections"].get(name) != got["sections"].get(name)]
            kinds = ["structure" if (want["sections"].get(name) or [0, 0])[1] != (got["sections"].get(name) or [0, 0])[1]
                     else "numbers" for name in sections]
            lines.append(f"sector {label}: expected {want['digest'][:12]}, got {got['digest'][:12]} "
                         f"(sections: {', '.join(f'{n} {k}' for n, k in zip(sections, kinds))})")
    return lines


@pytest.mark.db
def test_the_pinned_galaxy_is_the_pinned_one(make_database, monkeypatch):
    one = _one_worker_record(make_database, monkeypatch)
    expected = _check_or_write(golden_seed.GALAXY_FILE, one)
    if expected is None:
        return
    assert _sector_differences(expected, one) == []
    assert (expected["plan"], expected["region"]) == (one["plan"], one["region"])


@pytest.mark.db
@pytest.mark.skipif(golden_seed.writing() or not os.environ.get(VERSION_GATE),
                    reason="the end-of-version run only (scripts/test.py full), or set PLANETGEN_VERSION_GATE=1: "
                           "it builds the pinned galaxy a second time")
def test_the_pinned_galaxy_is_made_at_one_worker_and_at_four(make_database, monkeypatch):
    one = _one_worker_record(make_database, monkeypatch)
    four = _galaxy_record(make_database, monkeypatch, 4)
    assert _sector_differences(one, four) == [], "one worker and four make different sectors"
    assert one == four
