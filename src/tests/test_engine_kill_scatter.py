"""OPS.225, slice 1 (the real-MySQL kill tests, the scatter family): the stars and phenomena scattered with units killed at drawn commits are the ones an undisturbed run scatters (`tests/kill_harness.py`).

Commands: clear-bright-stars, scatter-mass, scatter-luminosity, refresh-block-lights, scatter-band, empty-phenomenon-scatter, scatter-phenomena, record-phenomenon-scatter. Replay a seed with `pytest -k seed_3`."""

import hashlib
import random

import pytest

from planetgen.db import store
from planetgen.generation import run_plan
from tests import inline_queue, kill_harness
from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]

SEEDS = ["control", *range(3)]
"""The first case kills nothing: it checks that two undisturbed runs scatter the same, which is what the others compare with."""
WORKERS = "2"
SKIPPED = ("id", "uid", "built_at", "epoch_unix", "next_update_due")
"""The columns of a stored row that the database or the clock gives: a reset keeps the Object ID counters (DB.20), so the same stars get new IDs."""


def _scatter_all(config):
    _seed_galaxy(config)
    with inline_queue.inline_queue():
        run_plan.scatter_bright_stars(_scatter_args(config, "--mass", "--workers", WORKERS))
        run_plan.add_bright_star_band(_scatter_args(config, "--band", "--down-to", "100", "--workers", WORKERS))
        run_plan.scatter_phenomena(_scatter_args(config, "--phenomena", "--workers", WORKERS))


def _digest(config):
    """The scattered rows by their own columns (no row ids, Object IDs, which a reset does not take back, or times of day), in an order of their own: (count, sha256) of the bright stars and of the phenomena."""
    conn = store.get_connection(config)
    try:
        found = {}
        for name, query in (("bright_stars", f"SELECT {', '.join(store.BRIGHT_STAR_COLUMNS)} FROM bright_stars"),
                            ("phenomenon_scatter", "SELECT * FROM phenomenon_scatter")):
            rows = sorted(repr(tuple(value for key, value in row.items() if key not in SKIPPED)) for row in conn.execute(query).fetchall())
            found[name] = (len(rows), hashlib.sha256("\n".join(rows).encode()).hexdigest())
        return found
    finally:
        conn.close()


def _wipe(config):
    from planetgen.cli import reset as reset_galaxy
    assert reset_galaxy.reset_database(config, assume_yes=True)


@pytest.fixture
def undisturbed(mysql_config):
    """What an undisturbed run scatters, and how many commits it took."""
    import _pytest.monkeypatch
    counter = kill_harness.Plan()
    patch = _pytest.monkeypatch.MonkeyPatch()
    try:
        kill_harness.install(patch, counter)
        _scatter_all(mysql_config)
    finally:
        patch.undo()
    return _digest(mysql_config), counter.seen


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_scatters_made_with_killed_units_are_the_scatters_an_undisturbed_run_makes(mysql_config, undisturbed, monkeypatch, seed):
    expected, commits = undisturbed
    assert commits > 0 and expected["bright_stars"][0] > 0 and expected["phenomenon_scatter"][0] > 0, "the undisturbed run scattered something to compare with"
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(commits), min(commits, rng.randint(1, 4)))}
    _wipe(mysql_config)
    plan = kill_harness.Plan(kills)
    kill_harness.install(monkeypatch, plan)
    _scatter_all(mysql_config)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (commits {commits}, plan {kills})"
    assert _digest(mysql_config) == expected, f"seed {seed}: kills {plan.fired} changed the scatter"
