"""OPS.225, slice 1 (the real-MySQL kill tests, the sky backfill): a galaxy run whose charting draws the sky (`scatter-backfill` units) with units killed at drawn commits leaves the stars, the cells' levels and the ledger an undisturbed run leaves (`tests/kill_harness.py`).

Commands: scatter-backfill (with the `generate-sectors` run around it). The toy tiers of `test_sky_sliced_draw` keep the sphere a few sectors across.
"""

import hashlib
import random
import sys

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import sky_tier_ledger as ledger_module
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import run_galaxy, sky_charting
from tests import inline_queue, kill_harness
from tests.bughunt_support import mysql_argv
from tests.test_sky_sliced_draw import PARTS, toy_tiers  # noqa: F401 (the fixture)

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]

SEEDS = ["control", *range(2)]
ARGV = ["--shell", "1", "--limit", "3", "--num-systems", "1", "--sky-limit", "7.5"]
SKIPPED = ("id", "uid", "built_at", "epoch_unix", "next_update_due")


def _generate(config):
    from tests.test_galaxy_gen import _seed_skeleton
    _seed_skeleton(config, shells=[(0, 0), (1, 1)])
    old_argv = sys.argv
    try:
        sys.argv = ["planetgen", "galaxy", *ARGV, *mysql_argv(config)]
        with inline_queue.inline_queue():
            generate_cli.main()
    finally:
        sys.argv = old_argv


def _digest(config):
    conn = store.get_connection(config)
    try:
        found = kill_harness.digest_of_galaxy(config)
        rows = sorted(repr(tuple(v for k, v in r.items() if k not in SKIPPED))
                      for r in conn.execute(f"SELECT {', '.join(store.BRIGHT_STAR_COLUMNS)} FROM bright_stars").fetchall())
        found["bright_stars"] = (len(rows), hashlib.sha256("\n".join(rows).encode()).hexdigest())
        cells = sorted((r["shell_index"], r["band_index"], r["wedge_index"], r["bright_level_sol"], r["bright_mass_sol"])
                       for r in conn.execute("SELECT shell_index, band_index, wedge_index, bright_level_sol, bright_mass_sol FROM sector_stats").fetchall())
        found["cells"] = (len(cells), hashlib.sha256(repr(cells).encode()).hexdigest())
        found["spheres"] = sorted((d.kind, d.low, d.high, d.center, d.radius_pc) for d in ledger_module.read(conn).draws)
        return found
    finally:
        conn.close()


def _wipe(config):
    from planetgen.cli import reset as reset_galaxy
    assert reset_galaxy.reset_database(config, assume_yes=True)


@pytest.fixture
def undisturbed(mysql_config, toy_tiers, monkeypatch):  # noqa: F811
    counter = kill_harness.Plan()
    import _pytest.monkeypatch
    patch = _pytest.monkeypatch.MonkeyPatch()
    try:
        kill_harness.install(patch, counter)
        _generate(mysql_config)
    finally:
        patch.undo()
    return _digest(mysql_config), counter


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_a_sky_drawn_with_killed_units_is_the_sky_an_undisturbed_run_draws(mysql_config, toy_tiers, undisturbed, monkeypatch, seed):  # noqa: F811
    expected, counter = undisturbed
    commits = counter.seen
    assert commits > 0 and expected["bright_stars"][0] > 0 and expected["spheres"], "the undisturbed run drew a sky to compare with"
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(commits), min(commits, rng.randint(1, 4)))}
    _wipe(mysql_config)
    plan = kill_harness.Plan(kills)
    kill_harness.install(monkeypatch, plan)
    _generate(mysql_config)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (commits {commits}, plan {kills})"
    assert _digest(mysql_config) == expected, f"seed {seed}: kills {plan.fired} changed the sky"
