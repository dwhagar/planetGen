"""OPS.225, slice 3 (the real-MySQL kill tests, sector making): a galaxy of sectors made with units killed at drawn commits ends as the galaxy of an undisturbed run (`tests/kill_harness.py`).

Commands: what a `planetgen galaxy` run does: generate-sectors, generate-system, settle-sector-paths and the stages around them. (A `sector` run of unplaced sectors draws its own place, so two undisturbed runs differ.) Replay a seed with `pytest -k seed_3`."""

import random
import sys

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.queue import runrecord, redisqueue
from tests import kill_harness
from tests.bughunt_support import mysql_argv

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]

SEEDS = ["control", *range(2)]
"""The first case kills nothing: two undisturbed runs make the same galaxy, which is what the others compare with."""
ARGV = ["--shell", "1", "--limit", "6", "--num-systems", "1"]


def _generate(config):
    from tests import inline_queue
    from tests.test_galaxy_gen import _seed_skeleton
    _seed_skeleton(config, shells=[(0, 0), (1, 1)])        # the galaxy's seed is stored, so two runs make the same galaxy
    old_argv = sys.argv
    try:
        sys.argv = ["planetgen", "galaxy", *ARGV, *mysql_argv(config)]
        with inline_queue.inline_queue():
            generate_cli.main()
    finally:
        sys.argv = old_argv


def _wipe(config):
    from planetgen.cli import reset as reset_galaxy
    assert reset_galaxy.reset_database(config, assume_yes=True)


@pytest.fixture
def undisturbed(mysql_config):
    """The galaxy an undisturbed run makes, and how many commits it took."""
    counter = kill_harness.Plan()
    import _pytest.monkeypatch
    patch = _pytest.monkeypatch.MonkeyPatch()
    try:
        kill_harness.install(patch, counter)
        _generate(mysql_config)
    finally:
        patch.undo()
    return kill_harness.digest_of_galaxy(mysql_config), counter.seen


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_sectors_made_with_killed_units_are_the_sectors_an_undisturbed_run_makes(mysql_config, undisturbed, monkeypatch, seed):
    expected, commits = undisturbed
    assert commits > 0 and expected["sectors"][0] > 0, "the undisturbed run made a galaxy to compare with"
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(commits), min(commits, rng.randint(1, 4)))}
    _wipe(mysql_config)
    plan = kill_harness.Plan(kills)
    kill_harness.install(monkeypatch, plan)
    _generate(mysql_config)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (commits {commits}, plan {kills})"
    found = kill_harness.digest_of_galaxy(mysql_config)
    assert found == expected, f"seed {seed}: kills {plan.fired} changed the galaxy"
