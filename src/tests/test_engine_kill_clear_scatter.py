"""OPS.225, slice 1 (the real-MySQL kill tests, the clear-scatter command): the scatter cleared a kind at a time with units killed at drawn commits leaves what an undisturbed clear leaves (`tests/kill_harness.py`).

Commands: clear-scatter (the phenomena, the luminosity pass and the mass pass, in the order a stored scatter can be taken back).
"""

import argparse
import random

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.generation import run_plan
from tests import inline_queue, kill_harness
from tests.test_engine_kill_scatter import SEEDS, _digest, _scatter_all, _wipe

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]

KINDS = ("phenomena", "luminosity", "mass")


def _clear_args(config, kind):
    parser = argparse.ArgumentParser(prefix_chars="-+")
    generate_cli.add_clear_scatter_arguments(parser)
    return parser.parse_args(["--mysql-host", config.host, "--mysql-port", str(config.port), "--mysql-user", config.user,
                              "--mysql-password", config.password, "--mysql-database", config.database, f"--{kind}"])


def _scatter_then_clear(config):
    """Scatters, then clears the kinds one at a time; what is left after each clear."""
    _scatter_all(config)
    left = []
    for kind in KINDS:
        with inline_queue.inline_queue():
            run_plan.clear_scatter(_clear_args(config, kind))
        left.append((kind, _digest(config)))
    return left


@pytest.fixture
def undisturbed(mysql_config):
    import _pytest.monkeypatch
    counter = kill_harness.Plan()
    patch = _pytest.monkeypatch.MonkeyPatch()
    try:
        kill_harness.install(patch, counter)
        left = _scatter_then_clear(mysql_config)
    finally:
        patch.undo()
    return left, counter.seen


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_a_scatter_cleared_with_killed_units_is_cleared_as_an_undisturbed_clear_clears_it(mysql_config, undisturbed, monkeypatch, seed):
    expected, commits = undisturbed
    assert commits > 0 and expected[0][1]["phenomenon_scatter"][0] == 0 and expected[-1][1]["bright_stars"][0] == 0, "the undisturbed run cleared the scatter"
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(commits), min(commits, rng.randint(1, 4)))}
    _wipe(mysql_config)
    plan = kill_harness.Plan(kills)
    kill_harness.install(monkeypatch, plan)
    found = _scatter_then_clear(mysql_config)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (commits {commits}, plan {kills})"
    assert found == expected, f"seed {seed}: kills {plan.fired} changed what the clears leave"
