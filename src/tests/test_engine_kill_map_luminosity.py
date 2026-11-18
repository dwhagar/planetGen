"""OPS.225, slice 3 (the real-MySQL kill tests, maintenance writes): the map luminosity filled with units killed at drawn commits is the one an undisturbed fill writes, and so are the block light parts rolled from it (`tests/kill_harness.py`).

Commands: fill-map-luminosity (and the block lights it ends with). The piece size is cut to 2,000 ids so a unit commits several times.
"""

import argparse
import hashlib
import random

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.generation import run_plan
from tests import inline_queue, kill_harness
from tests.test_engine_kill_scatter import SEEDS, _digest, _scatter_all, _wipe

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]


def _fill_args(config):
    parser = argparse.ArgumentParser(prefix_chars="-+")
    generate_cli.add_fill_arguments(parser)
    return parser.parse_args(["--mysql-host", config.host, "--mysql-port", str(config.port), "--mysql-user", config.user,
                              "--mysql-password", config.password, "--mysql-database", config.database])


def _parts(config):
    conn = store.get_connection(config)
    try:
        conn.commit()
        rows = sorted(repr(tuple(row.values())) for row in conn.execute("SELECT * FROM galaxy_block_light_parts").fetchall())
        return len(rows), hashlib.sha256("\n".join(rows).encode()).hexdigest()
    finally:
        conn.close()


def _scatter_blank_fill(config, monkeypatch):
    """Scatters, empties the map luminosity of every scattered row, then fills it."""
    monkeypatch.setattr(run_plan, "MAP_LUMINOSITY_PIECE_IDS", 2000)
    _scatter_all(config)
    conn = store.get_connection(config)
    try:
        with conn:
            conn.execute("UPDATE phenomenon_scatter SET map_luminosity_sol = NULL")
    finally:
        conn.close()
    with inline_queue.inline_queue():
        run_plan.fill_map_luminosity(_fill_args(config))
    return _digest(config), _parts(config)


@pytest.fixture
def undisturbed(mysql_config):
    import _pytest.monkeypatch
    counter = kill_harness.Plan()
    patch = _pytest.monkeypatch.MonkeyPatch()
    try:
        kill_harness.install(patch, counter)
        found = _scatter_blank_fill(mysql_config, patch)
    finally:
        patch.undo()
    return found, counter.seen


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_a_map_luminosity_filled_with_killed_units_is_the_one_an_undisturbed_fill_writes(mysql_config, undisturbed, monkeypatch, seed):
    expected, commits = undisturbed
    assert commits > 0 and expected[0]["phenomenon_scatter"][0] > 0, "the undisturbed run filled something to compare with"
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(commits), min(commits, rng.randint(1, 4)))}
    _wipe(mysql_config)
    plan = kill_harness.Plan(kills)
    kill_harness.install(monkeypatch, plan)
    found = _scatter_blank_fill(mysql_config, monkeypatch)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (commits {commits}, plan {kills})"
    assert found == expected, f"seed {seed}: kills {plan.fired} changed the map luminosity"
