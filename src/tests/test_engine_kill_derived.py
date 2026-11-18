"""OPS.225, slice 3 (the real-MySQL kill tests, maintenance writes): the habitability scores and the effect ledger refilled with units killed at drawn commits are the ones an undisturbed refill writes (`tests/kill_harness.py`).

Commands: rescore-habitability, refill-effect-ledger. Systems whose stamp is old are the ones a refill takes, so a retried unit goes on with the systems it left.
"""

import hashlib
import random

import pytest

from planetgen.db import store
from planetgen.generation import run_plan
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.galaxy import bring_to_epoch
from planetgen.physics import habitability_world
from planetgen.util import draw
from tests import inline_queue, kill_harness

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]

SEEDS = ["control", *range(3)]
SYSTEMS = 5


def _systems():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    made = []
    for seed in range(200):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if system.planets:
            made.append(system)
        if len(made) == SYSTEMS:
            return made, cfg
    pytest.fail("not enough systems with a planet")


def _stale_galaxy(config):
    """Systems saved with their scores and ledger blanked and their stamps old, as a galaxy made by older code is."""
    systems, cfg = _systems()
    conn = store.get_connection(config)
    try:
        with conn:
            sector_id = conn.execute("INSERT INTO sectors (name, edge_mpc, generator_epoch) VALUES ('Derived sector', 1000, 1)").lastrowid
            for system in systems:
                system_id = store.insert_star_system(conn, system, cfg)
                conn.execute("UPDATE star_systems SET sector_id = ? WHERE id = ?", (sector_id, system_id))
                store.refresh_habitability(conn, system_id)
                store.refresh_effect_ledger(conn, system_id)
            blank = ", ".join(f"{name} = NULL" for name in habitability_world.SCORE_FIELDS)
            for table in ("planets", "moons"):
                conn.execute(f"UPDATE {table} SET {blank}")
            conn.execute("DELETE FROM effect_ledger")
            conn.execute("UPDATE system_summary SET habitability_stamp = 1, ledger_stamp = 1")
    finally:
        conn.close()


def _refill(config):
    _stale_galaxy(config)
    with inline_queue.inline_queue():
        for family, run in (("habitability", run_plan.run_rescore_habitability), ("ledger", run_plan.run_refill_effect_ledger)):
            run(bring_to_epoch.derived_args(config, family))


def _digest(config):
    """The scores of every planet and moon, the ledger rows (no row ids) and the stamps, by the system's own columns."""
    conn = store.get_connection(config)
    try:
        conn.commit()
        found = {}
        fields = ", ".join(("flux_variation", *habitability_world.SCORE_FIELDS))
        for table in ("planets", "moons"):
            rows = sorted(repr(tuple(row.values())) for row in conn.execute(
                f"SELECT s.cell_key, s.name AS system_name, b.{fields.replace(', ', ', b.')} FROM {table} b JOIN star_systems s ON s.id = b.star_system_id").fetchall())
            found[table] = (len(rows), hashlib.sha256("\n".join(rows).encode()).hexdigest())
        ledger = sorted(repr(tuple(v for k, v in row.items() if k not in ("id",))) for row in conn.execute(
            "SELECT s.name AS system_name, e.* FROM effect_ledger e JOIN star_systems s ON s.id = e.star_system_id").fetchall())
        found["ledger"] = (len(ledger), hashlib.sha256("\n".join(ledger).encode()).hexdigest())
        found["stamps"] = [tuple(row.values()) for row in conn.execute(
            "SELECT habitability_stamp, ledger_stamp FROM system_summary ORDER BY star_system_id").fetchall()]
        return found
    finally:
        conn.close()


def _wipe(config):
    from planetgen.cli import reset as reset_galaxy
    assert reset_galaxy.reset_database(config, assume_yes=True)


@pytest.fixture
def undisturbed(mysql_config):
    import _pytest.monkeypatch
    counter = kill_harness.Plan()
    patch = _pytest.monkeypatch.MonkeyPatch()
    try:
        kill_harness.install(patch, counter)
        _refill(mysql_config)
    finally:
        patch.undo()
    return _digest(mysql_config), counter.seen


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_derived_values_refilled_with_killed_units_are_the_values_an_undisturbed_refill_writes(mysql_config, undisturbed, monkeypatch, seed):
    expected, commits = undisturbed
    assert commits > 0 and expected["planets"][0] > 0 and expected["moons"][0] > 0, "the undisturbed refill wrote something to compare with"
    assert all(stamp[0] > 1 and stamp[1] > 1 for stamp in expected["stamps"]), "every stamp is current"
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(commits), min(commits, rng.randint(1, 4)))}
    _wipe(mysql_config)
    plan = kill_harness.Plan(kills)
    kill_harness.install(monkeypatch, plan)
    _refill(mysql_config)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (commits {commits}, plan {kills})"
    assert _digest(mysql_config) == expected, f"seed {seed}: kills {plan.fired} changed the refill"
