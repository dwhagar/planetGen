"""OPS.225, slice 3 (the real-MySQL kill tests, maintenance writes): the species, civilizations, technology and polities made with units killed at drawn commits are the ones an undisturbed pass makes (`tests/kill_harness.py`).

Commands: scan-species, date-civilizations, refresh-technology, count-population and run-territories. Every world that reaches the civilization stage gets a civilization (`population_support`), so the passes have something to write.
"""

import hashlib
import random

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.generation import run_population
from tests import inline_queue, kill_harness
from tests.population_support import every_window_has_a_civilization
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.test_population import _civilized_system, _place

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]

SEEDS = ["control", *range(3)]


def _args(config, command):
    _parser, parsers = generate_cli.build_parser()
    return parsers[command].parse_args(["--mysql-host", config.host, "--mysql-port", str(config.port), "--mysql-user", config.user,
                                        "--mysql-password", config.password, "--mysql-database", config.database])


def _plain_system(config, number):
    """A system with no life, drawn from a stream of its own number so every build of the galaxy holds the same one (the fixture of test_population draws it unbound)."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M5V"
    cfg.INTELLIGENT_LIFE = False
    cfg.BINARY_SYSTEM = False
    with draw.bound(f"kill-population-plain:{number}"):
        system = StarSystem(system_config=cfg)
    return store.save_system(system, cfg, config=config)


def _stored_galaxy(config):
    """Two civilized systems and three plain ones placed in one sector (the `galaxy` fixture of test_population)."""
    capital_a, capital_b = _civilized_system(config), _civilized_system(config)
    plain = [_plain_system(config, number) for number in range(3)]
    conn = store.get_connection(config)
    try:
        with conn:
            sector_id = conn.execute("INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc, shell_index, band_index, wedge_index)"
                                     " VALUES ('Test Reach', 4000, 8000, 0, 0, 8000, 1, 0, 0)").lastrowid
            for system, offset in ((capital_a, -1.5), (capital_b, 1.5), (plain[0], -1.4), (plain[1], 1.4), (plain[2], 0.0)):
                _place(conn, sector_id, system, offset)
    finally:
        conn.close()


def _populate(config):
    with every_window_has_a_civilization(), inline_queue.inline_queue():
        run_population.run_population(_args(config, "population"))
        run_population.run_territories(_args(config, "territories"))


def _unpopulate(config):
    """The galaxy as it was before the passes: no species, polities or owners, and the scan's place forgotten. The draws of the passes read the stored ids, so the same galaxy is populated again, not a rebuilt one."""
    conn = store.get_connection(config)
    try:
        with conn:
            for table in ("system_owners", "polities", "species", "population_state"):
                conn.execute(f"DELETE FROM {table}")
    finally:
        conn.close()


def _digest(config):
    """The species and polities by their own columns and the system they are in, with no row ids, times or epochs."""
    conn = store.get_connection(config)
    try:
        conn.commit()
        found = {}
        species = sorted(repr(tuple(v for k, v in row.items())) for row in conn.execute(
            "SELECT s.name, y.cell_key, y.name AS system_name, s.life_class, s.life_stage, s.build, s.climate, s.size, s.civilization_age_years, s.era, s.spacefaring, s.energy_index,"
            " s.materials_index, s.information_index, s.medical_index, s.propulsion_index, s.defense_index, s.tech_level, s.tech_source"
            " FROM species s JOIN star_systems y ON y.id = s.star_system_id").fetchall())
        found["species"] = (len(species), hashlib.sha256("\n".join(species).encode()).hexdigest())
        polities = sorted(repr(tuple(v for k, v in row.items())) for row in conn.execute(
            "SELECT p.name, p.government, p.color, p.reach_ly, y.name AS capital FROM polities p JOIN star_systems y ON y.id = p.capital_system_id").fetchall())
        found["polities"] = (len(polities), hashlib.sha256("\n".join(polities).encode()).hexdigest())
        return found
    finally:
        conn.close()


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_a_population_made_with_killed_units_is_the_population_an_undisturbed_pass_makes(mysql_config, monkeypatch, seed):
    _stored_galaxy(mysql_config)
    counter = kill_harness.Plan()
    with monkeypatch.context() as patch:
        kill_harness.install(patch, counter)
        _populate(mysql_config)
    expected, commits = _digest(mysql_config), counter.seen
    assert commits > 0 and expected["species"][0] > 0, "the undisturbed pass made species to compare with"
    _unpopulate(mysql_config)
    assert _digest(mysql_config)["species"][0] == 0
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(commits), min(commits, rng.randint(1, 4)))}
    plan = kill_harness.Plan(kills)
    kill_harness.install(monkeypatch, plan)
    _populate(mysql_config)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (commits {commits}, plan {kills})"
    assert _digest(mysql_config) == expected, f"seed {seed}: kills {plan.fired} changed the population"
