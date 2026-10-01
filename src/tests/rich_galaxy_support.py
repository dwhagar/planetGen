# tests/rich_galaxy_support.py

"""
A small galaxy with something in nearly every table, built through the
real `generate.py` commands (TEST.11, TEST.14): a planned skeleton, a few
bright stars, a generated sector with the population pass, standalone
systems forced into each shape (close and wide binaries, comets, moons,
belts, life, a system file's slots), every phenomenon type (placed in that sector where the CLI
allows it), a polity holding a system, facilities, and two orbit updates.
"""

import os
import random
import sys

import pytest

import updateOrbits
from stellarObjects import _db
from tests.bughunt_support import mysql_argv, run_cli
from tests.conftest import _test_server_kwargs
from tests.db_schema_support import scratch_database
from tests.fuzz_support import deterministic_entropy

SYSTEM_SHAPES = (
    ["+planets", "+moons", "+asteroid_belt", "+comets", "+habitable_world", "+intelligent_life", "-binary_system"],
    ["+binary_system", "-wide_binary", "+planets", "+comets"],
    ["+binary_system", "+wide_binary", "+planets", "+moons"],
    ["+large_star", "+max_planets"],
)

SYSTEM_FILE = os.path.join(os.path.dirname(__file__), "..", "..", "examples", "systems", "arrakis_system.json")
"""A system file with orbit slots (`system_config_slots`)."""

PHENOMENA = ("black-hole", "neutron-star", "nebula", "supernova-remnant", "rogue-planet", "comet",
             "asteroid-field", "quasar")
PLACEABLE = ("nebula", "asteroid-field", "black-hole", "neutron-star")


@pytest.fixture(scope="module")
def rich_galaxy(_mysql_server_available):
    """A `MySQLConfig` for one seeded galaxy per test module; tests must
    leave it as they found it (roll back)."""
    with scratch_database(_db.MySQLConfig(database="", **_test_server_kwargs())) as config:
        build_rich_galaxy(config)
        yield config


def build_rich_galaxy(config, seed=20261001):
    """Builds the galaxy in `config`'s empty database (the same draws
    every time, `deterministic_entropy`) and returns the generated
    sector's id."""
    state = random.getstate()
    try:
        with deterministic_entropy(seed):
            return _build(config)
    finally:
        random.setstate(state)


def _build(config):
    target = mysql_argv(config)
    run_cli("plan", ["--quiet", "--no-bright-stars"] + target)
    run_cli("galaxy", ["--ring", "0", "--layer", "0", "--num-systems", "6", "+planets", "--population",
                       "--yes", "--quiet"] + target)
    conn = _db.get_connection(config)
    try:
        sector_id = conn.execute("SELECT MIN(id) AS id FROM sectors").fetchone()["id"]
    finally:
        conn.close()
    for shape in SYSTEM_SHAPES:
        run_cli("system", shape + ["--quiet"] + target)
    run_cli("system", ["--system-file", SYSTEM_FILE, "--quiet"] + target)
    for kind in PHENOMENA:
        extra = ["--sector-id", str(sector_id)] if kind in PLACEABLE else []
        run_cli("phenomenon", ["--type", kind, "--quiet"] + extra + target)
    for kind in ("black-hole", "neutron-star"):
        run_cli("phenomenon", ["--type", kind, "--anchor-system", "--quiet"] + target)
    _add_polity_and_facilities(config)
    for _ in range(2):  # the first only records a starting time
        _run_update_orbits(target)
    return sector_id


def _run_update_orbits(target):
    argv = sys.argv
    sys.argv = ["updateOrbits.py"] + target
    try:
        updateOrbits.main()
    finally:
        sys.argv = argv


def _add_polity_and_facilities(config):
    """A polity holding its species' home system (a real population pass
    makes civilizations too rarely to rely on), and one facility of each
    placement."""
    conn = _db.get_connection(config)
    try:
        def one(sql):
            return conn.execute(sql).fetchone()

        species = one("SELECT id, star_system_id FROM species ORDER BY id LIMIT 1")
        if species is None:
            planet = one("SELECT id, star_system_id FROM planets WHERE body_type = 't' ORDER BY id LIMIT 1")
            conn.execute(
                "INSERT INTO species (name, homeworld_planet_id, star_system_id, life_chemical, life_stage, build,"
                " climate, size, civilization_age_years, era, spacefaring) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                ("Velarix", planet["id"], planet["star_system_id"], "chlorophyll", "technological_civilization",
                 "medium", "temperate", "medium", 60000.0, "established", 1))
            species = one("SELECT id, star_system_id FROM species ORDER BY id LIMIT 1")
        conn.execute("INSERT INTO polities (name, species_id, capital_system_id, government, color, reach_ly)"
                     " VALUES (?, ?, ?, ?, ?, ?)",
                     ("Velarix Concord", species["id"], species["star_system_id"], "Concord", "#3366cc", 40.0))
        polity = one("SELECT id FROM polities ORDER BY id LIMIT 1")
        conn.execute("REPLACE INTO system_owners (star_system_id, polity_id, distance_ly) VALUES (?, ?, 0)",
                     (species["star_system_id"], polity["id"]))
        terrestrial = one("SELECT id FROM planets WHERE body_type = 't' ORDER BY id LIMIT 1")["id"]
        giant = one("SELECT id FROM planets WHERE body_type = 'g' ORDER BY id LIMIT 1")["id"]
        belt = one("SELECT id FROM asteroid_belts ORDER BY id LIMIT 1")["id"]
        _db.add_facility(conn, "New Hope", "colony", "terrestrial", "planet", terrestrial,
                         description="A colony.")
        _db.add_facility(conn, "High Yard", "starbase", "orbital", "planet", giant,
                         phase_deg=10.0)
        _db.add_facility(conn, "Rockpile", "mining-colony", "asteroid", "asteroid_belt", belt)
        conn.commit()
    finally:
        conn.close()
