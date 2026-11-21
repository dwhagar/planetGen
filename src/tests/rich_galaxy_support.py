# tests/rich_galaxy_support.py

"""
A small galaxy with something in nearly every table, built through the
real `planetgen` commands (TEST.11, TEST.14): a planned skeleton, a few
bright stars, a generated sector with the population pass, standalone
systems forced into each shape (close and wide binaries, comets, moons,
belts, life, a system file's slots), every phenomenon type (placed in that sector where the CLI
allows it), a polity holding a system, facilities, and two orbit updates.
"""

import contextlib
import os
import sys

from planetgen.util import draw
from planetgen.cli import orbits
from planetgen.generation.comet import Comet, comet_designation
from planetgen.generation.config import SystemConfig
from planetgen.physics.constants import SOLAR_MASS_TO_KG
from planetgen.db import store
from planetgen.population import model as population
from tests import rich_galaxy_cache
from tests.bughunt_support import forced_system_config, mysql_argv, run_cli
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


@contextlib.contextmanager
def scratch_rich_galaxy():
    """A `MySQLConfig` for one seeded galaxy in a scratch database, dropped on exit. The `rich_galaxy` fixture (conftest.py,
    one per test process, TEST.417: the build is four minutes and the modules that use it only read it) holds it; tests
    must leave it as they found it (roll back)."""
    with scratch_database(store.MySQLConfig(database="", **_test_server_kwargs())) as config:
        folder = rich_galaxy_cache.cache_folder()
        if folder is None:
            build_rich_galaxy(config)
        else:
            # TEST.788: with a shared folder, one process builds the galaxy of a code state and the rest load its dump (`tests/rich_galaxy_cache.py`)
            here = os.path.dirname(os.path.abspath(__file__))
            key = rich_galaxy_cache.key_of(config, os.path.join(here, "..", "planetgen"), here, os.environ.get("PLANETGEN_DATA_DIR"),
                                           extra_files=(SYSTEM_FILE, os.path.join(here, "..", "..", "config.json")))
            if rich_galaxy_cache.build_or_load(config, lambda: build_rich_galaxy(config), folder, key) == "loaded":
                _make_the_cluster_field(config)
        yield config


def _make_the_cluster_field(config):
    """
    TEST.807: the galaxy a process built has the galaxy's cluster field in its memory (and a file of its state folder), which the dump does not carry; a visitor's request for a
    page that needs it is refused (503, `CacheRebuilding`) until the warm-up has made it. The process that loads a dump makes it here, as the warm-up does, outside a request.
    """
    conn = store.get_connection(config)
    try:
        store.get_galaxy_shape(conn)
    finally:
        conn.close()


def build_rich_galaxy(config, seed=20261001):
    """Builds the galaxy in `config`'s empty database (the same draws
    every time, `deterministic_entropy`) and returns the generated
    sector's id."""
    state = draw.getstate()
    try:
        with deterministic_entropy(seed):
            return _build(config)
    finally:
        draw.setstate(state)


def _build(config):
    target = mysql_argv(config)
    run_cli("plan", ["--quiet", "--cell-star-ceiling", "1e12", "--nuclear-amplitude", "0"] + target)
    with forced_system_config(PLANETS=True):
        run_cli("sectors", ["--shell", "1", "--band", "0", "--limit", "3", "--num-systems", "6",
                            "--yes", "--quiet"] + target)
        run_cli("link", ["--quiet"] + target)
        run_cli("backfill", ["--quiet"] + target)
        run_cli("population", ["--quiet"] + target)
        run_cli("settle", ["--quiet"] + target)
    conn = store.get_connection(config)
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
    _add_parabolic_comet(config)
    _add_polity_and_facilities(config)
    for _ in range(2):  # the first only records a starting time
        _run_update_orbits(target)
    return sector_id


def _run_update_orbits(target):
    argv = sys.argv
    sys.argv = ["planetgen.cli.orbits"] + target
    try:
        orbits.main()
    finally:
        sys.argv = argv


def _add_parabolic_comet(config):
    """A parabolic comet on the first star system's first star, when the galaxy's comets hold none: a star holds comets only where its outer disk was scattered (GEN.1775) and 3 in 10
    comets are parabolic, so the few comets a seeded galaxy makes may all be elliptical, and `comets.parabolic_mean_anomaly` would be NULL in every row."""
    conn = store.get_connection(config)
    try:
        if conn.execute("SELECT 1 FROM comets WHERE parabolic_mean_anomaly IS NOT NULL LIMIT 1").fetchone() is not None:
            return
        host = conn.execute("SELECT s.id AS star_id, s.star_system_id, s.mass_kg, ss.name AS system_name FROM stars s JOIN star_systems ss ON ss.id = s.star_system_id"
                            " ORDER BY s.id LIMIT 1").fetchone()
        comet = Comet(SystemConfig(), host["mass_kg"] / SOLAR_MASS_TO_KG, orbit_type="parabolic")
        count = conn.execute("SELECT COUNT(*) AS n FROM comets WHERE star_system_id = ?", (host["star_system_id"],)).fetchone()["n"]
        comet.name = comet_designation(host["system_name"], count + 1, comet)
        store.insert_comet(conn, comet, host["star_system_id"], host["star_id"])
        conn.commit()
    finally:
        conn.close()


def _add_polity_and_facilities(config):
    """A polity holding its species' home system (a real population pass
    makes civilizations too rarely to rely on), and one facility of each
    placement."""
    conn = store.get_connection(config)
    try:
        def one(sql):
            return conn.execute(sql).fetchone()

        species = one("SELECT id, star_system_id FROM species ORDER BY id LIMIT 1")
        if species is None:
            planet = one("SELECT id, star_system_id FROM planets WHERE body_type = 't' ORDER BY id LIMIT 1")
            conn.execute(
                "INSERT INTO species (name, homeworld_planet_id, star_system_id, life_class, life_stage, build,"
                " climate, size, civilization_age_years, era, spacefaring) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                ("Velarix", planet["id"], planet["star_system_id"], "1a", "technological_civilization",
                 "medium", "temperate", "medium", 60000.0, "established", 1))
            species = one("SELECT id, star_system_id FROM species ORDER BY id LIMIT 1")
        population.refresh_tech_levels(conn)
        conn.execute("INSERT INTO polities (name, species_id, capital_system_id, government, color, reach_ly)"
                     " VALUES (?, ?, ?, ?, ?, ?)",
                     ("Velarix Concord", species["id"], species["star_system_id"], "Concord", "#3366cc", 40.0))
        polity = one("SELECT id FROM polities ORDER BY id LIMIT 1")
        conn.execute("REPLACE INTO system_owners (star_system_id, polity_id, distance_ly) VALUES (?, ?, 0)",
                     (species["star_system_id"], polity["id"]))
        terrestrial = one("SELECT id FROM planets WHERE body_type = 't' ORDER BY id LIMIT 1")["id"]
        giant = one("SELECT id FROM planets WHERE body_type = 'g' ORDER BY id LIMIT 1")["id"]
        belt = one("SELECT id FROM asteroid_belts ORDER BY id LIMIT 1")["id"]
        store.add_facility(conn, "New Hope", "colony", "terrestrial", "planet", terrestrial,
                         description="A colony.")
        store.add_facility(conn, "High Yard", "starbase", "orbital", "planet", giant,
                         phase_deg=10.0)
        store.add_facility(conn, "Rockpile", "mining-colony", "asteroid", "asteroid_belt", belt)
        conn.commit()
    finally:
        conn.close()
