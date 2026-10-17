# tests/test_gen337_point_fill.py

"""
GEN.337: the objects the free-object step (GEN.336) left out are written with the point columns of DB.55 when they are
made: a belt (position and vector at its star, a stated zero mass), the quasar (a stated zero vector and the orbit
epoch), a comet (the mass of its nucleus), a facility (the mass of what it sits on) and a scatter member (the epoch of
its tree).
"""

import pytest

from planetgen.db import scatter_systems, store
from planetgen.generation import multistar
from planetgen.generation.comet import Comet
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.system import StarSystem
from planetgen.physics import bulk_mass, constants
from planetgen.util import draw
from tests.test_facilities import _ids, _system_with_worlds
from tests.test_galactic_motion import SCATTER_RADIUS_PC, _insert_scatter_point, _seed_scatter_grid

pytestmark = pytest.mark.fresh_database

EPOCH = 1_700_000_000.0


def test_a_belt_is_at_its_star_with_no_vector_and_a_stated_zero_mass(mysql_config):
    system_id = _system_with_worlds(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        belt = conn.execute("SELECT * FROM asteroid_belts WHERE star_system_id = ?", (system_id,)).fetchone()
    finally:
        conn.close()
    for axis in "xyz":
        assert belt[f"position_{axis}_km"] == 0.0 and belt[f"velocity_{axis}_kms"] == 0.0
    assert belt["mass_kg"] == 0.0
    assert belt["epoch_unix"] is None and belt["next_update_due"] is None  # a belt follows its system's clock


def test_the_quasar_is_still_with_the_orbit_epoch(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("INSERT INTO orbit_simulation_state (id, last_updated_at) VALUES (1, FROM_UNIXTIME(?))"
                         " ON DUPLICATE KEY UPDATE last_updated_at = FROM_UNIXTIME(?)", (EPOCH, EPOCH))
            quasar_id = store.insert_quasar(conn, Quasar(SystemConfig()), placement={})
        row = conn.execute("SELECT * FROM quasars WHERE id = ?", (quasar_id,)).fetchone()
        epoch = store.get_orbit_epoch_unix(conn)
    finally:
        conn.close()
    assert (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]) == (0.0, 0.0, 0.0)
    assert epoch is not None and row["epoch_unix"] == epoch
    assert row["next_update_due"] is None


def test_a_comet_takes_the_mass_of_its_nucleus(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.COMETS = False
    system = StarSystem(system_config=cfg)
    comet = Comet(cfg, primary_mass_solar=system.star.mass / constants.SOLAR_MASS_TO_KG, orbit_type="elliptical")
    system.comets = [comet]
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        row = conn.execute("SELECT nucleus_diameter_km, mass_kg FROM comets WHERE star_system_id = ?",
                           (system_id,)).fetchone()
    finally:
        conn.close()
    assert row["mass_kg"] == pytest.approx(bulk_mass.comet_nucleus_mass_kg(row["nucleus_diameter_km"]), rel=1e-9)
    assert row["mass_kg"] > 0


def test_a_facility_takes_the_mass_of_what_it_sits_on(mysql_config):
    system_id = _system_with_worlds(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        ids = _ids(conn, system_id)
        with conn:
            colony = store.add_facility(conn, "New Hope", "colony", "terrestrial", "planet", ids["terrestrial"])
            outpost = store.add_facility(conn, "Sunwatch", "outpost", "orbital", "star", ids["star"],
                                         distance_km=constants.AU_TO_KM * 0.5)
            rockpile = store.add_facility(conn, "Rockpile", "mining-colony", "asteroid", "asteroid_belt", ids["belt"])
        mass = {name: conn.execute("SELECT mass_kg FROM facilities WHERE id = ?", (found,)).fetchone()["mass_kg"]
                for name, found in (("colony", colony), ("outpost", outpost), ("rockpile", rockpile))}
        planet = conn.execute("SELECT mass_kg FROM planets WHERE id = ?", (ids["terrestrial"],)).fetchone()["mass_kg"]
        star = conn.execute("SELECT mass_kg FROM stars WHERE id = ?", (ids["star"],)).fetchone()["mass_kg"]
    finally:
        conn.close()
    assert mass["colony"] == planet and mass["outpost"] == star
    assert mass["rockpile"] == 0.0


def test_a_scatter_member_holds_the_epoch_of_its_tree(mysql_config):
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 0.0, 0.0))
            tree = multistar.draw_member_tree(multistar.Placed("star", 1.0), 5e9, draw.Stream(3), n_stars=3)
            scatter_systems.store_member_tree(conn, "bright_stars", row_id, tree, EPOCH)
            rows = conn.execute("SELECT epoch_unix FROM scatter_members WHERE source = 1 AND source_id = ?",
                                (row_id,)).fetchall()
            batch_id = _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, 10.0, 0.0))
            scatter_systems.store_member_trees(conn, "bright_stars", [(batch_id, tree)], EPOCH + 5.0)
            batch = conn.execute("SELECT epoch_unix FROM scatter_members WHERE source = 1 AND source_id = ?",
                                 (batch_id,)).fetchall()
    finally:
        conn.close()
    assert len(rows) == 3 and {r["epoch_unix"] for r in rows} == {EPOCH}
    assert len(batch) == 3 and {r["epoch_unix"] for r in batch} == {EPOCH + 5.0}
