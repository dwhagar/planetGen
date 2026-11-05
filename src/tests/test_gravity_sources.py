# tests/test_gravity_sources.py

"""
MAP.386: the gravity map's point masses (`db.gravity_sources`): a system is one point with the mass of its stars, planets and moons, a remnant filed on a star is
not counted twice, and folding a system's bodies into one point changes the field at the zone size by less than a tenth of a percent. Needs MySQL for the reader.
"""

import math

import numpy as np
import pytest

from planetgen.db import gravity_sources, store
from planetgen.galaxy import gravity_kernel, spherical
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.system import StarSystem
from planetgen.physics import constants

ADDRESS = (1500, 2, 700)
EDGE_LY = 11.5
EDGE_PC = EDGE_LY * constants.LIGHTYEAR_M / constants.PARSEC_M
AU_PC = constants.AU_M / constants.PARSEC_M
KG_PER_MSUN = 1.98847e30


def _save(mysql_config):
    center_pc = spherical.sector_position_pc(*ADDRESS, EDGE_PC)
    sector = SpaceSector("Gravity Sources", edge_ly=EDGE_LY)
    for index in range(2):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = f"Gravstar {index}"
        sector.add_system(system, position=(float(index) - 0.5, 0.0, 0.0), system_config=cfg)
    sector.add_phenomenon(BlackHole(SystemConfig()), "black-hole", position=(0.0, 2.0, 0.0))
    sector.place_in_galaxy(tuple(c * constants.PARSEC_M / constants.LIGHTYEAR_M for c in center_pc))
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
        "galactic_radius_pc": math.sqrt(sum(c * c for c in center_pc)),
        "shell_index": ADDRESS[0], "band_index": ADDRESS[1], "wedge_index": ADDRESS[2],
    })
    return sector_id


def test_a_system_is_one_point_with_the_mass_of_its_stars_planets_and_moons(mysql_config):
    sector_id = _save(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        got = gravity_sources.system_point_masses(conn, [sector_id])
        systems = conn.execute("SELECT id FROM star_systems ORDER BY id").fetchall()
        expected = {}
        for row in systems:
            total = 0.0
            for table in ("stars", "planets", "moons"):
                total += conn.execute(f"SELECT COALESCE(SUM(mass_kg), 0) AS m FROM {table} WHERE star_system_id = ?", (row["id"],)).fetchone()["m"]
            expected[("star_systems", row["id"])] = total / KG_PER_MSUN
        holes = conn.execute("SELECT id, mass_solar, star_id FROM black_holes").fetchall()
        planet_kg = conn.execute("SELECT COALESCE(SUM(mass_kg), 0) AS m FROM planets").fetchone()["m"]
    finally:
        conn.close()
    assert planet_kg > 0.0       # the test is about folding something in
    systems_got = {key: mass for key, mass in zip(got.keys, got.masses_solar) if key[0] == "star_systems"}
    assert set(systems_got) == set(expected) and len(expected) == 2
    for key, mass in expected.items():
        assert systems_got[key] == pytest.approx(mass, rel=1e-4)
    stand_alone = [key for key in got.keys if key[0] == "black_holes"]
    assert len(stand_alone) == sum(1 for row in holes if row["star_id"] is None)
    assert got.positions_pc.shape == (len(got.keys), 3) and list(got.keys) == sorted(got.keys)
    conn = store.get_connection(mysql_config)
    try:
        assert gravity_sources.system_point_masses(conn, [sector_id]).keys == got.keys
    finally:
        conn.close()


def test_no_sector_or_an_unknown_one_is_empty(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        for ids in ([], [10 ** 9]):
            got = gravity_sources.system_point_masses(conn, ids)
            assert got.positions_pc.shape == (0, 3) and got.masses_solar.shape == (0,) and got.keys == ()
    finally:
        conn.close()


def test_folding_a_systems_bodies_into_one_point_changes_the_field_at_the_zone_size_by_under_a_tenth_of_a_percent():
    """The rule itself: a sun, a giant at 5 AU, a few planets and a moon of 0.3 percent of the system's mass in all, seen from a zone 0.25 pc away."""
    rng = np.random.default_rng(4)
    star = np.zeros(3)
    planets = [(5.2, 9.5e-4), (1.0, 3.0e-6), (30.0, 5.0e-5), (0.4, 1.7e-7), (9.5, 2.9e-4), (19.0, 4.4e-5)]
    positions = [star]
    masses = [1.0]
    for au, mass in planets:
        direction = rng.normal(size=3)
        positions.append(direction / np.linalg.norm(direction) * au * AU_PC)
        masses.append(mass)
    positions, masses = np.array(positions), np.array(masses)
    barycentre = (positions * masses[:, None]).sum(0) / masses.sum()
    zones = rng.normal(size=(200, 3))
    zones = zones / np.linalg.norm(zones, axis=1)[:, None] * 0.25
    explicit = gravity_kernel.fields(zones, positions, masses)
    folded = gravity_kernel.fields(zones, barycentre[None, :], np.array([masses.sum()]))
    star_only = gravity_kernel.fields(zones, star[None, :], np.array([1.0]))
    folded_error = np.linalg.norm(folded[:, 0:3] - explicit[:, 0:3], axis=1) / gravity_kernel.pull_size(explicit)
    star_only_error = np.linalg.norm(star_only[:, 0:3] - explicit[:, 0:3], axis=1) / gravity_kernel.pull_size(explicit)
    assert folded_error.max() < 1.0e-3
    assert star_only_error.max() > 5.0e-4        # leaving the planets out is the larger error, so they are folded in
    at_the_star = gravity_kernel.fields(zones, star[None, :], np.array([masses.sum()]))
    assert (np.linalg.norm(at_the_star[:, 0:3] - explicit[:, 0:3], axis=1) / gravity_kernel.pull_size(explicit)).max() < 1.0e-3     # the mass at the star, not the barycentre, is as good
