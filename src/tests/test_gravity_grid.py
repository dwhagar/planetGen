# tests/test_gravity_grid.py

"""
MAP.387: the gravity map's sector grid and its cache (`db.gravity_grid`): the smooth galaxy's stellar density in closed form, the local field against a
brute-force sum, the cache hit, a moved star changing the stamp only of the grids that reach it, and a corrupt file. Needs MySQL.
"""

import math

import numpy as np
import pytest

from planetgen.db import gravity_grid, gravity_sources, store
from planetgen.galaxy import gravity_cache, gravity_kernel, potential, spherical
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.system import StarSystem
from planetgen.physics import constants

EDGE_LY = 11.5
EDGE_PC = EDGE_LY * constants.LIGHTYEAR_M / constants.PARSEC_M
POT = potential.Potential(spheres=(potential.Hernquist("bulge", 1.0e10, 700.0),), disks=(potential.Disk("thin_disk", 4.0e10, 2600.0, 300.0),))


def test_the_stellar_density_is_the_poisson_density_of_the_stellar_terms():
    rng = np.random.default_rng(2)
    points = rng.normal(size=(6, 3)) * [3000.0, 3000.0, 400.0] + [500.0, 200.0, 50.0]
    h = 1.0
    laplacian = -6.0 * potential.potential(POT, *points.T)
    for axis in range(3):
        step = np.zeros(3)
        step[axis] = h
        laplacian = laplacian + potential.potential(POT, *(points + step).T) + potential.potential(POT, *(points - step).T)
    np.testing.assert_allclose(potential.stellar_density(POT, *points.T), laplacian / h ** 2 / (4.0 * math.pi * potential.G_PC_KMS2_PER_MSUN), rtol=2e-5)
    assert potential.stellar_density(POT, 100.0, 0.0, 0.0) > 0.0
    assert potential.stellar_density(potential.Potential(), 100.0, 0.0, 0.0) == 0.0       # no stellar term, no density


def _save(mysql_config, address, systems=2, black_hole=False, name="Grid"):
    center_pc = spherical.sector_position_pc(*address, EDGE_PC)
    sector = SpaceSector(name, edge_ly=EDGE_LY)
    for index in range(systems):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = f"{name} star {index}"
        sector.add_system(system, position=(float(index) - 0.5, 0.0, 0.0), system_config=cfg)
    if black_hole:
        sector.add_phenomenon(BlackHole(SystemConfig()), "black-hole", position=(0.0, 2.0, 0.0))
    sector.place_in_galaxy(tuple(c * constants.PARSEC_M / constants.LIGHTYEAR_M for c in center_pc))
    return store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
        "galactic_radius_pc": math.sqrt(sum(c * c for c in center_pc)),
        "shell_index": address[0], "band_index": address[1], "wedge_index": address[2]})


@pytest.fixture
def galaxy(mysql_config):
    near = _save(mysql_config, (1500, 2, 700), 2, True, "Near")
    beside = _save(mysql_config, (1500, 2, 701), 1, False, "Beside")
    far = _save(mysql_config, (1530, 2, 700), 1, False, "Far")
    return mysql_config, near, beside, far


def _conn(config):
    return store.get_connection(config)


def test_the_grid_is_the_generated_stars_less_the_smooth_galaxys_own(galaxy, tmp_path):
    config, near, beside, far = galaxy
    conn = _conn(config)
    try:
        result = gravity_grid.grid(conn, near, str(tmp_path), POT)
        target = gravity_grid._sector(conn, near)
        sectors = gravity_grid.sectors_in_reach(conn, target)
        assert {row["id"] for row in sectors} == {near, beside}
        points = gravity_grid.zone_points((target["center_x_pc"], target["center_y_pc"], target["center_z_pc"]), target["edge_mpc"] / 1000.0)
        stars = gravity_sources.system_point_masses(conn, [near, beside])
        own_positions, own_masses = gravity_grid.own_stars(sectors, POT)
    finally:
        conn.close()
    assert result.hit is False and result.values.shape == (4096, 10) and result.values.dtype == np.float32
    expected = gravity_kernel.fields(points, stars.positions_pc, stars.masses_solar) - gravity_kernel.fields(points, own_positions, own_masses, gravity_grid.OWN_SOFTENING_PC)
    scale = np.abs(expected).max(axis=0)
    assert (np.abs(result.values - expected).max(axis=0) <= 5e-3 * scale + 1e-12).all()
    assert len(stars.keys) == 4                    # two systems, one more, and the black hole filed alone
    assert np.isfinite(result.values).all()


def test_a_second_call_reads_the_file_and_a_moved_star_is_worked_out_again(galaxy, tmp_path):
    config, near, beside, far = galaxy
    conn = _conn(config)
    try:
        first = gravity_grid.grid(conn, near, str(tmp_path), POT)
        again = gravity_grid.grid(conn, near, str(tmp_path), POT)
        assert again.hit is True and again.stamp == first.stamp and again.values.tobytes() == first.values.tobytes()
        conn.execute("UPDATE star_systems SET position_x_mpc = position_x_mpc + 500 WHERE sector_id = ? ORDER BY id LIMIT 1", (beside,))
        conn.commit()
        moved = gravity_grid.grid(conn, near, str(tmp_path), POT)
        assert moved.hit is False and moved.stamp != first.stamp and moved.values.tobytes() != first.values.tobytes()
        assert len([n for n in __import__("os").listdir(tmp_path) if n.startswith(f"g{near}_")]) == 1       # the old file went
    finally:
        conn.close()


def test_a_star_moving_in_a_sector_out_of_reach_changes_nothing(galaxy, tmp_path):
    config, near, beside, far = galaxy
    conn = _conn(config)
    try:
        first = gravity_grid.grid(conn, near, str(tmp_path), POT)
        conn.execute("UPDATE star_systems SET position_x_mpc = position_x_mpc + 500 WHERE sector_id = ?", (far,))
        conn.commit()
        again = gravity_grid.grid(conn, near, str(tmp_path), POT)
        assert again.hit is True and again.stamp == first.stamp
    finally:
        conn.close()


def test_the_stamp_follows_the_potential_and_the_settings(galaxy, tmp_path):
    config, near, beside, far = galaxy
    conn = _conn(config)
    try:
        sectors = gravity_grid.sectors_in_reach(conn, gravity_grid._sector(conn, near))
        base = gravity_grid.source_stamp(conn, sectors, pot=POT)
        assert base == gravity_grid.source_stamp(conn, sectors, pot=POT)
        assert base != gravity_grid.source_stamp(conn, sectors, pot=potential.Potential())
        assert base != gravity_grid.source_stamp(conn, sectors, zones=8, pot=POT)
        assert base != gravity_grid.source_stamp(conn, sectors, reach_edges=3, pot=POT)
    finally:
        conn.close()


def test_a_corrupt_file_is_worked_out_again_and_an_unknown_sector_is_an_error(galaxy, tmp_path):
    config, near, beside, far = galaxy
    conn = _conn(config)
    try:
        first = gravity_grid.grid(conn, near, str(tmp_path), POT)
        path = tmp_path / gravity_cache.file_name(near, 16, first.stamp)
        path.write_bytes(path.read_bytes()[:100])
        again = gravity_grid.grid(conn, near, str(tmp_path), POT)
        assert again.hit is False and again.values.tobytes() == first.values.tobytes()
        with pytest.raises(ValueError):
            gravity_grid.grid(conn, 10 ** 9, str(tmp_path), POT)
    finally:
        conn.close()
