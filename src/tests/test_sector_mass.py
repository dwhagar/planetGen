"""
DB.108: the mass of each generated sector (`sector_stats.total_mass_kg`), of each block above it (`galaxy_block_stats.mass_sum_kg`)
and of each region (`region_stats.mass_sum_kg`). A sector's mass is the sum of the masses of the bodies saved under it, the saves
and deletes that `block_stats` hooks keep it, a block's and a region's is the sum of its sectors', a backfill fills a galaxy
saved before the column, a sector holding a body with no stored mass is left NULL and listed, and the order the sectors are
saved in does not change a bit of it.
"""

import math
import random

import pytest

from planetgen.db import block_stats, edits, region_stats, sector_mass, store
from planetgen.galaxy import spherical
from planetgen.galaxy.drill import STAGE_TOP, drill_chain_of
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen import tuning
from planetgen.physics import constants

pytestmark = pytest.mark.db

BODY_TABLES = ("stars", "planets", "moons", "asteroid_belts", "comets")


def _neighbours(start, size):
    """`size` sectors grown by face-neighbours from `start`, in a fixed order."""
    found = [start]
    frontier = [start]
    while frontier and len(found) < size:
        for neighbour in sorted(spherical.neighbor_addresses(*frontier.pop(0))):
            if neighbour not in found and len(found) < size:
                found.append(neighbour)
                frontier.append(neighbour)
    return found


def _addresses():
    """Three face-neighbours (one region) and a lone sector far away (another)."""
    low, high = spherical.band_range(40)
    near = _neighbours((40, (low + high) // 2, 0), 3)
    low, high = spherical.band_range(300)
    return near + [(300, (low + high) // 2, 5)]


def _make(address, systems, rogue=1, remnant=False, comet=False):
    """A sector at `address` holding `systems` star systems with their planets, `rogue` rogue planets and optionally a black hole and a neutron star."""
    random.seed(sum(address) + systems)
    sector = SpaceSector(f"Mass {address}", edge_ly=20.0)
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    for n in range(systems):
        sector.add_system(StarSystem(system_config=cfg), position=(1.0 + 3.0 * n, 1.0, 1.0), system_config=cfg)
    for _ in range(rogue):
        sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    if comet:
        sector.add_phenomenon(InterstellarComet(SystemConfig()), "comet")
    if remnant:
        sector.add_phenomenon(BlackHole(SystemConfig()), "black-hole")
        sector.add_phenomenon(NeutronStar(SystemConfig()), "neutron-star")
    return sector


def _save(config, address, sector):
    centre = spherical.sector_position_pc(*address, float(tuning.DEFAULT_SECTOR_EDGE_PC))
    return store.save_sector(sector, config=config, galaxy_position={
        "center_x_pc": float(centre[0]), "center_y_pc": float(centre[1]), "center_z_pc": float(centre[2]),
        "galactic_radius_pc": math.dist(centre, (0.0, 0.0, 0.0)),
        "shell_index": address[0], "band_index": address[1], "wedge_index": address[2]})


def _independent_mass(conn, sector_id):
    """The sector's mass read body by body from the tables, the long way."""
    masses = []
    for table in BODY_TABLES:
        masses += [r["m"] for r in conn.execute(
            f"SELECT b.mass_kg AS m FROM {table} b JOIN star_systems ss ON ss.id = b.star_system_id WHERE ss.sector_id = ?",
            (sector_id,)).fetchall()]
    for table in ("rogue_planets", "interstellar_comets", "asteroid_fields"):
        masses += [r["m"] for r in conn.execute(f"SELECT mass_kg AS m FROM {table} WHERE sector_id = ?", (sector_id,)).fetchall()]
    for table in ("black_holes", "neutron_stars"):
        masses += [r["m"] * constants.SOLAR_MASS_TO_KG for r in conn.execute(
            f"SELECT mass_solar AS m FROM {table} WHERE sector_id = ?", (sector_id,)).fetchall()]
    return math.fsum(masses)


def _stat(conn, address):
    return conn.execute("SELECT total_mass_kg FROM sector_stats WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                        address).fetchone()["total_mass_kg"]


def _id(conn, address):
    return conn.execute("SELECT id FROM sectors WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", address).fetchone()["id"]


def _rebuild(config):
    """The block and region rows of what is saved, as a galaxy run makes them when it ends."""
    conn = store.get_connection(config)
    try:
        block_stats.rebuild(conn)
        conn.commit()
    finally:
        conn.close()


@pytest.fixture
def galaxy(mysql_config):
    """Four saved sectors (three in one region), one with a black hole and a neutron star, one with a nebula, one empty of systems."""
    addresses = _addresses()
    _save(mysql_config, addresses[0], _make(addresses[0], 2, rogue=2, remnant=True))
    nebula_sector = _make(addresses[1], 1, rogue=0)
    nebula_sector.add_phenomenon(Nebula(SystemConfig(), name="Mass Veil"), "nebula")
    _save(mysql_config, addresses[1], nebula_sector)
    _save(mysql_config, addresses[2], _make(addresses[2], 0, rogue=0))
    _save(mysql_config, addresses[3], _make(addresses[3], 1, comet=True))
    _rebuild(mysql_config)
    return mysql_config


def test_a_sectors_mass_is_the_sum_of_its_bodies(galaxy):
    conn = store.get_connection(galaxy)
    try:
        for address in _addresses():
            sector_id = _id(conn, address)
            expected = _independent_mass(conn, sector_id)
            assert _stat(conn, address) == pytest.approx(expected, rel=1e-12)
            assert sector_mass.sector_mass_kg(conn, sector_id) == _stat(conn, address)
        assert _stat(conn, _addresses()[0]) > 0.0
    finally:
        conn.close()


def test_an_empty_sector_holds_zero_and_a_nebula_adds_nothing(galaxy):
    conn = store.get_connection(galaxy)
    try:
        empty, nebulous = _addresses()[2], _addresses()[1]
        assert _stat(conn, empty) == 0.0
        assert conn.execute("SELECT COUNT(*) AS n FROM nebulae WHERE sector_id = ?", (_id(conn, nebulous),)).fetchone()["n"] == 1
        assert _stat(conn, nebulous) == pytest.approx(_independent_mass(conn, _id(conn, nebulous)), rel=1e-12)
    finally:
        conn.close()


def test_blocks_and_regions_sum_their_sectors(galaxy):
    conn = store.get_connection(galaxy)
    try:
        masses = {address: _stat(conn, address) for address in _addresses()}
        chains = {address: {tuple(b) for b in drill_chain_of(*address)[:-1] if b.m <= STAGE_TOP} for address in masses}
        for block in set().union(*chains.values()):
            members = [m for address, m in masses.items() if block in chains[address]]
            row = conn.execute("SELECT mass_sum_kg FROM galaxy_block_stats WHERE m = ? AND shell_index = ? AND band_index = ?"
                               " AND wedge_index = ?", block).fetchone()
            assert row["mass_sum_kg"] == pytest.approx(math.fsum(members), rel=1e-12)
        cells = {}
        for address, mass in masses.items():
            cells.setdefault(tuple(region_stats.region_cell(address))[1:], []).append(mass)
        stored = {tuple(r["block"])[1:]: r for r in region_stats.regions(conn)}
        assert set(stored) == set(cells)
        for cell, members in cells.items():
            assert stored[cell]["mass_sum_kg"] == math.fsum(members)
    finally:
        conn.close()


def test_saving_and_deleting_a_sector_update_the_sums(galaxy):
    conn = store.get_connection(galaxy)
    try:
        near = _addresses()[:3]
        region = {tuple(r["block"])[1:]: r for r in region_stats.regions(conn)}[tuple(region_stats.region_cell(near[0]))[1:]]
        victim = near[0]
        lost = _stat(conn, victim)
        assert lost > 0.0
        edits.delete_sector_with_contents(conn, _id(conn, victim))
        conn.commit()
        assert _stat(conn, victim) is None
        after = {tuple(r["block"])[1:]: r for r in region_stats.regions(conn)}[tuple(region_stats.region_cell(near[0]))[1:]]
        assert after["mass_sum_kg"] == pytest.approx(region["mass_sum_kg"] - lost, rel=1e-9)
    finally:
        conn.close()
    again = _make(near[0], 1, rogue=0)
    _save(galaxy, near[0], again)
    conn = store.get_connection(galaxy)
    try:
        block_stats.refresh_chain(conn, near[0])     # what the end of a galaxy run does for the sectors it saved
        conn.commit()
        assert _stat(conn, victim) == pytest.approx(_independent_mass(conn, _id(conn, victim)), rel=1e-12)
        assert _stat(conn, victim) > 0.0
        region_row = {tuple(r["block"])[1:]: r for r in region_stats.regions(conn)}[tuple(region_stats.region_cell(near[0]))[1:]]
        assert region_row["mass_sum_kg"] == pytest.approx(math.fsum(_stat(conn, a) for a in near), rel=1e-12)
    finally:
        conn.close()


def test_a_backfill_fills_a_galaxy_saved_before_the_column(galaxy):
    conn = store.get_connection(galaxy)
    try:
        kept = {address: _stat(conn, address) for address in _addresses()}
        conn.execute("UPDATE sector_stats SET total_mass_kg = NULL")
        conn.execute("UPDATE galaxy_block_stats SET mass_sum_kg = 0")
        conn.execute("UPDATE region_stats SET mass_sum_kg = 0")
        skipped = []
        block_stats.rebuild(conn, on_skipped=skipped.extend)
        conn.commit()
        assert skipped == []
        assert {address: _stat(conn, address) for address in _addresses()} == kept
        assert all(r["mass_sum_kg"] > 0.0 for r in region_stats.regions(conn))
    finally:
        conn.close()


def test_a_sector_with_a_body_of_no_stored_mass_is_listed_and_left_null(galaxy):
    conn = store.get_connection(galaxy)
    try:
        address = _addresses()[3]
        conn.execute("UPDATE interstellar_comets SET mass_kg = NULL WHERE sector_id = ?", (_id(conn, address),))
        conn.commit()
        skipped = []
        block_stats.rebuild(conn, on_skipped=skipped.extend)
        conn.commit()
        assert skipped == [address]
        assert _stat(conn, address) is None
        others = [a for a in _addresses() if a != address]
        assert all(_stat(conn, a) is not None for a in others)
        far = {tuple(r["block"])[1:]: r for r in region_stats.regions(conn)}[tuple(region_stats.region_cell(address))[1:]]
        assert far["mass_sum_kg"] == 0.0              # the lone sector of its region counts as nothing
    finally:
        conn.close()


def test_the_order_the_sectors_are_saved_in_does_not_change_a_bit(mysql_config):
    addresses = _addresses()
    built = {a: _make(a, 1 + n % 2, rogue=1, remnant=n == 0) for n, a in enumerate(addresses)}
    for a in addresses:
        _save(mysql_config, a, built[a])
    conn = store.get_connection(mysql_config)
    try:
        block_stats.refresh_sectors(conn, addresses)
        conn.commit()
        forward = ({a: _stat(conn, a) for a in addresses},
                   sorted((tuple(r["block"]), r["mass_sum_kg"]) for r in region_stats.regions(conn)))
        block_stats.rebuild(conn)
        conn.commit()
        rebuilt = ({a: _stat(conn, a) for a in addresses},
                   sorted((tuple(r["block"]), r["mass_sum_kg"]) for r in region_stats.regions(conn)))
        assert rebuilt == forward
    finally:
        conn.close()
