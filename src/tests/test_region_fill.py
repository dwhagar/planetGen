"""
GEN.826: a normal fill charts one whole region. A random start is the anchor sector of a region nothing is generated in yet,
`--region` snaps a single-address or local-neighborhood run to the anchor of the region holding its center and fills the standard
sphere about it, and the Generate page's "This region" box asks for that. The snap changes which sectors a fill takes, never what
a sector generates.
"""

import argparse
import math
import random
import types

import pytest

from planetgen import tuning
from planetgen.db import region_stats, store
from planetgen.galaxy import named_regions, spherical
from planetgen.generation import run_galaxy
from planetgen.web import generate_page

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)


def _addresses(count=60, seed=826):
    """Sectors from the core out to far past the galaxy's edge, in a fixed order."""
    rng = random.Random(seed)
    found = []
    for shell in [0, 1, 2, 5, 9, 10, 40, 60, 200, 1000, 5000, 9000]:
        low, high = spherical.band_range(shell)
        for _ in range(count // 12):
            band = rng.randint(low, high)
            found.append((shell, band, rng.randrange(spherical.wedge_count(shell, band))))
    return found


def test_the_anchor_sector_is_in_its_region_and_the_snap_is_stable():
    for address in _addresses():
        region = named_regions.region_of(address)
        anchor = named_regions.snap_to_anchor(address, EDGE_PC)
        assert named_regions.region_of(anchor) == region
        assert named_regions.snap_to_anchor(anchor, EDGE_PC) == anchor
        assert named_regions.anchor_address(region, EDGE_PC) == anchor
        assert math.dist(spherical.sector_position_pc(*anchor, EDGE_PC), named_regions.anchor_pc(region, EDGE_PC)) < EDGE_PC


@pytest.fixture
def start_world(monkeypatch):
    """A galaxy reaching out to shell 100 holding nothing, the picks drawn as `_pick_random_start` draws them."""
    bounds = types.SimpleNamespace(
        contains=lambda shell, band: shell <= REACH_SHELL,
        random_address=lambda draw, max_shell=None: _random_address(draw, REACH_SHELL if max_shell is None else max_shell))
    batch = types.SimpleNamespace(
        bounds=bounds, skeleton=types.SimpleNamespace(expected_system_count_at_density_1=1.0),
        resolve=lambda args, address, position: types.SimpleNamespace(density=2.0, num_systems=None))
    charted = set()
    monkeypatch.setattr(run_galaxy.store, "get_sector_id_at", lambda conn, *address: None)
    monkeypatch.setattr(run_galaxy.region_stats, "region_row", lambda conn, cell: {"block": cell} if tuple(cell) in charted else None)
    args = types.SimpleNamespace(neighborhoods=1, neighborhood_gamma=0.0, max_shell=None, avoid_filled_space=False,
                                 min_start_density=None, density=None, num_systems=None)
    return batch, args, charted


REACH_SHELL = 100
"""int: The outermost shell of the fake galaxy the start tests draw in."""


def _random_address(draw, max_shell):
    shell = int(draw.random() * (max_shell + 1))
    low, high = spherical.band_range(shell)
    band = low + int(draw.random() * (high - low + 1))
    return shell, band, int(draw.random() * spherical.wedge_count(shell, band))


def test_a_random_start_is_the_anchor_of_a_region_nothing_is_generated_in(start_world):
    batch, args, charted = start_world
    for _ in range(200):
        address, _position, _args = run_galaxy._pick_random_start(args, EDGE_PC, batch, None, 12.0, [])
        anchor = named_regions.snap_to_anchor(address, EDGE_PC)
        assert address == anchor or anchor[0] > REACH_SHELL    # past the reach the drawn address starts the fill


def test_a_random_start_never_lands_in_a_charted_region(start_world):
    batch, args, charted = start_world
    seen = {tuple(named_regions.region_of(run_galaxy._pick_random_start(args, EDGE_PC, batch, None, 12.0, [])[0]))
            for _ in range(60)}
    assert len(seen) > 5
    charted.update(seen)
    for _ in range(200):
        address = run_galaxy._pick_random_start(args, EDGE_PC, batch, None, 12.0, [])[0]
        assert tuple(named_regions.region_of(address)) not in charted


def test_a_start_whose_anchor_is_beyond_max_shell_is_the_drawn_address(start_world):
    batch, args, _charted = start_world
    args.max_shell = 1                  # the anchor of the core cell is in shell 4
    for _ in range(100):
        address = run_galaxy._pick_random_start(args, EDGE_PC, batch, None, 12.0, [])[0]
        assert address[0] <= 1


def _namespace(**fields):
    base = dict(region=True, center_sector=None, shell=None, band=None, wedge=None, radius_pc=None, mysql_host=None)
    base.update(fields)
    return argparse.Namespace(**base)


def test_region_snaps_an_address_to_the_anchor_and_the_standard_radius():
    address = (40, 17, 3)
    args = _namespace(shell=40, band=17, wedge=3)
    run_galaxy._snap_to_region(args, EDGE_PC)
    assert (args.shell, args.band, args.wedge) == named_regions.snap_to_anchor(address, EDGE_PC)
    assert args.radius_pc == named_regions.REGION_RADIUS_PC == tuning.DEFAULT_GENERATE_RADIUS_PC == 12.0
    assert args.center_sector is None


def test_region_leaves_a_run_without_it_alone():
    args = _namespace(region=False, shell=40, band=17, wedge=3, radius_pc=5.0)
    run_galaxy._snap_to_region(args, EDGE_PC)
    assert (args.shell, args.band, args.wedge, args.radius_pc) == (40, 17, 3, 5.0)


def test_region_snaps_a_stored_sector_to_the_anchor_of_its_region(mysql_config):
    address = (40, 17, 3)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES ('Snap', 4000, ?, ?, ?)",
                     address)
        sector_id = conn.execute("SELECT id FROM sectors WHERE name = 'Snap'").fetchone()["id"]
        conn.commit()
    finally:
        conn.close()
    args = _namespace(center_sector=sector_id, mysql_host=mysql_config.host, mysql_port=mysql_config.port,
                      mysql_user=mysql_config.user, mysql_password=mysql_config.password, mysql_database=mysql_config.database)
    run_galaxy._snap_to_region(args, EDGE_PC)
    assert args.center_sector is None
    assert (args.shell, args.band, args.wedge) == named_regions.snap_to_anchor(address, EDGE_PC)
    assert args.radius_pc == named_regions.REGION_RADIUS_PC


def test_the_generate_page_asks_for_the_region_with_the_box_and_drops_the_radius():
    place = {"what": "place", "place": "address", "center_shell": "40", "center_band": "17", "center_wedge": "3"}
    argv, description = generate_page.center_argv({**place, "center_region": "1", "center_radius_pc": "5"})
    assert argv == ["--shell", "40", "--band", "17", "--wedge", "3", "--region"]
    assert description == "in the region of shell 40 band 17 wedge 3"
    argv, _description = generate_page.center_argv({**place, "center_radius_pc": "5"})
    assert argv == ["--shell", "40", "--band", "17", "--wedge", "3", "--radius-pc", "5.0"]
    sector = {"what": "place", "place": "sector", "center_sector": "9", "center_region": "1"}
    assert generate_page.center_argv(sector)[0] == ["--center-sector", "9", "--region"]


@pytest.mark.db
def test_a_finished_fill_of_one_anchor_sphere_is_charted_and_named(mysql_config, monkeypatch):
    """`--region` on an address fills the sphere about the region's anchor; once the block totals are made (a galaxy run's end) the
    region is stored with that fill's sectors and has its name."""
    from planetgen.db import block_stats
    from planetgen.names import naming_key
    from tests.test_galaxy_gen import (
        _address, _all_sectors, _fake_args_for_direct_generation, _mysql_argv, _run_cli, _seed_skeleton, _shells)

    monkeypatch.setattr(named_regions, "REGION_RADIUS_PC", 6.0)       # a smaller sphere than 12 pc keeps the fill short
    _seed_skeleton(mysql_config, shells=_shells(60, 60))
    address = (40, 17, 3)
    region = named_regions.region_of(address)
    anchor = named_regions.snap_to_anchor(address, EDGE_PC)
    run_galaxy.generate_and_save_sector_at(
        _fake_args_for_direct_generation(mysql_config), address, spherical.sector_position_pc(*address, EDGE_PC), EDGE_PC)
    sector_id = _all_sectors(mysql_config)[0]["id"]
    _run_cli(["--center-sector", str(sector_id), "--region", "--num-systems", "1"] + _mysql_argv(mysql_config), planets=False)
    made = {_address(row) for row in _all_sectors(mysql_config)}
    assert anchor in made and address in made and len(made) > 10
    assert all(named_regions.region_of(found) == region for found in made)
    conn = store.get_connection(mysql_config)
    try:
        assert region_stats.region_row(conn, region) is None          # the run's end makes the totals
        block_stats.rebuild(conn)
        conn.commit()
        row = region_stats.region_row(conn, region)
        assert row is not None and row["sector_count"] == len(made)
    finally:
        conn.close()
    assert named_regions.region_name(region, "DEADBEEF").endswith("Region")
