# tests/test_order_independence.py

"""
TEST.181: a sector's contents do not depend on the order sectors were generated in or on the worker count (GEN.57,
GEN.306, GEN.307). The same small set of sectors is generated in two orders and at one worker and four, and the
fingerprints of the set (GEN.58 over the canonical encoder of GEN.304) are compared:

- sectors that share an ID cell get their collision numbers by address, not by which saved first (GEN.307);
- a sector reached by the bright-star backfill holds the same stars whichever sector the run reached first;
- the same sectors at one worker and four are the same sectors.

Nearest-system links and the population data are rebuilt from content, so they are not in the comparison
(`db.fingerprint` leaves them out). Name decorations are compared as they are.
"""

import pytest

from planetgen.db import fingerprint, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.names import object_id as objectId
from planetgen.util import draw

from tests.test_galaxy_gen import _plan_wide_galaxy
from tests.test_object_id import _POSITION
from tests.test_parallel_generation import GALAXY_SEED, _run, make_database  # noqa: F401 -- a fixture

pytestmark = pytest.mark.db

SECTORS = ((1, 0, 0), (1, 0, 5), (1, 0, 11), (1, 1, 0))
"""tuple: `(shell, band, wedge)` of the set: three of the equator's twelve cells of shell 1 and the cell above them."""


def _fingerprint(config, addresses=SECTORS):
    """`(fingerprint, sections)` of the region: the sectors' digests and, per sector, the digests of its named sections
    (`fingerprint.sector_sections`), which say what differs when two do."""
    conn = store.get_connection(config)
    try:
        found = fingerprint.region_fingerprint(conn, addresses=addresses)
        sections = {f"{a[0]} {a[1]} {a[2]}": fingerprint.sector_sections(conn, address=a) for a in addresses}
        return found, sections
    finally:
        conn.close()


def _assert_same(first, second):
    """The two fingerprints agree, and a failure names the sectors, the sections and whether the shape or a number
    differs."""
    (first, first_sections), (second, second_sections) = first, second
    left, right = dict(first.sectors), dict(second.sectors)
    assert left.keys() == right.keys(), f"sectors only in one: {sorted(left.keys() ^ right.keys())}"
    differing = {}
    for label in left:
        if left[label] != right[label]:
            a, b = first_sections.get(label, {}), second_sections.get(label, {})
            differing[label] = {name: ("structure" if a[name].structural != b[name].structural else "numbers")
                                for name in sorted(a.keys() & b.keys()) if a[name].full != b[name].full}
            differing[label].update({name: "only in one" for name in sorted(a.keys() ^ b.keys())})
    assert differing == {}
    assert first.region == second.region


def _galaxy(make_database):
    config = make_database()
    _plan_wide_galaxy(config, galaxy_seed=GALAXY_SEED)
    return config


def _generate_in_order(config, order, monkeypatch, backfill="none"):
    """One `galaxy` run a sector, in `order` (a span of that one cell, 4 systems), then each run's backfill."""
    for shell, band, wedge in order:
        _run("galaxy", ["--shells", str(shell), "--bands", str(band), "--wedges", str(wedge), "--num-systems", "4",
                        "--workers", "1", "--backfill-from", backfill], config, monkeypatch)


def test_the_same_sectors_come_out_in_either_order(make_database, monkeypatch):
    results = []
    for order in (SECTORS, SECTORS[::-1]):
        config = _galaxy(make_database)
        _generate_in_order(config, order, monkeypatch)
        results.append(_fingerprint(config))
    assert len(results[0][0].sectors) == len(SECTORS)
    _assert_same(*results)


def test_a_sector_reached_by_the_backfill_comes_out_the_same_in_either_order(make_database, monkeypatch):
    """The first run generates one sector and backfills the bright stars around it (`--backfill-from edge`), so the
    sectors after it are reached by that backfill; they then come out the same in whichever order they are
    generated. (A sector generated before the backfill reaches its cell does not get those stars, so it differs from
    one the backfill reached first: that dependence is open, GEN.57's backfill clause.)"""
    first, rest = SECTORS[:1], SECTORS[1:]
    results = []
    for order in (rest, rest[::-1]):
        config = _galaxy(make_database)
        _generate_in_order(config, first, monkeypatch, backfill="edge")
        _generate_in_order(config, order, monkeypatch)
        results.append(_fingerprint(config))
    assert len(results[0][0].sectors) == len(SECTORS)
    _assert_same(*results)


@pytest.mark.parametrize("workers", [4])
def test_the_same_sectors_come_out_at_one_worker_and_several(make_database, monkeypatch, workers):
    results = []
    for count in (1, workers):
        config = _galaxy(make_database)
        _run("galaxy", ["--shells", "1", "--num-systems", "4", "--workers", str(count), "--backfill-from", "none"], config, monkeypatch)
        results.append(_fingerprint(config))
    assert len(results[0][0].sectors) == len(SECTORS)
    _assert_same(*results)


def _save_rogue_planets(config, order):
    """One rogue planet at one point in each of three sectors of shell 1, saved in `order`."""
    for wedge in order:
        sector = SpaceSector(f"Sector {wedge}", edge_ly=11.5)
        with draw.bound(5000 + wedge):   # each sector draws from a seed of its own, as a generated one does
            sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet", position=(1.0, 2.5, 1.0))
        position = {**_POSITION, "shell_index": 1, "band_index": 0, "wedge_index": wedge}
        store.save_sector(sector, config=config, galaxy_position=position)


def test_sectors_sharing_an_id_cell_come_out_the_same_in_either_order(make_database):
    """Three sectors whose rogue planets share one ID: the collision numbers follow the addresses (GEN.307), and the
    whole stored sector, not only the number, is the same whichever saved first."""
    results, collisions = [], []
    for order in ((0, 1, 2), (2, 1, 0)):
        config = make_database()
        _save_rogue_planets(config, order)
        addresses = [(1, 0, wedge) for wedge in (0, 1, 2)]
        results.append(_fingerprint(config, addresses))
        conn = store.get_connection(config)
        try:
            rows = conn.execute("SELECT s.wedge_index, r.name FROM rogue_planets r JOIN sectors s ON s.id = r.sector_id"
                                ).fetchall()
        finally:
            conn.close()
        collisions.append({row["wedge_index"]: objectId.collision(int(row["name"], 16)) for row in rows})
    assert collisions[0] == collisions[1] == {0: 0, 1: 1, 2: 2}
    _assert_same(*results)
