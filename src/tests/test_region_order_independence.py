# tests/test_region_order_independence.py

"""
TEST.501: charting regions in either order gives the same galaxy (docs/design/systems/region-object-store.md sections 5 and 8.2).

A region run (`planetgen sectors --around X Y Z --radius-pc R`, GEN.1128) charts the whole sectors a sphere touches. Two
overlapping regions A and B are charted in the two orders, each on a fresh database of the same galaxy seed:

- the sectors charted are the same set, and each sector's digest (GEN.58 over the canonical encoder) is the same, so nothing a
  sector holds depends on which region reached it first, on what a neighbouring region had already charted, or on the
  overlap being skipped the second time;
- the union of the objects the two runs made (the systems and the rogue planets, by object ID and place) is the same;
- charting a region that is already charted, in either order, changes nothing.

Not covered yet, because the pieces are not built: the "made here" and "arrived here" halves of a sector's digest (GEN.1129),
rogues landing in the sector of another donor (GEN.1111) and the top-up pass that closes a region run (GEN.1125, GEN.1127),
whose result the design says does depend on the order (section 8.5). They join this file when they land.
"""

import pytest

from planetgen.db import fingerprint, store

from tests.test_galaxy_gen import _plan_wide_galaxy
from tests.test_order_independence import _assert_same, _fingerprint
from tests.test_parallel_generation import GALAXY_SEED, _run, make_database  # noqa: F401 -- a fixture

pytestmark = [pytest.mark.db, pytest.mark.usefixtures("unit_density_plans", "redis_server")]

REGION_A = ("30", "0", "0")
REGION_B = ("33", "0", "0")
RADIUS = "3"
"""str: Two spheres of 3 pc whose centres are 3 pc apart: each about half inside the other."""


def _chart(config, around, monkeypatch):
    _run("sectors", ["--around", *around, "--radius-pc", RADIUS, "--unit-sectors", "4", "--num-systems", "1", "--workers", "1"],
         config, monkeypatch)


def _addresses(config):
    conn = store.get_connection(config)
    try:
        return sorted((row["shell_index"], row["band_index"], row["wedge_index"])
                      for row in conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors").fetchall())
    finally:
        conn.close()


def _objects(config):
    """`{(table, object ID): place}` of every system and rogue planet in the galaxy, whichever sector holds it (a rogue's place is
    inside its digest, so only its ID is listed here)."""
    found = {}
    conn = store.get_connection(config)
    try:
        for row in conn.execute("SELECT uid, galaxy_x_pc, galaxy_y_pc, galaxy_z_pc FROM star_systems").fetchall():
            found[("star_systems", bytes(row["uid"]).hex())] = tuple(round(float(row[axis]), 6) for axis in
                                                                      ("galaxy_x_pc", "galaxy_y_pc", "galaxy_z_pc"))
        for row in conn.execute("SELECT uid FROM rogue_planets").fetchall():
            found[("rogue_planets", bytes(row["uid"]).hex())] = None
    finally:
        conn.close()
    return found


def _galaxy(make_database):
    config = make_database()
    _plan_wide_galaxy(config, galaxy_seed=GALAXY_SEED)
    return config


def _charted(make_database, monkeypatch, order):
    config = _galaxy(make_database)
    for around in order:
        _chart(config, around, monkeypatch)
    addresses = _addresses(config)
    return config, addresses, _fingerprint(config, addresses), _objects(config)


def test_two_overlapping_regions_chart_the_same_galaxy_in_either_order(make_database, monkeypatch):
    first = _charted(make_database, monkeypatch, (REGION_A, REGION_B))
    second = _charted(make_database, monkeypatch, (REGION_B, REGION_A))
    assert first[1] == second[1], "the two orders did not chart the same sectors"
    assert len(first[1]) > 1
    assert first[3].keys() == second[3].keys(), "the two orders did not make the same objects"
    assert first[3] == second[3], "an object stands in a different place"
    _assert_same(first[2], second[2])


def test_the_overlap_is_charted_once_and_a_second_run_of_a_region_changes_nothing(make_database, monkeypatch):
    config, addresses, before, objects = _charted(make_database, monkeypatch, (REGION_A, REGION_B))
    for around in (REGION_A, REGION_B, REGION_A):
        _chart(config, around, monkeypatch)
    assert _addresses(config) == addresses
    assert _objects(config) == objects
    _assert_same(_fingerprint(config, addresses), before)


def test_a_region_charted_alone_is_the_same_inside_a_bigger_chart(make_database, monkeypatch):
    """The sectors of A come out the same whether A was charted alone or as the first half of A then B."""
    alone = _charted(make_database, monkeypatch, (REGION_A,))
    both = _charted(make_database, monkeypatch, (REGION_A, REGION_B))
    assert set(alone[1]) <= set(both[1])
    conn = store.get_connection(both[0])
    try:
        digests_both = dict(fingerprint.region_fingerprint(conn, addresses=alone[1]).sectors)
    finally:
        conn.close()
    assert digests_both == dict(alone[2][0].sectors), "a sector of A changed when B was charted after it"
