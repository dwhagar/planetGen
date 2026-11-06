# tests/test_order_independence.py

"""
TEST.181: a sector's contents do not depend on the order sectors were generated in or on the worker count (GEN.57,
GEN.306, GEN.307). The same small set of sectors is generated in two orders and at one worker and four, and the
fingerprints of the set (GEN.58 over the canonical encoder of GEN.304) are compared:

- sectors that draw the same point hold different Object ID names whichever saved first (GEN.307, GEN.1256);
- a sector reached by the bright-star backfill holds the same stars whichever sector the run reached first;
- the same sectors at one worker and four are the same sectors.

Nearest-system links and the population data are rebuilt from content, so they are not in the comparison
(`db.fingerprint` leaves them out). Name decorations are compared as they are.
"""

import pytest

from planetgen.db import fingerprint, models, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.galaxy import object_uid as objectUid
from planetgen.util import draw

from tests.test_galaxy_gen import _plan_wide_galaxy
from tests.test_object_id_names import _POSITION
from tests.test_parallel_generation import GALAXY_SEED, _run, make_database  # noqa: F401 -- a fixture

pytestmark = pytest.mark.db

SECTORS = ((1, 0, 0), (1, 0, 5), (1, 0, 11), (1, 1, 0))
"""tuple: `(shell, band, wedge)` of the set: three of the equator's twelve cells of shell 1 and the cell above them."""


def _canonical_scatter_serials(conn, addresses):
    """
    GEN.1302: gives the scatter-born objects of each sector the serials their brightness order would give them in one call, so the comparison
    does not see how the writes were batched. A scatter row's serial is its rank among the rows written in the same call, continuing the
    sector's counter (`db.scatter_uids`), so a sector whose dimmer stars were written by an earlier backfill (serials 0 and 1) and whose
    brighter ones by the first fill (2 and 3) holds the same four systems as one whose first fill wrote all four at once, with the
    serials in another order. That is by design: scatter serials depend on the run order (docs/design/architecture/object-ids.md,
    "Scatter rows"); what must not depend on it is everything else, which this test compares. The serials a sector's scatter stars hold are
    kept as a set and handed out again by `luminosity_w` (then `seed`), in every table that carries an Object ID.
    """
    tables = [name for name, table in models.metadata.tables.items() if "uid" in {column.name for column in table.columns}]
    for address in addresses:
        sector = objectUid.sector_bytes(*address)
        rows = conn.execute("SELECT uid, luminosity_w, seed FROM bright_stars WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                            address).fetchall()
        if not rows:
            continue
        held = sorted(objectUid.unpack(objectUid.from_bytes(row["uid"])).serial for row in rows)
        by_brightness = sorted(rows, key=lambda row: (-row["luminosity_w"], row["seed"]))
        remap = {objectUid.unpack(objectUid.from_bytes(row["uid"])).serial: serial for row, serial in zip(by_brightness, held)}
        if all(old == new for old, new in remap.items()):
            continue
        shift = 1 << 20   # two passes, so no row meets a unique key that another row still holds
        for phase in (0, 1):
            for table in tables:
                for row in conn.execute(f"SELECT uid FROM `{table}` WHERE uid IS NOT NULL AND LEFT(uid, ?) = ?",
                                        (len(sector), sector)).fetchall():
                    value = objectUid.from_bytes(row["uid"])
                    parts = objectUid.unpack(value)
                    if parts.serial_kind != objectUid.SERIAL_GENERATED:
                        continue
                    if phase:
                        if parts.serial < shift:
                            continue
                        target = parts.serial - shift
                    else:
                        if parts.serial not in remap:
                            continue
                        target = remap[parts.serial] + shift
                    moved = objectUid.pack(parts.shell, parts.band, parts.wedge, parts.serial_kind, target, parts.body)
                    conn.execute(f"UPDATE `{table}` SET uid = ? WHERE uid = ?", (objectUid.to_bytes(moved), row["uid"]))
    conn.commit()


def _fingerprint(config, addresses=SECTORS):
    """`(fingerprint, sections)` of the region: the sectors' digests and, per sector, the digests of its named sections
    (`fingerprint.sector_sections`), which say what differs when two do. The scatter serials are put in brightness order first
    (`_canonical_scatter_serials`)."""
    conn = store.get_connection(config)
    try:
        _canonical_scatter_serials(conn, addresses)
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


def _generate_in_order(config, order, monkeypatch, backfill=False):
    """One `sectors` run a sector, in `order` (a span of that one cell, 4 systems), then each run's backfill when
    `backfill` (the `backfill` command, for the sectors that run saved)."""
    for shell, band, wedge in order:
        _run("sectors", ["--shells", str(shell), "--bands", str(band), "--wedges", str(wedge), "--num-systems", "4",
                         "--workers", "1"], config, monkeypatch)
        if backfill:
            _run("backfill", [], config, monkeypatch)


def test_the_same_sectors_come_out_in_either_order(make_database, monkeypatch):
    results = []
    for order in (SECTORS, SECTORS[::-1]):
        config = _galaxy(make_database)
        _generate_in_order(config, order, monkeypatch)
        results.append(_fingerprint(config))
    assert len(results[0][0].sectors) == len(SECTORS)
    _assert_same(*results)


def test_a_sector_reached_by_the_backfill_comes_out_the_same_in_either_order(make_database, monkeypatch):
    """The first run generates one sector and backfills the bright stars around it (the `backfill` command), so the
    sectors after it are reached by that backfill; they then come out the same in whichever order they are
    generated. The case where the backfill reaches a sector that is already generated is the last test of the file
    (GEN.633)."""
    first, rest = SECTORS[:1], SECTORS[1:]
    results = []
    for order in (rest, rest[::-1]):
        config = _galaxy(make_database)
        _generate_in_order(config, first, monkeypatch, backfill=True)
        _generate_in_order(config, order, monkeypatch)
        results.append(_fingerprint(config))
    assert len(results[0][0].sectors) == len(SECTORS)
    _assert_same(*results)


@pytest.mark.parametrize("workers", [4])
def test_the_same_sectors_come_out_at_one_worker_and_several(make_database, monkeypatch, workers):
    results = []
    for count in (1, workers):
        config = _galaxy(make_database)
        _run("sectors", ["--shells", "1", "--num-systems", "4", "--workers", str(count)], config, monkeypatch)
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


def test_sectors_that_draw_one_point_come_out_the_same_in_either_order(make_database):
    """Three sectors whose rogue planets sit at one point: each planet's name is its own Object ID (GEN.1256), and the whole
    stored sector, not only the name, is the same whichever saved first."""
    results, names = [], []
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
        names.append({row["wedge_index"]: row["name"] for row in rows})
    assert names[0] == names[1] and len(set(names[0].values())) == 3
    _assert_same(*results)


def test_a_sector_comes_out_the_same_whether_the_backfill_reached_it_first_or_not(make_database, monkeypatch):
    """GEN.633: each run backfills the bright stars around its sector (the default `--backfill-from edge`), which used to
    reach only sectors not yet generated, so a sector generated before the backfill reached it held fewer systems. The
    first fill now takes a sector to the backfill's deepest floor itself, so the four sectors hold the same systems in
    either order."""
    results = []
    for order in (SECTORS, SECTORS[::-1]):
        config = _galaxy(make_database)
        _generate_in_order(config, order, monkeypatch, backfill="edge")
        results.append(_fingerprint(config))
    assert len(results[0][0].sectors) == len(SECTORS)
    _assert_same(*results)
