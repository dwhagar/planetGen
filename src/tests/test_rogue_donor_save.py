# tests/test_rogue_donor_save.py

"""GEN.1111 part 7: a donor system saved in a sector writes the rogue planets it threw out, each in the sector where it landed, once."""

import pytest

from planetgen.db import store
from planetgen.galaxy import object_uid, spherical
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import rogue_yield
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw

pytestmark = pytest.mark.fresh_database

EDGE_LY = 10.0
EDGE_PC = ly_to_pc(EDGE_LY)
DONOR_PC = (8000.0, 100.0, 20.0)
UID = bytes(range(1, 13))


def _donor():
    """A G-star system that threw something out (the first of a few seeds that did)."""
    for seed in range(40):
        config = SystemConfig()
        config.STAR_TYPE = "G2V"
        with draw.bound(f"rogue-donor-{seed}"):
            system = StarSystem(config)
        if rogue_yield.rogues_of(UID, system, DONOR_PC, EDGE_PC):
            return system
    pytest.skip("no system of forty threw a body out")


@pytest.fixture(scope="module")
def donor():
    return _donor()


def test_the_rogues_of_a_donor_are_a_function_of_its_uid_and_record(donor):
    first = rogue_yield.rogues_of(UID, donor, DONOR_PC, EDGE_PC)
    again = rogue_yield.rogues_of(UID, donor, DONOR_PC, EDGE_PC)
    assert [r.event.index for r in first] == [r.event.index for r in again]
    assert [r.planet.mass for r in first] == [r.planet.mass for r in again]
    assert [r.landing.address for r in first] == [r.landing.address for r in again]
    other = rogue_yield.rogues_of(bytes(reversed(UID)), donor, DONOR_PC, EDGE_PC)
    assert [r.planet.mass for r in other] != [r.planet.mass for r in first] or [r.landing.position_pc for r in other] != [r.landing.position_pc for r in first]
    assert len({r.event.index for r in first}) == len(first)
    for rogue in first:
        assert "Thrown out" in rogue.note and rogue.landing.separation_pc > 0.0


@pytest.fixture
def saved(donor, mysql_config, redis_server):
    address = spherical.sector_address_at(DONOR_PC, EDGE_PC)
    sector = SpaceSector("Donor Sector", edge_ly=EDGE_LY)
    sector.add_system(donor, position=(0.0, 0.0, 0.0))
    x, y, z = spherical.sector_position_pc(*address, EDGE_PC)
    position = {"center_x_pc": x, "center_y_pc": y, "center_z_pc": z, "galactic_radius_pc": (x * x + y * y + z * z) ** 0.5,
                "shell_index": address[0], "band_index": address[1], "wedge_index": address[2]}
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=position)
    conn = store.get_connection(mysql_config)
    try:
        yield conn, sector_id, address
    finally:
        conn.close()


def test_saving_a_donor_writes_its_rogues_in_the_sectors_where_they_landed(saved, donor):
    conn, sector_id, address = saved
    uid = conn.execute("SELECT uid FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["uid"]
    rows = conn.execute("SELECT id, uid, sector_id, event_index FROM rogue_planets WHERE donor_uid = ? ORDER BY event_index", (uid,)).fetchall()
    made = rogue_yield.rogues_of(bytes(uid), donor, _donor_position(conn, sector_id), EDGE_PC)
    assert rows and len(rows) == len(made)
    for row, rogue in zip(rows, made):
        assert row["event_index"] == rogue.event.index and row["uid"] is not None
        filed = conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors WHERE id = ?", (row["sector_id"],)).fetchone()
        assert (filed["shell_index"], filed["band_index"], filed["wedge_index"]) == tuple(rogue.landing.address)
        assert bytes(row["uid"])[:7] == object_uid.sector_bytes(*address)       # born in the donor's sector, whichever sector files it


def _donor_position(conn, sector_id):
    row = conn.execute("SELECT center_x_pc, center_y_pc, center_z_pc FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    return (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])


def test_saving_the_donor_again_gives_the_same_rogues_not_more(saved, donor):
    conn, sector_id, address = saved
    system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
    before = conn.execute("SELECT COUNT(*) AS n FROM rogue_planets").fetchone()["n"]
    placement = dict(zip(("center_x_pc", "center_y_pc", "center_z_pc"), _donor_position(conn, sector_id)))
    with conn:
        written = store.save_donor_rogues(conn, [(system_id, donor, placement)], sector_id, EDGE_PC)
    assert written == before == conn.execute("SELECT COUNT(*) AS n FROM rogue_planets").fetchone()["n"]


def test_a_landing_cell_that_was_never_charted_gets_a_row_with_no_ledger_entry(saved):
    conn, sector_id, address = saved
    rows = conn.execute("SELECT DISTINCT sector_id FROM rogue_planets").fetchall()
    for row in rows:
        cell = conn.execute("SELECT shell_index, band_index, wedge_index FROM sectors WHERE id = ?", (row["sector_id"],)).fetchone()
        cell = (cell["shell_index"], cell["band_index"], cell["wedge_index"])
        if cell != tuple(address):
            assert store.chart_ledger.state(conn, cell, 0) == "uncharted"


def test_deleting_the_donors_sector_takes_its_rogues_wherever_they_landed(saved):
    from planetgen.db import edits
    conn, sector_id, address = saved
    assert conn.execute("SELECT COUNT(*) AS n FROM rogue_planets").fetchone()["n"] > 0
    with conn:
        edits.delete_sector_with_contents(conn, sector_id)
    assert conn.execute("SELECT COUNT(*) AS n FROM rogue_planets").fetchone()["n"] == 0
    assert conn.execute("SELECT COUNT(*) AS n FROM moons WHERE rogue_planet_id IS NOT NULL").fetchone()["n"] == 0


def test_a_sector_that_holds_a_rogue_born_elsewhere_cannot_be_deleted(saved):
    from planetgen.db import edits
    conn, sector_id, address = saved
    other = conn.execute("SELECT DISTINCT sector_id FROM rogue_planets WHERE sector_id <> ?", (sector_id,)).fetchone()
    if other is None:
        pytest.skip("every rogue landed in its donor's own sector")
    with pytest.raises(edits.SectorHasResidents):
        with conn:
            edits.delete_sector_with_contents(conn, other["sector_id"])
