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
ADDRESS = spherical.sector_address_at(DONOR_PC, EDGE_PC)
CENTER_PC = spherical.sector_position_pc(*ADDRESS, EDGE_PC)


def _real_uid():
    """The ID the donor really gets when it is the only system saved in its sector: the sector's first generated serial."""
    return object_uid.to_bytes(object_uid.pack(*ADDRESS, object_uid.SERIAL_GENERATED, 0))


def _donor():
    """A G-star system that threw something out, by the ID and the place it will really be saved with (the first of a few seeds that did)."""
    for seed in range(40):
        config = SystemConfig()
        config.STAR_TYPE = "G2V"
        with draw.bound(f"rogue-donor-{seed}"):
            system = StarSystem(config)
        if rogue_yield.rogues_of(_real_uid(), system, CENTER_PC, EDGE_PC):
            return system
    pytest.skip("no system of forty threw a body out")


@pytest.fixture(scope="module")
def donor():
    return _donor()


def test_the_rogues_of_a_donor_are_a_function_of_its_uid_and_record(donor):
    first = rogue_yield.rogues_of(_real_uid(), donor, CENTER_PC, EDGE_PC)
    again = rogue_yield.rogues_of(_real_uid(), donor, CENTER_PC, EDGE_PC)
    assert [r.event.index for r in first] == [r.event.index for r in again]
    assert [r.planet.mass for r in first] == [r.planet.mass for r in again]
    assert [r.landing.address for r in first] == [r.landing.address for r in again]
    other = rogue_yield.rogues_of(bytes(reversed(_real_uid())), donor, CENTER_PC, EDGE_PC)
    assert [r.planet.mass for r in other] != [r.planet.mass for r in first] or [r.landing.position_pc for r in other] != [r.landing.position_pc for r in first]
    assert len({r.event.index for r in first}) == len(first)
    for rogue in first:
        assert "Thrown out" in rogue.note and rogue.landing.separation_pc > 0.0


@pytest.fixture
def saved(donor, mysql_config, redis_server):
    address = ADDRESS
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
    assert bytes(uid) == _real_uid()                # the donor was chosen by the ID it is saved with
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


def _leaf(conn, sector_id):
    row = conn.execute("SELECT leaf_digest, content_leaf FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    return (None if row["leaf_digest"] is None else bytes(row["leaf_digest"]),
            None if row["content_leaf"] is None else bytes(row["content_leaf"]))


def test_a_sector_with_a_rogue_in_it_has_a_leaf_that_holds_it_charted_or_not(saved):
    """TEST.602: the leaf the fingerprint tree reads is the content part joined to the rogues filed in the sector, and a sector made only for a landing has
    the rogues alone."""
    from planetgen.db import fingerprint_tree
    conn, sector_id, address = saved
    filed = [row["sector_id"] for row in conn.execute("SELECT DISTINCT sector_id FROM rogue_planets").fetchall()]
    assert filed
    for filed_in in filed:
        residents = fingerprint_tree.resident_digest(conn, filed_in)
        leaf, content = _leaf(conn, filed_in)
        assert residents is not None and leaf == fingerprint_tree.combined(content, residents)
        if filed_in != sector_id:               # a landing sector: never charted, no content part
            assert content is None and leaf == fingerprint_tree.combined(None, residents) and leaf is not None
        else:
            assert content is not None and leaf != content


def test_the_leaf_of_a_sector_with_no_rogue_is_its_content_part_and_saving_a_donor_again_leaves_every_leaf_as_it_was(saved, donor):
    conn, sector_id, address = saved
    before = {row["id"]: _leaf(conn, row["id"]) for row in conn.execute("SELECT id FROM sectors").fetchall()}
    system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["id"]
    placement = dict(zip(("center_x_pc", "center_y_pc", "center_z_pc"), _donor_position(conn, sector_id)))
    with conn:
        store.save_donor_rogues(conn, [(system_id, donor, placement)], sector_id, EDGE_PC)
    after = {row["id"]: _leaf(conn, row["id"]) for row in conn.execute("SELECT id FROM sectors").fetchall()}
    assert before == after


def test_the_landing_sectors_leaves_go_when_the_donors_sector_is_emptied(saved):
    from planetgen.db import edits
    conn, sector_id, address = saved
    landing = [row["sector_id"] for row in conn.execute("SELECT DISTINCT sector_id FROM rogue_planets WHERE sector_id <> ?", (sector_id,)).fetchall()]
    if not landing:
        pytest.skip("every rogue landed in its donor's own sector")
    with conn:
        edits.delete_sector_birth_contents(conn, sector_id)
    assert all(_leaf(conn, filed_in) == (None, None) for filed_in in landing)
    marked = conn.execute("SELECT COUNT(*) AS n FROM sector_tree_nodes WHERE dirty = 1").fetchone()["n"]
    assert marked > 0


def test_the_leaf_joins_the_content_part_and_the_residents():
    from planetgen.db import fingerprint_tree as tree
    content, residents = bytes([1]) * 32, bytes([2]) * 32
    assert tree.combined(content, None) == content and tree.combined(None, None) is None
    assert len(tree.combined(content, residents)) == 32 and tree.combined(content, residents) != content
    assert tree.combined(None, residents) != tree.combined(content, residents)
    assert tree.combined(content, residents) != tree.combined(content, bytes([3]) * 32)


CHILD_CHANGES = {
    "rogue_planet_reflection_spectrum": "value = CONCAT(value, 'x')",
    "rogue_planet_interior_layers": "mass_earth = mass_earth * 1.01",
    "rogue_planet_gases": "partial_pa = partial_pa * 2",
    "rogue_planet_life_class_scores": "score = IF(score = 9, 8, score + 1)",
}


def test_a_changed_row_of_a_rogue_planets_child_tables_changes_the_digest_of_the_sector_it_is_filed_in(saved):
    """TEST.647: the four child tables of a rogue planet are content of the fingerprint (the planet section) and of the leaf's residents digest, so a
    changed child row changes the digest of the sector the rogue is filed in, and the digest is the same again once the row is."""
    from planetgen.db import fingerprint, fingerprint_tree
    conn, sector_id, address = saved
    rogue = conn.execute("SELECT id FROM rogue_planets ORDER BY id LIMIT 1").fetchone()["id"]
    minimal = {"rogue_planet_reflection_spectrum": "(rogue_planet_id, spectrum_type, position, value) VALUES (?, 'visible', 0, 'blue')",
               "rogue_planet_interior_layers": "(rogue_planet_id, layer_index, kind, state, top_pressure_gpa, bottom_pressure_gpa, mass_earth, top_radius_fraction, "
                                               "bottom_radius_fraction, hydrogen, helium, iron, rock, water_liquid, ice_x, nh3_ice, ch4_ice) "
                                               "VALUES (?, 0, 'rock_core', 'molten', 0, 10, 0.1, 0.5, 1, 0, 0, 0.3, 0.7, 0, 0, 0, 0)",
               "rogue_planet_life_class_scores": "(rogue_planet_id, life_class, score) VALUES (?, 'A1', 3)"}
    seen = set()
    for table, change in CHILD_CHANGES.items():
        row = conn.execute(f"SELECT t.id AS id, r.sector_id AS sector_id FROM {table} t JOIN rogue_planets r ON r.id = t.rogue_planet_id "
                           f"ORDER BY t.id LIMIT 1").fetchone()
        if row is None:             # this donor's rogues hold no row of the table (a rock with no air has no gases, a small icy body no layers): give the first rogue one
            with conn:
                conn.execute(f"INSERT INTO {table} {minimal[table]}", (rogue,))
            row = conn.execute(f"SELECT t.id AS id, r.sector_id AS sector_id FROM {table} t JOIN rogue_planets r ON r.id = t.rogue_planet_id "
                               f"ORDER BY t.id LIMIT 1").fetchone()
        filed_in = row["sector_id"]
        digest, residents = fingerprint.sector_digest(conn, filed_in), fingerprint_tree.resident_digest(conn, filed_in)
        original = conn.execute(f"SELECT * FROM {table} WHERE id = ?", (row["id"],)).fetchone()
        with conn:
            conn.execute(f"UPDATE {table} SET {change} WHERE id = ?", (row["id"],))
        try:
            assert fingerprint.sector_digest(conn, filed_in) != digest, table
            assert fingerprint_tree.resident_digest(conn, filed_in) != residents, table
        finally:
            with conn:
                conn.execute(f"DELETE FROM {table} WHERE id = ?", (row["id"],))
                columns = list(original.keys())
                conn.execute(f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({', '.join('?' for _ in columns)})", tuple(original[c] for c in columns))
        assert fingerprint.sector_digest(conn, filed_in) == digest and fingerprint_tree.resident_digest(conn, filed_in) == residents, table
        seen.add(table)
    assert seen == set(CHILD_CHANGES)
