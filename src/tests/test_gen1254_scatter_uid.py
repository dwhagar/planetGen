# tests/test_gen1254_scatter_uid.py

"""GEN.1254: every scatter row has an Object ID (`uid`) from the moment it is written (`db/scatter_uids.py`). The serial is the row's rank in
its sector, brightest first, continuing from the sector's counter; stars and phenomena share a sector without sharing serials; the same
rows written again give the same IDs; a sector's own fill starts past the stars' serials; a system built round a scatter star keeps its ID."""

import pytest

from planetgen.db import scatter_uids, store
from planetgen.galaxy import object_uid as ou

SOLAR = 3.828e26


def _star(address, luminosity, seed, x=1000.0):
    values = dict(
        shell_index=address[0], band_index=address[1], wedge_index=address[2], population="old", star_type="G2V", yerkes_class="V",
        mass_kg=2e30, radius_km=7e5, temperature_k=5800.0, luminosity_w=luminosity * SOLAR, age_gy=4.0, lifespan_gy=10.0,
        initial_mass_sol=1.0, phase_end_age_gy=None, seed=seed, metallicity_fe_h=0.0, rotation_period_hours=600.0,
        nuclear_arrival_age_gy=None, position_x_mpc=x, position_y_mpc=2000.0, position_z_mpc=3000.0)
    return tuple(values[name] for name in store.BRIGHT_STAR_COLUMNS)


def _phenomenon(address, luminosity, seed, kind="black-hole"):
    values = {name: None for name in store.PHENOMENON_SCATTER_COLUMNS}
    values.update(shell_index=address[0], band_index=address[1], wedge_index=address[2], kind=kind, seed=seed, mass_solar=10.0,
                  map_luminosity_sol=luminosity, position_x_mpc=1000.0, position_y_mpc=2000.0, position_z_mpc=3000.0)
    return tuple(values[name] for name in store.PHENOMENON_SCATTER_COLUMNS)


def _uids(conn, table, address, order="id"):
    rows = conn.execute(
        f"SELECT uid FROM {table} WHERE shell_index = ? AND band_index = ? AND wedge_index = ? ORDER BY {order}", address).fetchall()
    return [ou.unpack(ou.from_bytes(row["uid"])) for row in rows]


@pytest.fixture
def conn(mysql_config):
    connection = store.get_connection(mysql_config)
    try:
        yield connection
    finally:
        connection.close()


def test_a_sectors_stars_are_numbered_brightest_first_in_the_sector_they_sit_in(conn):
    a, b = (5, 1, 7), (5, 2, 9)
    rows = [_star(a, 1.0, 11), _star(a, 9.0, 12), _star(b, 3.0, 13), _star(a, 4.0, 14)]
    with conn:
        assert store.insert_bright_stars(conn, rows) == 4
    ids = _uids(conn, "bright_stars", a, "luminosity_w DESC")
    assert [i.serial for i in ids] == [0, 1, 2]
    assert all((i.shell, i.band, i.wedge, i.serial_kind, i.body) == (*a, ou.SERIAL_GENERATED, 0) for i in ids)
    assert [i.serial for i in _uids(conn, "bright_stars", b)] == [0]
    assert conn.execute("SELECT COUNT(DISTINCT uid) AS n, COUNT(*) AS m FROM bright_stars").fetchone() == {"n": 4, "m": 4}


def test_a_later_write_for_the_sector_continues_the_numbers_and_never_repeats_one(conn):
    a = (5, 1, 7)
    with conn:
        store.insert_bright_stars(conn, [_star(a, 5.0, 1), _star(a, 6.0, 2)])
        store.insert_bright_stars(conn, [_star(a, 50.0, 3), _star(a, 0.5, 4)])
    by_seed = {row["seed"]: ou.unpack(ou.from_bytes(row["uid"])).serial
               for row in conn.execute("SELECT seed, uid FROM bright_stars").fetchall()}
    assert by_seed[2] == 0 and by_seed[1] == 1                 # the first write: brightest first
    assert by_seed[3] == 2 and by_seed[4] == 3                 # the second continues, brightest of its own rows first
    assert scatter_uids.next_star_serial(conn, *a) == 4 and scatter_uids.next_star_serial(conn, 9, 9, 9) == 0


def test_scattered_phenomena_take_the_upper_half_and_never_a_stars_serial(conn):
    a = (5, 1, 7)
    with conn:
        store.insert_bright_stars(conn, [_star(a, 5.0, 1)])
        store.insert_phenomenon_scatter(conn, [_phenomenon(a, 0.1, 21), _phenomenon(a, 2.0, 22, kind="neutron-star")])
    stars = {i.serial for i in _uids(conn, "bright_stars", a)}
    phenomena = _uids(conn, "phenomenon_scatter", a, "map_luminosity_sol DESC")
    assert [i.serial for i in phenomena] == [scatter_uids.PHENOMENON_SERIAL_BASE, scatter_uids.PHENOMENON_SERIAL_BASE + 1]
    assert stars == {0} and not stars & {i.serial for i in phenomena}


def test_the_same_rows_written_again_after_a_clear_get_the_same_ids(conn):
    a = (5, 1, 7)
    rows = [_star(a, 5.0, 1), _star(a, 6.0, 2), _star(a, 7.0, 3)]
    with conn:
        store.insert_bright_stars(conn, rows)
    first = {row["seed"]: bytes(row["uid"]) for row in conn.execute("SELECT seed, uid FROM bright_stars").fetchall()}
    store.clear_bright_stars(conn)
    with conn:
        store.insert_bright_stars(conn, list(reversed(rows)))        # inserted in another order: the rank does not depend on it
    again = {row["seed"]: bytes(row["uid"]) for row in conn.execute("SELECT seed, uid FROM bright_stars").fetchall()}
    assert first == again and len(set(first.values())) == 3


def test_a_row_written_with_no_id_by_other_means_stays_empty_until_the_writer_runs(conn):
    a = (5, 1, 7)
    columns = ", ".join(store.BRIGHT_STAR_COLUMNS)
    marks = ", ".join("?" for _ in store.BRIGHT_STAR_COLUMNS)
    with conn:
        conn.execute(f"INSERT INTO bright_stars ({columns}) VALUES ({marks})", _star(a, 1.0, 1))
        conn.execute(f"INSERT INTO bright_stars ({columns}) VALUES ({marks})", _star(a, 2.0, 2))
    assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE uid IS NULL").fetchone()["n"] == 2    # NULLs do not collide
    with conn:
        assert scatter_uids.assign(conn, "bright_stars") == 2
    assert [i.serial for i in _uids(conn, "bright_stars", a, "luminosity_w DESC")] == [0, 1]


def test_only_the_scatter_tables_are_numbered(conn):
    with pytest.raises(ValueError, match="not a scatter table"):
        scatter_uids.assign(conn, "stars")


def test_a_system_built_round_a_scatter_star_keeps_its_id_and_the_fill_starts_past_the_stars(mysql_config):
    from tests import test_bright_star_scatter as scatter_tests
    from planetgen.generation import run_galaxy, run_plan
    from planetgen.galaxy.spherical import sector_position_pc

    scatter_tests._seed_galaxy(mysql_config, e_value=scatter_tests.SMALL_E_VALUE)
    run_plan.scatter_bright_stars(scatter_tests._scatter_args(mysql_config, "--mass"))
    conn = store.get_connection(mysql_config)
    try:
        target = conn.execute(
            "SELECT shell_index, band_index, wedge_index, COUNT(*) AS n FROM bright_stars"
            " GROUP BY shell_index, band_index, wedge_index ORDER BY n DESC LIMIT 1").fetchone()
        address = (target["shell_index"], target["band_index"], target["wedge_index"])
        assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE uid IS NULL").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(DISTINCT uid) AS n FROM bright_stars").fetchone()["n"] == conn.execute(
            "SELECT COUNT(*) AS n FROM bright_stars").fetchone()["n"]
    finally:
        conn.close()

    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    sector_id, _name, _sector = run_galaxy.generate_and_save_sector_at(
        args, address, sector_position_pc(*address, scatter_tests.EDGE_PC), scatter_tests.EDGE_PC)

    conn = store.get_connection(mysql_config)
    try:
        held = conn.execute("SELECT MAX(uid) AS top FROM bright_stars WHERE shell_index = ? AND band_index = ? AND wedge_index = ?",
                            address).fetchone()["top"]
        stars_end = ou.unpack(ou.from_bytes(held)).serial + 1
        assert scatter_uids.next_star_serial(conn, *address) >= stars_end
        built = conn.execute(
            "SELECT b.uid AS star_uid, ss.uid AS system_uid FROM bright_stars b JOIN star_systems ss ON ss.id = b.star_system_id"
            " WHERE ss.sector_id = ?", (sector_id,)).fetchall()
        assert built and all(bytes(row["star_uid"]) == bytes(row["system_uid"]) for row in built)
        others = conn.execute(
            "SELECT uid FROM star_systems WHERE sector_id = ? AND uid NOT IN (SELECT uid FROM bright_stars WHERE uid IS NOT NULL)",
            (sector_id,)).fetchall()
        floor = scatter_uids.next_star_serial(conn, *address)
        assert all(ou.unpack(ou.from_bytes(row["uid"])).serial >= floor for row in others)
        assert conn.execute("SELECT COUNT(*) AS n, COUNT(DISTINCT uid) AS m FROM star_systems WHERE sector_id = ?",
                            (sector_id,)).fetchone()["n"] == conn.execute(
            "SELECT COUNT(DISTINCT uid) AS m FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchone()["m"]
    finally:
        conn.close()
