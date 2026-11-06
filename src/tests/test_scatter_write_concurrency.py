# tests/test_scatter_write_concurrency.py

"""DB.113: several scatter workers writing batches at once. Each batch gives its new rows their Object IDs in the same transaction
(`scatter_uids.assign` reads every row of the table that has no ID yet), so under the default isolation level the workers lock each
other's uncommitted rows and MySQL kills one with error 1213 (a deadlock). `scatter_companions.write_committed` writes the batch READ
COMMITTED and starts it over after 1213 or 1205, so every worker's batches land and the IDs are the ones a quiet run gives.

The writers are real connections in real threads against the test database. To make the clash certain instead of likely, each writer waits at a
barrier between its inserts and its ID step the first time (all of them then hold uncommitted rows, which is the state a 15-worker run reaches
by chance); the barrier is the only thing the test adds to the code under test.
"""

import threading

import pymysql
import pytest

from planetgen.db import scatter_uids, store
from planetgen.galaxy import object_uid as ou
from planetgen.generation import scatter_companions

from tests.bright_star_rows import star_row

WRITERS = 8
BATCH = 40
BATCHES = 3
TIMEOUT_S = 120.0


def _stars(writer, batch):
    """One batch of one writer: its own sector, so IDs rank within the writer's rows; luminosity rises with the seed."""
    address = (3, writer, 5)
    return [star_row(address, (1000 + n, 2000, 3000), population="old", star_type="G2V", yerkes_class="V", mass_kg=2e30, radius_km=7e5,
                     temperature_k=5800.0, luminosity_w=1e26 * (1 + n), age_gy=4.0, lifespan_gy=10.0, initial_mass_sol=1.0,
                     metallicity_fe_h=0.0, seed=writer * 10000 + batch * 100 + n)
            for n in range(BATCH)]


def _phenomena(writer, batch):
    values = {name: None for name in store.PHENOMENON_SCATTER_COLUMNS}
    rows = []
    for n in range(BATCH):
        values.update(shell_index=3, band_index=writer, wedge_index=5, kind="black-hole", seed=writer * 10000 + batch * 100 + n,
                      mass_solar=10.0, map_luminosity_sol=1.0 + n, position_x_mpc=1000.0, position_y_mpc=2000.0, position_z_mpc=3000.0)
        rows.append(tuple(values[name] for name in store.PHENOMENON_SCATTER_COLUMNS))
    return rows


def _gated_assign(monkeypatch, parties):
    """`scatter_uids.assign` that waits, once per thread, until every writer has inserted its batch (and so holds its rows)."""
    barrier = threading.Barrier(parties, timeout=TIMEOUT_S)
    seen = threading.local()
    real = scatter_uids.assign

    def assign(conn, table):
        if not getattr(seen, "waited", False):
            seen.waited = True
            barrier.wait()
        return real(conn, table)

    monkeypatch.setattr(scatter_uids, "assign", assign)


def _run(mysql_config, table, make_rows, insert, write):
    """`WRITERS` threads, each writing `BATCHES` batches through `write(conn, insert, rows)`; returns the errors by writer."""
    errors = {}

    def work(writer):
        conn = store.get_connection(mysql_config)
        try:
            for batch in range(BATCHES):
                write(conn, insert, make_rows(writer, batch))
        except pymysql.err.OperationalError as exc:
            errors[writer] = exc.args[0]
            conn.rollback()
        except Exception as exc:            # anything else is a failure of the test or the code, not the clash: show it
            errors[writer] = repr(exc)
        finally:
            conn.close()

    threads = [threading.Thread(target=work, args=(writer,)) for writer in range(WRITERS)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(TIMEOUT_S)
    assert not any(thread.is_alive() for thread in threads), "a writer never finished"
    return errors


def _plain(conn, insert, rows):
    """The write the units made before DB.113: the default isolation level, one commit, no retry."""
    insert(conn, rows, 0.0)
    conn.commit()


def _fixed(conn, insert, rows):
    scatter_companions.write_committed(conn, insert, rows, 0.0)


@pytest.fixture
def conn(mysql_config):
    connection = store.get_connection(mysql_config)
    try:
        yield connection
    finally:
        connection.close()


def _assert_numbered(conn, table, order):
    """Every row has an ID, no two share one, and each writer's sector is numbered 0..n-1 in the order the rows rank by."""
    assert conn.execute(f"SELECT COUNT(*) AS n FROM {table} WHERE uid IS NULL").fetchone()["n"] == 0
    total = WRITERS * BATCHES * BATCH
    assert conn.execute(f"SELECT COUNT(*) AS n, COUNT(DISTINCT uid) AS d FROM {table}").fetchone() == {"n": total, "d": total}
    for writer in range(WRITERS):
        rows = conn.execute(f"SELECT uid FROM {table} WHERE shell_index = 3 AND band_index = ? AND wedge_index = 5 ORDER BY {order}",
                            (writer,)).fetchall()
        serials = [ou.unpack(ou.from_bytes(row["uid"])).serial for row in rows]
        assert len(serials) == BATCHES * BATCH
        first = min(serials)                                      # 0 for stars; phenomena start in the upper half of the serials
        assert sorted(serials) == list(range(first, first + len(serials)))    # no gap, no repeat: a rolled-back batch gave nothing away


def test_the_unfixed_write_deadlocks_when_the_writers_meet_in_the_id_step(mysql_config, monkeypatch):
    """The test's own proof: with the old write at least one of the eight writers is killed with 1213 (or 1205)."""
    _gated_assign(monkeypatch, WRITERS)
    errors = _run(mysql_config, "bright_stars", _stars, scatter_companions.insert_bright_stars, _plain)
    assert errors, "no writer was killed: the test does not reach the clash it is meant to guard against"
    assert set(errors.values()) <= set(store.RETRYABLE_ERRORS)


def test_eight_writers_of_bright_stars_all_land_with_the_ids_a_quiet_run_gives(mysql_config, conn, monkeypatch):
    _gated_assign(monkeypatch, WRITERS)
    errors = _run(mysql_config, "bright_stars", _stars, scatter_companions.insert_bright_stars, _fixed)
    assert errors == {}
    _assert_numbered(conn, "bright_stars", "luminosity_w DESC")
    # a quiet run ranks a batch's rows brightest first, continuing the sector's counter: so the first batch's brightest star is serial 0
    first = conn.execute("SELECT uid FROM bright_stars WHERE seed = ?", (3 * 10000 + 0 * 100 + BATCH - 1,)).fetchone()
    assert ou.unpack(ou.from_bytes(first["uid"])).serial == 0


def test_eight_writers_of_phenomena_all_land_with_distinct_ids(mysql_config, conn, monkeypatch):
    _gated_assign(monkeypatch, WRITERS)
    errors = _run(mysql_config, "phenomenon_scatter", _phenomena, scatter_companions.insert_phenomenon_scatter, _fixed)
    assert errors == {}
    _assert_numbered(conn, "phenomenon_scatter", "map_luminosity_sol DESC")
