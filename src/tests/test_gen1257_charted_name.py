# tests/test_gen1257_charted_name.py

"""
GEN.1257: when charting gives an uncharted star its own word-salad name, the name its Object ID spelled until then (a local name, GEN.1255)
is kept as the first retired name of its system, reason `charted`: the "Formerly" line shows it, and a search for it finds the system.
This file runs the single-star route (`POST /api/uncharted-systems/<id>/generate`); the sector fill's own path is checked with the fill that
builds scatter systems, in `test_gen1254_scatter_uid.py` (a charted sector costs ten minutes cold, so it is built once).
"""

import pytest

from planetgen.api import lookup
from planetgen.db import designation_rows, store
from planetgen.galaxy import designations, object_uid as ou
from planetgen.names import naming_key
from planetgen.web import formerly

pytestmark = pytest.mark.db

KEY = "00AB12CD"


@pytest.fixture(autouse=True)
def key():
    naming_key.set_resolver(lambda: KEY)
    yield
    naming_key.set_resolver(None)


def test_the_uncharted_route_keeps_the_name_too_and_keeping_it_twice_leaves_one_row(mysql_config):
    from tests import test_bright_star_scatter as scatter_tests
    from planetgen.api import routes
    from planetgen.generation import run_plan

    scatter_tests._seed_galaxy(mysql_config, e_value=scatter_tests.SMALL_E_VALUE)
    run_plan.scatter_bright_stars(scatter_tests._scatter_args(mysql_config, "--mass"))
    conn = store.get_connection(mysql_config)
    try:
        star = conn.execute("SELECT id, uid FROM bright_stars ORDER BY id LIMIT 1").fetchone()
    finally:
        conn.close()
    routes.generate_uncharted_system_job(mysql_config, star["id"])
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT star_system_id FROM bright_stars WHERE id = ?", (star["id"],)).fetchone()["star_system_id"]
        old_name = naming_key.local_name(ou.from_bytes(star["uid"]), KEY)
        rows = conn.execute("SELECT designation_name, reason FROM designations WHERE object_kind = 'star_systems' AND object_row_id = ?",
                            (system_id,)).fetchall()
        assert [(row["designation_name"], row["reason"]) for row in rows] == [(old_name, "charted")]
        status, body = lookup.lookup(conn, old_name)
        assert status == 200 and body["status"] == designations.REDIRECT
        with conn:
            assert store.keep_charted_name(conn, system_id, star["uid"]) == old_name       # again: still one row
        assert conn.execute("SELECT COUNT(*) AS n FROM designations WHERE reason = 'charted' AND object_row_id = ?",
                            (system_id,)).fetchone()["n"] == 1
        naming_key.set_resolver(None)
        before = conn.execute("SELECT COUNT(*) AS n FROM designations").fetchone()["n"]
        with conn:
            assert store.keep_charted_name(conn, system_id, star["uid"]) is None        # no key in force: nothing can be spelled
            assert store.keep_charted_name(conn, system_id, None) is None               # not a scatter star
        assert conn.execute("SELECT COUNT(*) AS n FROM designations").fetchone()["n"] == before
    finally:
        conn.close()


def test_charted_is_a_reason_and_every_reason_has_its_words():
    assert "charted" in designations.REASONS
    assert designations.check_reason("charted") == "charted"
    assert set(designations.REASONS) <= set(formerly.REASON_TEXT)
