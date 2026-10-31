# tests/test_api_retired_ids.py

"""
GEN.1139: an ID an object has left. A read of it is a 301 to the same URL under the object's current ID, a write against it is a 409 that changes
nothing, and a deleted object's ID answers 410. A page asked for one is redirected the same way. A retired ID is not any other object's, and an ID
nobody held is still a 404.
"""

import pytest

from planetgen.db import designation_history, store
from planetgen.galaxy import object_uid
from tests.publicids import pid
from tests.test_admin_edits import _saved_nebula, admin, client  # noqa: F401 -- fixtures

pytestmark = pytest.mark.db


def _renumbered(mysql_config):
    """`(old printed ID, new printed ID, row id)` of a saved nebula after `renumber` gave it a free ID."""
    _sector_id, nebula_id = _saved_nebula(mysql_config)
    old = pid("nebula", nebula_id, mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        before = object_uid.from_bytes(conn.execute("SELECT uid FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()["uid"])
        free = next(object_uid.with_body(before, body) for body in range(1, 4096)
                    if not designation_history.is_live_uid(conn, object_uid.with_body(before, body)))
        with conn:
            designation_history.renumber(conn, "nebulae", nebula_id, free)
        return old, object_uid.format_id(free), nebula_id
    finally:
        conn.close()


def test_a_read_of_a_retired_id_redirects_to_the_current_one(admin, mysql_config):
    old, new, _row = _renumbered(mysql_config)
    response = admin.get(f"/api/phenomena/nebula/{old}?x=1")
    assert response.status_code == 301
    assert response.headers["Location"] == f"/api/phenomena/nebula/{new}?x=1"
    assert response.headers["X-Designation-Formerly"] == old
    assert response.get_json() == {"moved_to": new, "formerly": old}
    followed = admin.get(response.headers["Location"])
    assert followed.status_code == 200 and followed.get_json()["id"] == new


def test_a_write_against_a_retired_id_answers_409_and_changes_nothing(admin, mysql_config):
    old, new, _row = _renumbered(mysql_config)
    response = admin.delete(f"/api/phenomena/nebula/{old}")
    assert response.status_code == 409
    assert response.get_json()["current_id"] == new and response.get_json()["formerly"] == old
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM nebulae").fetchone()["n"] == 1
    finally:
        conn.close()
    assert admin.delete(f"/api/phenomena/nebula/{new}").status_code == 200


def test_a_deleted_objects_id_answers_410(admin, mysql_config):
    old, _new, row = _renumbered(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE designations SET gone_epoch = 5 WHERE object_kind = 'nebulae' AND object_row_id = ?", (row,))
    finally:
        conn.close()
    for method in (admin.get, admin.delete):
        response = method(f"/api/phenomena/nebula/{old}")
        assert response.status_code == 410 and response.get_json()["gone"] is True


def test_an_id_nobody_held_is_still_a_404(admin, mysql_config):
    _renumbered(mysql_config)
    assert admin.get("/api/phenomena/nebula/FFFFFFFFFFFF-FFFFFFF-FFF").status_code == 404


def test_a_retired_id_does_not_name_another_kind(admin, mysql_config):
    old, _new, _row = _renumbered(mysql_config)
    assert admin.get(f"/api/phenomena/rogue_planet/{old}").status_code == 404


def test_a_page_asked_for_a_retired_id_redirects_to_the_current_page(admin, mysql_config):
    old, new, row = _renumbered(mysql_config)
    response = admin.get(f"/phenomenon/nebula/{old}?code=json")
    assert response.status_code == 301
    assert response.headers["Location"] == f"/phenomenon/nebula/{new}?code=json"
    assert admin.get(response.headers["Location"]).status_code == 200
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE designations SET gone_epoch = 5 WHERE object_kind = 'nebulae' AND object_row_id = ?", (row,))
    finally:
        conn.close()
    assert admin.get(f"/phenomenon/nebula/{old}").status_code == 410
