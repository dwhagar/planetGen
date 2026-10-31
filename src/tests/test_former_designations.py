# tests/test_former_designations.py

"""
GEN.1140: what an object used to be called. The answers of a system, a sector and a phenomenon carry `former_designations` (newest first); the
object's page prints a "Formerly" line; a search by an old name finds the object; a reference written with a retired ID resolves to the object.
"""

import pytest

from planetgen.api import ids
from planetgen.db import designation_history, store
from planetgen.galaxy import object_uid
from tests.publicids import pid
from tests.test_admin_edits import _saved_nebula, _saved_system, admin, client  # noqa: F401 -- fixtures

pytestmark = pytest.mark.db


def _conn(config):
    return store.get_connection(config)


def _rename(config, system_id, name):
    conn = _conn(config)
    try:
        with conn:
            old = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"]
            store.rename_star_system(conn, system_id, name)
        return old
    finally:
        conn.close()


def test_a_systems_answer_lists_its_old_names_and_its_bodies(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    old = _rename(mysql_config, system_id, "Zzyzx Prime")
    detail = admin.get(f"/api/systems/{pid('system', system_id, mysql_config)}").get_json()
    assert detail["former_designations"] == [{"designation": old, "kind": "name", "reason": "rename", "retired_epoch": detail["former_designations"][0]["retired_epoch"]}]
    followed = [p for p in detail["planets"] if p["former_designations"]]
    assert followed and all(p["former_designations"][0]["reason"] == "parent change" for p in followed)
    assert all("former_designations" in row for row in detail["stars"] + detail["planets"] + detail["belts"] + detail["comets"])


def test_an_object_that_never_changed_has_an_empty_list(admin, mysql_config):
    sector_id, system_id = _saved_system(mysql_config)
    assert admin.get(f"/api/systems/{pid('system', system_id, mysql_config)}").get_json()["former_designations"] == []
    assert admin.get(f"/api/sectors/{pid('sector', sector_id, mysql_config)}").get_json()["former_designations"] == []


def test_a_phenomenon_lists_its_retired_ids_printed(admin, mysql_config):
    _sector_id, nebula_id = _saved_nebula(mysql_config)
    old = pid("nebula", nebula_id, mysql_config)
    conn = _conn(mysql_config)
    try:
        before = object_uid.from_bytes(conn.execute("SELECT uid FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()["uid"])
        with conn:
            designation_history.renumber(conn, "nebulae", nebula_id, object_uid.with_body(before, 77))
        new = object_uid.format_id(object_uid.with_body(before, 77))
    finally:
        conn.close()
    detail = admin.get(f"/api/phenomena/nebula/{new}").get_json()
    assert [(line["kind"], line["designation"], line["reason"]) for line in detail["former_designations"]] == [("uid", old, "admin renumber")]
    # a reference written with the retired ID names the object
    conn = _conn(mysql_config)
    try:
        assert ids.resolve_ref(conn, f"nebula:{old}") == ids.resolve_ref(conn, f"nebula:{new}") is not None
    finally:
        conn.close()


def test_the_object_pages_print_a_formerly_line(admin, mysql_config):
    sector_id, system_id = _saved_system(mysql_config)
    old = _rename(mysql_config, system_id, "Zzyzx Prime")
    page = admin.get(f"/system/{pid('system', system_id, mysql_config)}")
    assert page.status_code == 200
    text = page.get_data(as_text=True)
    assert 'class="formerly"' in text and old in text and "rename, epoch" in text
    plain = admin.get(f"/sector/{pid('sector', sector_id, mysql_config)}").get_data(as_text=True)
    assert 'class="formerly"' not in plain


def test_a_search_by_an_old_name_finds_the_object(admin, mysql_config):
    _sector_id, system_id = _saved_system(mysql_config)
    old = _rename(mysql_config, system_id, "Zzyzx Prime")
    word = old.split()[0]
    rows = admin.get(f"/api/search?system_q={word}&panels=systems").get_json()["results"]["systems"]["rows"]
    assert [row["name"] for row in rows if row["name"] == "Zzyzx Prime"] == ["Zzyzx Prime"]
    assert admin.get("/api/search?system_q=Zzyzx&panels=systems").get_json()["results"]["systems"]["total"] == 1
