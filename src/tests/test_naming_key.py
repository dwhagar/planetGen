# tests/test_naming_key.py

"""
The naming key (GEN.70, `planetgen/names/naming_key.py`): drawn from the
galaxy seed, stored per galaxy database in the control database, changed
by an admin through the API, and the codec names it gives.

Tests that take `real_app` or `mysql_config` need MySQL and Redis and skip
without them.
"""

import pytest

from tests.publicids import pid

from planetgen.names import naming_key
from tests.test_login_backoff import _admin_client, real_app  # noqa: F401 -- the fixture




def test_the_key_is_drawn_from_the_seed_and_is_eight_hex_digits():
    seed = bytes(range(16))
    assert naming_key.draw_key(seed) == naming_key.draw_key(seed)
    assert naming_key.draw_key(seed) != naming_key.draw_key(bytes(16))
    assert len(naming_key.draw_key(seed)) == naming_key.KEY_DIGITS
    assert naming_key.parse_key(" 0a1b2c3d ") == "0A1B2C3D"


@pytest.mark.parametrize("text", ["", "1234567", "123456789", "0A1B2C3G", "0x1B2C3D"])
def test_a_key_must_be_eight_hex_digits(text):
    with pytest.raises(ValueError):
        naming_key.parse_key(text)


def test_only_the_region_and_constellation_kinds_have_a_domain():
    assert naming_key.codec_domain("00AB12CD", "region") == "00AB12CD:region"
    with pytest.raises(ValueError):
        naming_key.codec_domain("00AB12CD", "nebula")


def test_a_key_is_stored_kept_and_changed(real_app):  # noqa: F811
    from planetgen.db import store

    conn = store.get_control_connection(real_app.config["CONTROL_MYSQL_CONFIG"])
    try:
        assert naming_key.get(conn, "galaxy_a") is None
        seed = bytes(range(16))
        key = naming_key.draw(conn, "galaxy_a", seed)
        assert key == naming_key.draw_key(seed)
        assert naming_key.change(conn, "galaxy_a", "deadbeef", "admin") == "DEADBEEF"
        # Planning again over the same seed keeps the admin's key; a new seed draws a new one.
        assert naming_key.draw(conn, "galaxy_a", seed) == "DEADBEEF"
        assert naming_key.draw(conn, "galaxy_a", bytes(16), replace=True) == naming_key.draw_key(bytes(16))
        assert naming_key.get(conn, "galaxy_a")["changed_by"] is None
        with pytest.raises(LookupError):
            naming_key.change(conn, "never_planned", "00000000", "admin")
    finally:
        conn.close()


def test_an_admin_reads_and_changes_the_key_through_the_api(real_app):  # noqa: F811
    from planetgen.db import store

    admin = _admin_client(real_app)
    database = real_app.config["MYSQL_CONFIG"].database
    body = admin.get("/api/admin/naming-key").get_json()
    assert body["key"] is None and body["database"] == database
    assert admin.post("/api/admin/naming-key", json={"key": "00112233"}).status_code == 409

    conn = store.get_control_connection(real_app.config["CONTROL_MYSQL_CONFIG"])
    try:
        naming_key.draw(conn, database, bytes(16))
    finally:
        conn.close()
    assert admin.post("/api/admin/naming-key", json={"key": "nope"}).status_code == 400
    changed = admin.post("/api/admin/naming-key", json={"key": "a1b2c3d4"}).get_json()
    assert changed["key"] == "A1B2C3D4" and changed["changed_by"] == "admin"
    drawn = admin.post("/api/admin/naming-key", json={"draw": True}).get_json()
    assert len(drawn["key"]) == 8 and drawn["key"] != "A1B2C3D4"
    assert admin.get("/api/admin/naming-key").get_json()["key"] == drawn["key"]


def test_planning_draws_a_key_for_a_new_seed_and_keeps_an_admins_over_the_same_seed(real_app):  # noqa: F811
    from planetgen.db import store
    from planetgen.generation import run_plan

    from planetgen.admin import auth as adminAuth

    config = real_app.config["MYSQL_CONFIG"]
    control = store.control_mysql_config(config)  # where planning looks for it
    adminAuth.bootstrap_control_schema(control)
    seed = bytes(range(16))
    run_plan._draw_naming_key(config, seed, new_seed=True)
    conn = store.get_control_connection(control)
    try:
        assert naming_key.key_of(conn, config.database) == naming_key.draw_key(seed)
        naming_key.change(conn, config.database, "00000042", "admin")
        run_plan._draw_naming_key(config, seed, new_seed=False)
        conn.commit()  # a new snapshot: planning wrote on its own connection
        assert naming_key.key_of(conn, config.database) == "00000042"
        other = bytes(16)
        run_plan._draw_naming_key(config, other, new_seed=True)
        conn.commit()
        assert naming_key.key_of(conn, config.database) == naming_key.draw_key(other)
    finally:
        conn.close()


def _rogue_planet_with_id(mysql_config):
    from planetgen.db import store
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.phenomena.rogue import RoguePlanet

    planet = RoguePlanet(SystemConfig())
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            row_id = store.insert_rogue_planet(conn, planet, placement={
                "center_x_pc": 10.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 10.0})
            store.assign_uids(conn, phenomenon=("rogue_planets", row_id))
        stored = conn.execute("SELECT name FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()["name"]
    finally:
        conn.close()
    return stored, row_id


def test_the_api_shows_the_codec_name_and_a_key_change_renames_it(real_app, mysql_config):  # noqa: F811
    from planetgen.api import naming
    from planetgen.db import store

    stored, row_id = _rogue_planet_with_id(mysql_config)
    uid = naming_key.stored_uid_of(stored)
    assert uid is not None
    admin = _admin_client(real_app)
    detail = f"/api/phenomena/rogue_planet/{pid('rogue_planet', row_id)}"
    # No key drawn yet: the ID is the name.
    assert admin.get(detail).get_json()["name"] == stored

    conn = store.get_control_connection(real_app.config["CONTROL_MYSQL_CONFIG"])
    try:
        naming_key.draw(conn, real_app.config["MYSQL_CONFIG"].database, bytes(16))
    finally:
        conn.close()
    naming.forget()
    key = naming_key.draw_key(bytes(16))
    first = admin.get(detail).get_json()["name"]
    assert first == naming_key.local_name(uid, key) != stored
    listed = [item["name"] for item in admin.get("/api/phenomena").get_json()["items"]]
    assert first in listed and stored not in listed

    admin.post("/api/admin/naming-key", json={"key": "FFFFFFFF"})
    second = admin.get(detail).get_json()["name"]
    assert second == naming_key.local_name(uid, "FFFFFFFF") and second != first
    shown = admin.get(detail).get_json()
    assert shown["name_sector"] == naming_key.sector_text(uid)
    assert naming_key.uid_of_local_name(second, shown["name_sector"], "FFFFFFFF") == uid


def test_a_key_change_changes_the_tile_stamp(real_app, mysql_config):  # noqa: F811
    from planetgen.api import naming
    from planetgen.db import store

    admin = _admin_client(real_app)
    conn = store.get_control_connection(real_app.config["CONTROL_MYSQL_CONFIG"])
    try:
        naming_key.draw(conn, real_app.config["MYSQL_CONFIG"].database, bytes(16))
    finally:
        conn.close()
    naming.forget()
    before = admin.get("/api/galaxy/stamp").get_json()["stamp"]
    admin.post("/api/admin/naming-key", json={"key": "FFFFFFFF"})
    assert admin.get("/api/galaxy/stamp").get_json()["stamp"] != before


def _printed(serial, body=0):
    from planetgen.galaxy import object_uid

    return object_uid.format_id(object_uid.pack(3, 4, 5, object_uid.SERIAL_GENERATED, serial, body))


def test_only_object_id_names_are_renamed_and_the_key_is_asked_once():
    from planetgen.galaxy import object_uid

    nebula, core = _printed(7), _printed(8)
    asked = []

    def key():
        asked.append(1)
        return "00AB12CD"

    payload = {"items": [{"name": nebula}, {"name": core}, {"name": "Sol"}, {"name": "0123456789ABCDEF012"},
                         {"nested": {"name": nebula, "other": nebula}}], "total": 5}
    out = naming_key.rename_names(payload, key)
    first = naming_key.local_name(object_uid.parse_id(nebula), "00AB12CD")
    assert out["items"][0] == {"name": first, "name_sector": naming_key.sector_text(object_uid.parse_id(nebula))}
    assert out["items"][1]["name"] == naming_key.local_name(object_uid.parse_id(core), "00AB12CD")
    assert out["items"][2] == {"name": "Sol"} and out["items"][3] == {"name": "0123456789ABCDEF012"}
    assert out["items"][4]["nested"]["name"] == first and out["items"][4]["nested"]["other"] == nebula
    assert payload["items"][0]["name"] == nebula  # the input is not changed
    assert len(asked) == 1
    plain = {"items": [{"name": "Sol"}]}
    assert naming_key.rename_names(plain, lambda: 1 / 0) is plain  # no ID, no key lookup
    assert naming_key.rename_names({"name": nebula}, lambda: None) == {"name": nebula}  # no key drawn: the ID stays


def test_a_typed_local_name_with_its_sector_finds_its_id():
    from planetgen.galaxy import object_uid

    stored = _printed(11, 2)
    uid = object_uid.parse_id(stored)
    typed = naming_key.local_name_with_sector(uid, "00AB12CD")
    assert naming_key.uid_of_typed_name(typed.upper(), "00AB12CD") == uid
    assert naming_key.uid_of_typed_name(typed.replace(", sector", ","), "00AB12CD") == uid
    assert naming_key.uid_of_typed_name("Sol", "00AB12CD") is None
    assert naming_key.uid_of_typed_name("Sol, sector 12", "00AB12CD") is None
    assert naming_key.display_name(stored, "00AB12CD") == naming_key.local_name(uid, "00AB12CD")
    assert naming_key.display_name("Sol", "00AB12CD") == "Sol"
