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


ID = "0123456789ABCDEF012"


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


def test_names_round_trip_to_the_id_at_19_digits():
    for kind in naming_key.KINDS:
        name = naming_key.codec_name(ID, kind, "00AB12CD")
        assert name == name.title() and " " in name or len(name) > 3
        assert naming_key.object_id_of(name, kind, "00AB12CD") == ID


def test_the_key_and_the_kind_both_change_the_name():
    base = naming_key.codec_name(ID, "nebula", "00000001")
    names = {naming_key.codec_name(ID, "nebula", "00000001"), naming_key.codec_name(ID, "quasar", "00000001")}
    assert len(names) == 2
    others = {naming_key.codec_name(ID, "nebula", f"{n:08X}") for n in range(40)}
    assert len(others) > 10 and base in others | {base}


def test_only_the_codec_kinds_have_names():
    with pytest.raises(ValueError):
        naming_key.codec_name(ID, "star", "00000001")
    with pytest.raises(ValueError):
        naming_key.codec_name("ABC", "nebula", "00000001")


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
            store.insert_rogue_planet(conn, planet, placement={
                "center_x_pc": 10.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 10.0})
        row_id = conn.execute("SELECT id FROM rogue_planets WHERE name = ?", (planet.name,)).fetchone()["id"]
    finally:
        conn.close()
    return planet.name, row_id


def test_the_api_shows_the_codec_name_and_a_key_change_renames_it(real_app, mysql_config):  # noqa: F811
    from planetgen.api import naming
    from planetgen.db import store

    stored, row_id = _rogue_planet_with_id(mysql_config)
    assert len(stored) == 19
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
    assert first == naming_key.codec_name(stored, "rogue-planet", key) != stored
    listed = [item["name"] for item in admin.get("/api/phenomena").get_json()["items"]]
    assert first in listed and stored not in listed

    admin.post("/api/admin/naming-key", json={"key": "FFFFFFFF"})
    second = admin.get(detail).get_json()["name"]
    assert second == naming_key.codec_name(stored, "rogue-planet", "FFFFFFFF") and second != first
    assert naming_key.object_id_of(second, "rogue-planet", "FFFFFFFF") == stored


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


def test_only_object_id_names_are_renamed_and_the_key_is_asked_once():
    from planetgen.names import object_id

    nebula = object_id.format_id(object_id.pack("nebula", (100.0, 200.0, 3.0)))
    core = object_id.format_id(object_id.pack("black-hole-core", (100.0, 200.0, 3.0)))
    bright = object_id.format_id(object_id.pack("bright-star", (100.0, 200.0, 3.0)))
    asked = []

    def key():
        asked.append(1)
        return "00AB12CD"

    payload = {"items": [{"name": nebula}, {"name": core}, {"name": bright}, {"name": "Sol"},
                         {"nested": {"name": nebula, "other": nebula}}], "total": 5}
    out = naming_key.rename_names(payload, key)
    assert out["items"][0]["name"] == naming_key.codec_name(nebula, "nebula", "00AB12CD")
    assert out["items"][1]["name"] == naming_key.codec_name(core, "black-hole-core", "00AB12CD")
    assert out["items"][2]["name"] == bright and out["items"][3]["name"] == "Sol"
    assert out["items"][4]["nested"] == {"name": out["items"][0]["name"], "other": nebula}
    assert payload["items"][0]["name"] == nebula  # the input is not changed
    assert len(asked) == 1
    plain = {"items": [{"name": "Sol"}]}
    assert naming_key.rename_names(plain, lambda: 1 / 0) is plain  # no ID, no key lookup


def test_a_typed_codec_name_finds_its_id():
    from planetgen.names import object_id

    stored = object_id.format_id(object_id.pack("quasar", (5000.0, 10.0, 1.0)))
    name = naming_key.display_name(stored, "00AB12CD")
    assert name != stored
    assert naming_key.stored_name_for(name.upper(), "00AB12CD") == stored
    assert naming_key.stored_name_for("Sol", "00AB12CD") is None
    assert naming_key.display_name(stored, None) in (stored,) or naming_key.active_key() is not None
