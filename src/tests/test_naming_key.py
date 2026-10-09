# tests/test_naming_key.py

"""
The naming key (GEN.70, `planetgen/names/naming_key.py`): drawn from the
galaxy seed, stored per galaxy database in the control database, changed
by an admin through the API, and the codec names it gives.

Tests that take `real_app` or `mysql_config` need MySQL and Redis and skip
without them.
"""

import pytest

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
