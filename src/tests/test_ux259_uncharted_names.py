# tests/test_ux259_uncharted_names.py

"""UX.259: an uncharted scatter object is shown by its local name and birth sector, with its real Object ID, never as
"Uncharted <class> star <number>". With no naming key in force the Object ID shows and the older label stays as a fallback."""

import pytest

from planetgen.names import naming_key
from tests.test_uncharted_phenomenon import _scatter, db_app  # noqa: F401
from tests.test_api import client  # noqa: F401

KEY = "00AB12CD"


@pytest.fixture
def key_in_force():
    naming_key.set_resolver(lambda: KEY)
    try:
        yield KEY
    finally:
        naming_key.set_resolver(None)


def _uid_of(mysql_config, scatter_id):
    from planetgen.db import store
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute("SELECT uid FROM phenomenon_scatter WHERE id = ?", (scatter_id,)).fetchone()["uid"]
    finally:
        conn.close()


def test_scatter_identity_is_empty_without_an_object_id():
    assert naming_key.scatter_identity(None) == {}


def test_the_phenomenon_page_carries_the_local_name_sector_and_object_id(db_app, mysql_config, key_in_force):
    scatter_id = _scatter(mysql_config, "black-hole", "stellar", mass=12.5)
    identity = naming_key.scatter_identity(_uid_of(mysql_config, scatter_id))
    html = db_app.test_client().get(f"/uncharted-phenomenon/{scatter_id}").get_data(as_text=True)
    assert identity["local_name"] in html and identity["sector_text"] in html
    assert identity["uid"] in html and "Object ID" in html
    assert "(uncharted)" not in html


def test_without_a_naming_key_the_object_id_shows_and_the_old_label_is_the_name(db_app, mysql_config):
    scatter_id = _scatter(mysql_config, "black-hole", "stellar", mass=12.5)
    html = db_app.test_client().get(f"/uncharted-phenomenon/{scatter_id}").get_data(as_text=True)
    assert "Black Hole (uncharted)" in html
    assert naming_key.scatter_identity(_uid_of(mysql_config, scatter_id))["uid"] in html
