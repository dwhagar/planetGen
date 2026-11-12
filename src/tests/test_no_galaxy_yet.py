# tests/test_no_galaxy_yet.py

"""
MAP.407: a database that holds nothing but `migrate`'s bookkeeping (or nothing at all) is the galaxy's own schema before its first plan, and the
pages say there is no galaxy yet, not "no such database" or "no page at this address". A schema that holds something else (API.32) stays "no such
database", and everything that handles a `NotFoundError` still handles this one.
"""

import uuid

import pytest

from planetgen.api.config import Config
from planetgen.db import store as _db
from planetgen.web.app import create_app
from planetgen.web.lib import apiclient

from tests.test_api import client  # noqa: F401 -- fixture


@pytest.fixture
def schemas(mysql_config):
    """Throwaway schemas on the test server (prefixed so `list_databases` offers them), dropped afterwards: `make(*tables)` gives the name."""
    admin = _db.get_connection(
        _db.MySQLConfig(host=mysql_config.host, port=mysql_config.port, user=mysql_config.user, password=mysql_config.password, database=""),
        ensure_schema=False)
    made = []

    def make(*tables):
        name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
        admin.execute(f"CREATE DATABASE `{name}`")
        made.append(name)
        for table in tables:
            admin.execute(f"CREATE TABLE `{name}`.`{table}` (id INT PRIMARY KEY)")
        return name

    yield make
    for name in made:
        admin.execute(f"DROP DATABASE IF EXISTS `{name}`")
    admin.close()


@pytest.mark.parametrize("tables", [(), ("schema_migrations",), ("alembic_version",), ("schema_migrations", "alembic_version")])
def test_a_schema_with_only_the_bookkeeping_has_no_galaxy_yet(mysql_config, schemas, tables):
    name = schemas(*tables)
    with pytest.raises(_db.NoGalaxyError, match="There is no galaxy in database .* yet") as found:
        _db.resolve_database(mysql_config, name)
    assert isinstance(found.value, ValueError)


def test_a_schema_that_holds_something_else_is_still_no_such_database(mysql_config, schemas):
    name = schemas("admins")
    with pytest.raises(ValueError, match="No such database") as found:
        _db.resolve_database(mysql_config, name)
    assert not isinstance(found.value, _db.NoGalaxyError)
    name = schemas("schema_migrations", "admins")
    with pytest.raises(ValueError, match="No such database") as found:
        _db.resolve_database(mysql_config, name)
    assert not isinstance(found.value, _db.NoGalaxyError)


def test_an_unlisted_name_is_still_no_such_database(mysql_config):
    with pytest.raises(ValueError, match="No such database") as found:
        _db.resolve_database(mysql_config, "not_a_planetgen_schema")
    assert not isinstance(found.value, _db.NoGalaxyError)


def test_a_schema_with_a_galaxy_is_selected(mysql_config):
    assert _db.resolve_database(mysql_config, mysql_config.database).database == mysql_config.database


def test_the_api_answers_404_and_says_there_is_no_galaxy_yet(client, schemas):  # noqa: F811
    name = schemas("alembic_version")
    response = client.get(f"/api/sectors?db={name}")
    body = response.get_json()
    assert response.status_code == 404 and body["noGalaxy"] is True and "no galaxy" in body["error"]
    stray = schemas("admins")
    other = client.get(f"/api/sectors?db={stray}")
    assert other.status_code == 404 and "noGalaxy" not in other.get_json()


def test_a_no_galaxy_error_is_a_not_found_error():
    assert issubclass(apiclient.NoGalaxyError, apiclient.NotFoundError)


@pytest.fixture
def web(mysql_config, schemas):
    name = schemas("alembic_version")

    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        WEB_DATABASE = name
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    application = create_app(RealConfig)
    application.testing = True
    return application.test_client(), name


@pytest.mark.parametrize("path", ["/", "/sectors", "/galaxy", "/search", "/nearby"])
def test_every_page_says_there_is_no_galaxy_yet(web, path):
    client, name = web
    response = client.get(path)
    html = response.get_data(as_text=True)
    assert response.status_code == 404, path
    assert "<title>No galaxy yet" in html and "There is no galaxy in database" in html and "Generate page" in html
    assert "no page at this address" not in html.lower() and "no such database" not in html.lower()


def test_a_database_is_resolved_without_listing_every_schema(mysql_config, monkeypatch):
    """A request with `?db=` resolves its schema on every call: it asks the server about that one name, not for the size of each of its
    schemas (`list_databases` scanned `information_schema.tables` once per schema; 0.5 s a request on a server with 147 test schemas)."""
    def listed(*args, **kwargs):
        raise AssertionError("resolve_database listed every schema")
    monkeypatch.setattr(_db, "list_databases", listed)
    assert _db.resolve_database(mysql_config, mysql_config.database).database == mysql_config.database
    with pytest.raises(ValueError, match="No such database"):
        _db.resolve_database(mysql_config, "not_a_planetgen_schema")


def test_a_name_is_matched_as_written_and_inside_the_prefix(mysql_config):
    """The listing's rules hold: a case-folded name and a name outside the prefix are not databases to select."""
    with pytest.raises(ValueError, match="No such database"):
        _db.resolve_database(mysql_config, mysql_config.database.upper())
    with pytest.raises(ValueError, match="No such database"):
        _db.resolve_database(mysql_config, mysql_config.database, prefix="no_such_prefix_")
    assert _db.resolve_database(mysql_config, mysql_config.database, prefix=mysql_config.database[:8]).database == mysql_config.database


@pytest.mark.parametrize("name", ["\U0001F4A5", "caf\u00e9", "a'; DROP DATABASE x; --", "\x00", "x" * 300])
def test_a_hostile_name_is_no_such_database_and_never_a_server_error(mysql_config, name):
    """The name is compared in Python: a character the schema column's collation cannot mix with the name's was a server error (503) before."""
    with pytest.raises(ValueError, match="No such database"):
        _db.resolve_database(mysql_config, name)
