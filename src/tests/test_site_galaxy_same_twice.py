# tests/test_site_galaxy_same_twice.py

"""
GEN.1883: the sample galaxy the site tests walk (`test_web_a11y._build_site_galaxy`, loaded by about twenty files through the dump cache of TEST.807) is the same galaxy every time it is built.

`plan` took a new galaxy seed from the operating system on each build, so two builds of one commit made two galaxies: the heaviest planet was 006 in one and 002 in the next, and
`test_irregular_moon_pages` (which reads the heaviest planet) failed by chance. The build is made twice here, in two scratch databases, and the rows of the generated tables are held
to be the same. Skipped without the MySQL test server.
"""

import hashlib
import uuid

import pymysql

from planetgen.db import store
from planetgen.db.store import MySQLConfig
from tests import dbpool
from tests.conftest import _test_server_kwargs
from tests.test_web_a11y import _build_site_galaxy

TABLES = ("galaxy_shape", "sectors", "star_systems", "stars", "planets", "moons", "irregular_sets", "irregular_moons", "asteroid_belts", "planet_features")
"""tuple: The tables of what the build generates, each read in the order of its first column."""

NOT_GENERATED = ("_at", "time", "updated", "created", "stamp")
"""tuple: Parts of a column's name that mark a time of the build (the machine's clock, not a draw)."""


def _digest(config):
    """`{table: (rows, digest of every column that is not a time)}` of the generated tables of `config`'s database."""
    conn = store.get_connection(config)
    found = {}
    try:
        for table in TABLES:
            try:
                rows = conn.execute(f"SELECT * FROM {table} ORDER BY 1").fetchall()
            except Exception:       # a table the schema does not have is not compared
                continue
            h = hashlib.sha256()
            for row in rows:
                h.update(repr(sorted((k, v) for k, v in dict(row).items() if not any(s in k for s in NOT_GENERATED))).encode())
            found[table] = (len(rows), h.hexdigest())
    finally:
        conn.close()
    return found


def test_the_site_galaxy_built_twice_holds_the_same_rows(_mysql_server_available):
    kwargs = _test_server_kwargs()
    names, digests = [], []
    try:
        for _ in range(2):
            name = f"planetgen_twice_{uuid.uuid4().hex[:12]}"
            admin = pymysql.connect(**kwargs)
            try:
                with admin.cursor() as cur:
                    cur.execute(f"CREATE DATABASE `{name}`")
                admin.commit()
            finally:
                admin.close()
            names.append(name)
            config = MySQLConfig(database=name, **kwargs)
            _build_site_galaxy(config)
            digests.append(_digest(config))
            store.close_pool(config)
    finally:
        for name in names:
            dbpool.drop_database(kwargs, name)
    first, second = digests
    assert first["planets"][0] > 0, "the sample galaxy holds no planet"
    assert first == second, "the two builds differ in " + ", ".join(t for t in first if first[t] != second.get(t))
