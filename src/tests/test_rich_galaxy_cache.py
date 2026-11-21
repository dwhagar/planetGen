# tests/test_rich_galaxy_cache.py

"""
TEST.788: the rich test galaxy is built once for every process that tests the same code, and the rest load its dump (`tests/rich_galaxy_cache.py`).
"""

import os
import time

import pymysql
import pytest

from planetgen.db import store
from tests import rich_galaxy_cache as cache
from tests.conftest import _test_server_kwargs
from tests.db_schema_support import scratch_database


def _server():
    return store.MySQLConfig(database="", **_test_server_kwargs())


def _fill(config):
    """Every kind of value the galaxy's tables hold that a dump could lose: doubles to the last bit, a JSON document, text with quotes, a blob, NULLs, a trigger, a foreign key."""
    connection = pymysql.connect(host=config.host, port=config.port, user=config.user, password=config.password, database=config.database)
    try:
        with connection.cursor() as cursor:
            cursor.execute("CREATE TABLE parent (id INT AUTO_INCREMENT PRIMARY KEY, name VARCHAR(40) NOT NULL) ENGINE=InnoDB")
            cursor.execute("CREATE TABLE child (id INT AUTO_INCREMENT PRIMARY KEY, parent_id INT NOT NULL, x DOUBLE, y DOUBLE NULL, f FLOAT, doc JSON, note TEXT, raw BLOB, "
                           "seen INT DEFAULT 0, FOREIGN KEY (parent_id) REFERENCES parent(id)) ENGINE=InnoDB")
            cursor.execute("CREATE TRIGGER bump BEFORE UPDATE ON child FOR EACH ROW SET NEW.seen = OLD.seen + 1")
            cursor.executemany("INSERT INTO parent (name) VALUES (%s)", [("Sol",), ("Alpha 'Centauri\"",)])
            rows = [(1, 0.1 + 0.2, None, 1.1, '{"a": [1, 2.5, null]}', "it's \"quoted\"\nnew line ☃", b"\x00\xff\x10", 0),
                    (2, 5e-324, 1.7976931348623157e308, 3.4e38, '{}', None, b"", 0),
                    (2, -0.0, 123456789.12345678, 0.0, '[1, 2]', "", None, 3)]
            cursor.executemany("INSERT INTO child (parent_id, x, y, f, doc, note, raw, seen) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)", rows)
        connection.commit()
    finally:
        connection.close()


def _checksums(config):
    connection = pymysql.connect(host=config.host, port=config.port, user=config.user, password=config.password, database=config.database)
    try:
        with connection.cursor() as cursor:
            out = {}
            for table in ("parent", "child"):
                cursor.execute(f"CHECKSUM TABLE `{table}`")
                out[table] = cursor.fetchone()[1]
            cursor.execute("SELECT x, y, f, note, raw FROM child ORDER BY id")
            out["rows"] = cursor.fetchall()
            cursor.execute("SHOW CREATE TABLE child")
            out["create"] = cursor.fetchone()[1].split(" AUTO_INCREMENT=")[0]
            cursor.execute("SELECT TRIGGER_NAME FROM information_schema.TRIGGERS WHERE TRIGGER_SCHEMA = DATABASE()")
            out["triggers"] = cursor.fetchall()
            return out
    finally:
        connection.close()


def test_a_dump_loaded_into_a_new_database_is_the_same_rows_to_the_last_bit(tmp_path, _mysql_server_available):
    with scratch_database(_server()) as source, scratch_database(_server()) as target:
        _fill(source)
        path = str(tmp_path / "galaxy.sql.gz")
        cache.dump(source, path)
        assert os.path.getsize(path) > 0 and not [name for name in os.listdir(tmp_path) if ".part" in name]
        cache.load(target, path)
        assert _checksums(source) == _checksums(target)
        assert [repr(row[0]) for row in _checksums(target)["rows"]] == ["0.30000000000000004", "5e-324", "-0.0"]         # the minus zero is a minus zero


def test_a_minus_zero_in_a_dump_stays_a_minus_zero_and_a_string_that_looks_like_one_is_not_touched():
    line = b"INSERT INTO `t` VALUES (1,-0,'a,-0,b','it\\'s -0',-0),(2,-0.5,-0,NULL);\n"
    assert cache.keep_negative_zero(line) == b"INSERT INTO `t` VALUES (1,-0e0,'a,-0,b','it\\'s -0',-0e0),(2,-0.5,-0e0,NULL);\n"
    other = b"CREATE TABLE `t` (`a` int DEFAULT -0);\n"
    assert cache.keep_negative_zero(other) == other


def test_the_first_process_builds_and_keeps_and_the_next_loads_without_building(tmp_path, _mysql_server_available):
    built = []
    with scratch_database(_server()) as first, scratch_database(_server()) as second, scratch_database(_server()) as third:
        assert cache.build_or_load(first, lambda: (built.append("first"), _fill(first)), str(tmp_path), "key1") == "built and kept"
        assert cache.build_or_load(second, lambda: built.append("second"), str(tmp_path), "key1") == "loaded"
        assert built == ["first"] and _checksums(first) == _checksums(second)
        assert cache.build_or_load(third, lambda: (built.append("third"), _fill(third)), str(tmp_path), "key2") == "built and kept"          # another key is another galaxy
        assert built == ["first", "third"]
    assert sorted(name for name in os.listdir(tmp_path)) == ["rich-key1.sql.gz", "rich-key2.sql.gz"]            # no lock or temporary file is left


def test_a_process_that_finds_the_build_taken_waits_for_its_dump_and_loads_it(tmp_path, _mysql_server_available):
    with scratch_database(_server()) as maker, scratch_database(_server()) as waiter:
        _fill(maker)
        os.mkdir(tmp_path / "rich-key.lock")                    # another process is building
        slept = []

        def sleep(seconds):
            slept.append(seconds)
            cache.dump(maker, str(tmp_path / "rich-key.sql.gz"))        # it finishes while this one waits
        result = cache.build_or_load(waiter, lambda: pytest.fail("built what another process was building"), str(tmp_path), "key", sleep=sleep)
        assert result == "loaded" and slept == [cache.POLL_SECONDS] and _checksums(maker) == _checksums(waiter)


def test_a_lock_left_by_a_dead_build_is_taken_over_and_a_wait_that_ends_builds_without_caching(tmp_path, _mysql_server_available):
    lock = tmp_path / "rich-old.lock"
    os.mkdir(lock)
    old = time.time() - cache.LOCK_STALE_SECONDS - 60
    os.utime(lock, (old, old))
    with scratch_database(_server()) as one, scratch_database(_server()) as two:
        assert cache.build_or_load(one, lambda: _fill(one), str(tmp_path), "old") == "built and kept"
        os.mkdir(tmp_path / "rich-busy.lock")                   # a live build that never finishes
        clock = iter(range(0, 10 ** 6, 100))
        result = cache.build_or_load(two, lambda: _fill(two), str(tmp_path), "busy", sleep=lambda seconds: None, now=lambda: next(clock), wait_seconds=250)
        assert result == "built" and not (tmp_path / "rich-busy.sql.gz").exists()


def test_the_key_changes_with_a_file_of_the_package_or_the_data_and_not_otherwise(tmp_path, _mysql_server_available):
    package, tests, data = tmp_path / "planetgen", tmp_path / "tests", tmp_path / "data"
    for folder in (package, tests, data, package / "__pycache__"):
        folder.mkdir(exist_ok=True)
    (package / "a.py").write_text("x = 1\n")
    (data / "table.csv").write_text("1,2\n")
    for name in cache.SUPPORT_FILES:
        (tests / name).write_text(name)
    config = _server()

    def key():
        return cache.key_of(config, str(package), str(tests), str(data))
    first = key()
    (package / "__pycache__" / "a.cpython-313.pyc").write_bytes(b"compiled")           # compiled files and caches are not the code
    assert key() == first
    (package / "a.py").write_text("x = 2\n")
    changed = key()
    assert changed != first
    (data / "table.csv").write_text("1,2,3\n")
    assert key() not in (first, changed)
    (tests / "rich_galaxy_support.py").write_text("changed")
    third = key()
    assert len({first, changed, third}) == 3
    system_file = tmp_path / "arrakis.json"
    system_file.write_text("{}")
    with_file = cache.key_of(config, str(package), str(tests), str(data), extra_files=(str(system_file), str(tmp_path / "config.json")), environ={})        # a file that is not there is left out
    system_file.write_text('{"a": 1}')
    assert cache.key_of(config, str(package), str(tests), str(data), extra_files=(str(system_file),), environ={}) != with_file
    plain = cache.key_of(config, str(package), str(tests), str(data), environ={"PLANETGEN_MYSQL_HOST": "one", "PLANETGEN_REDIS_URL": "x", "PLANETGEN_TEST_MYSQL_PORT": "1"})
    assert plain == cache.key_of(config, str(package), str(tests), str(data), environ={"PLANETGEN_MYSQL_HOST": "two", "PLANETGEN_REDIS_URL": "y"})      # where the servers are is not what the build makes
    assert plain != cache.key_of(config, str(package), str(tests), str(data), environ={"PLANETGEN_SOME_SETTING": "on"})                                    # a setting the generator reads is


def test_the_cache_is_on_by_default_in_the_machines_folder_and_off_only_when_set_empty():
    """TEST.807: unset keeps the dumps in the machine's folder; a folder shares them; set but empty switches the cache off."""
    assert cache.cache_folder({}) == cache.default_folder()
    assert cache.cache_folder({cache.ENVIRONMENT_VARIABLE: "/shared/dumps"}) == "/shared/dumps"
    assert cache.cache_folder({cache.ENVIRONMENT_VARIABLE: "  "}) is None


def test_prune_removes_old_dumps_and_dead_locks_and_keeps_new_ones_and_other_files(tmp_path):
    day = 86400
    old, fresh, other = tmp_path / "rich-a.sql.gz", tmp_path / "rich-b.sql.gz", tmp_path / "notes.txt"
    lock = tmp_path / "rich-c.lock"
    for path in (old, fresh, other):
        path.write_text("x")
    lock.mkdir()
    now = old.stat().st_mtime + 10 * day
    os.utime(fresh, (now, now))
    os.utime(other, (now - 10 * day,) * 2)
    cache.prune(str(tmp_path), now=lambda: now)
    assert not old.exists() and not lock.exists()
    assert fresh.exists() and other.exists()
