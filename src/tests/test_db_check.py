# tests/test_db_check.py

"""
The database check (DB.8, `planetgen.db.check`): a healthy saved sector
passes, and each kind of damage is found and named; nothing is written.
"""

import pytest

from planetgen.cli import generate
from planetgen.db import check, store
from tests.test_object_uids import PLACE, _sector


@pytest.fixture
def saved(mysql_config):
    store.save_sector(_sector("Check One", 3, "Check Wanderer"), config=mysql_config, galaxy_position=PLACE)
    return mysql_config


def _run(config, **kwargs):
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return check.run_checks(conn, config, **kwargs)
    finally:
        conn.close()


def _sql(config, *statements):
    conn = store.get_connection(config)
    try:
        for statement in statements:
            conn.execute(statement)
        conn.commit()
    finally:
        conn.close()


def _status(report):
    return {result.name: result.status for result in report.results}


def _problem_text(report, name):
    [result] = [r for r in report.results if r.name == name]
    return "\n".join(line for problem in result.problems for line in problem.lines())


def _heal_key(config):
    _sql(config, "UPDATE sectors SET version_key = '0008000000030F030D1000' WHERE version_key IS NULL")


def test_a_healthy_database_passes(saved):
    _heal_key(saved)
    report = _run(saved)
    assert report.exit_code() == check.EXIT_OK, "\n".join(report.lines())
    assert set(_status(report).values()) <= {check.PASS, check.WARN}
    assert report.lines()[-1] == "The database passed every check."


def test_a_row_without_its_parent_is_found(saved):
    _sql(saved, "SET foreign_key_checks = 0", "DELETE FROM star_systems WHERE id = (SELECT MIN(id) FROM (SELECT id FROM star_systems) x)")
    report = _run(saved)
    assert report.damaged and report.exit_code() == check.EXIT_DAMAGED
    assert _status(report)["orphans"] == check.FAIL
    assert "no parent in star_systems" in _problem_text(report, "orphans")


def test_a_pointer_to_a_missing_object_is_found(saved):
    _sql(saved, "INSERT INTO sector_paths (sector_id, object_table, object_id, exited, duration_years)"
                " SELECT id, 'rogue_planets', 999999, 0, 1 FROM sectors LIMIT 1")
    assert "rogue_planets id 999999" in _problem_text(_run(saved), "orphans")


def test_an_id_block_behind_the_rows_is_found(saved):
    _sql(saved, "UPDATE id_blocks SET next_id = 1 WHERE table_name = 'stars'")
    assert "id_blocks.next_id for stars" in _problem_text(_run(saved), "ids")


def test_an_impossible_value_is_found(saved):
    _sql(saved, "UPDATE planets SET radius_km = -5 WHERE id = (SELECT MIN(id) FROM (SELECT id FROM planets) x)",
         "UPDATE stars SET galactic_orbital_phase_deg = 400 WHERE id = (SELECT MIN(id) FROM (SELECT id FROM stars) x)")
    text = _problem_text(_run(saved), "values")
    assert "radius_km = -5" in text and "galactic_orbital_phase_deg = 400" in text


def test_scoping_to_another_sector_skips_the_damage(saved):
    _sql(saved, "UPDATE planets SET radius_km = -5")
    elsewhere = check.Scope(addresses=((1, 0, 0),))
    assert _status(_run(saved, scope=elsewhere))["values"] == check.PASS
    here = check.Scope(addresses=((PLACE["ring_index"], PLACE["layer_index"], PLACE["ring_slot_index"]),))
    assert _status(_run(saved, scope=here))["values"] == check.FAIL


def test_a_stale_system_count_is_a_warning_not_damage(saved):
    _heal_key(saved)
    _sql(saved, "INSERT INTO sector_stats (ring_index, layer_index, ring_slot_index, actual_systems)"
                f" VALUES ({PLACE['ring_index']}, {PLACE['layer_index']}, {PLACE['ring_slot_index']}, 99)"
                " ON DUPLICATE KEY UPDATE actual_systems = 99")
    report = _run(saved)
    assert _status(report)["counts"] == check.WARN
    assert report.exit_code() == check.EXIT_OK


def test_a_malformed_version_key_is_found(saved):
    _sql(saved, "UPDATE sectors SET version_key = 'nonsense'")
    assert "'nonsense'" in _problem_text(_run(saved), "version keys")


def test_a_schema_that_disagrees_with_the_code_is_found(saved):
    _sql(saved, "DELETE FROM schema_migrations", "ALTER TABLE planets ADD COLUMN stray_column INT")
    text = _problem_text(_run(saved), "revision")
    assert "schema_migrations has no version" in text and "stray_column" in text


def test_a_check_that_cannot_run_is_not_called_damage(saved, monkeypatch):
    def broken(ctx):
        raise KeyError("boom")
    monkeypatch.setattr(check, "CHECKS", (("ids", broken),))
    report = _run(saved)
    assert report.exit_code() == check.EXIT_UNCHECKED and not report.damaged
    assert report.results[0].status == check.ERROR


def test_the_check_writes_nothing(saved):
    _heal_key(saved)
    before = [(row["table_name"], row["table_rows"]) for row in _rows(saved, "SELECT table_name, table_rows FROM information_schema.tables WHERE table_schema = DATABASE()")]
    _run(saved)
    after = [(row["table_name"], row["table_rows"]) for row in _rows(saved, "SELECT table_name, table_rows FROM information_schema.tables WHERE table_schema = DATABASE()")]
    assert before == after


def _rows(config, sql):
    conn = store.get_connection(config)
    try:
        return conn.execute(sql).fetchall()
    finally:
        conn.close()


def test_the_command_exits_one_on_damage_and_zero_when_clean(saved, capsys, monkeypatch):
    _heal_key(saved)
    argv = ["check-db", "--mysql-host", saved.host, "--mysql-port", str(saved.port), "--mysql-user", saved.user,
            "--mysql-password", saved.password, "--mysql-database", saved.database, "--no-check-table"]
    monkeypatch.setattr("sys.argv", ["generate"] + argv)
    generate.main()
    _sql(saved, "UPDATE planets SET radius_km = -5")
    with pytest.raises(SystemExit) as raised:
        generate.main()
    assert raised.value.code == check.EXIT_DAMAGED
