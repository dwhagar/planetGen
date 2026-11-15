# tests/test_check_pieces.py

"""
OPS.204: the checks run on the engine's loop. Each check lists its pieces (a table, a key, a sector, a star system, a check, a dataset) and the unit's bar counts them; a check that
cannot list its pieces is "could not check"; the math check and the data check give the same text as before.
"""

import pytest

from planetgen.db import check, store
from planetgen.cli import datasets
from planetgen.physics import mathcheck
from planetgen.queue import commands, wrapper
from tests.test_db_check import _heal_key, _sql, saved  # noqa: F401 -- `saved` is a fixture


def _args(config, **more):
    return {"database": wrapper.database_args(config), "shells": [], "addresses": [], "wide_binaries": False, **more}


def _counted(kind, args):
    reports = []
    result = commands.run_alone(kind, args, lambda done, total: reports.append((done, total)))
    return result, reports


@pytest.mark.db
def test_a_table_check_counts_its_tables(saved):
    result, reports = _counted("check-tables", _args(saved))
    conn = store.get_connection(saved, ensure_schema=False)
    try:
        tables = len(check._Context(conn, saved, check.Scope()).tables())
    finally:
        conn.close()
    assert reports[0] == (0, tables) and reports[-1] == (tables, tables) and len(reports) == tables + 1
    assert check.CheckResult.from_dict(result).status == check.PASS


@pytest.mark.db
def test_the_systems_check_counts_systems_and_names_the_one_that_fails(saved):
    _heal_key(saved)
    _sql(saved, "UPDATE stars SET mass_kg = -1", "UPDATE planets SET distance_km = 1e30")
    result, reports = _counted("check-systems", _args(saved))
    systems = reports[-1][1]
    assert systems >= 1 and reports[-1] == (systems, systems) and len(reports) == systems + 1
    found = check.CheckResult.from_dict(result)
    assert found.status == check.FAIL and any("star system" in line for problem in found.problems for line in problem.lines())


@pytest.mark.db
def test_the_storage_check_counts_sectors_and_names_a_damaged_one(saved):
    from planetgen.db import storage_checksum
    conn = store.get_connection(saved)
    try:
        for row in conn.execute("SELECT id FROM sectors").fetchall():
            storage_checksum.seal(conn, row["id"])
        conn.commit()
    finally:
        conn.close()
    result, reports = _counted("check-storage", _args(saved))
    sectors = reports[-1][1]
    assert sectors >= 1 and reports[-1] == (sectors, sectors) and check.CheckResult.from_dict(result).status == check.PASS
    _sql(saved, "UPDATE stars SET temperature_k = temperature_k + 1")        # written behind the seal's back
    found = check.CheckResult.from_dict(commands.run_alone("check-storage", _args(saved)))
    assert found.status == check.FAIL and any("sector id" in line for problem in found.problems for line in problem.lines())


@pytest.mark.db
def test_a_check_that_cannot_list_its_pieces_is_could_not_check(saved, monkeypatch):
    def broken(ctx):
        raise KeyError("boom")
    monkeypatch.setitem(check.BY_NAME, "ids", check.Check("ids", "ids", broken, broken))
    result, reports = _counted("check-ids", _args(saved))
    found = check.CheckResult.from_dict(result)
    assert found.status == check.ERROR and "boom" in found.error and reports[-1] == (1, 1)


@pytest.mark.db
def test_a_check_that_fails_in_one_piece_is_could_not_check_not_damage(saved, monkeypatch):
    def listing(ctx):
        return ["a", "b"]

    def piece(ctx, name):
        if name == "b":
            raise KeyError("no such table")
        return []
    monkeypatch.setitem(check.BY_NAME, "ids", check.Check("ids", "ids", listing, piece))
    found = check.CheckResult.from_dict(commands.run_alone("check-ids", _args(saved)))
    assert found.status == check.ERROR and "no such table" in found.error


def test_the_math_check_counts_its_checks_and_gives_the_report_it_always_gave():
    result, reports = _counted("check-math", {})
    total = len(mathcheck.all_checks())
    assert reports[0] == (0, total) and reports[-1] == (total, total)
    assert result["failed"] == 0 and result["report"].splitlines()[-1].startswith(f"Math check passed: all {total} checks in ")
    verbose = commands.run_alone("check-math", {"verbose": True})
    assert len(verbose["report"].splitlines()) == total + 1


def test_the_data_check_counts_its_datasets_and_says_each_one():
    keys = [dataset.key for dataset in datasets.catalogue.ALL][:2]
    result, reports = _counted("check-datasets", {"keys": keys})
    assert reports[0] == (0, len(keys)) and reports[-1] == (len(keys), len(keys))
    assert [line.split(":")[0] for line in result["lines"] if not line.startswith(" ")][:len(keys)] == keys
    lines, absent = datasets.status_lines(datasets._chosen(keys), False)
    assert result["lines"] == lines and result["absent"] == absent
