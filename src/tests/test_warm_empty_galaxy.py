# tests/test_warm_empty_galaxy.py

"""
A warm-up or a whole-map build of a database that holds no galaxy, e.g. right after a New galaxy wipe, does nothing and succeeds
(OPS.114): one information line, exit status 0, and update.sh says the warm-up was skipped, not that it failed. It runs normally once a
galaxy exists.
"""

import sys
from types import SimpleNamespace

import pytest

from planetgen.cli import tiles as tiles_cli
from planetgen.cli import warm as warm_cli
from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector


def _place_a_sector(mysql_config):
    position = {"center_x_pc": 5.0, "center_y_pc": 5.0, "center_z_pc": 5.0, "galactic_radius_pc": 8.7,
                "shell_index": 1, "band_index": 0, "wedge_index": 0}
    store.save_sector(SpaceSector("Placed Sector", edge_ly=10.0), config=mysql_config, galaxy_position=position)


def test_a_wiped_database_is_empty_and_one_with_a_sector_is_not(mysql_config):
    from planetgen.cli import reset as reset_cli
    conn = store.get_connection(mysql_config)
    try:
        assert store.galaxy_is_empty(conn)
    finally:
        conn.close()
    _place_a_sector(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        assert not store.galaxy_is_empty(conn)
    finally:
        conn.close()
    reset_cli.wipe(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        assert store.galaxy_is_empty(conn)
    finally:
        conn.close()


def test_a_missing_table_counts_as_empty_and_other_errors_are_raised():
    import pymysql

    class Conn:
        def __init__(self, error):
            self.error = error

        def execute(self, _sql):
            raise self.error

    assert store.galaxy_is_empty(Conn(pymysql.err.ProgrammingError(1146, "no such table")))
    with pytest.raises(pymysql.err.ProgrammingError):
        store.galaxy_is_empty(Conn(pymysql.err.ProgrammingError(1064, "syntax")))


def _args(database="galaxy"):
    return SimpleNamespace(database=database, on_change=None, follow=False, families=[])


def test_an_empty_galaxy_is_not_warmed_and_the_run_succeeds(monkeypatch, tmp_path, capsys):
    outcome = tmp_path / "outcome"
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    monkeypatch.setenv(warm_cli.OUTCOME_FILE_ENV, str(outcome))
    monkeypatch.setattr(warm_cli, "galaxy_is_empty", lambda database: True)
    monkeypatch.setattr(warm_cli, "_warm_run", lambda *a, **k: pytest.fail("nothing to warm"))
    warm_cli._warm_once(_args())      # returns, raises nothing: exit 0
    seen = capsys.readouterr()
    assert warm_cli.EMPTY_GALAXY_LINE in seen.out + seen.err
    assert "error" not in (seen.out + seen.err).lower()
    assert outcome.read_text().strip() == "empty"


def test_a_galaxy_that_exists_is_warmed_as_before(monkeypatch, tmp_path):
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path))
    monkeypatch.delenv(warm_cli.OUTCOME_FILE_ENV, raising=False)
    ran = []
    monkeypatch.setattr(warm_cli, "galaxy_is_empty", lambda database: False)
    monkeypatch.setattr(warm_cli, "_warm_run", lambda args, reason, record: ran.append(args.database))
    warm_cli._warm_once(_args())
    assert ran == ["galaxy"]


def test_the_whole_map_build_of_an_empty_galaxy_does_nothing_and_succeeds(monkeypatch, capsys):
    monkeypatch.setattr(warm_cli, "galaxy_is_empty", lambda database: True)
    from planetgen.web import tilebuild
    monkeypatch.setattr(tilebuild, "build_all", lambda *a, **k: pytest.fail("nothing to build"))
    status = tiles_cli._build(SimpleNamespace(workers=0, force=False, no_swap=False),
                              SimpleNamespace(database="galaxy"))
    assert status == 0
    seen = capsys.readouterr()
    assert warm_cli.EMPTY_GALAXY_LINE in seen.out + seen.err


def test_update_sh_reports_a_skipped_warm_up_not_a_failure():
    text = open("scripts/deploy-common.sh", encoding="utf-8").read()
    assert "PLANETGEN_WARM_OUTCOME_FILE=$outcome_file" in text
    assert "the galaxy is empty" in text and "skipped, the galaxy is empty" in text
