"""OPS.195 and OPS.200: a long statement of a command is named in the beat (`planetgen/queue/longsql.py`, `db/statementmark.py`)."""

import threading
import time

import pytest

from planetgen.db import statementmark, store
from planetgen.generation import run_plan
from planetgen.queue import longsql, unitbeat


def test_longsql_names_the_statement_while_it_runs_and_not_after(mysql_config):
    conn = store.get_connection(mysql_config)
    seen = []

    def watch():
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline and not seen:
            mark = statementmark.running()
            if mark is not None:
                seen.append(mark)
            time.sleep(0.02)

    watcher = threading.Thread(target=watch)
    try:
        watcher.start()
        longsql.run(conn, "SELECT SLEEP(0.6)", "Sleeping a moment")
        watcher.join()
        assert seen and seen[0][0] == "Sleeping a moment", "the statement is named for the beat while it runs"
        assert statementmark.running() is None, "and not after"
        assert longsql.stage_events_readable(mysql_config) in (True, False)
    finally:
        conn.close()


def test_the_beat_carries_the_work_done_and_estimated_where_the_server_gives_it():
    """Where MySQL's stage events can be read a statement's `progress` gives `(done, estimated)`; a reading that fails or gives nothing adds nothing and never stops the beat."""
    with statementmark.labelled("Building an index", lambda: (30, 120)), statementmark.statement("ALTER TABLE x ADD INDEX i (a)"):
        found = unitbeat.evidence(100.0)
    assert found["sql"] == "Building an index" and (found["sql_done"], found["sql_estimated"]) == (30, 120)
    for progress in (lambda: None, lambda: 1 / 0):
        with statementmark.labelled("Building an index", progress), statementmark.statement("ALTER TABLE x ADD INDEX i (a)"):
            found = unitbeat.evidence(100.0)
        assert "sql_done" not in found and found["sql"] == "Building an index"


class _Recorder:
    """Steps that record what each step is labelled, and the label the statement saw while it ran."""

    def __init__(self):
        self.labels = []

    def sql(self, label):
        self.labels.append(label)
        return statementmark.labelled(label)

    def nothing(self, text):
        pytest.fail(f"the clear found nothing to clear ({text}): the test must give it something")


@pytest.mark.parametrize("clear, plan", [(store.clear_bright_stars, store.clear_bright_stars_plan), (store.clear_phenomenon_scatter, store.clear_phenomenon_scatter_plan),
                                         (store.clear_nebula_field, store.clear_nebula_field_plan)])
def test_a_clear_takes_exactly_the_steps_of_its_plan(mysql_config, monkeypatch, clear, plan):
    """OPS.196: the plan is the one list of a clear's statements, and the pieces of the unit that runs it are its labels. A new galaxy has nothing to clear (OPS.254, OPS.257), so the test says it has."""
    for name in ("bright_stars_hold_nothing", "phenomenon_scatter_holds_nothing", "nebula_field_holds_nothing"):
        monkeypatch.setattr(store, name, lambda conn: False)
    steps = _Recorder()
    conn = store.get_connection(mysql_config)
    try:
        clear(conn, steps)
    finally:
        conn.close()
    labels = [label for label, _ in plan()]
    assert len(set(labels)) == len(labels) and steps.labels == labels, steps.labels


def test_a_pass_clear_and_the_index_drop_take_no_more_steps_than_they_say(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        steps = _Recorder()
        store.drop_phenomenon_scatter_indexes(conn, steps)
        assert 1 <= len(steps.labels) <= 1 + len(store.PHENOMENON_SCATTER_INDEXES)
        from planetgen.galaxy import scatter_limits
        limits = scatter_limits.ScatterLimits.of(store.get_galaxy_shape(conn), 0.0, 1.0) if store.get_galaxy_shape(conn) is not None else None
        if limits is not None:
            steps = _Recorder()
            store.clear_bright_star_pass(conn, "luminosity", limits, steps)
            assert steps.labels == [label for label, _ in store.clear_bright_star_pass_plan('luminosity', limits, store.pass_masses(conn, limits))]
    finally:
        conn.close()


def test_the_unit_counts_one_step_to_a_statement():
    seen = []
    phases = run_plan._UnitPhases(lambda done, total: seen.append((done, total)), 3, loud=False)
    for text in ("a", "b"):
        with phases.sql(text):
            pass
    phases.finish()
    assert seen == [(0, 3), (1, 3), (3, 3)]


def test_the_reset_names_each_table_it_empties(mysql_config, monkeypatch):
    """OPS.206: one `longsql.run` to a table that holds something (OPS.254 leaves the empty ones alone), labelled with the table and its place in the list, and the unit's count moves a table at a time."""
    from planetgen.cli import reset
    from tests.test_api import _place_sector
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    labels, counts = [], []
    real = longsql.run
    monkeypatch.setattr(longsql, "run", lambda conn, statement, label, *a, **k: labels.append((statement, label)) or real(conn, statement, label, *a, **k))
    wiped = reset.wipe(mysql_config, lambda done, total: counts.append((done, total)))
    assert 0 < len(labels) <= len(wiped.tables) and all(s.startswith("TRUNCATE TABLE ") for s, _ in labels)
    assert labels[0][1].startswith(f"Emptying {labels[0][0].split()[-1]} (") and f" of {len(wiped.tables)}" in labels[0][1]
    assert counts[-1] == (len(wiped.tables), len(wiped.tables))


def test_the_reset_unit_has_a_piece_to_a_table_and_one_to_tidy_after(mysql_config):
    """OPS.196: the unit's count is the pieces, a table at a time and the tidying last, and what it returns is what `wipe` returns."""
    from planetgen.cli import reset
    from planetgen.queue import commands, wrapper
    from tests.test_api import _place_sector
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    reports = []
    result = commands.run_alone("reset-galaxy", {"database": wrapper.database_args(mysql_config)}, lambda done, total: reports.append((done, total)))
    assert reports[0][0] == 0 and reports[-1] == (result["tables"] + 1, result["tables"] + 1) and result["with_rows"] > 0
    assert reset.reset_galaxy_pieces({"database": wrapper.database_args(mysql_config)}, 0)[-1] == reset.TIDY_PIECE
