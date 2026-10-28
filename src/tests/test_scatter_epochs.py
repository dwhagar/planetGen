# tests/test_scatter_epochs.py

"""
GEN.764: every stage of a stored galaxy records the generator epoch it was made under. A scatter's draw writes its row in
`scatter_epochs` and a clear of it deletes the row; a species carries the epoch of the pass that made it; and
`galaxy/epoch_record.py` reads them together with the plan's and the sectors', and says what is behind the running epoch.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy import epoch_record
from planetgen.generation import run_plan
from planetgen.util import generator_epoch

from tests.test_bright_star_scatter import _clear_args, _scatter_args, _seed_galaxy, light_steps  # noqa: F401 -- a fixture


def _epochs(config):
    conn = store.get_connection(config)
    try:
        return store.scatter_epochs(conn)
    finally:
        conn.close()


def test_each_scatter_records_its_epoch_and_a_clear_forgets_it(mysql_config, monkeypatch):
    current = generator_epoch.generator_epoch()
    _seed_galaxy(mysql_config)
    run_plan.run_scatter(_scatter_args(mysql_config, "--mass"))
    assert _epochs(mysql_config) == {"mass": current}
    run_plan.run_scatter(_scatter_args(mysql_config, "--luminosity"))
    run_plan.run_scatter(_scatter_args(mysql_config, "--phenomena"))
    run_plan.run_scatter(_scatter_args(mysql_config, "--nebula-field"))
    assert _epochs(mysql_config) == {"mass": current, "luminosity": current, "phenomena": current,
                                     "nebula-field": current}
    run_plan.run_clear_scatter(_clear_args(mysql_config, "--nebula-field"))
    run_plan.run_clear_scatter(_clear_args(mysql_config, "--luminosity"))
    assert _epochs(mysql_config) == {"mass": current, "phenomena": current}
    # drawn again under a later epoch, the record follows
    monkeypatch.setattr(generator_epoch, "generator_epoch", lambda: current + 1)
    run_plan.run_scatter(_scatter_args(mysql_config, "--luminosity"))
    assert _epochs(mysql_config) == {"mass": current, "luminosity": current + 1, "phenomena": current}


def test_a_whole_scatter_clear_forgets_the_star_kinds(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.run_scatter(_scatter_args(mysql_config, "--mass"))
    run_plan.run_scatter(_scatter_args(mysql_config, "--luminosity"))
    conn = store.get_connection(mysql_config)
    try:
        store.record_scatter_epoch(conn, "band")
        conn.commit()
        store.clear_bright_stars(conn)
        assert store.scatter_epochs(conn) == {}
    finally:
        conn.close()


def test_a_kind_that_does_not_exist_is_refused(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with pytest.raises(ValueError, match="kind must be"):
            store.record_scatter_epoch(conn, "backfill")
    finally:
        conn.close()


def test_the_record_says_what_is_behind_the_running_epoch(mysql_config):
    current = generator_epoch.generator_epoch()
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE galaxy_shape SET generator_epoch = ? WHERE id = 1", (current,))
        store.record_scatter_epoch(conn, ["mass", "phenomena"])
        conn.commit()
        found = epoch_record.recorded(conn)
        assert found["plan"] == current and found["scatter"] == {"mass": current, "phenomena": current}
        assert found["fill"] == {} and found["population"] == {}
        assert epoch_record.behind(conn) == {"stale": {"scatter": {}, "fill": {}, "population": {}}, "stages": ()}
        # the plan at epoch 59 and a scatter at 95: epochs 60 to 96 changed the plan (64) and the fill, and 96 only the fill
        conn.execute("UPDATE galaxy_shape SET generator_epoch = 59 WHERE id = 1")
        conn.execute("UPDATE scatter_epochs SET generator_epoch = 95 WHERE kind = 'mass'")
        conn.commit()
        found = epoch_record.behind(conn, current=96)
        assert found["stale"]["plan"] == 59 and found["stale"]["scatter"] == {} and found["stages"] == generator_epoch.stages()
        # a scatter at 95 needs no redo for epoch 96, which changed the fill alone; a fill at 95 does
        conn.execute("UPDATE galaxy_shape SET generator_epoch = 96 WHERE id = 1")
        conn.commit()
        found = epoch_record.behind(conn, current=96)
        assert found["stale"] == {"scatter": {}, "fill": {}, "population": {}} and found["stages"] == ()
    finally:
        conn.close()
