# tests/test_scatter_params_at_open.py

"""
PERF.212: the bright-star and phenomenon scatters take their seed, mass limit and clearing from the database when their stage opens,
not when the run is built, so a run can plan the galaxy and scatter it. The followers are built before the galaxy has a skeleton and
the stored output is what the verb gives after it.
"""

import pytest

from planetgen.db import store
from planetgen.generation import run_plan
from planetgen.util import log
from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy
from tests.test_phenomenon_scatter import _scatter_args as _phenomena_args

pytestmark = pytest.mark.db


@pytest.fixture(autouse=True)
def _the_console_log_writes_to_stdout():
    log.configure(log.NORMAL)


def _rows(rows):
    """The rows without the time they were written."""
    return [{key: value for key, value in dict(row).items() if key != "created_at"} for row in rows]


def _stars(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute("SELECT * FROM bright_stars ORDER BY id").fetchall()
    finally:
        conn.close()


def _phenomena(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute("SELECT * FROM phenomenon_scatter ORDER BY id").fetchall()
    finally:
        conn.close()


def test_a_star_pass_built_before_the_skeleton_exists_scatters_what_the_verb_does_after_it(mysql_config):
    args = _scatter_args(mysql_config, "--mass")
    planned = run_plan.planned_bright_stars(args, ("mass",))     # no skeleton yet: nothing is read to build it
    _seed_galaxy(mysql_config)
    from planetgen.generation import followed
    followed.run_planned(args, planned)
    built_first = _rows(_stars(mysql_config))
    assert built_first
    store.clear_bright_stars(store.get_connection(mysql_config))
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"), ("mass",))
    assert _rows(_stars(mysql_config)) == built_first


def test_the_phenomena_built_before_the_skeleton_exists_scatter_what_the_verb_does_after_it(mysql_config):
    args = _phenomena_args(mysql_config)
    planned = run_plan.planned_phenomena(args)
    _seed_galaxy(mysql_config)
    from planetgen.generation import followed
    followed.run_planned(args, planned)
    built_first = _rows(_phenomena(mysql_config))
    assert built_first
    run_plan.scatter_phenomena(_phenomena_args(mysql_config))
    assert _rows(_phenomena(mysql_config)) == built_first


def test_a_scatter_with_no_skeleton_still_refuses_in_words(mysql_config):
    with pytest.raises(RuntimeError, match="skeleton has never been built"):
        run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"), ("mass",))
    with pytest.raises(RuntimeError, match="skeleton has never been built"):
        run_plan.scatter_phenomena(_phenomena_args(mysql_config))


def test_a_luminosity_pass_with_no_mass_pass_stored_still_refuses_in_words(mysql_config):
    _seed_galaxy(mysql_config)
    with pytest.raises(RuntimeError, match="No massive-star pass is stored yet"):
        run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--luminosity"), ("luminosity",))


def test_both_passes_work_at_the_mass_limit_stored_before_the_scatter_cleared_it(mysql_config):
    """The mass pass empties the scatter's record when it opens, so the luminosity pass must not read the limit again: both take what
    was read when the first opened, the stored limit when no `--phenomenon-min-mass` is given."""
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--phenomenon-min-mass", "8"), ("mass",))
    conn = store.get_connection(mysql_config)
    try:
        assert store.bright_star_mass_limit(conn) == 8.0
    finally:
        conn.close()
    args = _scatter_args(mysql_config, "--mass")
    args.phenomenon_min_mass = None
    run_plan.scatter_bright_stars(args)
    conn = store.get_connection(mysql_config)
    try:
        assert store.bright_star_mass_limit(conn) == 8.0
    finally:
        conn.close()
