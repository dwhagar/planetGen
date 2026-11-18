# tests/test_worker_warm.py

"""PERF.42: the warm set (the bright-star tables) is the table a cold call builds, and a unit run after it finds them all. OPS.235: the queue worker and its release check are gone; the
template of the supervisor (`planetgen.queue.template`) warms these functions."""

import pytest

from planetgen.generation import star_population


def test_a_warm_table_is_the_table_a_cold_call_builds():
    from planetgen.queue import units
    star_population._bright_table.cache_clear()
    cold = star_population._bright_table(1000.0, star_population.POPULATIONS[0], None, None, -0.5)
    star_population._bright_table.cache_clear()
    units.warm_star_tables()
    assert star_population._bright_table(1000.0, star_population.POPULATIONS[0], None, None, -0.5) == cold


def test_every_warm_function_a_command_names_is_there():
    from planetgen.queue import commands
    named = [kind for kind, spec in commands.REGISTRY.items() if spec.warm is not None]
    assert {"scatter-mass", "scatter-luminosity", "scatter-band"} <= set(named)
    for kind in named:
        assert callable(commands.warm_function(kind))
    assert commands.warm_function("plan-skeleton") is None


@pytest.mark.fresh_database
def test_a_warmed_scatter_unit_builds_no_sampling_table_of_its_own(mysql_config, monkeypatch):
    """The worker's warm function fills `star_population._bright_table` with the very tables the unit looks up, so a unit
    run after it finds them all."""
    from planetgen.db import store
    from planetgen.generation import run_plan
    from planetgen.queue import commands
    from tests.test_bright_star_scatter import _scatter_args, _seed_galaxy
    seen = []
    real = run_plan.scatter_shells_pieces

    def recording(unit_args, cursor):
        seen.append(unit_args)
        return real(unit_args, cursor)

    monkeypatch.setattr(run_plan, "scatter_shells_pieces", recording)
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass", "--phenomenon-min-mass", "8",
                                                "--bright-star-min-luminosity", "2500", "--workers", "1"), ("mass",))
    assert seen
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM bright_stars")
        conn.commit()
    finally:
        conn.close()
    star_population._bright_table.cache_clear()
    run_plan.scatter_shells_warm(seen[0])
    warmed = star_population._bright_table.cache_info().misses
    assert warmed > 0
    commands.run_alone("scatter-mass", seen[0])
    assert star_population._bright_table.cache_info().misses == warmed
