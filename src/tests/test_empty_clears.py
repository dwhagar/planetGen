"""
OPS.257: a run on a wiped galaxy clears no scatter. A clearing stage whose target holds nothing (no row, a counter never moved, no recorded limits or epoch) is not planned, by whatever route
the run was started: the stage is listed as skipped with its reason and its follower is left out, in the CLI `scatter` and `galaxy`, the web's New galaxy steps (the same commands, each planned after the
wipe step has run) and `run-job`. A scatter that holds something is cleared as before.
"""

from planetgen.db import store
from planetgen.generation import job, run_plan, stages
from tests.test_phenomenon_scatter import _scatter_args, _seed_galaxy
from tests.test_run_job import FLOOR, _job_args
from tests.test_scatter_skips import _galaxy_args


def _holds(mysql_config, work):
    conn = store.get_connection(mysql_config)
    try:
        return work(conn)
    finally:
        conn.close()


def test_a_galaxy_with_nothing_stored_holds_nothing_to_clear(mysql_config):
    _seed_galaxy(mysql_config)
    assert _holds(mysql_config, store.phenomenon_scatter_holds_nothing)
    assert _holds(mysql_config, store.nebula_field_holds_nothing)


def test_a_stored_scatter_is_something_to_clear(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.run_scatter(_scatter_args(mysql_config, "--phenomena", "--phenomenon-min-mass", "12"))
    assert not _holds(mysql_config, store.phenomenon_scatter_holds_nothing)
    assert not _holds(mysql_config, store.nebula_field_holds_nothing)


def test_the_record_of_a_drawn_scatter_is_something_to_clear_even_with_no_row(mysql_config):
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        store.record_scatter_epoch(conn, ("phenomena", "nebula-field"), None)
        conn.commit()
    finally:
        conn.close()
    assert not _holds(mysql_config, store.phenomenon_scatter_holds_nothing)
    assert not _holds(mysql_config, store.nebula_field_holds_nothing)


def test_the_clears_of_a_new_galaxy_are_not_planned(mysql_config):
    _seed_galaxy(mysql_config)
    args = _scatter_args(mysql_config, "--phenomena", "--phenomenon-min-mass", "12")
    skips = run_plan.empty_clear_skips(args, ("phenomena",))
    assert set(skips) == {"phenomena-clear", "nebula-clear"}
    args.stage_skips = skips
    listed = {stage.key: stage.skip for stage in stages.scatter_stages(args)}
    assert listed["phenomena-clear"] and listed["nebula-clear"] and listed["phenomena"] is None
    followers, _book = run_plan.phenomena_followers(args)
    assert [f.key for f in followers if f.key.endswith("clear")] == []


def test_a_run_with_no_skeleton_or_an_earlier_scatter_in_its_job_plans_every_clear(mysql_config):
    args = _scatter_args(mysql_config, "--phenomena", "--phenomenon-min-mass", "12")
    assert run_plan.empty_clear_skips(args, ("phenomena",)) == {}
    _seed_galaxy(mysql_config)
    assert run_plan.empty_clear_skips(args, ("phenomena",), fresh=False) == {}


def test_the_bright_star_clear_stays_for_the_limits_it_stores(mysql_config):
    _seed_galaxy(mysql_config)
    assert run_plan.empty_clear_skips(_scatter_args(mysql_config, "--mass"), ("mass", "luminosity")) == {}


def test_a_cli_scatter_of_the_phenomena_on_a_new_galaxy_clears_nothing(mysql_config, monkeypatch):
    _seed_galaxy(mysql_config)
    for name in ("clear_phenomenon_scatter", "clear_nebula_field"):
        monkeypatch.setattr(store, name, lambda *a, **k: (_ for _ in ()).throw(AssertionError("a clear of nothing ran")))
    summary = run_plan.run_scatter(_scatter_args(mysql_config, "--phenomena", "--phenomenon-min-mass", "12"))
    assert summary is None or isinstance(summary, dict)


def test_a_cli_galaxy_names_the_clears_it_does_not_plan(mysql_config):
    _seed_galaxy(mysql_config)
    args = _galaxy_args(mysql_config, "--bright-star-min-luminosity", "18000", "--phenomenon-min-mass", "12")
    args.stage_skips = {**run_plan.empty_clear_skips(args, ("phenomena",)), **run_plan.stored_scatter_skips(args)}
    assert {"phenomena-clear", "nebula-clear"} <= set(args.stage_skips)
    assert "mass-clear" not in args.stage_skips
    listed = {stage.key: stage.skip for stage in stages.galaxy_stages(args)}
    assert listed["nebula-clear"] and listed["phenomena-clear"] and listed["mass-clear"] is None


def test_a_job_plans_no_clear_of_nothing_and_its_estimates_follow_its_stages(mysql_config, monkeypatch):
    from planetgen.generation import stats
    monkeypatch.setattr(stats, "stage_seconds", lambda conn, key, settings: 1.0)
    _seed_galaxy(mysql_config)
    args = _job_args(mysql_config, [("scatter", ["--phenomena"]), ("scatter", ["--mass", *FLOOR])])
    conn = store.get_connection(mysql_config)
    try:
        parts = job.build_parts(args.job, args)
        keys = [stage.key for stage in job.stage_list(parts)]
        assert "phenomena-clear" not in keys and "nebula-clear" not in keys and "phenomena" in keys
        for part in parts:
            assert [f.key for f in part.followers] == [stage.key for stage in part.stages]
        assert len(job.stage_estimates(conn, args.job)) == len(keys)
    finally:
        conn.close()
