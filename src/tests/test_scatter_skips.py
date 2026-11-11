"""
PERF.256: a `galaxy` rerun does not plan a scatter that is stored and current. A scatter stands when its rows are stored, the epoch it was drawn under is not one
that changed the scatter stage, and the settings the request would draw it with are the stored ones; the block light stage stands too when every scatter does,
the totals hold every unbuilt object and the run generates no sector. The stage list shows each such stage as skipped with its reason, and its estimate is no time.
"""

import argparse

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.db import block_lights, store
from planetgen.generation import run_galaxy, run_plan, stages
from planetgen.util import generator_epoch
from tests.test_phenomenon_scatter import _scatter_args, _seed_galaxy

FLOOR = "9000"
"""str: The luminosity floor `_scatter_args` draws at."""


def _galaxy_args(mysql_config, *extra):
    """The arguments of `planetgen galaxy` against `mysql_config` (the scatters' options are the galaxy's own)."""
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_galaxy_arguments(parser)
    return parser.parse_args(["--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port), "--mysql-user", mysql_config.user,
                              "--mysql-password", mysql_config.password, "--mysql-database", mysql_config.database, *extra])


def _stored(mysql_config, *kinds):
    """The scatters of `kinds` drawn on the seeded galaxy at the tests' threshold."""
    _seed_galaxy(mysql_config)
    for kind in kinds:
        run_plan.run_scatter(_scatter_args(mysql_config, f"--{kind}", "--phenomenon-min-mass", "12"))


def _skips(mysql_config, *extra):
    args = _galaxy_args(mysql_config, "--bright-star-min-luminosity", FLOOR, "--phenomenon-min-mass", "12", *extra)
    return run_plan.stored_scatter_skips(args)


def test_a_galaxy_with_nothing_stored_skips_nothing(mysql_config):
    _seed_galaxy(mysql_config)
    assert _skips(mysql_config) == {}


def test_a_galaxy_with_no_skeleton_skips_nothing(mysql_config):
    assert _skips(mysql_config) == {}


def test_the_star_passes_stand_together_with_the_settings_they_were_drawn_with(mysql_config):
    _stored(mysql_config, "mass", "luminosity")
    found = _skips(mysql_config)
    assert set(found) == {"mass-clear", "mass", "luminosity"} and all("current epoch" in reason for reason in found.values())


def test_the_star_passes_are_drawn_again_when_only_one_of_them_is_stored(mysql_config):
    _stored(mysql_config, "mass")
    assert _skips(mysql_config) == {}


@pytest.mark.parametrize("flags", (("--bright-star-min-luminosity", "5000"), ("--phenomenon-min-mass", "16")))
def test_the_star_passes_are_drawn_again_at_another_floor_or_mass_limit(mysql_config, flags):
    _stored(mysql_config, "mass", "luminosity")
    args = _galaxy_args(mysql_config, "--bright-star-min-luminosity", FLOOR, "--phenomenon-min-mass", "12")
    for name, value in zip(flags[::2], flags[1::2]):
        setattr(args, name.lstrip("-").replace("-", "_"), float(value))
    assert run_plan.stored_scatter_skips(args) == {}


def test_a_scatter_drawn_under_an_epoch_that_changed_the_scatter_stage_is_drawn_again(mysql_config, monkeypatch):
    _stored(mysql_config, "mass", "luminosity")
    assert _skips(mysql_config)
    monkeypatch.setattr(generator_epoch, "stages_between", lambda first, last: ("scatter",))
    assert _skips(mysql_config) == {}


def test_a_scatter_drawn_under_an_older_epoch_that_left_the_scatter_stage_alone_stands(mysql_config, monkeypatch):
    _stored(mysql_config, "mass", "luminosity")
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE scatter_epochs SET generator_epoch = 1")
        conn.commit()
    finally:
        conn.close()
    monkeypatch.setattr(generator_epoch, "stages_between", lambda first, last: ("fill", "population"))
    assert set(_skips(mysql_config)) == {"mass-clear", "mass", "luminosity"}


def test_the_phenomena_stand_with_their_indexes_and_the_nebula_field_they_were_drawn_with(mysql_config):
    _stored(mysql_config, "phenomena")
    assert set(_skips(mysql_config)) == {"phenomena-clear", "phenomena-unindex", "phenomena", "phenomena-index", "nebula-clear", "nebula-draw"}, "a phenomena scatter ends with the field, drawn and perhaps empty"
    conn = store.get_connection(mysql_config)
    try:
        store.forget_scatter_epoch(conn, "nebula-field")
        conn.commit()
    finally:
        conn.close()
    assert set(_skips(mysql_config)) == {"phenomena-clear", "phenomena-unindex", "phenomena", "phenomena-index"}
    conn = store.get_connection(mysql_config)
    try:
        store.record_scatter_epoch(conn, "nebula-field")
        store.forget_scatter_epoch(conn, "phenomena")
        conn.commit()
    finally:
        conn.close()
    assert _skips(mysql_config) == {}, "the field is drawn again whenever the phenomena are"


def test_the_phenomena_are_drawn_again_at_another_compact_mass(mysql_config):
    _stored(mysql_config, "phenomena")
    args = _galaxy_args(mysql_config, "--bright-star-min-luminosity", FLOOR, "--phenomenon-min-mass", "12", "--compact-min-mass", "30")
    assert "phenomena" not in run_plan.stored_scatter_skips(args)


def test_the_block_light_totals_stand_only_when_they_hold_every_unbuilt_object(mysql_config):
    _stored(mysql_config, "mass", "luminosity")
    conn = store.get_connection(mysql_config)
    try:
        assert block_lights.totals_current(conn), "the scatter's last stage made them"
        conn.execute("DELETE FROM galaxy_block_lights")
        conn.commit()
        assert not block_lights.totals_current(conn), "a stage that did not end left none"
        block_lights.rebuild(conn)
        conn.commit()
        assert block_lights.totals_current(conn)
        conn.execute("DELETE FROM galaxy_block_lights WHERE m = (SELECT MIN(m) FROM (SELECT m FROM galaxy_block_lights) AS levels)")
        conn.commit()
        assert not block_lights.totals_current(conn), "a level is missing"
    finally:
        conn.close()


def test_an_empty_galaxy_has_current_totals(mysql_config):
    _seed_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        assert block_lights.totals_current(conn)
    finally:
        conn.close()


def test_the_block_light_stage_is_skipped_only_when_every_scatter_stands(mysql_config):
    args = _galaxy_args(mysql_config)
    every = {key: "stands" for keys in run_plan.STORED_SCATTER_STAGES.values() for key in keys}
    _seed_galaxy(mysql_config)
    assert run_plan.block_lights_skip(args, every) is not None, "no object at all: the totals are whole"
    assert run_plan.block_lights_skip(args, {"mass": "stands", "luminosity": "stands"}) is None


def test_a_stage_that_stands_is_listed_as_skipped_with_its_reason_and_costs_no_time(mysql_config):
    args = _galaxy_args(mysql_config)
    args.stage_skips = {"mass": "stored", "luminosity": "stored", "block-lights": "whole"}
    listed = {stage.key: stage.skip for stage in stages.galaxy_stages(args)}
    assert listed["mass"] == listed["luminosity"] == "stored" and listed["block-lights"] == "whole" and listed["phenomena"] is None
    assert not any(stage.skip for stage in stages.galaxy_stages(_galaxy_args(mysql_config))), "no skips, no stage skipped"
    shared = stages.BLOCK_LIGHTS_STAGE
    assert shared.skip is None, "the shared stage is not changed"


def test_a_run_that_generates_no_sector_adds_the_block_light_stage_to_what_stands(mysql_config, monkeypatch):
    args = _galaxy_args(mysql_config)
    args.stage_skips = {"mass": "stored"}
    monkeypatch.setattr(run_galaxy.run_plan, "block_lights_skip", lambda found, skips: "whole")
    none = argparse.Namespace(key="sectors", count=0)
    run_galaxy._skip_block_lights(args, [none])
    assert args.stage_skips["block-lights"] == "whole"
    args.stage_skips = {"mass": "stored"}
    run_galaxy._skip_block_lights(args, [argparse.Namespace(key="sectors", count=3)])
    run_galaxy._skip_block_lights(args, [argparse.Namespace(key="neighborhood", count=0)])
    assert "block-lights" not in args.stage_skips, "a sector to make, or a stage that plans when it opens, keeps the stage"
    args.stage_skips = {}
    run_galaxy._skip_block_lights(args, [none])
    assert args.stage_skips == {}, "no scatter stands"


def _galaxy_run(mysql_config, monkeypatch, sectors=0):
    """What `planetgen galaxy --shell 3` plans: `(keys of the followers, {stage key: skip} of the stage list)`, the sectors stage holding `sectors` sectors to make."""
    from tests.test_bright_star_scatter import EDGE_PC
    from planetgen.generation import run_common
    monkeypatch.setattr(run_common, "_edge_pc", lambda: EDGE_PC)
    monkeypatch.setattr(run_galaxy, "_plan_sectors", lambda args, edge_pc: [run_galaxy.SectorsFollower(args, "sectors", "Sectors", edge_pc, count=sectors)])
    monkeypatch.setattr(run_galaxy, "_check_plan", lambda args, edge_pc: None)
    monkeypatch.setattr(run_galaxy, "_snap_to_region", lambda args, edge_pc: None)
    monkeypatch.setattr(run_galaxy, "_start_clock", lambda args: None)
    seen = {}
    monkeypatch.setattr(run_galaxy.followed, "run_followed", lambda args, followers, title, seed=None: seen.update(keys=[f.key for f in followers]))
    monkeypatch.setattr(run_galaxy.stages, "begin", lambda listed, args, title: seen.update(listed={stage.key: stage.skip for stage in listed}))
    monkeypatch.setattr(run_galaxy.stages, "finish", lambda reason=None: None)
    args = _galaxy_args(mysql_config, "--shell", "3", "--bright-star-min-luminosity", FLOOR, "--phenomenon-min-mass", "12")
    run_galaxy.run_galaxy(args)
    return seen["keys"], seen["listed"]


def test_a_second_identical_galaxy_run_plans_no_scatter_and_no_block_light_stage(mysql_config, monkeypatch):
    _stored(mysql_config, "mass", "luminosity", "phenomena", "nebula-field")
    keys, listed = _galaxy_run(mysql_config, monkeypatch)
    assert keys == ["sectors", "link", "backfill", "population", "territories", "settle", "run-blocks"]
    assert [key for key, reason in listed.items() if reason] == ["mass-clear", "mass", "luminosity", "phenomena-clear", "phenomena-unindex", "phenomena", "phenomena-index", "nebula-clear", "nebula-draw", "block-lights"]
    assert list(listed) == ["sectors", "link", "mass-clear", "mass", "luminosity", "phenomena-clear", "phenomena-unindex", "phenomena", "phenomena-index", "nebula-clear", "nebula-draw", "backfill", "population", "territories", "settle",
                            "run-blocks", "block-lights"], "the stage list keeps every stage, numbered"


def test_a_galaxy_run_that_makes_a_sector_keeps_the_block_light_stage_and_a_run_with_a_scatter_missing_plans_it(mysql_config, monkeypatch):
    _stored(mysql_config, "mass", "luminosity", "phenomena", "nebula-field")
    keys, _listed = _galaxy_run(mysql_config, monkeypatch, sectors=2)
    assert keys == ["sectors", "link", "backfill", "population", "territories", "settle", "run-blocks", "block-lights"]
    conn = store.get_connection(mysql_config)
    try:
        run_plan.store.clear_phenomenon_scatter(conn)
        conn.commit()
    finally:
        conn.close()
    keys, listed = _galaxy_run(mysql_config, monkeypatch)
    assert keys == ["sectors", "link", "phenomena-clear", "phenomena-unindex", "phenomena", "phenomena-index", "nebula-clear", "nebula-draw", "backfill", "population", "territories", "settle", "run-blocks", "block-lights"]
    assert listed["mass"] and listed["luminosity"] and not listed["phenomena"] and not listed["block-lights"]
