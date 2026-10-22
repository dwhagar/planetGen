"""PERF.123 (slice B): `planetgen sectors` generates sectors in the galaxy's modes and nothing else; `link`, `backfill` and
`settle` then do the stages after them for the sectors saved since the run's start time."""

import sys

import pytest

from planetgen.cli import generate
from planetgen.generation import stages
from planetgen.queue import commands, wrapper
from tests.bughunt_support import mysql_argv
from tests.test_galaxy_gen import _all_sectors, _plan_wide_galaxy

pytestmark = [pytest.mark.db, pytest.mark.usefixtures("unit_density_plans")]


def _args(*argv):
    parser, command_parsers = generate.build_parser()
    args = parser.parse_args(list(argv))
    generate.validate_shared_generation_args(args, command_parsers[args.command])
    generate.validate_sectors_args(args, command_parsers[args.command])
    return args


def _main(monkeypatch, mysql_config, *argv):
    monkeypatch.setattr(sys, "argv", ["generate", *argv, *mysql_argv(mysql_config)])
    generate.main()


def test_the_stages_of_a_sectors_run_are_the_sectors_alone():
    assert [stage.key for stage in stages.for_command("sectors", _args("sectors"))] == ["start", "neighborhood"]
    assert [stage.key for stage in stages.for_command("sectors", _args("sectors", "--shell", "1"))] == ["sectors"]
    assert [step.kind for step in wrapper.steps_for("sectors", _args("sectors", "--shell", "1"))] == ["generate-sectors"]
    assert [step.kind for step in wrapper.steps_for("sectors", _args("sectors"))] == [
        "generate-start-sector", "generate-neighborhood"]
    assert all(step.kind in commands.REGISTRY for step in wrapper.steps_for("sectors", _args("sectors")))


def test_the_sector_options_are_the_same_on_sectors_and_galaxy_and_the_later_stages_are_only_on_galaxy(capsys):
    sectors_help = generate.build_parser()[1]["sectors"].format_help()
    galaxy_help = generate.build_parser()[1]["galaxy"].format_help()
    for option in ("--shell", "--band", "--wedge", "--shells", "--course", "--center-sector",
                   "--radius-pc", "--neighborhoods", "--estimate-only", "--limit"):
        assert option in sectors_help and option in galaxy_help, option
    for option in ("--then-scatter", "--no-settle", "--backfill-from"):
        assert option in galaxy_help and option not in sectors_help, option
    with pytest.raises(SystemExit):
        generate.build_parser()[0].parse_args(["sectors", "--then-scatter"])


def test_the_modes_are_refused_the_same_way_on_sectors():
    with pytest.raises(SystemExit):
        _args("sectors", "--course", "a", "b", "--shell", "1")
    with pytest.raises(SystemExit):
        _args("sectors", "--band", "0")


def test_sectors_saves_the_sectors_and_says_how_to_link_them_and_the_later_stages_do_the_rest(
        mysql_config, monkeypatch, redis_server, capsys):
    _plan_wide_galaxy(mysql_config)
    _main(monkeypatch, mysql_config, "sectors", "--shell", "1", "--band", "0", "--num-systems", "1")
    out = capsys.readouterr().out
    assert len(_all_sectors(mysql_config)) == 12
    assert "Run started at " in out and "planetgen link --since" in out
    assert _count(mysql_config, "sector_paths") == 0, "the paths wait for 'planetgen settle'"
    _main(monkeypatch, mysql_config, "settle", "--since", "2000-01-01 00:00:00")
    assert _count(mysql_config, "sector_paths") > 0


def _count(mysql_config, table):
    from planetgen.db import store
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
    finally:
        conn.close()


def _planned_params(mysql_config, *argv):
    from planetgen.generation import run_common, run_galaxy
    args = _args("sectors", *argv)
    return args, {"database": wrapper.database_args(mysql_config), "args": run_galaxy._sector_args_json(args),
                  "edge_pc": run_common._edge_pc(), "items": [], "planned": True}


def test_a_span_is_cut_into_the_sectors_the_outline_holds_less_those_that_exist_and_up_to_the_limit(mysql_config):
    from planetgen.generation import run_galaxy
    _plan_wide_galaxy(mysql_config)
    args, params = _planned_params(mysql_config, "--shells", "1:2", "--bands", "0", "--num-systems", "1")
    bounds = run_galaxy._BatchDensity(mysql_config).bounds
    wanted = [index for address in args.span.addresses(bounds) for index in address]
    assert wanted and run_galaxy.planned_items(params) == wanted
    first = run_galaxy.expand_generate_sectors(params)
    assert [index for unit, _cost in first.units for index in unit["items"]] == wanted
    _args_limited, limited = _planned_params(mysql_config, "--shells", "1:2", "--bands", "0", "--limit", "3",
                                             "--num-systems", "1")
    assert run_galaxy.planned_items(limited) == wanted[:9]


def test_a_shell_run_hands_its_stage_the_arguments_that_name_the_sectors_not_the_sectors(
        mysql_config, monkeypatch, redis_server):
    _plan_wide_galaxy(mysql_config)
    seen = []
    real = wrapper.step_for

    def spy(args, key, title, params, *rest, **kwargs):
        seen.append((key, params))
        return real(args, key, title, params, *rest, **kwargs)

    monkeypatch.setattr(wrapper, "step_for", spy)
    _main(monkeypatch, mysql_config, "sectors", "--shell", "1", "--band", "0", "--limit", "2", "--num-systems", "1")
    [(key, params)] = seen
    assert key == "sectors" and params["items"] == [] and params["planned"] is True
    assert len(_all_sectors(mysql_config)) == 2
