"""GEN.1128 (second slice): `planetgen sectors --around X Y Z --radius-pc R [--grow-shells K] [--unit-sectors N]` runs the plan of
`galaxy/region_run.py` through the sector stage: one unit of work per region unit, nearest first, and a rerun after a stop builds only the
sectors not saved yet."""

import sys

import pytest

from planetgen.cli import generate
from planetgen.galaxy import region_run, spherical
from planetgen.generation import run_galaxy, stages
from planetgen.queue import wrapper
from tests.bughunt_support import mysql_argv
from tests.test_galaxy_gen import EDGE_PC, _all_sectors, _plan_wide_galaxy

pytestmark = [pytest.mark.db, pytest.mark.usefixtures("unit_density_plans")]

AROUND = ("30", "0", "0")
RADIUS = "6"


def _args(*argv):
    parser, command_parsers = generate.build_parser()
    args = parser.parse_args(["sectors", *argv])
    generate.validate_shared_generation_args(args, command_parsers[args.command])
    generate.validate_sectors_args(args, command_parsers[args.command])
    return args


def _main(monkeypatch, mysql_config, *argv):
    monkeypatch.setattr(sys, "argv", ["generate", *argv, *mysql_argv(mysql_config)])
    generate.main()


def _expansion(args, mysql_config):
    [follower] = run_galaxy._plan_sectors(args, EDGE_PC, deferred=True)
    return run_galaxy.expand_generate_sectors({"database": wrapper.database_args(mysql_config), **follower.params()})


def _expected(args):
    return region_run.plan(tuple(args.around), args.radius_pc, EDGE_PC, block_sectors=args.unit_sectors or region_run.BLOCK_SECTORS,
                           grow_shells=args.grow_shells or 0)


def test_a_region_run_is_one_sectors_stage_on_the_generate_sectors_command():
    args = _args("--around", *AROUND, "--radius-pc", RADIUS)
    assert [stage.key for stage in stages.for_command("sectors", args)] == ["sectors"]
    assert [step.kind for step in wrapper.steps_for("sectors", args)] == ["generate-sectors"]


def test_the_units_of_the_cut_are_the_region_units_nearest_first(mysql_config):
    _plan_wide_galaxy(mysql_config)
    args = _args("--around", *AROUND, "--radius-pc", RADIUS, "--unit-sectors", "5", *mysql_argv(mysql_config))
    expansion = _expansion(args, mysql_config)
    expected = _expected(args)
    assert len(expansion.units) == len(expected.units) > 1
    for (unit, cost), planned in zip(expansion.units, expected.units):
        assert unit["items"] == region_run.unit_items(planned, spherical.wedge_count)
        assert cost == planned.count <= 5 and unit["again"] is True
    assert expansion.plan["sectors"] == expected.sectors and expansion.plan["skipped"] == 0


def test_growing_shells_makes_the_region_larger_and_a_unit_size_changes_only_the_cut(mysql_config):
    _plan_wide_galaxy(mysql_config)
    base = _expansion(_args("--around", *AROUND, "--radius-pc", RADIUS, *mysql_argv(mysql_config)), mysql_config)
    grown = _expansion(_args("--around", *AROUND, "--radius-pc", RADIUS, "--grow-shells", "1", *mysql_argv(mysql_config)), mysql_config)
    small = _expansion(_args("--around", *AROUND, "--radius-pc", RADIUS, "--unit-sectors", "2", *mysql_argv(mysql_config)), mysql_config)
    assert grown.plan["sectors"] > base.plan["sectors"]
    assert small.plan["sectors"] == base.plan["sectors"] and len(small.units) > len(base.units)


def test_a_rerun_after_a_stop_builds_only_what_is_missing(mysql_config, monkeypatch, redis_server):
    _plan_wide_galaxy(mysql_config)
    argv = ("--around", *AROUND, "--radius-pc", "3", "--unit-sectors", "4", "--num-systems", "1")
    args = _args(*argv, *mysql_argv(mysql_config))
    total = _expected(args).sectors
    _main(monkeypatch, mysql_config, "sectors", *argv)
    assert len(_all_sectors(mysql_config)) == total > 0
    assert _expansion(args, mysql_config).units == []
    _main(monkeypatch, mysql_config, "sectors", *argv)
    assert len(_all_sectors(mysql_config)) == total


@pytest.mark.parametrize("argv,message", [
    (["--around", "0", "0", "0"], "--around requires --radius-pc."),
    (["--around", "0", "0", "0", "--radius-pc", "0"], "--radius-pc must be a positive number."),
    (["--around", "0", "0", "0", "--radius-pc", "5", "--grow-shells", "-1"], "--grow-shells must be >= 0."),
    (["--around", "0", "0", "0", "--radius-pc", "5", "--unit-sectors", "0"], "--unit-sectors must be >= 1."),
    (["--around", "0", "0", "0", "--radius-pc", "5", "--limit", "3"], "--around can't be combined with --limit."),
    (["--shell", "1", "--grow-shells", "1"], "--grow-shells requires --around."),
    (["--radius-pc", "5", "--unit-sectors", "4"], "--unit-sectors requires --around."),
])
def test_a_region_run_refuses_what_cannot_go_with_it(argv, message, capsys):
    with pytest.raises(SystemExit):
        _args(*argv)
    assert message in capsys.readouterr().err


def test_a_region_run_cannot_be_given_two_modes(capsys):
    with pytest.raises(SystemExit):
        _args("--around", "0", "0", "0", "--radius-pc", "5", "--shell", "1")
    assert "not allowed with argument" in capsys.readouterr().err
