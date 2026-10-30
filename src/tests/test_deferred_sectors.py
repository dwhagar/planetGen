"""PERF.214: the sectors planner runs when the sectors stage opens, in a worker, and not in the process that starts the run
(`run_galaxy._plan_sectors(..., deferred=True)`). The stages are the ones the starter's plan gives, they are told only the arguments
that name their sectors, and the stage works them out when it opens (`planned_items`); the neighborhood after a random start reads the
starts the stage before it planned (`wrapper.plan_detail`). The sectors, and so the stored galaxy, are the ones the starter's plan makes."""

import sys

import pytest

from planetgen.cli import generate
from planetgen.db import store
from planetgen.generation import run_common, run_galaxy, stages
from planetgen.queue import wrapper
from planetgen.util import draw
from tests.bughunt_support import mysql_argv
from tests.test_galaxy_gen import EDGE_PC, WIDE_GALAXY_HIGH_BAND, _address, _all_sectors, _plan_wide_galaxy

pytestmark = [pytest.mark.db, pytest.mark.usefixtures("unit_density_plans")]

MODES = (
    ("random start", ["--radius-pc", "4"]),
    ("random starts", ["--neighborhoods", "2", "--radius-pc", "4"]),
    ("a shell", ["--shell", "1", "--band", "0"]),
    ("a span", ["--shells", "1:2", "--bands", "0"]),
    ("a single sector", ["--shell", "1", "--band", "0", "--wedge", "0"]),
    ("a single sector and its neighborhood", ["--shell", "1", "--band", "0", "--wedge", "0", "--radius-pc", "4"]),
    ("a neighborhood", ["--center-sector", "1", "--radius-pc", "4"]),
)


def _args(*argv):
    parser, command_parsers = generate.build_parser()
    args = parser.parse_args(["sectors", *argv])
    generate.validate_shared_generation_args(args, command_parsers[args.command])
    generate.validate_sectors_args(args, command_parsers[args.command])
    return args


def _params(follower, mysql_config):
    return {"database": wrapper.database_args(mysql_config), **follower.params()}


@pytest.mark.parametrize("name,argv", MODES, ids=[mode[0] for mode in MODES])
def test_the_deferred_stages_are_the_stages_of_the_run_and_read_nothing(name, argv, monkeypatch):
    args = _args(*argv)

    def refuse(*_args, **_kwargs):
        raise AssertionError("a deferred plan read the database")

    monkeypatch.setattr(store, "get_connection", refuse)
    followers = run_galaxy._plan_sectors(args, EDGE_PC, deferred=True)
    assert [f.key for f in followers] == [stage.key for stage in stages.sectors_stages(args)]
    assert all(f.items == [] and f.planned for f in followers)


def test_a_deferred_plan_cannot_be_estimated():
    with pytest.raises(ValueError):
        run_galaxy._plan_sectors(_args("--shell", "1", "--band", "0", "--estimate-only"), EDGE_PC, deferred=True)


def test_a_wedge_out_of_range_is_refused_before_the_run_starts():
    with pytest.raises(SystemExit):
        run_galaxy._plan_sectors(_args("--shell", "1", "--band", "0", "--wedge", "9999"), EDGE_PC, deferred=True)


def test_the_neighborhood_after_a_random_start_reads_the_starts_the_stage_before_planned():
    start, neighborhood = run_galaxy._plan_sectors(_args("--radius-pc", "4"), EDGE_PC, deferred=True)
    assert start.params()["planned"] == "random-start" and isinstance(start.params()["run_seed"], int)
    assert neighborhood.params()["centers"] == wrapper.plan_detail("start", "centers")


def test_a_single_sector_run_gives_its_neighborhood_the_sector_itself():
    _sector, neighborhood = run_galaxy._plan_sectors(
        _args("--shell", "1", "--band", "0", "--wedge", "0", "--radius-pc", "4"), EDGE_PC, deferred=True)
    assert neighborhood.params()["centers"] == [[1, 0, 0]]


def _cut(follower, mysql_config, **more):
    """The units' sectors of `follower`'s stage, flat, and the plan summary it records."""
    expansion = run_galaxy.expand_generate_sectors({**_params(follower, mysql_config), **more})
    return [index for unit, _cost in expansion.units for index in unit["items"]], expansion.plan


def test_the_random_starts_are_the_starts_of_the_starters_plan_and_the_draws_after_them_are_the_same(mysql_config):
    _plan_wide_galaxy(mysql_config)
    args = _args("--neighborhoods", "2", "--radius-pc", "8", "--max-shell", "60", "--num-systems", "1", *mysql_argv(mysql_config))
    progress = run_common._generation_progress(disable=True)
    draw.set_run_seed(4321)
    with progress as shown:
        from planetgen.util import log
        log.set_console(shown.console)
        try:
            starter = run_galaxy._plan_random_start(args, EDGE_PC, shown)
        finally:
            log.reset_console()
    after_starter = draw.getrandbits(32)
    draw.set_run_seed(4321)
    deferred = run_galaxy._plan_sectors(args, EDGE_PC, deferred=True)
    assert draw.getrandbits(32) == after_starter
    starts, plan = _cut(deferred[0], mysql_config)
    assert starts == starter[0].items and 1 <= len(starts) // 3 <= 2
    assert plan["centers"] == [list(starts[place:place + 3]) for place in range(0, len(starts), 3)]
    batch_density = run_galaxy._BatchDensity(mysql_config)
    for place in range(0, len(starts), 3):          # the stage before saves the starts, then the neighborhoods are cut
        address = tuple(starts[place:place + 3])
        position = run_galaxy.sector_position_pc(*address, EDGE_PC)
        run_galaxy.generate_and_save_sector_at(batch_density.resolve(args, address, position), address, position, EDGE_PC)
    around, _plan = _cut(deferred[1], mysql_config, centers=plan["centers"])
    expected, _plan = _cut(starter[1], mysql_config)
    assert around == expected and around


def test_a_shell_cut_when_the_stage_opens_is_the_sectors_the_outline_holds_less_those_that_exist(mysql_config):
    _plan_wide_galaxy(mysql_config)
    args = _args("--shells", "1:2", "--bands", "0", "--num-systems", "1", *mysql_argv(mysql_config))
    bounds = run_galaxy._BatchDensity(mysql_config).bounds
    wanted = [index for address in args.span.addresses(bounds) for index in address]
    [follower] = run_galaxy._plan_sectors(args, EDGE_PC, deferred=True)
    assert _cut(follower, mysql_config)[0] == wanted and wanted


def test_a_single_sector_that_exists_is_not_cut_again_and_one_beyond_the_reach_is_refused(mysql_config):
    _plan_wide_galaxy(mysql_config)
    args = _args("--shell", "1", "--band", "0", "--wedge", "0", *mysql_argv(mysql_config))
    [follower] = run_galaxy._plan_sectors(args, EDGE_PC, deferred=True)
    assert _cut(follower, mysql_config)[0] == [1, 0, 0]
    run_galaxy.generate_and_save_sector_at(
        run_galaxy._BatchDensity(mysql_config).resolve(args, (1, 0, 0), run_galaxy.sector_position_pc(1, 0, 0, EDGE_PC)),
        (1, 0, 0), run_galaxy.sector_position_pc(1, 0, 0, EDGE_PC), EDGE_PC)
    assert _cut(follower, mysql_config)[0] == []


# --- whole runs: the stored galaxy is the one the starter's plan makes ---------------------------------------------------------------

@pytest.fixture
def deferred(monkeypatch):
    real = run_galaxy._plan_sectors
    monkeypatch.setattr(run_galaxy, "_plan_sectors", lambda args, edge_pc: real(args, edge_pc, deferred=True))


def _main(monkeypatch, mysql_config, *argv):
    monkeypatch.setattr(sys, "argv", ["generate", *argv, *mysql_argv(mysql_config)])
    generate.main()


def _addresses(mysql_config):
    return sorted(_address(row) for row in _all_sectors(mysql_config))


def test_a_shell_run_planned_when_its_stage_opens_saves_the_sectors_the_starters_plan_does(
        mysql_config, fresh_mysql_config, monkeypatch, redis_server):
    _plan_wide_galaxy(mysql_config)
    _plan_wide_galaxy(fresh_mysql_config)
    argv = ["sectors", "--shell", "1", "--band", "0", "--num-systems", "1"]
    _main(monkeypatch, fresh_mysql_config, *argv)
    real = run_galaxy._plan_sectors
    monkeypatch.setattr(run_galaxy, "_plan_sectors", lambda args, edge_pc: real(args, edge_pc, deferred=True))
    _main(monkeypatch, mysql_config, *argv)
    assert len(_addresses(mysql_config)) == 12
    assert _addresses(mysql_config) == _addresses(fresh_mysql_config)


def test_a_random_start_planned_when_its_stage_opens_is_a_start_and_the_sectors_around_it(mysql_config, monkeypatch, deferred, redis_server):
    import math
    _plan_wide_galaxy(mysql_config)
    radius_pc = 8.0
    _main(monkeypatch, mysql_config, "sectors", "--max-shell", "1", "--radius-pc", str(radius_pc), "--num-systems", "1")
    sectors = _all_sectors(mysql_config)
    assert len(sectors) > 1
    start = min(sectors, key=lambda row: row["id"])
    center = (start["center_x_pc"], start["center_y_pc"], start["center_z_pc"])
    assert start["shell_index"] <= 1
    assert all(math.dist(center, (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])) <= radius_pc + EDGE_PC for row in sectors)
    assert all(abs(row["band_index"]) <= WIDE_GALAXY_HIGH_BAND for row in sectors)


def test_a_single_sector_and_its_neighborhood_planned_when_the_stages_open(mysql_config, monkeypatch, deferred, redis_server):
    # Shell 30, as the single-address test of test_galaxy_gen does: the toy shape's dense bulge makes a neighborhood of thousands of systems.
    from planetgen.galaxy.spherical import enumerate_sectors_within_radius
    _plan_wide_galaxy(mysql_config)
    _main(monkeypatch, mysql_config, "sectors", "--shell", "30", "--band", "0", "--wedge", "4", "--radius-pc", "5")
    addresses = set(_addresses(mysql_config))
    assert (30, 0, 4) in addresses and len(addresses) > 1
    center = run_galaxy.sector_position_pc(30, 0, 4, EDGE_PC)
    assert addresses <= {(s, b, w) for s, b, w, *_ in enumerate_sectors_within_radius(center, 5.0, EDGE_PC)}
