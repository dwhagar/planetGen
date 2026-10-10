"""
GEN.196: `planetgen plan --redo-scatters` redoes only the scatters named. A star pass alone clears and
rewrites its own stars and keeps the other's; the phenomena redo keeps the stars and always has its nucleus;
the stage list names every stage with the ones not chosen skipped.
"""

import argparse

import pytest

from planetgen import tuning
from planetgen.cli import generate as generate_cli
from planetgen.db import store
from planetgen.generation import run_plan, stages
from planetgen.generation import phenomenon_scatter as scatter
from tests.test_bright_star_scatter import EDGE_PC, EXTENTS, E_VALUE, SHAPE, THRESHOLD, light_rings  # noqa: F401


def _args(mysql_config, *extra):
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_plan_arguments(parser)
    args = parser.parse_args([
        "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
        "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password,
        "--mysql-database", mysql_config.database, "--workers", "1",
        "--bright-star-min-luminosity", str(THRESHOLD), *extra,
    ])
    generate_cli.validate_plan_args(args, parser)
    return args


def _seed(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=8,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config)
    store.replace_galaxy_layers(EXTENTS, config=mysql_config)


def _star_ids(mysql_config, where):
    conn = store.get_connection(mysql_config)
    try:
        return {row["id"] for row in conn.execute(f"SELECT id FROM bright_stars WHERE {where}").fetchall()}
    finally:
        conn.close()


def _settings(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.bright_star_scatter_settings(conn)[0], store.bright_star_mass_limit(conn)
    finally:
        conn.close()


def _first_scatter(mysql_config, mass="12"):
    _seed(mysql_config)
    run_plan.scatter_bright_stars(_args(mysql_config, "--phenomenon-min-mass", mass))
    return _star_ids(mysql_config, f"initial_mass_sol >= {mass}"), _star_ids(mysql_config, f"initial_mass_sol < {mass}")


def test_redoing_the_luminosity_pass_keeps_the_mass_pass_and_rewrites_the_rest(mysql_config):
    heavy, light = _first_scatter(mysql_config)
    assert heavy and light
    run_plan.run_plan(_args(mysql_config, "--redo-scatters", "luminosity", "--bright-star-min-luminosity", "5000"))
    assert heavy <= _star_ids(mysql_config, "initial_mass_sol >= 12")
    assert not light & _star_ids(mysql_config, "initial_mass_sol < 12"), "the lighter stars are new rows"
    assert _settings(mysql_config) == (5000.0, 12.0)


def test_redoing_the_mass_pass_at_the_same_limit_keeps_the_luminosity_pass(mysql_config):
    heavy, light = _first_scatter(mysql_config)
    run_plan.run_plan(_args(mysql_config, "--redo-scatters", "mass"))
    assert light <= _star_ids(mysql_config, "initial_mass_sol < 12")
    assert not heavy & _star_ids(mysql_config, "initial_mass_sol >= 12")
    assert _settings(mysql_config) == (THRESHOLD, 12.0)


def test_a_new_mass_limit_redoes_the_luminosity_pass_with_it(mysql_config):
    _first_scatter(mysql_config)
    run_plan.run_plan(_args(mysql_config, "--redo-scatters", "mass", "--phenomenon-min-mass", "16"))
    assert _settings(mysql_config) == (THRESHOLD, 16.0)
    # Both passes ran at 16: the mass pass holds only stars of 16 or more, the luminosity pass only lighter ones.
    assert _star_ids(mysql_config, "initial_mass_sol >= 16") and _star_ids(mysql_config, "initial_mass_sol < 16")


def test_redoing_the_phenomena_leaves_the_stars_alone_and_keeps_the_nucleus(mysql_config):
    heavy, light = _first_scatter(mysql_config)
    run_plan.run_plan(_args(mysql_config, "--redo-scatters", "phenomena", "--compact-min-mass", "6"))
    assert _star_ids(mysql_config, "1 = 1") == heavy | light
    conn = store.get_connection(mysql_config)
    try:
        assert store.phenomenon_scatter_settings(conn)[1] == 6.0
        nuclei = conn.execute("SELECT COUNT(*) AS n FROM phenomenon_scatter WHERE kind = 'quasar'"
                              " OR (kind = 'black-hole' AND subtype = 'supermassive')").fetchone()["n"]
    finally:
        conn.close()
    assert nuclei == 1


def test_a_single_pass_with_no_scatter_yet_runs_both(mysql_config):
    _seed(mysql_config)
    run_plan.run_plan(_args(mysql_config, "--redo-scatters", "luminosity", "--phenomenon-min-mass", "12"))
    assert _star_ids(mysql_config, "initial_mass_sol >= 12") and _star_ids(mysql_config, "initial_mass_sol < 12")
    assert _settings(mysql_config) == (THRESHOLD, 12.0)


def _plan(*flags):
    _parser, parsers = generate_cli.build_parser()
    return parsers["plan"].parse_args(list(flags))


def test_the_redo_stages_list_every_scatter_with_the_ones_not_chosen_skipped():
    found = stages.plan_stages(_plan("--redo-scatters", "luminosity"))
    assert [stage.key for stage in found] == ["phenomena", "mass", "luminosity"]
    assert [stage.skip is not None for stage in found] == [True, True, False]
    assert all("not chosen to be redone" in stage.skip for stage in found[:2])
    assert [s.skip for s in stages.plan_stages(_plan("--redo-scatters", "mass", "luminosity", "phenomena"))] == [None] * 3
    only_mass = stages.plan_stages(_plan("--redo-scatters", "mass"))
    assert "redone if the mass limit changes" in only_mass[2].skip
    assert stages.plan_stages(_plan("--redo-scatters", "mass", "--phenomenon-min-mass", "14"))[2].skip is None


def test_the_redo_option_refuses_the_other_plan_modes():
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_plan_arguments(parser)
    for other in ("--phenomena-only", "--bright-stars-only", "--no-bright-stars"):
        with pytest.raises(SystemExit):
            generate_cli.validate_plan_args(parser.parse_args(["--redo-scatters", "mass", other]), parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["--redo-scatters", "everything"])
    assert tuple(parser.parse_args(["--redo-scatters", "mass", "phenomena"]).redo_scatters) == ("mass", "phenomena")
    assert tuning.REDO_SCATTERS == ("mass", "luminosity", "phenomena")
