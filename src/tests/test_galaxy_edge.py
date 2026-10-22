# tests/test_galaxy_edge.py

"""The Edge dial (GEN.605): how far out the plan walks the outline, in scale lengths (`12L`), parsecs, kiloparsecs or
megaparsecs. It replaces the plan's `--max-shell`, and Auto walks to the outline's own edge."""

import pytest

from planetgen.db import store
from planetgen.galaxy import dials, shape_draw
from tests.bughunt_support import mysql_argv, run_cli

EDGE = dials.lever("edge")


@pytest.mark.parametrize("text, pc", [("160pc", 160.0), ("2kpc", 2000.0), ("0.1Mpc", 100000.0)])
def test_an_absolute_edge_is_converted_to_parsecs(text, pc):
    assert EDGE.parse(text) == pc


def test_a_relative_edge_is_kept_until_the_galaxy_is_known():
    value = EDGE.parse("12L")
    assert isinstance(value, dials.RelativeLength) and (value.value, value.unit) == (12.0, "L")
    assert shape_draw.edge_limit_pc(value, 2600.0) == 12.0 * 2600.0
    assert shape_draw.edge_limit_pc(None, 2600.0) is None
    assert EDGE.format(value) == "12L" and EDGE.format(160.0) == "160pc"


@pytest.mark.parametrize("text", ["2L", "21L", "3pc", "300000pc", "1.1Mpc", "300"])
def test_an_edge_outside_its_range_or_without_a_unit_is_refused(text):
    with pytest.raises(Exception):
        EDGE.parse(text)


def test_the_old_option_is_gone():
    from planetgen.cli import generate as generate_cli
    _parser, parsers = generate_cli.build_parser()
    assert "--max-shell" not in {option for action in parsers["plan"]._actions for option in action.option_strings}


@pytest.mark.db
def test_a_plan_stops_at_the_edge_and_stores_it_as_a_pin(mysql_config, capsys):
    run_cli("plan", ["--edge=40pc", "--nuclear-amplitude", "0", "--seed",
                     "00112233445566778899AABBCCDDEEFF"] + mysql_argv(mysql_config))
    out = capsys.readouterr().out
    assert "WARNING: the galactic plane still qualified at the --edge (40pc)" in out
    conn = store.get_connection(mysql_config)
    try:
        stored = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    assert stored.outer_shell_index == 10 and stored.pins["edge"] == 40.0


@pytest.mark.db
def test_an_edge_in_scale_lengths_follows_the_disc(mysql_config, capsys):
    run_cli("plan", ["--edge=3L", "--disk-scale-length=100pc", "--nuclear-amplitude", "0",
                     "--seed", "00112233445566778899AABBCCDDEEFF"] + mysql_argv(mysql_config))
    conn = store.get_connection(mysql_config)
    try:
        stored = store.get_galaxy_shape(conn)
    finally:
        conn.close()
    assert stored.outer_shell_index <= 75 and stored.pins["edge"] == 300.0     # 3 L of 100 pc, in 4 pc shells
