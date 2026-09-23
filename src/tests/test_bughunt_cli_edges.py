# tests/test_bughunt_cli_edges.py

"""
Tier 1 bug-hunt coverage: every documented incompatible-option combo from
`README.md`'s "Note on Incompatible Options", plus numeric boundary
values for every option that takes one -- driven through the real CLI
entry point (`generate.main()`, via `sys.argv`, same convention
`test_galaxy_gen.py` already uses) so a validation regression (an
incompatible combo that stops being rejected, or a numeric boundary that
starts crashing deep in generation instead of failing cleanly at argument
-parsing time) is caught at the actual user-facing surface, not just in
the validator functions' own unit tests (`test_appconfig.py` covers
`SystemConfig` itself; this file is the CLI layer above it).

Every "should be rejected" case asserts `SystemExit` -- `parser.error()`
and `build_system_config`'s own explicit `raise SystemExit(1)` both exit
the process this way; that's the correct, already-intended outcome, not a
bug. What this file actually hunts for is a combo that *stops* exiting
(silently accepted, then either crashes unpredictably deep in generation
or produces a system violating the very constraint the option combo was
supposed to prevent) -- so `pytest.raises(SystemExit)` is the Tier 1
assertion throughout, deliberately not narrowed to a specific exit code
(argparse's own `parser.error` uses 2; `build_system_config`'s uses 1).
"""

import pytest

from tests.bughunt_support import mysql_argv, run_cli


# --- README "Note on Incompatible Options" ------------------------------

@pytest.mark.parametrize("argv", [
    ["-planets", "+moons"],
    ["-planets", "+max_planets"],
    ["-planets", "+habitable_world"],
    ["--star-type", "G2V", "+large_star"],
    ["+intelligent_life", "-habitable_world"],
    ["-intelligent_life", "-habitable_world"],
    ["--num-orbits", "-1"],
    ["--num-orbits", "3", "-planets"],
    ["--flavor-chance-system", "-0.1"],
    ["--flavor-chance-system", "1.1"],
    ["--flavor-chance-planet", "-0.1"],
    ["--flavor-chance-planet", "1.1"],
])
def test_system_incompatible_combo_is_rejected(argv, mysql_config):
    with pytest.raises(SystemExit):
        run_cli("system", argv + mysql_argv(mysql_config))


def test_system_habitable_and_belt_forbidding_large_star_is_rejected(mysql_config):
    """The one incompatible combo build_system_config itself rejects
    (SystemExit(1), not parser.error) rather than add_system_arguments."""
    with pytest.raises(SystemExit):
        run_cli("system", ["+habitable_world", "+asteroid_belt", "-large_star"] + mysql_argv(mysql_config))


def test_system_habitable_and_belt_with_large_star_forced_is_accepted(mysql_config):
    """Sanity check the *complement*: the same combo with +large_star (or
    no large_star opinion) is legitimate and must still generate."""
    run_cli("system", ["+habitable_world", "+asteroid_belt", "+large_star"] + mysql_argv(mysql_config))


@pytest.mark.parametrize("argv", [
    ["--num-systems", "0"],
    ["--num-systems", "-1"],
    ["--min-habitable", "-1"],
])
def test_sector_incompatible_or_invalid_value_is_rejected(argv, mysql_config):
    with pytest.raises(SystemExit):
        run_cli("sector", argv + mysql_argv(mysql_config))


def test_sector_min_habitable_exceeding_num_systems_is_rejected(mysql_config):
    with pytest.raises(SystemExit):
        run_cli("sector", ["--num-systems", "3", "--min-habitable", "4"] + mysql_argv(mysql_config))


def test_sector_min_habitable_with_forbidden_habitable_world_is_rejected(mysql_config):
    with pytest.raises(SystemExit):
        run_cli("sector", ["--min-habitable", "1", "-habitable_world"] + mysql_argv(mysql_config))


def test_sector_star_type_with_large_star_is_rejected(mysql_config):
    with pytest.raises(SystemExit):
        run_cli("sector", ["--star-type", "G2V", "+large_star"] + mysql_argv(mysql_config))


# --- Numeric boundary values that SHOULD be accepted --------------------

def test_system_num_orbits_zero_is_accepted(mysql_config):
    """--num-orbits 0 is explicitly documented as valid ("zero or a
    positive integer") -- a system with no orbital slots at all."""
    run_cli("system", ["--num-orbits", "0"] + mysql_argv(mysql_config))


def test_system_flavor_chance_boundary_values_are_accepted(mysql_config):
    run_cli("system", ["--flavor-chance-system", "0.0", "--flavor-chance-planet", "1.0"] + mysql_argv(mysql_config))


def test_sector_min_habitable_equal_to_num_systems_is_accepted(mysql_config):
    """min_habitable == num_systems (the exact boundary, not exceeding
    it) is legitimate -- every system in the sector is habitable."""
    run_cli("sector", ["--num-systems", "2", "--min-habitable", "2"] + mysql_argv(mysql_config))


def test_system_num_orbits_large_is_accepted_or_clean_error(mysql_config):
    """A very large explicit --num-orbits isn't itself documented as
    invalid -- the generator either honors it (however slowly) or fails
    cleanly (e.g. no room in any zone), but must not crash with a raw,
    unexplained exception."""
    try:
        run_cli("system", ["--num-orbits", "1000", "+large_star"] + mysql_argv(mysql_config))
    except SystemExit:
        pass  # a clean, argparse/log.error-mediated rejection is acceptable


# --- --version and no-argv smoke checks ----------------------------------

def test_version_flag_exits_cleanly(capsys):
    import sys

    old_argv = sys.argv
    try:
        sys.argv = ["generate.py", "--version"]
        with pytest.raises(SystemExit) as exc_info:
            import generate
            generate.main()
        assert exc_info.value.code in (0, None)
    finally:
        sys.argv = old_argv
    out = capsys.readouterr().out
    assert out.strip(), "--version printed nothing"


def test_no_subcommand_exits_with_usage(capsys):
    import sys

    old_argv = sys.argv
    try:
        sys.argv = ["generate.py"]
        with pytest.raises(SystemExit):
            import generate
            generate.main()
    finally:
        sys.argv = old_argv
