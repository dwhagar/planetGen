# tests/test_star_count_options.py

"""
GEN.267: the star count settings. `--stars N` forces the count of one system and `--max-stars N` raises or lowers the cap
(`SystemConfig.STARS`, `MAX_STARS`); the multi_star prevalence moves the chance of five or more stars. Each is carried in
the config for the member-tree host (GEN.630) to read, and what the cap or the pair options rule out is refused with the
reason. No database: these parse options and build configs.
"""

import json

import pytest

from planetgen import tuning
from planetgen.cli import generate as generate_cli
from planetgen.generation import multistar, prevalence, run_system
from planetgen.generation.config import SystemConfig
from planetgen.physics import imf
from planetgen.web import generate_page


def _parse(*argv):
    parser, command_parsers = generate_cli.build_parser()
    args = parser.parse_args(list(argv))
    command_parser = command_parsers[args.command]
    if args.command == "system":
        generate_cli.validate_system_args(args, command_parser)
    else:
        generate_cli.validate_shared_generation_args(args, command_parser)
    return args


def _refused(capsys, *argv):
    with pytest.raises(SystemExit):
        _parse(*argv)
    return capsys.readouterr().err


def test_a_config_leaves_the_count_to_chance_under_the_default_cap():
    config = SystemConfig()
    assert config.STARS is None
    assert config.MAX_STARS == tuning.MAX_STARS_DEFAULT == multistar.MAX_MEMBERS == 12


def test_the_count_and_cap_round_trip_through_the_config_dict():
    config = SystemConfig()
    config.STARS, config.MAX_STARS = 5, 20
    config.PREVALENCE = {"multi_star": 150.0}
    again = SystemConfig.from_dict(json.loads(json.dumps(config.to_dict())))
    assert (again.STARS, again.MAX_STARS, again.PREVALENCE) == (5, 20, {"multi_star": 150.0})


def test_the_system_command_carries_the_count_and_cap_into_the_config():
    config = run_system.build_system_config(_parse("system", "--stars", "4", "--max-stars", "9"))
    assert (config.STARS, config.MAX_STARS) == (4, 9)
    plain = run_system.build_system_config(_parse("system"))
    assert (plain.STARS, plain.MAX_STARS) == (None, 12)


def test_a_system_file_sets_the_count_and_the_command_line_wins(tmp_path):
    path = tmp_path / "system.json"
    path.write_text(json.dumps({"stars": 3, "max_stars": 6}))
    config = run_system.build_system_config(_parse("system", "-f", str(path)))
    assert (config.STARS, config.MAX_STARS) == (3, 6)
    config = run_system.build_system_config(_parse("system", "-f", str(path), "--stars", "5"))
    assert (config.STARS, config.MAX_STARS) == (5, 6)


@pytest.mark.parametrize("argv, reason", [
    (["--stars", "0"], "whole number of at least 1"),
    (["--stars", "-2"], "whole number of at least 1"),
    (["--stars", "13"], "more than the cap of 12 stars; raise the cap with --max-stars"),
    (["--stars", "9", "--max-stars", "8"], "more than the cap of 8 stars"),
    (["--stars", "65", "--max-stars", "64"], "more than any system may hold (64)"),
    (["--max-stars", "0"], "max_stars 0 must be a whole number from 1 to 64"),
    (["--max-stars", "65"], "max_stars 65 must be a whole number from 1 to 64"),
    (["--stars", "1", "+binary_system"], "single star, which +binary_system contradicts"),
    (["--stars", "1", "+wide_binary"], "no companion for +wide_binary"),
    (["--stars", "2", "-binary_system"], "more than one star, which -binary_system contradicts"),
    (["--stars", "7", "-binary_system"], "more than one star, which -binary_system contradicts"),
])
def test_a_count_the_cap_or_the_pair_options_rule_out_is_refused_with_the_reason(argv, reason, capsys):
    assert reason in _refused(capsys, "system", *argv)


@pytest.mark.parametrize("argv", [
    ["--stars", "1"], ["--stars", "1", "-binary_system"], ["--stars", "1", "-wide_binary"], ["--stars", "2", "+binary_system"],
    ["--stars", "3", "+binary_system"], ["--stars", "12"], ["--stars", "20", "--max-stars", "20"], ["--stars", "64", "--max-stars", "64"],
    ["--max-stars", "1"],
])
def test_a_count_that_fits_is_accepted(argv):
    _parse("system", *argv)


def test_a_system_file_that_contradicts_itself_is_refused_with_the_reason(tmp_path, capsys):
    path = tmp_path / "system.json"
    path.write_text(json.dumps({"stars": 20}))
    assert "more than the cap of 12 stars" in _refused(capsys, "system", "-f", str(path))
    path.write_text(json.dumps({"stars": 20, "max_stars": 30}))
    _parse("system", "-f", str(path))
    # The command line's cap replaces the file's, so the file's count is judged against it.
    assert "more than the cap of 8 stars" in _refused(capsys, "system", "-f", str(path), "--max-stars", "8")


def test_a_config_built_from_a_file_the_parser_did_not_see_is_refused_by_the_run():
    config = SystemConfig()
    config.STARS, config.BINARY_SYSTEM = 1, True
    assert "+binary_system" in run_system.single_system_conflict(config)
    config = SystemConfig()
    config.STARS, config.MAX_STARS = 3, 2
    assert "cap of 2" in run_system.single_system_conflict(config)
    assert run_system.single_system_conflict(SystemConfig()) is None


@pytest.mark.parametrize("command", ["sector", "galaxy"])
def test_a_sector_or_galaxy_run_refuses_a_forced_count_but_takes_the_cap(command, capsys):
    with pytest.raises(SystemExit):
        _parse(command, "--stars", "3")
    assert "--stars forces one system" in capsys.readouterr().err
    assert _parse(command, "--max-stars", "20").max_stars == 20
    assert "max_stars 99 must be" in _refused(capsys, command, "--max-stars", "99")


def test_a_sector_run_gives_every_system_the_cap():
    args = _parse("sector", "--max-stars", "20", "--prevalence", "multi_star=50")
    args.system_file = args.num_orbits = args.name = None  # what the sector command sets before it builds configs
    config = run_system.build_system_config(args)
    assert config.MAX_STARS == 20 and config.PREVALENCE == {"multi_star": 50.0}


def test_the_multi_star_prevalence_scales_the_five_or_more_columns():
    config = SystemConfig()
    assert prevalence.multi_star_scale(config) == 1.0
    config.PREVALENCE = {"multi_star": 50}
    assert prevalence.multi_star_scale(config) == 1.5
    config.PREVALENCE = {"multi_star": -100}
    assert prevalence.multi_star_scale(config) == 0.0
    plain = multistar.star_count_weights(1.0)
    more = multistar.star_count_weights(1.0, 1.5)
    assert more[:4] == plain[:4] and all(m > p for m, p in zip(more[4:], plain[4:]))
    assert all(w == 0 for w in multistar.star_count_weights(1.0, 0.0)[4:])


def test_the_usual_share_of_five_or_more_stars_is_the_tables_over_the_mass_function():
    """`tuning.PREVALENCE_DRAW_SHARES["multi_star"]` is not measured from generated systems yet (GEN.630): it is the
    `star-count` table's five-or-more columns averaged over the primary mass function (the placed star's, GEN.928)."""
    import numpy as np
    low, high = tuning.IMF_BREAKS_SOL[0], tuning.IMF_BREAKS_SOL[-1]
    grid = np.geomspace(low, high, 400)
    total = five = 0.0
    for a, b in zip(grid[:-1], grid[1:]):
        number = imf.integral(0.0, a, b, imf.PRIMARY)
        weights = multistar.star_count_weights((a * b) ** 0.5)
        total += number
        five += number * sum(weights[4:]) / sum(weights)
    assert tuning.PREVALENCE_DRAW_SHARES["multi_star"] == pytest.approx(five / total, rel=0.03)
    assert prevalence.USUAL_SHARES["multi_star"] == tuning.PREVALENCE_DRAW_SHARES["multi_star"]


def test_the_prevalence_option_takes_multi_star():
    assert generate_cli.prevalence_setting("multi-star=+50") == ("multi_star", 50.0)
    with pytest.raises(Exception):
        generate_cli.prevalence_setting("multi_star=-101")


def test_the_generate_page_offers_the_share_and_the_cap():
    names = {name: usual for name, _label, usual in generate_page.PREVALENCE_FIELDS}
    assert names["prevalence_multi_star"] == "0.46"
    argv = generate_page.prevalence_argv({"prevalence_multi_star": "0.92"})
    assert argv[:1] == ["--prevalence"] and argv[1].startswith("multi_star=") and float(argv[1].split("=")[1]) == pytest.approx(100, abs=3)
    assert generate_page.prevalence_argv({}) == []
    assert generate_page.prevalence_argv({"max_stars": "12"}) == []
    assert generate_page.prevalence_argv({"max_stars": "24"}) == ["--max-stars", "24"]


@pytest.mark.parametrize("value", ["0", "65", "2.5", "many"])
def test_the_generate_page_refuses_a_cap_outside_one_to_the_ceiling(value):
    with pytest.raises(generate_page.FormError):
        generate_page.prevalence_argv({"max_stars": value})
