# tests/test_star_mix.py

"""ADM.45: the star mix (single, close and wide shares of systems) totals 100%
and turns into the binary and wide-pair prevalences."""

import argparse

import pytest

from planetgen.cli import generate as generate_cli
from planetgen.generation import prevalence


def test_the_usual_mix_gives_no_change():
    mix = prevalence.usual_star_mix()
    assert sum(mix.values()) == pytest.approx(1.0)
    got = prevalence.star_mix_prevalences(*(mix[kind] * 100 for kind in prevalence.STAR_MIX))
    assert got["binary_system"] == pytest.approx(0.0, abs=1e-9) and got["wide_binary"] == pytest.approx(0.0, abs=1e-9)


def test_a_mix_becomes_binary_and_wide_prevalences():
    got = prevalence.star_mix_prevalences(50, 30, 20)
    # Half the systems are binary (0.303 usual); two fifths of those wide (0.499 usual).
    assert got["binary_system"] == pytest.approx((0.5 / 0.303 - 1) * 100)
    assert got["wide_binary"] == pytest.approx((0.4 / 0.499 - 1) * 100)


def test_all_single_leaves_nothing_binary():
    got = prevalence.star_mix_prevalences(100, 0, 0)
    assert got["binary_system"] == -100.0


@pytest.mark.parametrize("shares", [(50, 30, 30), (50, 30, 10), (110, -5, -5)])
def test_a_mix_that_is_not_100_is_refused(shares):
    with pytest.raises(ValueError):
        prevalence.star_mix_prevalences(*shares)


def _parse(*argv):
    parser = argparse.ArgumentParser(prefix_chars='-+')
    generate_cli.add_shared_generation_options(parser)
    args = parser.parse_args(list(argv))
    return args, parser


def test_the_cli_takes_a_star_mix_and_rejects_a_bad_total():
    args, parser = _parse("--star-mix", "50", "30", "20")
    generate_cli.validate_shared_generation_args(args, parser)
    assert {feature for feature, _ in args.prevalence} == {"binary_system", "wide_binary"}
    bad, parser = _parse("--star-mix", "50", "30", "30")
    with pytest.raises(SystemExit):
        generate_cli.validate_shared_generation_args(bad, parser)
    clash, parser = _parse("--star-mix", "50", "30", "20", "--prevalence", "binary_system=10")
    with pytest.raises(SystemExit):
        generate_cli.validate_shared_generation_args(clash, parser)
