# tests/test_math_check.py

"""
The math check gate (TEST.63 to TEST.67): every check in
`stellarObjects.mathCheck` as its own test, marked `mathcheck` so
`conftest.py` runs them first and stops the suite if one fails, plus a few
tests that the gate itself can tell a broken function from a working one.
"""

import math
import random

import pytest

from stellarObjects import mathCheck

pytestmark = pytest.mark.mathcheck

CHECKS = mathCheck.all_checks()


@pytest.mark.parametrize("check", CHECKS, ids=[c.name for c in CHECKS])
def test_math_check(check):
    result = mathCheck.run_check(check)
    assert result.passed, result.describe()


def test_check_names_are_unique_and_documented():
    names = [c.name for c in CHECKS]
    assert len(names) == len(set(names))
    for check in CHECKS:
        assert check.group in ("reference", "invariant", "distribution")
        assert check.source and check.function
        assert check.mode in ("rel", "abs", "max", "min")


def test_every_group_is_present():
    groups = {c.group for c in CHECKS}
    assert groups == {"reference", "invariant", "distribution"}


def test_whole_run_is_fast():
    results = mathCheck.run_all()
    assert sum(r.seconds for r in results) < 5.0, mathCheck.format_report(results, verbose=True)


def test_run_leaves_global_random_alone():
    random.seed(12345)
    expected = random.random()
    random.seed(12345)
    mathCheck.run_all()
    assert random.random() == expected


def test_chi_square_p_value_matches_known_values():
    # chi-square = 3.841 with 1 degree of freedom is p = 0.05; 2 observed
    # bins of 60/40 against 50/50 gives chi-square 4.0, p = 0.0455.
    assert mathCheck._regularized_gamma_q(0.5, 3.841459 / 2) == pytest.approx(0.05, rel=1e-4)
    assert mathCheck.chi_square_p_value([60, 40], [1, 1]) == pytest.approx(0.0455003, rel=1e-4)
    # 10 degrees of freedom at chi-square 18.307 is p = 0.05.
    assert mathCheck._regularized_gamma_q(5.0, 18.307 / 2) == pytest.approx(0.05, rel=1e-3)
    assert mathCheck.chi_square_p_value([100, 100, 100], [1, 1, 1]) == pytest.approx(1.0)


def test_chi_square_merges_rare_bins():
    # The 1-in-10,000 bin would expect 0.1 draws; merged, the test stays valid.
    assert mathCheck.chi_square_p_value([500, 499, 1], [0.5, 0.4999, 0.0001]) > 0.5


def test_a_broken_sampler_fails(monkeypatch):
    """A sampler whose every value is in range but whose shares are wrong
    (ages piling up young) is caught by the distribution check."""
    from stellarObjects import stellarEvolution

    def skewed(age_bias=None, rng=random, population=None):
        return 10.0 * rng.random() ** 2

    monkeypatch.setattr(stellarEvolution, "sample_star_age_gy", skewed)
    check = next(c for c in CHECKS if c.name == "star_ages_uniform")
    assert not mathCheck.run_check(check).passed


def test_a_wrong_constant_fails(monkeypatch):
    from stellarObjects import physical_constants

    monkeypatch.setattr(physical_constants, "SPEED_OF_LIGHT_M_S", 2.998e8)
    check = next(c for c in CHECKS if c.name == "speed_of_light_consistent")
    assert not mathCheck.run_check(check).passed


def test_an_exception_is_a_failure_not_a_crash():
    def boom():
        raise ZeroDivisionError("division by zero")

    check = mathCheck.Check("boom", "invariant", "boom()", boom, 0.0, 0.0, "test", mode="max")
    result = mathCheck.run_check(check)
    assert not result.passed
    assert "ZeroDivisionError" in result.describe()
    assert "FAILED" in mathCheck.format_report([result])


def test_nan_never_passes():
    check = mathCheck.Check("nan", "invariant", "nan", lambda: math.nan, 0.0, 0.0, "test", mode="max")
    assert not mathCheck.run_check(check).passed


def test_command_line_exit_codes(capsys, monkeypatch):
    assert mathCheck.main([]) == 0
    assert "Math check passed" in capsys.readouterr().out
    bad = mathCheck.Check("bad", "invariant", "bad", lambda: 1.0, 0.0, 0.0, "test", mode="max")
    monkeypatch.setattr(mathCheck, "all_checks", lambda: [bad])
    assert mathCheck.main([]) == 1
    assert "FAIL bad" in capsys.readouterr().out
