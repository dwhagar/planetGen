# tests/test_math_check.py

"""
The math check gate (TEST.63 to TEST.67): every check in
`planetgen.physics.mathcheck` as its own test, marked `mathcheck` so
`conftest.py` runs them first and stops the suite if one fails, plus a few
tests that the gate itself can tell a broken function from a working one.
"""

import math
import time

import pytest

from planetgen.util import draw
from planetgen.physics import mathcheck as mathCheck

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
    # The process's CPU time, not the clock: on a loaded machine the clock also counts the time the process waited
    # for a core, which the run does not control (TEST.170).
    started = time.process_time()
    results = mathCheck.run_all()
    cpu_seconds = time.process_time() - started
    assert cpu_seconds < 5.0, f"{cpu_seconds:.2f} s of CPU\n" + mathCheck.format_report(results, verbose=True)


def test_run_leaves_global_random_alone():
    draw.set_run_seed(12345)
    expected = draw.random()
    draw.set_run_seed(12345)
    mathCheck.run_all()
    assert draw.random() == expected


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
    from planetgen.physics import stellar_evolution

    def skewed(age_bias=None, rng=draw, population=None):
        return 10.0 * rng.random() ** 2

    monkeypatch.setattr(stellar_evolution, "sample_star_age_gy", skewed)
    check = next(c for c in CHECKS if c.name == "star_ages_uniform")
    assert not mathCheck.run_check(check).passed


def test_a_wrong_constant_fails(monkeypatch):
    from planetgen.physics import constants

    monkeypatch.setattr(constants, "SPEED_OF_LIGHT_M_S", 2.998e8)
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


def test_a_broken_wedge_rule_fails(monkeypatch):
    """A wedge count that is not 6*m with m made of 2, 3, 5 and 7 (an 11 in it) is caught."""
    from planetgen.galaxy import spherical

    real = spherical.wedge_count
    monkeypatch.setattr(spherical, "wedge_count", lambda shell, band: real(shell, band) * 11 // 6 * 6 if shell == 9 and band == 0 else real(shell, band))
    check = next(c for c in CHECKS if c.name == "wedge_count_rule")
    assert not mathCheck.run_check(check).passed


def test_a_sector_found_in_the_wrong_cell_fails(monkeypatch):
    """A lookup that puts a point one wedge over from its true cell is caught."""
    from planetgen.galaxy import spherical

    real = spherical.sector_address_at

    def shifted(position, edge_pc):
        shell, band, wedge = real(position, edge_pc)
        return (shell, band, wedge + 1) if shell == 5 and spherical.wedge_count(shell, band) > 1 else (shell, band, wedge)

    monkeypatch.setattr(spherical, "sector_address_at", shifted)
    for name in ("sector_address_round_trip", "point_lies_in_its_sector"):
        check = next(c for c in CHECKS if c.name == name)
        assert not mathCheck.run_check(check).passed, name


def test_shells_that_stop_tiling_fail(monkeypatch):
    """A cell volume that loses a sliver of every cell is caught by the shell volume check."""
    from planetgen.galaxy import spherical

    real = spherical.cell_volume
    monkeypatch.setattr(spherical, "cell_volume", lambda *a, **k: real(*a, **k) * (1 - 1e-9))
    check = next(c for c in CHECKS if c.name == "shell_cells_fill_shell")
    assert not mathCheck.run_check(check).passed
