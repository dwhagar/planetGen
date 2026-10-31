# tests/test_scatter.py

"""The scatter step (GEN.1089): giant-free and giant-rich discs; the Solar System is a test target, never an input."""

import math

from planetgen.physics import scatter

STAR = 1.0


def _shares(giant_mass=0.0, giant_au=30.0, radius=3.9, migration=0.0, inner=15.0, edge=35.0, instability=0.0):
    return scatter.scatter_shares(STAR, giant_mass, giant_au, radius, migration, inner, edge, instability)


def test_a_disc_with_no_giant_stays_cold():
    shares = _shares(0.0)
    assert shares["cold_classical"] == 1.0
    assert sum(shares.values()) == 1.0
    assert scatter.scattered_share(shares) == 0.0
    assert "no giant" in scatter.history_line(shares, 30.0, 0.0)


def test_the_shares_always_sum_to_one_and_stay_in_range():
    for mass, au, migration, instability in ((17.0, 30.0, 12.0, 1.0), (318.0, 5.0, 0.0, 0.0), (5.0, 10.0, 3.0, 0.5), (900.0, 60.0, 30.0, 1.0)):
        shares = _shares(mass, au, 3.9 if mass < 100 else 11.0, migration, 3.0, 80.0, instability)
        assert math.isclose(sum(shares.values()), 1.0, abs_tol=1e-9)
        assert all(0.0 <= value <= 1.0 for value in shares.values())
        assert set(shares) == set(scatter.POPULATIONS)


def test_a_giant_that_never_moved_scatters_only_its_zone():
    shares = _shares(17.1, 30.0, 3.9, 0.0, 15.0, 50.0)
    assert shares["cold_classical"] > 0.8
    assert shares["resonant"] == 0.0 and shares["hot_classical"] == 0.0
    assert 0.0 < scatter.scattered_share(shares) < 0.2
    assert "never moved" in scatter.history_line(shares, 30.0, 0.0)


def test_a_migrating_giant_scatters_most_of_the_disc_it_crossed():
    # a Neptune-sized giant that moved from 18 to 30 AU through a disc of 15 to 35 AU (test target; the line of the design: "scattered 82 percent")
    shares = _shares(17.1, 30.0, 3.9, 12.0, 15.0, 35.0, instability=1.0)
    assert 0.7 < scatter.scattered_share(shares) < 0.95
    assert shares["resonant"] > 0.0 and shares["hot_classical"] > 0.0 and shares["cold_classical"] < 0.2
    assert "moved from 18 to 30 AU" in scatter.history_line(shares, 30.0, 12.0)
    # more migration, less cold disc
    assert _shares(17.1, 30.0, 3.9, 6.0, 15.0, 35.0)["cold_classical"] > shares["cold_classical"]


def test_a_heavier_giant_has_a_wider_zone_and_ejects_more():
    assert scatter.zone_half_width_au(5.0, 318.0, STAR) > scatter.zone_half_width_au(5.0, 17.0, STAR)
    light = _shares(17.1, 30.0, 3.9, 12.0, 15.0, 35.0)
    heavy = _shares(318.0, 30.0, 11.2, 12.0, 15.0, 35.0)
    assert heavy["ejected"] / scatter.scattered_share(heavy) >= light["ejected"] / scatter.scattered_share(light) - 1e-9
    assert scatter.scattered_share(heavy) > scatter.scattered_share(light)


def test_the_resonance_positions_and_the_safronov_number():
    assert math.isclose(scatter.resonance_au(30.0, 3, 2) / 30.0, 1.31, abs_tol=0.01)
    assert math.isclose(scatter.resonance_au(30.0, 2, 1) / 30.0, 1.587, abs_tol=0.005)
    # a Jupiter-mass planet at 5 AU and a Neptune-mass one at 30 AU are both above 1: scattered bodies are thrown out
    assert scatter.safronov_number(318.0, 11.2, 5.2, STAR) > 1.0
    assert scatter.safronov_number(17.1, 3.9, 30.0, STAR) > 1.0
    # a super-Earth with no envelope at 1 AU keeps its scattered bodies
    assert scatter.safronov_number(5.0, 1.5, 1.0, STAR) < 1.0


def test_the_disc_fraction():
    assert math.isclose(scatter.disc_fraction(0.0, 100.0, 10.0, 40.0), 1.0)
    assert scatter.disc_fraction(50.0, 60.0, 10.0, 40.0) == 0.0
    # more mass near the star for a falling surface density
    assert scatter.disc_fraction(10.0, 20.0, 10.0, 40.0) > scatter.disc_fraction(30.0, 40.0, 10.0, 40.0)
    assert math.isclose(scatter.disc_fraction(10.0, 25.0, 10.0, 40.0, slope=2.0), math.log(2.5) / math.log(4.0))


def test_the_split_of_the_swept_mass_sums_to_one():
    assert math.isclose(sum(share for _, share in scatter.SWEPT_SPLIT), 1.0)
    assert {key for key, _ in scatter.SWEPT_SPLIT} <= set(scatter.POPULATIONS)
