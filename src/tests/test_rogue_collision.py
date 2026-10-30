# tests/test_rogue_collision.py

"""The collision of two rogue planets (GEN.1099, `physics/rogue_collision.py`): the rate of about one in the life of the Galaxy, the shares of the collisions by speed (design section 3.1),
the disruption speed of a pair, the largest remnant and the order of the regimes, the speed of the pieces and the shock."""

import math

import pytest

from planetgen.physics import constants, rogue_collision as rc

EARTH = constants.EARTH_MASS_TO_KG
EARTH_RADIUS_M = 6.371e6
KM = 1000.0


def test_two_earths_have_the_focusing_and_escape_speed_of_the_design():
    v_esc = rc.escape_speed_touching_ms(EARTH, EARTH, EARTH_RADIUS_M, EARTH_RADIUS_M)
    assert v_esc / KM == pytest.approx(11.2, abs=0.1)
    sigma_40 = rc.cross_section_m2(EARTH, EARTH, EARTH_RADIUS_M, EARTH_RADIUS_M, 40 * KM)
    assert sigma_40 / (math.pi * (2 * EARTH_RADIUS_M) ** 2) == pytest.approx(1.078, abs=0.005)       # "at 40 km/s focusing is 1.08"


@pytest.mark.parametrize("kms,rate", [(10, 9.6e-24), (20, 1.1e-23), (40, 1.8e-23), (80, 3.5e-23)])
def test_the_rate_per_rogue_follows_the_designs_table(kms, rate):
    found = rc.rate_per_rogue_per_year(0.78, EARTH, EARTH, EARTH_RADIUS_M, EARTH_RADIUS_M, kms * KM)
    assert found == pytest.approx(rate, rel=0.35)           # the design's table is rounded and uses 0.78 per pc3 for Earth-size pairs


def test_the_galaxy_has_about_one_collision_in_its_life():
    disc = rc.galaxy_rate_per_year(6.5, 3.4e10, EARTH, EARTH, EARTH_RADIUS_M, EARTH_RADIUS_M, 40 * KM)
    bulge = rc.galaxy_rate_per_year(6.5, 5.1e10, EARTH, EARTH, EARTH_RADIUS_M, EARTH_RADIUS_M, 120 * KM)
    per_ten_gyr = (disc + bulge) * 1e10
    assert 0.1 < per_ten_gyr < 10                            # design 2.3: 0.9 for Earth-size pairs, range 0.1 to 10
    assert (disc + bulge) == pytest.approx(8.9e-11, rel=0.1)


@pytest.mark.parametrize("dispersion,expected", [
    (25, (0.03, 0.17, 0.60, 0.20)), (50, (0.002, 0.016, 0.15, 0.83)),
])
def test_the_collisions_fall_in_the_speed_bins_of_the_design(dispersion, expected):
    v_esc = rc.escape_speed_touching_ms(EARTH, EARTH, EARTH_RADIUS_M, EARTH_RADIUS_M)
    found = rc.speed_fractions(dispersion * KM, v_esc, [15 * KM, 30 * KM, 60 * KM])
    assert sum(found) == pytest.approx(1.0)
    for got, want in zip(found, expected):
        assert got == pytest.approx(want, abs=0.006)


def test_the_disruption_speed_scales_with_the_pair():
    assert rc.combined_radius_m(2 * EARTH) / KM == pytest.approx(14180, rel=0.005)
    assert rc.disruption_speed_ms(EARTH, EARTH) / KM == pytest.approx(22.6, abs=0.1)
    assert rc.disruption_speed_ms(0.3 * EARTH, 0.3 * EARTH) / KM == pytest.approx(15.2, abs=0.1)
    assert rc.disruption_speed_ms(0.1 * EARTH, 0.1 * EARTH) / KM == pytest.approx(10.5, abs=0.1)
    assert rc.disruption_energy_j_kg(EARTH, EARTH) == pytest.approx(6.4e7, rel=0.02)
    assert rc.combined_escape_speed_ms(2 * EARTH) / KM == pytest.approx(10.6, abs=0.1)


def test_a_small_projectile_needs_more_speed_than_an_equal_one_of_the_same_total_mass():
    equal = rc.disruption_speed_ms(0.5 * EARTH, 0.5 * EARTH)
    small = rc.disruption_speed_ms(0.9 * EARTH, 0.1 * EARTH)
    assert small > 2 * equal
    assert rc.disruption_speed_ms(0.1 * EARTH, 0.9 * EARTH) == rc.disruption_speed_ms(0.9 * EARTH, 0.1 * EARTH)       # the target is the larger


def test_the_rubble_window_runs_from_the_threshold_to_a_tenth_of_the_mass():
    low, high = rc.rubble_speed_window_ms(EARTH, EARTH)
    assert (low / KM, high / KM) == pytest.approx((22.6, 30.3), abs=0.15)
    assert rc.largest_remnant_fraction(rc.specific_impact_energy_j_kg(EARTH, EARTH, low) / rc.disruption_energy_j_kg(EARTH, EARTH)) == pytest.approx(0.5)
    assert rc.largest_remnant_fraction(rc.specific_impact_energy_j_kg(EARTH, EARTH, high) / rc.disruption_energy_j_kg(EARTH, EARTH)) == pytest.approx(0.1, abs=0.001)


def test_the_largest_remnant_follows_the_universal_law_and_the_power_law():
    assert rc.largest_remnant_fraction(0.0) == 1.0
    assert rc.largest_remnant_fraction(1.0) == 0.5
    assert rc.largest_remnant_fraction(1.8) == pytest.approx(0.1)
    assert rc.largest_remnant_fraction(1.8 + 1e-9) == pytest.approx(0.1, abs=1e-6)       # continuous at the join
    assert rc.largest_remnant_fraction(1.8 * 4) == pytest.approx(0.1 * 4 ** -1.5)
    values = [rc.largest_remnant_fraction(r) for r in (0.1, 0.5, 1.0, 1.5, 1.8, 3.0, 10.0, 100.0)]
    assert values == sorted(values, reverse=True) and values[-1] > 0.0
    assert rc.second_remnant_fraction(0.5) == pytest.approx(0.15 * 0.5 / (2 * 2.85))
    assert rc.second_remnant_fraction(1.0) == 0.0


def test_the_regimes_come_in_order_of_energy():
    names = [rc.regime(r) for r in (0.05, 0.5, 1.0, 1.5, 1.8, 5.0)]
    assert names == ["merger", "partial_erosion", "catastrophic", "catastrophic", "super_catastrophic", "super_catastrophic"]
    assert [rc.REGIMES.index(n) for n in names if n != "hit_and_run"] == sorted(rc.REGIMES.index(n) for n in names)
    assert rc.regime(0.5, grazing=True) == "hit_and_run" and rc.regime(2.0, grazing=True) == "super_catastrophic"
    with pytest.raises(rc.RogueCollisionError):
        rc.regime(-1.0)


def test_an_oblique_impact_needs_more_energy():
    r1, r2 = EARTH_RADIUS_M, EARTH_RADIUS_M
    head_on = rc.oblique_disruption_energy_j_kg(EARTH, EARTH, r1, r2, 0.0)
    assert head_on == pytest.approx(rc.disruption_energy_j_kg(EARTH, EARTH))
    assert rc.interacting_share(r1, r2, 0.0) == 1.0
    energies = [rc.oblique_disruption_energy_j_kg(EARTH, EARTH, r1, r2, b) for b in (0.0, 0.35, 0.7, 0.9)]
    assert energies == sorted(energies) and energies[-1] > 2 * energies[0]
    assert rc.critical_impact_parameter(r1, r2) == 0.5
    assert rc.impact_angle_cdf(0.7071) == pytest.approx(0.5, abs=0.001) and rc.impact_angle_cdf(1.0) == 1.0
    with pytest.raises(rc.RogueCollisionError):
        rc.impact_angle_cdf(1.5)


def test_the_pieces_fly_at_a_few_escape_speeds_and_about_a_tenth_to_a_seventh_leaves():
    shares = [rc.unbound_fraction(m) for m in (0.5, 0.122, 0.044, 0.005)]
    assert shares[1] == pytest.approx(0.109, abs=0.003) and shares[-1] == pytest.approx(0.135, abs=0.003)
    assert all(0.0 < s < 0.14 for s in shares)
    assert rc.unbound_fraction(0.1, 0.0) == pytest.approx(0.9)                  # at a threshold of 0 every fragment counts
    rms = rc.unbound_rms_speed_ms(0.1, 2 * EARTH)
    assert 8 * KM < rms < 18 * KM                                               # inside the design's stand-in bracket of 5 to 15 km/s, a little over at its top
    with pytest.raises(rc.RogueCollisionError):
        rc.unbound_fraction(1.0)


@pytest.mark.parametrize("kms,gpa", [(5, 70), (10, 200), (20, 600), (40, 2100), (60, 4300)])
def test_the_peak_shock_pressure_of_two_rocks_follows_the_designs_table(kms, gpa):
    assert rc.peak_shock_pressure_pa(kms * KM) / 1e9 == pytest.approx(gpa, rel=0.1)


def test_arguments_must_be_positive():
    with pytest.raises(rc.RogueCollisionError):
        rc.cross_section_m2(EARTH, EARTH, EARTH_RADIUS_M, EARTH_RADIUS_M, 0.0)
    with pytest.raises(rc.RogueCollisionError):
        rc.disruption_speed_ms(0.0, EARTH)
    with pytest.raises(rc.RogueCollisionError):
        rc.speed_fractions(1.0, 1.0, [3.0, 2.0])
