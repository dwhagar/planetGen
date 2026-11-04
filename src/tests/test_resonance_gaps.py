# tests/test_resonance_gaps.py

"""GEN.1062: the resonance zones, chaotic zones, Yarkovsky cut and stirring of a belt (`physics/resonance_gaps.py`; belt-evolution.md sections 2.3 and 2.4)."""

import math

import pytest

from planetgen.physics import constants, resonance_gaps as rg

SUN = constants.SOLAR_MASS_TO_KG
EARTH = constants.EARTH_MASS_TO_KG
JUPITER = rg.Giant(5.203, 317.8 * EARTH, 0.048)
SATURN = rg.Giant(9.537, 95.16 * EARTH, 0.054)
URANUS = rg.Giant(19.19, 14.54 * EARTH, 0.047)
NEPTUNE = rg.Giant(30.07, 17.15 * EARTH, 0.009)


def by_label(belt):
    return {(r.label, r.giant): r for r in belt.resonances}


def test_jupiter_makes_the_kirkwood_gaps_where_they_are():
    """A test and not an input: 3:1 at 2.50, 5:2 at 2.82, 7:3 at 2.96 and 2:1 at 3.27 AU, inside a belt of 2.2 to 3.4 AU."""
    belt = rg.belt_map([JUPITER], SUN, 2.2, 3.4, 4500.0)
    found = by_label(belt)
    for label, au in (("3:1", 2.50), ("5:2", 2.82), ("7:3", 2.96), ("2:1", 3.27)):
        assert found[(label, 0)].a_au == pytest.approx(au, abs=0.01)
        assert found[(label, 0)].side == rg.INTERIOR
    assert {r.label for r in belt.resonances} == {"3:1", "5:2", "7:3", "2:1"}
    assert [r.a_au for r in belt.resonances] == sorted(r.a_au for r in belt.resonances)


def test_the_orders_and_widths_follow_the_resonances():
    belt = rg.belt_map([JUPITER], SUN, 2.2, 3.4, 4500.0)
    found = by_label(belt)
    assert [found[(label, 0)].order for label in ("2:1", "3:1", "5:2", "7:3")] == [1, 2, 3, 4]
    assert found[("2:1", 0)].width_au > found[("3:1", 0)].width_au > found[("5:2", 0)].width_au > found[("7:3", 0)].width_au > 0.0
    assert 0.05 < found[("2:1", 0)].width_au < 0.6        # a libration width, wider than the cleared part of the 2:1 Kirkwood gap
    for r in belt.resonances:
        assert r.inner_au < r.a_au < r.outer_au


def test_the_first_order_width_is_the_two_laws_of_mustill_and_wyatt():
    mu = 1.0e-3
    low = rg.resonance_width(mu, 2, 1, 1e-4)
    assert low == pytest.approx(3.2 * mu ** (2 / 3) * 2 ** (1 / 3))
    high = rg.resonance_width(mu, 2, 1, 0.2)
    assert high == pytest.approx(7.78 * math.sqrt(0.2 * mu * 2))
    assert rg.resonance_width(mu, 3, 2, 0.1) > rg.resonance_width(mu, 3, 1, 0.1)
    with pytest.raises(ValueError):
        rg.resonance_width(mu, 1, 2, 0.1)


def test_a_gas_giants_zone_is_the_classical_law_at_low_eccentricity_and_the_extended_one_above():
    mu = 9.5e-4
    half, regime = rg.chaotic_zone(mu, 0.001)
    assert (half, regime) == (pytest.approx(1.3 * mu ** (2 / 7)), "classical")
    critical = 0.21 * mu ** (3 / 7)
    half, regime = rg.chaotic_zone(mu, 2 * critical)
    assert regime == "extended" and half == pytest.approx(1.8 * (2 * critical * mu) ** 0.2)
    assert half > 1.3 * mu ** (2 / 7)
    half, regime = rg.chaotic_zone(mu, 0.9)
    assert regime == "crossing" and half >= 0.9
    # the two laws meet at the critical eccentricity [R3226]
    assert 1.8 * (critical * mu) ** 0.2 == pytest.approx(1.3 * mu ** (2 / 7), rel=0.05)


def test_the_empty_zone_inside_neptunes_chaotic_zone():
    """A test and not an input: at zero eccentricity Neptune clears 27.7 to 32.4 AU; with the belt's small free eccentricity the zone is a little wider."""
    half, regime = rg.chaotic_zone(NEPTUNE.mass_kg / SUN, 0.0)
    assert regime == "classical"
    assert NEPTUNE.a_au * (1 - half) == pytest.approx(27.7, abs=0.1) and NEPTUNE.a_au * (1 + half) == pytest.approx(32.4, abs=0.1)
    belt = rg.belt_map([NEPTUNE], SUN, 25.0, 50.0, 4500.0)
    zone = belt.zones[0]
    assert zone.regime in ("classical", "extended")
    assert zone.inner_au < NEPTUNE.a_au < zone.outer_au
    assert zone.inner_au == pytest.approx(27.7, abs=1.0) and zone.outer_au == pytest.approx(32.4, abs=1.0)
    assert belt.survivors[-1][1] == 50.0
    assert all(low >= zone.outer_au - 1e-9 or high <= zone.inner_au + 1e-9 for low, high in belt.survivors)
    assert not any(zone.inner_au < r.a_au < zone.outer_au for r in belt.resonances)


def test_neptunes_exterior_resonances_are_the_kuiper_belts():
    belt = rg.belt_map([NEPTUNE], SUN, 35.0, 55.0, 4500.0)
    found = by_label(belt)
    assert found[("3:2", 0)].a_au == pytest.approx(39.4, abs=0.1) and found[("3:2", 0)].side == rg.EXTERIOR
    assert found[("2:1", 0)].a_au == pytest.approx(47.7, abs=0.1)
    assert found[("5:3", 0)].a_au == pytest.approx(42.3, abs=0.1)


def test_a_belt_far_from_the_giants_has_no_zone():
    belt = rg.belt_map([JUPITER], SUN, 20.0, 60.0, 100.0)
    assert belt.zones == () and belt.survivors == ((20.0, 60.0),)


def test_the_survivors_are_the_belt_outside_every_zone():
    belt = rg.belt_map([JUPITER], SUN, 3.0, 9.0, 1000.0)       # across Jupiter's own orbit
    assert len(belt.zones) == 1
    zone = belt.zones[0]
    assert belt.survivors == ((3.0, zone.inner_au), (zone.outer_au, 9.0))


def test_the_yarkovsky_cut_grows_with_the_age_and_the_solar_systems_is_a_few_km():
    young = rg.belt_map([JUPITER], SUN, 2.2, 3.4, 100.0)
    old = rg.belt_map([JUPITER], SUN, 2.2, 3.4, 4500.0)
    assert old.cut_diameter_km == pytest.approx(45.0 * young.cut_diameter_km, rel=1e-9)
    assert 1.0 < old.cut_diameter_km < 40.0           # bodies above 40 km do not move in 4.5 Gyr
    assert old.gaps_open and not rg.belt_map([JUPITER], SUN, 2.2, 3.4, 1.0).gaps_open


def test_the_drift_is_the_measured_bennu_scaled():
    assert rg.drift_au_per_myr(0.493, 1.1) == pytest.approx(18.9e-4, rel=0.02)
    assert rg.drift_au_per_myr(1.0, 2.5) == pytest.approx(1.8e-4, rel=0.1)
    assert rg.drift_au_per_myr(40.0, 2.5) * 4500.0 < 0.05
    assert rg.cut_diameter_km(2.5, 4500.0, rg.drift_au_per_myr(5.0, 2.5) * 4500.0) == pytest.approx(5.0)


def test_an_eccentric_giant_stirs_the_belt():
    quiet = rg.belt_map([rg.Giant(5.203, 317.8 * EARTH, 0.0)], SUN, 2.9, 3.0, 1000.0)
    stirred = rg.belt_map([JUPITER], SUN, 2.9, 3.0, 1000.0)
    assert quiet.forced_eccentricity == 0.0
    assert stirred.forced_eccentricity == pytest.approx(1.25 * (2.95 / 5.203) * 0.048 / (1 - 0.048 ** 2), rel=1e-9)
    assert stirred.eccentricity == pytest.approx(max(stirred.forced_eccentricity, (3.0 - 2.9) / 5.9))
    assert stirred.eccentricity > quiet.eccentricity or stirred.forced_eccentricity <= quiet.eccentricity
    assert rg.forced_eccentricity([], 3.0) == 0.0


def test_bad_input_is_refused():
    with pytest.raises(ValueError):
        rg.belt_map([JUPITER], SUN, 3.0, 2.0, 100.0)
    with pytest.raises(ValueError):
        rg.Giant(-1.0, 1.0)
    with pytest.raises(ValueError):
        rg.Giant(1.0, 1.0, 1.0)
