# tests/test_super_puff.py

"""Super-puff radius and its fade with age (GEN.1042): Kepler-51 b and d, WASP-107 b, Millholland's tidal inflation."""

import math

from planetgen.physics import bulk_structure, constants, envelope_escape, super_puff

SUN_L = 0.72        # Kepler-51, from 0.87 solar radii at about 5700 K


def _envelope_for_radius(mass_earth, target_radius_earth, flux_earth, age_gy):
    """The envelope fraction the law gives the radius, by bisection; None when the radius lies outside what 0.5 of the mass can give."""
    low, high = 0.002, envelope_escape.MAX_ENVELOPE
    if not (super_puff.radius_earth(mass_earth, low, flux_earth, age_gy) <= target_radius_earth <= super_puff.radius_earth(mass_earth, high, flux_earth, age_gy)):
        return None
    for _ in range(60):
        middle = 0.5 * (low + high)
        if super_puff.radius_earth(mass_earth, middle, flux_earth, age_gy) < target_radius_earth:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def _radius_at_age(mass_earth, envelope, flux_earth, host_mass_solar, luminosity_solar, age_gy):
    """The radius after the XUV has taken its share of the envelope (the loss `envelope_capture.kept` takes), under the age term of the law."""
    left = envelope_escape.envelope_left_earth(mass_earth * (1.0 - envelope), envelope, flux_earth, host_mass_solar, luminosity_solar, age_gy)
    if left <= 0.0:
        return bulk_structure.core_radius_earth(mass_earth * (1.0 - envelope))
    total = mass_earth * (1.0 - envelope) + left
    return super_puff.radius_earth(total, left / total, flux_earth, age_gy)


def test_kepler_51_d_is_a_thick_envelope_at_the_density_measured():
    # 5.6 Earth masses, 9.32 Earth radii, 3.5 Earth fluxes (350 K), 0.5 Gyr [R2673]: "an envelope of over 30 percent"
    envelope = _envelope_for_radius(5.6, 9.32, 3.5, 0.5)
    assert 0.3 < envelope < 0.5
    assert abs(super_puff.density_g_cm3(5.6, 9.32) - 0.038) < 0.009
    assert abs(super_puff.radius_earth(5.6, envelope, 3.5, 0.5) - 9.32) < 0.01


def test_kepler_51_b_is_inside_the_envelope_range_of_lee_and_chiang():
    # 2.1 Earth masses and 7.1 Earth radii at 0.2 AU; "gas-to-core ratio 7 to 28 percent"
    envelope = _envelope_for_radius(2.1, 7.1, SUN_L / 0.2 ** 2, 0.5)
    assert 0.07 < envelope < 0.28


def test_wasp_107_b_leaves_the_envelope_law_for_the_giant_tables():
    # 30.5 Earth masses, 0.96 Jupiter radii (10.7 Earth radii), 0.134 g/cc, core under 4.6 Earth masses [R3516]
    assert not super_puff.covers(30.5, 0.88)
    radius = super_puff.radius_earth(30.5, 1.0 - 3.7 / 30.5, 70.0, 3.4)
    assert abs(radius - 0.96 * constants.JUPITER_RADIUS_KM / constants.EARTH_RADIUS_KM) < 0.05 * radius
    assert abs(super_puff.density_g_cm3(30.5, radius) - 0.134) < 0.03
    assert not super_puff.is_super_puff(30.5, radius, 3.4, 5.7)


def test_inside_the_law_the_radius_is_lopez_and_fortneys():
    assert super_puff.radius_earth(5.0, 0.2, 10.0, 2.0) == bulk_structure.envelope_radius_earth(5.0, 0.2, 10.0, 2.0)
    assert _envelope_for_radius(5.0, 50.0, 10.0, 2.0) is None
    assert _envelope_for_radius(40.0, 8.0, 10.0, 2.0) is None


def test_the_gate_is_the_design_row():
    assert super_puff.is_super_puff(5.6, 9.32, 0.5, 130.0) is False       # 130 d is past the period range
    assert super_puff.is_super_puff(5.6, 9.32, 0.5, 85.0)
    assert not super_puff.is_super_puff(5.6, 9.32, 2.0, 85.0)             # too old
    assert not super_puff.is_super_puff(12.0, 9.32, 0.5, 45.0)            # too heavy
    assert not super_puff.is_super_puff(5.6, 3.0, 0.5, 45.0)              # dense, a normal sub-Neptune


def test_tides_inflate_by_millholland_eq_14_only_from_20_earth_fluxes():
    assert super_puff.tidal_inflation_ratio(5.6, 0.3, 3.5, 1.0e4) == 1.0
    ratio = super_puff.tidal_inflation_ratio(5.6, 0.3, 30.0, 1.0e4)
    assert 1.2 < ratio < 1.8
    assert super_puff.tidal_inflation_ratio(2.0, 0.4, 500.0, 1.0e2) == super_puff.TIDAL_MAX_RATIO
    # a stronger tide (smaller Q) inflates more, a heavier planet less
    assert super_puff.tidal_inflation_ratio(5.6, 0.3, 30.0, 1.0e3) > ratio
    assert super_puff.tidal_inflation_ratio(20.0, 0.3, 30.0, 1.0e4) < ratio


def test_tidal_power_goes_as_the_distance_to_minus_7_5():
    star = constants.SOLAR_MASS_TO_KG
    radius = 3.0 * constants.EARTH_RADIUS_KM * 1e3
    near = super_puff.tidal_luminosity_w(star, radius, 0.1 * constants.AU_TO_M, 1e4)
    far = super_puff.tidal_luminosity_w(star, radius, 0.2 * constants.AU_TO_M, 1e4)
    assert math.isclose(near / far, 2 ** 7.5, rel_tol=1e-9)
    ratio = super_puff.tidal_luminosity_ratio(star, radius, 0.1 * constants.AU_TO_M, 100.0, 1e4)
    assert 0.0 < ratio < 1.0


def test_a_cool_puffy_planet_contracts_with_age_and_keeps_its_envelope():
    young = _radius_at_age(5.6, 0.41, 3.5, 0.96, SUN_L, 0.5)
    old = _radius_at_age(5.6, 0.41, 3.5, 0.96, SUN_L, 5.0)
    assert old < young
    # the age term of the law alone: the envelope thickness falls as (age / 5)^-0.18, so a tenth of the age is 0.66 of it
    assert 0.6 < old / young < 0.95
    assert young > 8.5


def test_a_hot_puffy_planet_loses_its_envelope_and_shrinks_to_a_core():
    bare = _radius_at_age(5.0, 0.3, 80.0, 1.0, 1.0, 8.0)
    young = _radius_at_age(5.0, 0.3, 80.0, 1.0, 1.0, 0.5)
    assert bare < 0.5 * young
    assert bare >= bulk_structure.core_radius_earth(5.0 * 0.7) - 1e-9


def test_the_rows_are_in_the_design_ranges():
    assert super_puff.PUFF_ENVELOPE[0] < super_puff.PUFF_ENVELOPE[1] <= 0.5
    assert super_puff.OCCURRENCE_SHARE < 0.01
    assert super_puff.PERIOD_RANGE_D == (10.0, 100.0)
