"""The heavy-gas loss term (GEN.1558 prebuild): pure, bounded, monotone, and a window of one efficiency that separates the Earth from Mars."""

import math

import pytest

from planetgen.physics import constants, heavy_gas_loss as hgl

ME = constants.EARTH_MASS_TO_KG
SUN = constants.SOLAR_LUMINOSITY


def lifted(mass_earth, radius_km, xuv, age=4.5):
    return hgl.lifted_mass_kg(mass_earth * ME, radius_km * 1e3, xuv, age, 1.0, SUN)


def air_mass(mass_earth, radius_km, bar):
    gravity = constants.G * mass_earth * ME / (radius_km * 1e3) ** 2
    return 4.0 * math.pi * (radius_km * 1e3) ** 2 * bar * 1e5 / gravity


def test_no_xuv_or_no_age_lifts_nothing():
    assert lifted(1.0, 6371.0, 0.0) == 0.0
    assert hgl.lifted_mass_kg(ME, 6.371e6, 1.0, 0.0, 1.0, SUN) == 0.0


def test_the_bound_grows_with_the_flux_and_falls_with_the_mass():
    assert lifted(1.0, 6371.0, 2.0) > lifted(1.0, 6371.0, 1.0)
    assert lifted(0.01, 1700.0, 1.0) > lifted(1.0, 6371.0, 1.0)


def test_the_share_kept_is_bounded_and_monotone():
    assert hgl.kept_share(1e18, 1e21, 0.0) == 1.0
    assert hgl.kept_share(0.0, 1e21, 0.01) == 0.0
    shares = [hgl.kept_share(1e18, 1e21, e) for e in (1e-4, 1e-3, 1e-2, 1e-1)]
    assert all(0.0 <= s <= 1.0 for s in shares) and shares == sorted(shares, reverse=True)
    assert hgl.kept_share(1e18, 1e25, 1.0) == pytest.approx(math.exp(-700.0))


def test_a_field_keeps_more_than_no_field():
    assert hgl.kept_share(1e18, 1e21, 0.01, magnetopause_rp=10.0) > hgl.kept_share(1e18, 1e21, 0.01, magnetopause_rp=None)
    assert hgl.kept_share(1e18, 1e21, 0.01, magnetopause_rp=None) == hgl.kept_share(1e18, 1e21, 0.01, magnetopause_rp=0.0)


def test_negative_efficiency_is_refused():
    with pytest.raises(ValueError):
        hgl.kept_share(1e18, 1e21, -0.1)


def test_one_efficiency_separates_the_earth_from_mars_and_the_moon_sized_body():
    """The window of the module: at 0.007 the Earth (closed at ten radii) keeps most of its air, an open Mars and an open Moon-sized body lose most of theirs."""
    efficiency = 0.007
    earth = hgl.kept_share(air_mass(1.0, 6371.0, 1.0), lifted(1.0, 6371.0, 1.0), efficiency, 10.0)
    mars = hgl.kept_share(air_mass(0.107, 3389.5, 5.0), lifted(0.107, 3389.5, 0.43), efficiency, None)
    moon = hgl.kept_share(air_mass(0.0123, 1737.4, 1.0), lifted(0.0123, 1737.4, 1.0), efficiency, None)
    assert earth > 0.6 and mars < 0.4 and moon < 1e-6
