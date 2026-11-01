"""
The dynamo field module (GEN.966, `planetgen.physics.dynamo`): the dipole moment and its share from the power of a core and the spin, tested on
Earth (the unit), Jupiter and Saturn (the giants), Ganymede (a moon), Mercury (the thin shell the scaling overshoots) and a dead Mars.
"""

import math

import pytest

from planetgen.physics import dynamo, magnetism


def _spin(hours):
    return 2.0 * math.pi / (hours * 3600.0)


def _earth(**kw):
    args = dict(power_w=dynamo.EARTH_DYNAMO_POWER_W, outer_m=dynamo.EARTH_OUTER_CORE_M, inner_m=dynamo.EARTH_INNER_CORE_M,
                density_kg_m3=dynamo.EARTH_OUTER_CORE_DENSITY_KG_M3, spin_rad_s=dynamo.EARTH_SPIN_RAD_S)
    args.update(kw)
    return dynamo.field(**args)


def test_earth_gives_its_own_moment_rossby_number_and_surface_field():
    earth = _earth()
    assert earth.moment_a_m2 == pytest.approx(7.8e22, rel=0.01)
    assert dynamo.EARTH_DYNAMO_POWER_W == pytest.approx(1.15e12, rel=0.01)                  # inside the published 0.2 to 3 TW
    assert earth.rossby == pytest.approx(0.09, rel=1e-9)
    assert earth.dipole_fraction == pytest.approx(0.80, abs=0.01)
    assert earth.surface_equatorial_t(6.371e6) == pytest.approx(3.0e-5, rel=0.05)      # 30 microtesla at the equator
    assert earth.dipolar_field_t == pytest.approx(1.86e-4, rel=0.01)                      # the dipole at the core-mantle boundary: 0.2 mT
    assert earth.top_equatorial_t == pytest.approx(earth.dipolar_field_t)


def test_a_dead_mars_has_no_field_however_it_turns():
    assert _earth(power_w=0.0).moment_a_m2 == 0.0
    assert dynamo.field(0.0, 1.7e6, 0.0, 6000.0, _spin(24.6)).moment_a_m2 == 0.0
    assert dynamo.field(-5.0e9, 1.7e6, 0.0, 6000.0, _spin(24.6)).moment_a_m2 == 0.0
    assert dynamo.field(1.0e11, 0.0, 0.0, 6000.0, _spin(24.6)).moment_a_m2 == 0.0         # no region
    assert dynamo.field(1.0e11, 1.0e6, 1.0e6, 6000.0, _spin(24.6)).moment_a_m2 == 0.0     # no thickness


def test_jupiter_is_driven_by_a_tenth_of_its_heat_and_gives_its_moment_and_four_gauss():
    jupiter = dynamo.field(dynamo.giant_power_w(3.35e17), 0.84 * 7.1492e7, 0.0, 2000.0, _spin(9.925))
    assert jupiter.moment_a_m2 == pytest.approx(1.5e27, rel=0.15)                         # 1.5e27 A m2
    assert jupiter.dipole_fraction > 0.99 and jupiter.rossby < 0.05                      # a fast, strongly driven dynamo is dipolar
    assert jupiter.surface_equatorial_t(7.1492e7) == pytest.approx(4.2e-4, rel=0.15)      # 4.2 gauss


def test_saturn_comes_out_within_a_factor_of_seven_of_its_weak_field():
    """The same share of the heat gives Saturn six times its moment: a stable helium-rain layer weakens Saturn's dynamo, which the scaling does not model."""
    saturn = dynamo.field(dynamo.giant_power_w(8.6e16), 0.6 * 6.0268e7, 0.0, 1500.0, _spin(10.656))
    assert 4.6e25 < saturn.moment_a_m2 < 7.0 * 4.6e25


def test_ganymedes_slow_turn_makes_a_weak_multipolar_field_of_its_moment_to_a_factor_of_two():
    ganymede = dynamo.field(1.0e10, 8.0e5, 0.0, 5500.0, _spin(7.155 * 24.0))
    assert ganymede.dipole_fraction < 0.2 and ganymede.rossby > dynamo.DIPOLE_ROSSBY      # tidally locked: slow for its power
    assert 1.3e20 / 2.0 < ganymede.moment_a_m2 < 1.3e20 * 2.0


def test_mercurys_thin_shell_overshoots_by_the_stable_layer_the_scaling_does_not_hold():
    shell = dict(outer_m=2.0e6, inner_m=1.4e6, density_kg_m3=7000.0, spin_rad_s=_spin(58.65 * 24.0))
    scaled = dynamo.field(1.0e10, **shell)
    assert 10.0 * 3.0e19 < scaled.moment_a_m2 < 30.0 * 3.0e19                             # 3e19 observed: 17 times less
    layered = dynamo.field(1.0e10, attenuation=0.06, **shell)
    assert layered.moment_a_m2 == pytest.approx(0.06 * scaled.moment_a_m2)
    assert 1.0e19 < layered.moment_a_m2 < 1.0e20
    assert dynamo.field(1.0e10, attenuation=2.0, **shell).moment_a_m2 == pytest.approx(scaled.moment_a_m2)       # a layer never adds field
    assert dynamo.field(1.0e10, attenuation=-1.0, **shell).moment_a_m2 == 0.0


def test_the_moment_follows_the_cube_root_of_the_power_while_the_dynamo_stays_dipolar():
    full, eighth = _earth(), _earth(power_w=dynamo.EARTH_DYNAMO_POWER_W / 8.0)
    assert eighth.dipole_fraction > 0.95
    assert eighth.moment_a_m2 / full.moment_a_m2 == pytest.approx(0.5, rel=1e-9)


def test_a_slower_spin_makes_the_field_multipolar_and_lighter():
    day, week = _earth(), _earth(spin_rad_s=dynamo.EARTH_SPIN_RAD_S / 7.0)
    assert week.rossby == pytest.approx(day.rossby * 7.0 ** (3.0 * dynamo.FLUX_RAYLEIGH_EXPONENT))        # the cube of the spin, to the 0.41
    assert week.dipole_fraction < 0.2 < day.dipole_fraction
    assert week.moment_a_m2 < 0.3 * day.moment_a_m2
    assert math.isinf(_earth(spin_rad_s=0.0).rossby) and _earth(spin_rad_s=0.0).dipole_fraction == pytest.approx(dynamo.MIN_DIPOLE_FRACTION)


def test_a_larger_core_makes_a_stronger_field():
    assert _earth(outer_m=1.5 * dynamo.EARTH_OUTER_CORE_M, inner_m=1.5 * dynamo.EARTH_INNER_CORE_M).moment_a_m2 > 2.0 * _earth().moment_a_m2


def test_the_share_of_the_dipole_follows_the_rotation_at_earths_power():
    """GEN.86's rule (the Rossby number from the rotation period alone at Earth's power) is this rule: a locked planet keeps about 0.15."""
    for hours in (6.0, 24.0, 48.0, 100.0, 1000.0):
        rossby = dynamo.local_rossby(dynamo.EARTH_DYNAMO_POWER_W, dynamo.EARTH_OUTER_CORE_M, dynamo.EARTH_INNER_CORE_M,
                                     dynamo.EARTH_OUTER_CORE_DENSITY_KG_M3, _spin(hours))
        assert rossby == pytest.approx(0.09 * (hours / 23.934) ** 1.23, rel=0.002)
        assert dynamo.dipole_fraction(rossby) == pytest.approx(dynamo.dipole_fraction(0.09 * (hours / 24.0) ** 1.23), abs=0.01)
    assert dynamo.dipole_fraction(0.0) == 1.0
    assert dynamo.dipole_fraction(10.0) == pytest.approx(0.15, abs=0.01)


def test_the_lifetime_is_the_run_of_ages_with_power():
    ages = [0.1, 0.5, 1.0, 2.0, 4.5]
    mars = dynamo.lifetime(ages, [1.0e11, 5.0e10, 0.0, 0.0, 0.0])
    assert (mars.started_gy, mars.ended_gy, mars.running) == (0.1, 1.0, False)
    earth = dynamo.lifetime(ages, [1.0e12, 2.0e12, 2.4e12, 2.4e12, 2.4e12])
    assert (earth.started_gy, earth.ended_gy, earth.running) == (0.1, None, True)
    late = dynamo.lifetime(ages, [0.0, 0.0, 1.0e10, 2.0e10, 2.0e10])
    assert (late.started_gy, late.running) == (1.0, True)
    never = dynamo.lifetime(ages, [0.0] * 5)
    assert (never.started_gy, never.ended_gy, never.running) == (None, None, False)


def test_a_giant_is_driven_by_a_share_of_its_heat_and_no_heat_gives_none():
    assert dynamo.giant_power_w(1.0e17) == pytest.approx(1.0e16)
    assert dynamo.giant_power_w(0.0) == 0.0 and dynamo.giant_power_w(-1.0) == 0.0
