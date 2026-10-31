# tests/test_moon_orbit_evolution.py

"""
GEN.958: the orbits of a planet's moons under its tide (`planetgen.physics.moon_orbits`): the Moon's recession, Phobos's and Triton's fall, the Roche
fall into a ring, and the conservation of angular momentum between the planet's spin and the moons' orbits.
"""

import math

import pytest

from planetgen.physics import constants, moon_orbits, spin_torques

YEAR = constants.SECONDS_PER_YEAR
EARTH = spin_torques.Planet(5.972e24, 6.371e6)
MOON = moon_orbits.Moon(7.342e22, 3.844e8, 3344.0)
EARTH_SPIN = 2.0 * math.pi / 86164.0


def _earth_torques(k2_over_q):
    value = k2_over_q
    return spin_torques.Torques(EARTH, (), spin_torques.constant_dissipation(value * 100.0, 100.0), None, 1.5, moon_orbits.COROTATION_SOFTNESS)


def test_the_moon_recedes_at_the_rate_lunar_laser_ranging_measures():
    """3.82 cm a year (Dickey et al. 1994, recalled): the Earth's effective k2/Q of 0.025 and the standard 3/2 torque give it to 10 percent."""
    rate = moon_orbits.drift_rate_m_s(_earth_torques(0.025), EARTH_SPIN, EARTH.mass_kg, MOON)
    assert rate * YEAR == pytest.approx(0.0382, rel=0.1)
    assert rate > 0.0


def test_the_integration_recedes_the_moon_by_the_rate_over_a_million_years():
    result = moon_orbits.evolve(EARTH, None, spin_torques.constant_dissipation(2.5, 100.0), [MOON], EARTH_SPIN, 1.0e6 * YEAR, factor=1.5)
    (fate,) = result.moons
    assert fate.fate == "stays"
    assert (fate.distance_m - MOON.distance_m) / 1.0e6 == pytest.approx(0.0382, rel=0.15)     # m a year
    assert result.spin_rate < EARTH_SPIN


def test_angular_momentum_is_conserved_between_the_spin_and_the_orbit():
    """With no star: `I omega + L_orbit` before and after 1 Gy for a Moon that starts close and a fast Earth."""
    close = moon_orbits.Moon(7.342e22, 1.0e8, 3344.0)
    spin0 = 2.0 * math.pi / (6.0 * 3600.0)
    result = moon_orbits.evolve(EARTH, None, spin_torques.constant_dissipation(2.5, 100.0), [close], spin0, 1.0e9 * YEAR, factor=1.5)
    (fate,) = result.moons
    assert fate.fate == "stays" and fate.distance_m > close.distance_m
    before = EARTH.inertia * spin0 + moon_orbits.orbital_angular_momentum(EARTH.mass_kg, close)
    after = EARTH.inertia * result.spin_rate + moon_orbits.orbital_angular_momentum(
        EARTH.mass_kg, moon_orbits.Moon(close.mass_kg, fate.distance_m, close.density_kg_m3))
    assert after == pytest.approx(before, rel=1.0e-5)
    assert result.planet_heat_j > 0.0


def test_the_heat_in_the_planet_is_the_work_the_tide_does():
    """The dissipated energy is the loss of spin and orbital energy: `I omega^2 / 2 - G M m / 2a` summed, before minus after."""
    close = moon_orbits.Moon(7.342e22, 1.0e8, 3344.0)
    spin0 = 2.0 * math.pi / (6.0 * 3600.0)
    result = moon_orbits.evolve(EARTH, None, spin_torques.constant_dissipation(2.5, 100.0), [close], spin0, 1.0e9 * YEAR, factor=1.5)
    (fate,) = result.moons

    def energy(spin, a):
        return 0.5 * EARTH.inertia * spin ** 2 - constants.G * EARTH.mass_kg * close.mass_kg / (2.0 * a)

    lost = energy(spin0, close.distance_m) - energy(result.spin_rate, fate.distance_m)
    assert result.planet_heat_j == pytest.approx(lost, rel=0.02)


def test_phobos_spirals_in_at_a_rate_of_the_order_measured():
    """Phobos: Mars turns slower than Phobos goes round, so it falls, about 1.8 cm a year (Bills et al. 2005, recalled; Q of 88 and k2 of 0.14)."""
    mars = spin_torques.Planet(6.417e23, 3.3895e6)
    phobos = moon_orbits.Moon(1.066e16, 9.376e6, 1876.0)
    torques = spin_torques.Torques(mars, (), spin_torques.constant_dissipation(0.14, 88.0), None, 1.5, moon_orbits.COROTATION_SOFTNESS)
    rate = moon_orbits.drift_rate_m_s(torques, 2.0 * math.pi / (24.62 * 3600.0), mars.mass_kg, phobos)
    assert rate < 0.0
    assert -rate * YEAR == pytest.approx(0.018, rel=1.5)


def test_triton_goes_in_because_it_is_retrograde_and_never_recedes():
    """Triton at Neptune: retrograde, so the tide takes its orbital angular momentum; the decay time to the Roche limit is some Gy at Q of 1e4 (Chyba et al. 1989)."""
    neptune = spin_torques.Planet(1.024e26, 2.4622e7)
    triton = moon_orbits.Moon(2.139e22, 3.548e8, 2061.0, retrograde=True)
    torques = spin_torques.Torques(neptune, (), spin_torques.constant_dissipation(0.41, 1.0e4), None, 1.5, moon_orbits.COROTATION_SOFTNESS)
    rate = moon_orbits.drift_rate_m_s(torques, 2.0 * math.pi / (16.1 * 3600.0), neptune.mass_kg, triton)
    assert rate < 0.0
    roche = moon_orbits.roche_limit_m(neptune.radius_m, neptune.mean_density, triton.density_kg_m3)
    fall_gy = (2.0 / 13.0) * (triton.distance_m - roche) / -rate / YEAR / 1.0e9
    assert 1.0 < fall_gy < 30.0


def test_a_moon_of_a_locked_planet_falls_into_a_ring_and_the_spin_takes_its_momentum():
    """A planet that keeps one face to its star turns far slower than a close moon goes round: the moon falls, its mass is ring material, and the
    ring's angular momentum spins the planet up."""
    slow = 2.0 * math.pi / (30.0 * 86400.0)
    near = moon_orbits.Moon(7.0e21, 4.0e7, 3000.0)
    result = moon_orbits.evolve(EARTH, None, spin_torques.constant_dissipation(2.5, 100.0), [near], slow, 4.5e9 * YEAR, factor=1.5)
    (fate,) = result.moons
    assert fate.fate == "fell"
    assert fate.time_s < 4.5e9 * YEAR
    assert fate.distance_m == pytest.approx(moon_orbits.roche_limit_m(EARTH.radius_m, EARTH.mean_density, near.density_kg_m3))
    assert result.ring_mass_kg == near.mass_kg
    assert result.spin_rate > slow


def test_a_moon_inside_the_roche_limit_from_the_start_falls_at_once():
    inside = moon_orbits.Moon(1.0e20, 1.0e7, 3000.0)
    result = moon_orbits.evolve(EARTH, None, spin_torques.constant_dissipation(2.5, 100.0), [inside], EARTH_SPIN, 1.0e9 * YEAR, factor=1.5)
    assert result.moons[0].fate == "fell" and result.moons[0].time_s == 0.0


def test_a_moon_pushed_past_the_hill_sphere_escapes():
    hill = 1.0e8
    far = moon_orbits.Moon(7.342e22, 3.9e7, 3344.0)
    spin0 = 2.0 * math.pi / (3.0 * 3600.0)
    result = moon_orbits.evolve(EARTH, None, spin_torques.constant_dissipation(25.0, 100.0), [far], spin0, 1.0e10 * YEAR, hill_radius_m=hill, factor=1.5)
    assert result.moons[0].fate == "escaped"
    assert result.moons[0].distance_m >= moon_orbits.HILL_ESCAPE_FRACTION * hill * 0.999


def test_the_star_locking_the_planet_pulls_a_moon_inside_the_synchronous_radius_in():
    """A close-in planet is locked by its star in tens of millions of years and then turns once an orbit (4 days); a moon inside the radius where it
    would go round that fast (about 1e8 m) is outside the planet's turn, and falls. The same moon round a fast Earth drifts out."""
    a = 0.05 * constants.AU_TO_M
    sun = 1.989e30
    star = spin_torques.Raiser(sun, a, spin_torques.orbital_rate(sun, a))
    moon = moon_orbits.Moon(7.342e22, 5.0e7, 3344.0)
    locked = moon_orbits.evolve(EARTH, star, spin_torques.constant_dissipation(0.5, 100.0), [moon], EARTH_SPIN, 4.5e9 * YEAR)
    assert locked.moons[0].fate == "fell"
    free = moon_orbits.evolve(EARTH, None, spin_torques.constant_dissipation(2.5, 100.0), [moon], EARTH_SPIN, 1.0e8 * YEAR, factor=1.5)
    assert free.moons[0].distance_m > moon.distance_m


def test_several_moons_are_followed_together():
    moons = [moon_orbits.Moon(1.0e20, 2.0e7, 3000.0), MOON, moon_orbits.Moon(1.0e19, 1.2e9, 2000.0)]
    result = moon_orbits.evolve(EARTH, None, spin_torques.constant_dissipation(2.5, 100.0), moons, 2.0 * math.pi / (6.0 * 3600.0), 1.0e9 * YEAR, factor=1.5)
    assert [f.moon for f in result.moons] == moons
    assert result.moons[0].fate == "stays" or result.moons[0].fate == "fell"
    assert result.moons[2].distance_m >= moons[2].distance_m


def test_a_planet_held_near_a_balance_of_the_tides_integrates_in_seconds():
    """A planet of 0.2 Earth masses with a thermal tide and three moons, one inside its Roche limit, whose spin hovers at a balance of the thermal and rheology
    tides: at a relative tolerance of 1e-7 LSODA took minutes; at the module's tolerance it takes about a second, and the answer (the spin it ends at) agrees."""
    import time
    from planetgen.physics import tidal_rheology
    planet = spin_torques.Planet(1.1944335735582758e+24, 3931635.019609854, 0.33)
    star = spin_torques.Raiser(1.9471535274991518e+30, 68448422524.914894, 6.365870419227959e-07)
    thermal = spin_torques.ThermalTide(star, 17.65749914006499, 7.491435544865652e-06)

    def dissipation(frequency):
        return tidal_rheology.response(planet.mass_kg, planet.radius_m, max(frequency, 1.0e-20), 1746.3423804466247, 640405843720.1215, 265994.4516871315, 0.0).k2_over_q

    moons = [moon_orbits.Moon(6.082467392073513e+21, 8958795.337131515, 2355.2), moon_orbits.Moon(9.227560476330121e+20, 128270066.15492842, 2310.6),
             moon_orbits.Moon(2.257114046278461e+22, 169499247.46339804, 2061.5)]
    start = time.time()
    result = moon_orbits.evolve(planet, star, dissipation, moons, 1.1462362329176932e-05, 3.3536155063031206e+17, thermal, 4691.9546812192075, 403252317.70733774)
    assert time.time() - start < 30.0
    assert result.spin_rate == pytest.approx(8.2504e-6, rel=1.0e-3)
    assert all(f.fate == "fell" for f in result.moons) and result.ring_mass_kg == pytest.approx(sum(m.mass_kg for m in moons))
