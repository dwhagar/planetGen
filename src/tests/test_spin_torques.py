# tests/test_spin_torques.py

"""
GEN.956: the torques on a planet's spin (`planetgen.physics.spin_torques`): the star's tide, a moon's, the thermal tide of the air (Leconte et al.
2015), checked against the built despinning, the Earth-Moon day, Venus's balance and the closed form of the asynchronous equilibria.
"""

import math
import random
import types

import pytest

from planetgen.physics import constants, interior_heat, spin_torques, tidal_heating, tidal_rheology, tidal_spin

from tests.tide_support import rocky

EARTH_MASS, EARTH_RADIUS = 5.972e24, 6.371e6
MOON_MASS, MOON_A = 7.342e22, 3.844e8
SUN = 1.989e30
AU = constants.AU_TO_M
YEAR = constants.SECONDS_PER_YEAR


def _body(mass, radius):
    return rocky(mass, radius / 1000.0)


def _star_only(body, host_mass, a, factor=spin_torques.GRAVITATIONAL_TORQUE_FACTOR):
    k2, q = interior_heat.tidal_response(body, spin_torques.orbital_rate(host_mass, a))
    planet = spin_torques.Planet(body.mass, body.radius * 1000.0)
    star = spin_torques.Raiser(host_mass, a, spin_torques.orbital_rate(host_mass, a))
    return spin_torques.Torques(planet, (star,), spin_torques.constant_dissipation(k2, q), None, factor)


def _constant_q_star_despin(body, host, a, spin, age, e):
    """The closed form of a constant-Q despinning by one star: the spin falls (or rises) linearly at `F K k2 / (I Q)` to the first stop, the highest resonance of
    the window that holds `e` below the spin, else synchronous (the despinning of GEN.731, kept here as the check of the module)."""
    k2, q = interior_heat.tidal_response(body, spin_torques.orbital_rate(host, a))
    n = spin_torques.orbital_rate(host, a)
    radius = body.radius * 1000.0
    inertia = tidal_spin.moment_factor(body) * body.mass * radius ** 2
    rate = spin_torques.GRAVITATIONAL_TORQUE_FACTOR * constants.G * host ** 2 * radius ** 5 / a ** 6 * k2 / q / inertia
    ratio = 1.0
    for target, low, high in tidal_spin.RESONANCES:
        if target < spin / n and low < e < high:
            ratio = target
            break
    target_rate = ratio * n
    step = rate * age
    if abs(spin - target_rate) <= step:
        return target_rate, ratio
    return spin - math.copysign(step, spin - target_rate), None


def test_one_star_with_constant_q_is_the_closed_form_despinning():
    """The module is the linear constant-Q despinning when only the star's tide acts (the torque factor is the built one)."""
    rng = random.Random(7)
    for _ in range(60):
        body = _body(rng.uniform(0.05, 5.0) * EARTH_MASS, rng.uniform(0.4, 1.8) * EARTH_RADIUS)
        host = rng.uniform(0.15, 1.5) * SUN
        a = 10.0 ** rng.uniform(-2.0, 0.3) * AU
        n = spin_torques.orbital_rate(host, a)
        spin = n * 10.0 ** rng.uniform(-0.2, 3.0) * rng.choice((1.0, 1.0, -1.0))
        if spin < 0.0:
            continue
        age = rng.uniform(0.1, 10.0) * 1.0e9 * YEAR
        e = rng.choice((0.0, 0.0, 0.1, 0.206, 0.4))
        built_rate, built_ratio = _constant_q_star_despin(body, host, a, spin, age, e)
        found = spin_torques.evolve(_star_only(body, host, a), spin, age, e)
        assert found.rate == pytest.approx(built_rate, rel=1.0e-6, abs=1.0e-12 * n)
        assert found.ratio == built_ratio


def test_mercury_is_caught_in_the_three_halves_resonance():
    """Mercury: the planet's e of 0.206 is inside the 3:2 window, and a spin coming down from far above is caught there, as `tidal_spin` has it."""
    body = _body(3.30e23, 2.44e6)
    a = 0.387 * AU
    n = spin_torques.orbital_rate(SUN, a)
    result = spin_torques.evolve(_star_only(body, SUN, a), 40.0 * n, 20.0e9 * YEAR, 0.206)
    assert (result.ratio, result.kind) == (1.5, "resonance")
    assert result.rate == pytest.approx(1.5 * n)
    circular = spin_torques.evolve(_star_only(body, SUN, a), 40.0 * n, 20.0e9 * YEAR, 0.0)
    assert (circular.ratio, circular.kind) == (1.0, "synchronous")


def test_a_spin_that_misses_the_resonance_goes_on_to_the_synchronous_state():
    """GEN.1172: capture in the 3:2 window is a chance; a body whose uniform is above it slips by and is locked (e of 0.206 is outside the 1:1 window of the closed form, but the flip is a stop of its own)."""
    body = _body(3.30e23, 2.44e6)
    a = 0.387 * AU
    n = spin_torques.orbital_rate(SUN, a)
    chance = tidal_spin.RESONANCE_CAPTURE_PROBABILITY[1.5]
    caught = spin_torques.evolve(_star_only(body, SUN, a), 40.0 * n, 20.0e9 * YEAR, 0.206, (0.0, chance * 0.99))
    missed = spin_torques.evolve(_star_only(body, SUN, a), 40.0 * n, 20.0e9 * YEAR, 0.206, (0.0, chance * 1.01))
    assert (caught.ratio, caught.kind) == (1.5, "resonance")
    assert (missed.ratio, missed.kind) == (1.0, "synchronous")


def test_the_two_to_one_window_is_passed_without_capture_to_the_three_halves():
    """An orbit with e of 0.38 holds both the 2:1 and the 3:2 windows: the first uniform decides the 2:1, a pass lets the second decide the 3:2."""
    body = _body(3.30e23, 2.44e6)
    a = 0.387 * AU
    n = spin_torques.orbital_rate(SUN, a)
    pass_two, take_three = 0.5, 0.0
    first = spin_torques.evolve(_star_only(body, SUN, a), 40.0 * n, 20.0e9 * YEAR, 0.38, (0.0, 0.99))
    assert (first.ratio, first.kind) == (2.0, "resonance")
    second = spin_torques.evolve(_star_only(body, SUN, a), 40.0 * n, 20.0e9 * YEAR, 0.38, (pass_two, take_three))
    assert (second.ratio, second.kind) == (1.5, "resonance")
    neither = spin_torques.evolve(_star_only(body, SUN, a), 40.0 * n, 20.0e9 * YEAR, 0.38, (pass_two, pass_two))
    assert (neither.ratio, neither.kind) == (1.0, "synchronous")


def test_the_share_caught_is_the_chance_over_uniform_units():
    body = _body(3.30e23, 2.44e6)
    a = 0.387 * AU
    n = spin_torques.orbital_rate(SUN, a)
    units = [(i + 0.5) / 200.0 for i in range(200)]
    caught = sum(spin_torques.evolve(_star_only(body, SUN, a), 40.0 * n, 20.0e9 * YEAR, 0.206, (0.9, u)).kind == "resonance" for u in units)
    assert caught / 200.0 == pytest.approx(tidal_spin.RESONANCE_CAPTURE_PROBABILITY[1.5], abs=0.02)


def test_a_far_planet_keeps_its_spin():
    body = _body(EARTH_MASS, EARTH_RADIUS)
    a = AU
    n = spin_torques.orbital_rate(SUN, a)
    spin = 2.0 * math.pi / (5.0 * 3600.0)
    result = spin_torques.evolve(_star_only(body, SUN, a), spin, 4.5e9 * YEAR, 0.0167)
    assert result.kind == "free" and result.ratio is None
    assert 0.3 * spin < result.rate < spin     # the Andrade Q of tens takes between a tenth and two thirds of a five-hour spin off in 4.5 Gy


def test_the_moon_slows_the_earth_by_the_measured_two_milliseconds_a_century():
    """Earth-Moon: the day lengthens 2.3 ms per century (the Moon's share about 2.0); the lunar tide alone with the Earth's effective k2/Q of
    0.3/12 and the standard 3/2 torque gives it to 30 percent."""
    planet = spin_torques.Planet(EARTH_MASS, EARTH_RADIUS)
    moon = spin_torques.Raiser(MOON_MASS, MOON_A, spin_torques.orbital_rate(EARTH_MASS, MOON_A))
    torques = spin_torques.Torques(planet, (moon,), spin_torques.constant_dissipation(0.30, 12.0), None, 1.5)
    spin = 2.0 * math.pi / 86164.0
    rate = torques.total(spin) / planet.inertia
    lengthening_s_per_century = -rate / spin * 86164.0 * 100.0 * YEAR
    assert lengthening_s_per_century == pytest.approx(2.0e-3, rel=0.3)


def test_a_high_q_barely_moves_the_earth_in_the_age_of_the_system():
    """At a k2 of 0.3 and a Q of 100 (the rock Q the code had before GEN.955) the Moon changes the Earth's day by well under an hour in 4.5 Gy (the spin
    report's 0.1 h); at the Q of tens that the Andrade mantle gives, a Moon that stays where it is would lock the Earth, which is why the moon's
    recession (GEN.959) goes with the moon's torque."""
    k2, q = 0.3, 100.0
    planet = spin_torques.Planet(EARTH_MASS, EARTH_RADIUS)
    moon = spin_torques.Raiser(MOON_MASS, MOON_A, spin_torques.orbital_rate(EARTH_MASS, MOON_A))
    sun = spin_torques.Raiser(SUN, AU, spin_torques.orbital_rate(SUN, AU))
    spin = 2.0 * math.pi / (5.0 * 3600.0)
    with_moon = spin_torques.evolve(spin_torques.Torques(planet, (sun, moon), spin_torques.constant_dissipation(k2, q)), spin, 4.5e9 * YEAR)
    without = spin_torques.evolve(spin_torques.Torques(planet, (sun,), spin_torques.constant_dissipation(k2, q)), spin, 4.5e9 * YEAR)
    day = lambda w: 2.0 * math.pi / w / 3600.0
    assert 0.0 < day(with_moon.rate) - day(without.rate) < 1.0


def test_a_moon_inside_corotation_spins_the_planet_up():
    """A moon that orbits faster than the planet turns drags the planet's spin up (Phobos, Triton) and tends to catch it at its own orbital rate."""
    planet = spin_torques.Planet(EARTH_MASS, EARTH_RADIUS)
    close = spin_torques.Raiser(1.0e22, 2.0e7, spin_torques.orbital_rate(EARTH_MASS, 2.0e7))
    far_star = spin_torques.Raiser(SUN, 5.0 * AU, spin_torques.orbital_rate(SUN, 5.0 * AU))
    torques = spin_torques.Torques(planet, (far_star, close), spin_torques.constant_dissipation(0.3, 100.0))
    assert torques.total(0.1 * close.orbital_rate) > 0.0 > torques.total(10.0 * close.orbital_rate)
    result = spin_torques.evolve(torques, 0.0, 1.0e12 * YEAR)
    assert (result.kind, result.rate) == ("corotation", pytest.approx(close.orbital_rate))


def test_thermal_amplitude_and_relaxation_reproduce_leconte_table_1():
    bar = 1.0e5
    assert spin_torques.thermal_pressure_amplitude_pa(1.0 * bar, 1366.0) == pytest.approx(1180.0, rel=0.02)
    assert spin_torques.thermal_pressure_amplitude_pa(10.0 * bar, 1366.0) == pytest.approx(4050.0, rel=0.02)
    assert spin_torques.thermal_pressure_amplitude_pa(92.0 * bar, 2610.0) == pytest.approx(201.0, rel=0.02)
    assert spin_torques.thermal_relaxation_rate(1.0 * bar, 1366.0) == pytest.approx(2.30e-6, rel=0.02)
    assert spin_torques.thermal_relaxation_rate(10.0 * bar, 1366.0) == pytest.approx(1.46e-6, rel=0.02)
    assert spin_torques.thermal_relaxation_rate(92.0 * bar, 2610.0) == pytest.approx(3.77e-7, rel=0.02)
    # the other flux of the table
    assert spin_torques.thermal_pressure_amplitude_pa(1.0 * bar, 450.0) == pytest.approx(890.0, rel=0.1)
    assert spin_torques.thermal_pressure_amplitude_pa(10.0 * bar, 450.0) == pytest.approx(2960.0, rel=0.1)
    assert spin_torques.thermal_relaxation_rate(1.0 * bar, 450.0) == pytest.approx(1.18e-6, rel=0.25)
    assert spin_torques.thermal_relaxation_rate(10.0 * bar, 450.0) == pytest.approx(7.17e-7, rel=0.25)


def test_thin_air_has_a_proportionally_smaller_thermal_tide():
    assert spin_torques.thermal_pressure_amplitude_pa(1.0e3, 1366.0) == pytest.approx(11.8, rel=0.05)
    assert spin_torques.thermal_pressure_amplitude_pa(1.0e3, 1366.0) < spin_torques.thermal_pressure_amplitude_pa(1.0e5, 1366.0)


def _venus(factor=1.5, q=100.0, k2=0.25):
    planet = spin_torques.Planet(4.867e24, 6.052e6)
    a = 0.723 * AU
    star = spin_torques.Raiser(SUN, a, spin_torques.orbital_rate(SUN, a))
    pressure, flux = 92.0e5, 2610.0
    thermal = spin_torques.ThermalTide(star, spin_torques.thermal_pressure_amplitude_pa(pressure, flux),
                                       spin_torques.thermal_relaxation_rate(pressure, flux))
    return spin_torques.Torques(planet, (star,), spin_torques.constant_dissipation(k2, q), thermal, factor), star, thermal


def test_venus_thermal_torque_is_the_size_of_the_gravitational_one():
    """Leconte et al. 2015 (supplement section 3): Venus's thermal torque is roughly equal to the gravitational one for Q of about 100."""
    torques, star, thermal = _venus()
    peak_thermal = torques.atmospheric(star.orbital_rate + thermal.relaxation_rate)
    gravitational = abs(torques.gravitational(star.orbital_rate + thermal.relaxation_rate, star))
    assert 0.3 < peak_thermal / gravitational < 3.0
    # it opposes the gravitational torque, on both sides of synchronous rotation
    for spin in (star.orbital_rate + thermal.relaxation_rate, star.orbital_rate - thermal.relaxation_rate):
        assert torques.atmospheric(spin) * torques.gravitational(spin, star) < 0.0


def test_the_thermal_torque_is_zero_at_synchronous_rotation_and_vanishes_at_high_frequency():
    torques, star, thermal = _venus()
    assert torques.atmospheric(star.orbital_rate) == 0.0
    peak = abs(torques.atmospheric(star.orbital_rate + thermal.relaxation_rate))
    assert abs(torques.atmospheric(star.orbital_rate + 1.0e3 * thermal.relaxation_rate)) < 1.0e-2 * peak


def test_the_asynchronous_equilibria_are_leconte_equation_13():
    """With constant Q the balance of the Debye thermal torque and the sign-flipping gravitational one is `x = (1 +- sqrt(1 - 4 c^2)) / (2 c)`,
    the outer one stable, in units of omega0 from the orbital rate, on both sides (the mirror pair)."""
    planet = spin_torques.Planet(EARTH_MASS, EARTH_RADIUS)
    a = 0.1 * AU
    star = spin_torques.Raiser(0.3 * SUN, a, spin_torques.orbital_rate(0.3 * SUN, a))
    k2, q = 0.3, 3000.0
    thermal = spin_torques.ThermalTide(star, 3000.0, 1.0e-6)
    torques = spin_torques.Torques(planet, (star,), spin_torques.constant_dissipation(k2, q), thermal, 1.5)
    k_g = constants.G * star.mass_kg ** 2 * planet.radius_m ** 5 / a ** 6
    k_a = 3.0 * star.mass_kg * planet.radius_m ** 3 / (5.0 * planet.mean_density * a ** 3)
    c = k_g * (k2 / q) / (k_a * thermal.amplitude_pa)
    assert c < 0.5
    expected_outer = (1.0 + math.sqrt(1.0 - 4.0 * c * c)) / (2.0 * c)
    expected_inner = (1.0 - math.sqrt(1.0 - 4.0 * c * c)) / (2.0 * c)
    found = spin_torques.equilibria(torques)
    stable = sorted(rate for rate, kind in found if kind == "asynchronous")
    unstable = sorted(rate for rate, kind in found if kind == "unstable")
    n, w0 = star.orbital_rate, thermal.relaxation_rate
    assert stable == pytest.approx([n - expected_outer * w0, n + expected_outer * w0], rel=1.0e-6)
    assert unstable == pytest.approx([n - expected_inner * w0, n + expected_inner * w0], rel=1.0e-6)
    assert (n, "synchronous") in found


def test_a_weak_thermal_tide_leaves_only_the_synchronous_state():
    planet = spin_torques.Planet(EARTH_MASS, EARTH_RADIUS)
    a = 0.1 * AU
    star = spin_torques.Raiser(0.3 * SUN, a, spin_torques.orbital_rate(0.3 * SUN, a))
    thermal = spin_torques.ThermalTide(star, 10.0, 1.0e-6)
    torques = spin_torques.Torques(planet, (star,), spin_torques.constant_dissipation(0.3, 100.0), thermal, 1.5)
    assert [kind for rate, kind in spin_torques.equilibria(torques)] == ["synchronous"]


def test_a_strong_thermal_tide_holds_a_fast_planet_away_from_synchronous_rotation():
    planet = spin_torques.Planet(EARTH_MASS, EARTH_RADIUS)
    a = 0.1 * AU
    star = spin_torques.Raiser(0.3 * SUN, a, spin_torques.orbital_rate(0.3 * SUN, a))
    thermal = spin_torques.ThermalTide(star, 3000.0, 1.0e-6)
    torques = spin_torques.Torques(planet, (star,), spin_torques.constant_dissipation(0.3, 3000.0), thermal, 1.5)
    n = star.orbital_rate
    result = spin_torques.evolve(torques, 20.0 * n, 10.0e9 * YEAR)
    assert result.kind == "asynchronous" and result.ratio is None
    assert abs(result.rate - n) > thermal.relaxation_rate
    retrograde = spin_torques.evolve(torques, -20.0 * n, 10.0e9 * YEAR)
    assert retrograde.kind == "asynchronous" and retrograde.rate < n


def test_the_rheology_makes_the_torque_depend_on_the_tide_frequency():
    """With `tidal_rheology` the gravitational torque is largest where the mantle's relaxation time matches the tide and falls to the cold floor
    (`COLD_Q_CEILING`) at both ends; a constant-Q body has no such dependence."""
    planet = spin_torques.Planet(EARTH_MASS, EARTH_RADIUS)
    a = 0.3 * AU
    star = spin_torques.Raiser(SUN, a, spin_torques.orbital_rate(SUN, a))
    solidus, liquidus = 1700.0, 2150.0
    dissipation = spin_torques.rheology_dissipation(planet, 1700.0, 1.0e12, 2.5e5, solidus, liquidus)
    relaxation = tidal_rheology.rigidity_pa(0.0) / tidal_rheology.viscosity_pa_s(1700.0, 1.0e12, 2.5e5)
    floor = tidal_rheology.rock_response(EARTH_MASS, EARTH_RADIUS, 1.0, 1000.0, 1.0e12, 2.5e5, solidus, liquidus).k2_over_q
    assert dissipation(relaxation) > 3.0 * dissipation(1.0e4 * relaxation)
    assert dissipation(1.0e8 * relaxation) == pytest.approx(floor, rel=0.5)
    torques = spin_torques.Torques(planet, (star,), dissipation)
    n = star.orbital_rate
    assert torques.total(n * 1.5) < 0.0 < torques.total(n * 0.5)
    result = spin_torques.evolve(torques, 3.0 * n, 10.0e9 * YEAR)
    assert math.isfinite(result.rate)


def test_a_balance_is_reached_in_finite_time():
    """Where the thermal and the rheology's gravitational torques cancel away from synchronous rotation the torque is zero, and the flow into it
    (whose integral diverges as a log) must end, not divide by zero: the planet sits within a billionth of the balance and reports it."""
    planet = spin_torques.Planet(3.2e24, 5.32e6)
    a = 4.7e11
    star = spin_torques.Raiser(2.24e30, a, spin_torques.orbital_rate(2.24e30, a))
    thermal = spin_torques.ThermalTide(star, 72.0, 7.9e-7)
    n = star.orbital_rate
    kinds = set()
    for temperature in (900.0, 1200.0, 1500.0, 1800.0):
        dissipation = spin_torques.rheology_dissipation(planet, temperature, 1.0e12, 2.5e5, 1700.0, 2150.0)
        torques = spin_torques.Torques(planet, (star,), dissipation, thermal)
        for start in (0.5 * n, 2.0 * n, 25.0 * n, 300.0 * n):
            result = spin_torques.evolve(torques, start, 10.0e9 * YEAR)
            assert math.isfinite(result.rate)
            kinds.add(result.kind)
            if result.kind == "asynchronous":
                assert abs(torques.total(result.rate)) < 1.0e-6 * abs(torques.total(start))
    assert kinds


def test_softness_smooths_the_sign_flip_and_leaves_the_far_torque_alone():
    planet = spin_torques.Planet(EARTH_MASS, EARTH_RADIUS)
    raiser = spin_torques.Raiser(MOON_MASS, MOON_A, spin_torques.orbital_rate(EARTH_MASS, MOON_A))
    sharp = spin_torques.Torques(planet, (raiser,), spin_torques.constant_dissipation(0.3, 100.0), None, 1.5)
    soft = spin_torques.Torques(planet, (raiser,), spin_torques.constant_dissipation(0.3, 100.0), None, 1.5, 0.01)
    n = raiser.orbital_rate
    assert soft.total(n) == 0.0
    assert abs(soft.total(n * 1.0001)) < 0.1 * abs(sharp.total(n * 1.0001))
    assert soft.total(n * 2.0) == pytest.approx(sharp.total(n * 2.0), rel=1.0e-6)
    assert soft.total(n * 0.5) == pytest.approx(sharp.total(n * 0.5), rel=1.0e-6)
