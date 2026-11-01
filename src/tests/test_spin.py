"""
Spin vectors and axial tilts (GEN.104, `planetgen.physics.spin`): the
draws follow docs/design/orbits/orbital-updates.md section 6, and every generated
star, planet, moon, comet, remnant, rogue planet and interstellar comet
carries a unit spin axis whose angle from its reference normal is its tilt.
"""

import math
import statistics

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen.physics import chemically_peculiar, constants, spin
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy

SUN_RADIUS_KM = 695_700.0


def _norm(v):
    return math.sqrt(sum(c * c for c in v))


def _axis(body):
    return (body.spin_axis_x, body.spin_axis_y, body.spin_axis_z)


def _assert_spin(body, normal):
    axis = _axis(body)
    assert _norm(axis) == pytest.approx(1.0, abs=1e-9)
    assert 0.0 <= body.axial_tilt_deg <= 180.0
    assert spin.obliquity_between(axis, normal) == pytest.approx(body.axial_tilt_deg, abs=1e-6)


def test_orbit_normal_matches_the_orbit_frame():
    assert spin.orbit_normal(0.0, 123.0) == pytest.approx((0.0, 0.0, 1.0))
    # Inclined 90 degrees with the node on +x: the orbit runs through z
    # and y, so its normal lies along -y.
    assert spin.orbit_normal(90.0, 0.0) == pytest.approx((0.0, -1.0, 0.0), abs=1e-12)


@pytest.mark.parametrize("obliquity", [0.0, 23.4, 90.0, 177.0])
def test_tilted_axis_keeps_its_obliquity(obliquity):
    normal = spin.orbit_normal(30.0, 40.0)
    for precession in (0.0, 1.0, 4.0):
        axis = spin.tilted_axis(normal, obliquity, precession)
        assert _norm(axis) == pytest.approx(1.0)
        assert spin.obliquity_between(axis, normal) == pytest.approx(obliquity, abs=1e-6)


def test_tilt_distributions():
    draw.set_run_seed(104)
    rayleigh = [spin.rayleigh_tilt_deg() for _ in range(4000)]
    assert statistics.fmean(rayleigh) == pytest.approx(
        spin.STELLAR_TILT_SIGMA_DEG * math.sqrt(math.pi / 2), rel=0.05)
    impact = [spin.accretion_spin(5.972e24, 6371.0)[1] for _ in range(2000)]
    assert all(0.0 <= t <= 180.0 for t in impact)
    assert sum(t > 90.0 for t in impact) > 0  # some retrograde spinners
    yorp = [spin.yorp_tilt_deg() for _ in range(2000)]
    assert all(min(abs(t - 10.0), abs(t - 170.0)) < 45.0 for t in yorp)
    assert 0.3 < sum(t > 90.0 for t in yorp) / len(yorp) < 0.7


def test_the_sun_spins_in_about_a_month():
    draw.set_run_seed(1)
    periods = [spin.star_rotation_period_hours(constants.SOLAR_MASS_TO_KG, SUN_RADIUS_KM, 5772, 4.6, "V")
               for _ in range(500)]
    median_days = statistics.median(periods) / 24.0
    assert 20.0 < median_days < 30.0


def test_hot_stars_stay_under_breakup():
    draw.set_run_seed(2)
    mass = 8 * constants.SOLAR_MASS_TO_KG
    radius = 4 * SUN_RADIUS_KM
    cap = spin.BREAKUP_FRACTION * spin.breakup_speed_kms(mass, radius)
    for _ in range(2000):
        period = spin.star_rotation_period_hours(mass, radius, 22000, 0.02, "V")
        speed = 2 * math.pi * radius / (period * 3600)
        assert speed <= cap * (1 + 1e-9)


def _median_speed_kms(temperature_k, mass_sol, age_gy=3.0, radius_sol=1.4, draws=300):
    # each draw is a star of its own (GEN.1027 keys an Ap or Am dwarf's spin to its mass, temperature and age); the peculiar ones, whose slow
    # spin is `chemically_peculiar`'s and not the cascade's speed law, are left out
    radius = radius_sol * SUN_RADIUS_KM
    masses = [mass_sol * constants.SOLAR_MASS_TO_KG * (1.0 + 1.0e-7 * k) for k in range(draws)]
    speeds = [2 * math.pi * radius / (spin.star_rotation_period_hours(mass, radius, temperature_k, age_gy, "V") * 3600)
              for mass in masses if chemically_peculiar.of_star("V", mass, radius, temperature_k, age_gy) is None]
    return statistics.median(speeds)


def test_a_dwarf_cooler_than_the_kraft_break_spins_in_days_not_hours():
    # GEN.1163: a 6200 K to 6450 K dwarf was given the hot-star speed (180 km/s); the stars there turn at a few km/s [R0378]
    draw.set_run_seed(11)
    for temperature, mass in ((6000, 1.1), (6200, 1.15), (6300, 1.2), (6450, 1.3)):
        speed = _median_speed_kms(temperature, mass)
        assert 2.0 < speed < 12.0, (temperature, speed)


def test_the_speed_climbs_through_the_kraft_break_to_the_hot_star_speed():
    draw.set_run_seed(12)
    speeds = [_median_speed_kms(t, m) for t, m in ((6450, 1.3), (6500, 1.35), (6550, 1.38), (6600, 1.4), (6650, 1.45), (7000, 1.55), (7500, 1.7), (8500, 2.0))]
    assert all(later > earlier for earlier, later in zip(speeds, speeds[1:]))
    assert speeds[0] < 12.0
    assert speeds[4] == pytest.approx(spin.hot_star_median_speed_kms(1.45), rel=0.2)
    assert speeds[-1] == pytest.approx(spin.HOT_STAR_SPEED_KMS, rel=0.2)


def test_the_hot_star_median_runs_from_the_break_speed_to_the_full_speed_over_mass():
    assert spin.hot_star_median_speed_kms(1.0) == pytest.approx(spin.HOT_STAR_BREAK_SPEED_KMS)
    assert spin.hot_star_median_speed_kms(spin.KRAFT_BREAK_TOP_MASS_SOL) == pytest.approx(spin.HOT_STAR_BREAK_SPEED_KMS)
    assert spin.hot_star_median_speed_kms(spin.HOT_STAR_FULL_SPEED_MASS_SOL) == pytest.approx(spin.HOT_STAR_SPEED_KMS)
    assert spin.hot_star_median_speed_kms(spin.HOT_STAR_PEAK_MASS_SOL) == pytest.approx(spin.HOT_STAR_PEAK_SPEED_KMS)
    assert spin.hot_star_median_speed_kms(spin.HOT_STAR_FALL_MASS_SOL) == pytest.approx(spin.HOT_STAR_SPEED_KMS)
    assert spin.hot_star_median_speed_kms(20.0) == pytest.approx(spin.HOT_STAR_SPEED_KMS)
    rising = [spin.hot_star_median_speed_kms(m) for m in (1.8, 2.5, 3.0, 3.5)]      # GEN.1170: Zorec and Royer's rise to the 3.15 to 3.85 solar mass bin
    assert all(b > a for a, b in zip(rising, rising[1:]))
    falling = [spin.hot_star_median_speed_kms(m) for m in (3.5, 5.0, 6.5, 8.0)]
    assert all(b < a for a, b in zip(falling, falling[1:]))
    masses = [1.41, 1.5, 1.6, 1.7, 1.8]
    values = [spin.hot_star_median_speed_kms(m) for m in masses]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_the_kraft_weight_is_zero_below_one_above_and_smooth_between():
    assert spin.kraft_weight(5772) == 0.0
    assert spin.kraft_weight(spin.KRAFT_BREAK_LOW_K) == 0.0
    assert spin.kraft_weight(6550) == pytest.approx(0.5)
    assert spin.kraft_weight(spin.KRAFT_BREAK_HIGH_K) == 1.0
    assert spin.kraft_weight(22000) == 1.0
    grid = [spin.kraft_weight(t) for t in range(6400, 6700, 10)]
    assert all(b >= a for a, b in zip(grid, grid[1:]))


def test_the_gyrochronology_period_is_given_only_below_the_break():
    sun = constants.SOLAR_MASS_TO_KG
    assert spin.nominal_rotation_period_hours(sun, 5772, 4.6) == pytest.approx(spin.gyro_median_days(spin.b_minus_v(5772) - spin.GYRO_C, 4.6) * 24.0)
    assert spin.nominal_rotation_period_hours(1.25 * sun, 6300, 3.0) is not None
    assert spin.nominal_rotation_period_hours(1.4 * sun, 6550, 3.0) is None
    assert spin.nominal_rotation_period_hours(2.0 * sun, 8500, 1.0) is None


def test_small_bodies_respect_the_spin_barrier():
    draw.set_run_seed(3)
    assert min(spin.small_body_period_hours() for _ in range(5000)) >= spin.SPIN_BARRIER_HOURS


def test_black_hole_spin_follows_the_beta_draw():
    draw.set_run_seed(4)
    spins = [spin.black_hole_spin() for _ in range(4000)]
    assert all(0.0 <= a <= spin.BLACK_HOLE_MAX_SPIN for a in spins)
    alpha, beta = spin.BLACK_HOLE_SPIN_BETA
    assert statistics.fmean(spins) == pytest.approx(alpha / (alpha + beta), abs=0.02)
    assert spin.black_hole_horizon_period_s(10 * constants.SOLAR_MASS_TO_KG, 0.0) == math.inf


@pytest.mark.parametrize("seed", range(6))
def test_every_generated_system_body_spins(seed):
    with deterministic_entropy(seed):
        system = StarSystem(system_config=SystemConfig())
    for star in system.stars:
        _assert_spin(star, spin.GALACTIC_POLE)
        assert star.rotation_period_hours > 0
    for planet in system.planets:
        if not hasattr(planet, "moons"):
            continue  # an asteroid belt
        for body in [planet, *planet.moons]:
            _assert_spin(body, spin.orbit_normal(body.orbital_inclination_deg, body.orbital_ascending_node_deg))
            orbit_hours = body.period * constants.SECONDS_PER_YEAR / 3600
            if body.rotation_period_hours == pytest.approx(orbit_hours):
                assert body.axial_tilt_deg <= min(body.orbital_inclination_deg, 180.0 - body.orbital_inclination_deg) + 1e-9  # locked: Cassini state 1 (GEN.1165)
    for comet in system.comets:
        _assert_spin(comet, spin.orbit_normal(comet.inclination_deg, comet.ascending_node_deg))
        assert comet.rotation_period_hours >= spin.SPIN_BARRIER_HOURS


def test_close_in_planets_lock_to_their_star():
    locked = 0
    for seed in range(20):
        with deterministic_entropy(seed):
            system = StarSystem(system_config=SystemConfig())
        for planet in system.planets:
            if not hasattr(planet, "moons"):
                continue
            orbit_hours = planet.period * constants.SECONDS_PER_YEAR / 3600
            locked += planet.rotation_period_hours == pytest.approx(orbit_hours)
    assert locked > 0


@pytest.mark.parametrize("cls", [BlackHole, NeutronStar, RoguePlanet, InterstellarComet])
def test_free_objects_spin_about_the_galactic_pole(cls):
    with deterministic_entropy(7):
        body = cls(SystemConfig())
    _assert_spin(body, spin.GALACTIC_POLE)
    restored = cls.from_dict(body.to_dict(), SystemConfig())
    assert _axis(restored) == _axis(body)
    assert restored.axial_tilt_deg == body.axial_tilt_deg


def test_rocky_rogue_planets_are_born_at_the_accretion_spin_and_giants_keep_the_gas_giant_day():
    """GEN.842: a rocky rogue planet turns at the rate its giant impacts left (`spin.accretion_spin`, GEN.731), so its day is not held to the
    old 8 to 48 hours, and never turns faster than the breakup fraction; a rogue giant keeps the gas-giant day of 8 to 20 hours."""
    rocky, giants = [], []
    for seed in range(300):
        with deterministic_entropy(seed):
            body = RoguePlanet(SystemConfig())
        if body.planet_type == 't':
            critical = spin.critical_spin_rad_s(body.mass_kg, body.radius_km)
            assert body.rotation_period_hours >= 2.0 * math.pi / (spin.BREAKUP_FRACTION * critical) / 3600.0 * (1.0 - 1e-9)
            rocky.append(body.rotation_period_hours)
        else:
            giants.append(body.rotation_period_hours)
    low, high = constants.ROTATION_PERIOD_RANGE_HOURS['t']
    assert rocky and any(not low <= hours <= high for hours in rocky), "the old log-uniform law never leaves 8 to 48 hours"
    glow, ghigh = constants.ROTATION_PERIOD_RANGE_HOURS['g']
    assert all(glow <= hours <= ghigh for hours in giants)
