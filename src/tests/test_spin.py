"""
Spin vectors and axial tilts (GEN.104, `planetgen.physics.spin`): the
draws follow docs/design/orbital-updates.md section 6, and every generated
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
from planetgen.physics import constants, spin
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
    impact = [spin.impact_tilt_deg() for _ in range(2000)]
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
                assert body.axial_tilt_deg == pytest.approx(0.0, abs=1e-9)  # locked: upright
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
