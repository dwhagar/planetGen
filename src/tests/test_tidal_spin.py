# tests/test_tidal_spin.py

"""
GEN.731: a planet's birth spin from the giant impacts that built it (`spin.accretion_spin`), the tide that wears it down with the body's own
Love number and quality factor, and the spin-orbit resonances it can be caught in (`tidal_spin`, `planets.generate_spin`): Mercury, Venus,
Earth, a hot Jupiter and the Moon.
"""

import math
import statistics
import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants as pc
from planetgen.physics import spin, tidal_heating, tidal_spin
from planetgen.physics.planets import generate_spin
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy

SUN = 1.989e30
EARTH_KG = 5.972e24
GYR_S = 1.0e9 * pc.SECONDS_PER_YEAR
BIRTH = 2.0 * math.pi / (5.0 * 3600.0)
"""float: A birth spin of five hours, rad/s (Kokubo and Genda's planets turn at a median near three)."""


def _body(mass_kg, radius_km, body_type="t", **fields):
    return types.SimpleNamespace(mass=mass_kg, radius=radius_km, body_type=body_type, **fields)


MERCURY = _body(3.30e23, 2439.7)
VENUS = _body(4.867e24, 6051.8)
EARTH = _body(EARTH_KG, 6371.0)
MOON = _body(7.346e22, 1737.4)
HOT_JUPITER = _body(1.898e27, 1.2 * 71492.0, "g")
JUPITER = _body(1.898e27, 71492.0, "g")


def _evolve(body, host_kg, au, eccentricity, age_gyr=4.5, birth=BIRTH):
    return tidal_spin.evolve(body, host_kg, au * pc.AU_TO_M, birth, age_gyr * GYR_S, eccentricity)


def _hours(rate):
    return 2.0 * math.pi / rate / 3600.0


def test_the_accretion_spin_follows_the_published_law():
    """Kokubo and Genda 2010: Gaussian components of 1.3 rad/h at 3000 kg/m3 (rms 2.35 rad/h, 0.71 of the breakup rate), an isotropic axis."""
    draw.set_run_seed(731)
    mass = 1.0e25
    radius_km = (3.0 * mass / (4.0 * math.pi * 3000.0)) ** (1.0 / 3.0) / 1000.0  # the paper's 3000 kg/m3
    samples = [spin.accretion_spin(mass, radius_km) for _ in range(6000)]
    critical = spin.critical_spin_rad_s(mass, radius_km)
    rates = [rate / critical for rate, _ in samples]
    assert all(rate <= spin.BREAKUP_FRACTION + 1e-12 for rate in rates)
    # the uncut rms is 0.41 sqrt(3) = 0.71 of breakup; the cut at 0.85 trims it a little
    assert statistics.fmean(rate ** 2 for rate in rates) ** 0.5 == pytest.approx(0.69, abs=0.05)
    # isotropic: cos(tilt) is uniform in -1..1
    cosines = [math.cos(math.radians(tilt)) for _, tilt in samples]
    assert statistics.fmean(cosines) == pytest.approx(0.0, abs=0.03)
    assert statistics.pstdev(cosines) == pytest.approx(1.0 / math.sqrt(3.0), abs=0.03)


def test_the_paper_planet_turns_in_about_three_hours():
    draw.set_run_seed(5)
    mass = 5.0e24
    radius_km = (3.0 * mass / (4.0 * math.pi * 3000.0)) ** (1.0 / 3.0) / 1000.0
    hours = [_hours(spin.accretion_spin(mass, radius_km)[0]) for _ in range(4000)]
    assert 2.0 < statistics.median(hours) < 4.0


def test_the_rheology_is_each_bodys_own():
    k2, q = tidal_heating.response(EARTH)
    assert (q, k2) == (100.0, pytest.approx(tidal_heating.rigid_k2(EARTH_KG, 6.371e6)))
    assert tidal_heating.response(MERCURY)[0] < tidal_heating.response(EARTH)[0]  # a small body is stiffer
    assert tidal_heating.response(MOON)[0] == pytest.approx(0.021, rel=0.2)
    assert tidal_heating.response(JUPITER) == (tidal_heating.GIANT_K2, tidal_heating.GIANT_TIDAL_Q)
    icy = _body(1.0e22, 1000.0, water_mass_fraction=0.5)
    assert tidal_heating.response(icy) == (tidal_heating.ICY_K2, 100.0)


def test_mercury_is_caught_in_the_three_halves_resonance():
    """At its eccentricity of 0.2056 and its age Mercury is brought down from a fast spin into 3:2: 58.6 days, three turns per two orbits."""
    rate, ratio = _evolve(MERCURY, SUN, 0.387, 0.2056)
    assert ratio == 1.5
    period_h = 87.969 * 24.0
    assert _hours(rate) == pytest.approx(period_h / 1.5, rel=2e-3)
    assert spin.spin_state(_hours(rate), 87.969 / 365.25) == ("resonant", 1.5)
    # the same planet on a nearly round orbit is locked 1:1, and one not eccentric enough for 2:1 passes it
    assert _evolve(MERCURY, SUN, 0.387, 0.02)[1] == 1.0
    assert _evolve(MERCURY, SUN, 0.387, 0.40)[1] == 2.0


def test_the_resonance_needs_the_time_to_reach_it():
    """Younger than the tide takes, Mercury is still on its way down and has been caught by nothing."""
    young_rate, young_ratio = _evolve(MERCURY, SUN, 0.387, 0.2056, age_gyr=0.5)
    assert young_ratio is None
    assert young_rate < BIRTH and _hours(young_rate) < 87.969 * 24.0 / 1.5
    assert tidal_spin.locking_time_s(MERCURY, SUN, 0.387 * pc.AU_TO_M, BIRTH, 0.2056) / GYR_S == pytest.approx(3.8, rel=0.25)


def test_venus_and_earth_are_not_locked_by_the_sun():
    """Their locking times against the sun's tide alone are tens to hundreds of Gyr: they keep their spin (Venus's 243 days is the
    atmospheric tide, not in this model)."""
    for body, au, e in ((VENUS, 0.723, 0.0068), (EARTH, 1.0, 0.0167)):
        rate, ratio = _evolve(body, SUN, au, e)
        assert ratio is None
        assert tidal_spin.locking_time_s(body, SUN, au * pc.AU_TO_M, BIRTH, e) > 10.0 * GYR_S
        assert 0.8 * BIRTH < rate <= BIRTH  # worn down by a few percent to a tenth, no more
        assert spin.spin_state(_hours(rate), au ** 1.5)[0] == "free"


def test_a_hot_jupiter_is_locked_within_a_few_million_years():
    rate, ratio = _evolve(HOT_JUPITER, SUN, 0.05, 0.0, age_gyr=0.1, birth=2.0 * math.pi / (10.0 * 3600.0))
    assert ratio == 1.0
    assert _hours(rate) == pytest.approx(0.05 ** 1.5 * pc.SECONDS_PER_YEAR / 3600.0, rel=2e-2)
    time = tidal_spin.locking_time_s(HOT_JUPITER, SUN, 0.05 * pc.AU_TO_M, 2.0 * math.pi / (10.0 * 3600.0))
    assert 1.0e5 * pc.SECONDS_PER_YEAR < time < 1.0e8 * pc.SECONDS_PER_YEAR
    # and a Jupiter at 5.2 AU keeps its spin
    assert _evolve(JUPITER, SUN, 5.2, 0.048, birth=2.0 * math.pi / (10.0 * 3600.0))[1] is None


def test_the_moon_is_locked_to_the_earth_in_tens_of_millions_of_years():
    time = tidal_spin.locking_time_s(MOON, EARTH_KG, 3.84399e8, 2.0 * math.pi / (10.0 * 3600.0), 0.0549)
    assert 1.0e7 * pc.SECONDS_PER_YEAR < time < 3.0e8 * pc.SECONDS_PER_YEAR
    assert _evolve(MOON, EARTH_KG, 3.84399e8 / pc.AU_TO_M, 0.0549)[1] == 1.0


def test_a_slow_spinner_is_spun_up_to_the_orbit_and_never_past_it():
    """A hot Jupiter born turning slower than it orbits is carried up to the orbital rate, not through it."""
    slow = 2.0 * math.pi / (400.0 * 3600.0)
    rate, ratio = _evolve(HOT_JUPITER, SUN, 0.03, 0.0, birth=slow)
    assert ratio == 1.0
    orbital = math.sqrt(pc.G * SUN / (0.03 * pc.AU_TO_M) ** 3)
    assert rate == pytest.approx(orbital)


def test_spin_state_reads_the_two_stored_periods():
    year_h = pc.SECONDS_PER_YEAR / 3600.0
    assert spin.spin_state(year_h, 1.0) == ("locked", 1.0)
    assert spin.spin_state(year_h * 2.0 / 3.0, 1.0) == ("resonant", 1.5)
    assert spin.spin_state(year_h / 2.0, 1.0) == ("resonant", 2.0)
    assert spin.spin_state(24.0, 1.0)[0] == "free"
    assert spin.spin_state(year_h * 1.5, 1.0)[0] == "free"  # a spin slower than the orbit is not a resonance
    assert spin.spin_state(None, 1.0) == (None, None) and spin.spin_state(24.0, 0.0) == (None, None)


def _planet(body, au, period_years, eccentricity, age_gyr, is_moon=False):
    planet = types.SimpleNamespace(
        is_moon=is_moon, distance=au, period=period_years, eccentricity=eccentricity, orbital_inclination_deg=1.0,
        orbital_ascending_node_deg=30.0, star=types.SimpleNamespace(age=age_gyr), **vars(body))
    return planet


def test_generate_spin_stores_the_state_the_tide_leaves(monkeypatch):
    """`planets.generate_spin` on Mercury's orbit: the stored rotation period is its resonance's, and the axis is upright."""
    import planetgen.physics.planets as planets
    for name in ("magnetism", "radiation", "habitability_world"):
        monkeypatch.setattr(getattr(planets, name), "generate_field" if name == "magnetism" else "generate_dose" if name == "radiation"
                            else "generate", lambda *a, **k: None)
    monkeypatch.setattr(planets, "internal_heat_flux_w_m2", lambda planet: 0.0)
    monkeypatch.setattr(spin, "accretion_spin", lambda mass, radius: (BIRTH, 40.0))  # a spin of five hours, tilted 40 degrees
    draw.set_run_seed(3)
    mercury = _planet(MERCURY, 0.387, 87.969 / 365.25, 0.2056, 4.5)
    generate_spin(mercury, SUN)
    assert mercury.axial_tilt_deg == 0.0
    assert spin.spin_state(mercury.rotation_period_hours, mercury.period) == ("resonant", 1.5)
    assert mercury.axial_tilt_deg == 0.0
    venus = _planet(VENUS, 0.723, 224.7 / 365.25, 0.0068, 4.5)
    generate_spin(venus, SUN)
    assert spin.spin_state(venus.rotation_period_hours, venus.period)[0] == "free"
    assert venus.axial_tilt_deg == pytest.approx(40.0, abs=1e-6)  # the impacts' tilt, kept: the tide has not finished
    assert _hours(spin.critical_spin_rad_s(VENUS.mass, VENUS.radius)) / spin.BREAKUP_FRACTION <= venus.rotation_period_hours
    hot = _planet(HOT_JUPITER, 0.05, 0.05 ** 1.5, 0.0, 4.5)
    generate_spin(hot, SUN)
    assert spin.spin_state(hot.rotation_period_hours, hot.period) == ("locked", 1.0)


def test_every_generated_planet_is_in_a_state_the_tide_allows():
    seen = {"locked": 0, "resonant": 0, "free": 0}
    for seed in range(12):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets + system.secondary_planets:
            if not hasattr(planet, "moons") or planet.body_type == "a":
                continue
            state, turns = spin.spin_state(planet.rotation_period_hours, planet.period)
            seen[state] += 1
            if state != "free":
                assert planet.axial_tilt_deg == 0.0
            host = planet.primary_mass_kg
            locking = tidal_spin.locking_time_s(planet, host, planet.distance * pc.AU_TO_M, 2.0 * math.pi / (planet.rotation_period_hours * 3600.0),
                                                 planet.eccentricity)
            if state == "free":
                assert locking > 0.0
    assert seen["free"] > 0 and seen["locked"] + seen["resonant"] > 0
