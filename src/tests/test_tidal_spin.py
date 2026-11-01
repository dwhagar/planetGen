# tests/test_tidal_spin.py

"""
GEN.731: a planet's birth spin from the giant impacts that built it (`spin.accretion_spin`), the tide that wears it down with the body's own
Love number and quality factor, and the spin-orbit resonances it can be caught in (`body_spin`, `planets.generate_spin`): Mercury, Venus,
Earth, a hot Jupiter and the Moon.
"""

import math
import statistics
import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import icy_tide, constants as pc
from planetgen.physics import body_spin, impact_spin, interior_heat, spin, spin_torques, tidal_heating, tidal_spin
from planetgen.physics.planets import generate_spin
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy
from tests.tide_support import locking_time_s, rocky

SUN = 1.989e30
EARTH_KG = 5.972e24
GYR_S = 1.0e9 * pc.SECONDS_PER_YEAR
BIRTH = 2.0 * math.pi / (5.0 * 3600.0)
"""float: A birth spin of five hours, rad/s (Kokubo and Genda's planets turn at a median near three)."""


def _body(mass_kg, radius_km, body_type="t", **fields):
    """A giant, or a rocky body with the content of the bulk silicate Earth and a mantle of the middle of the drawn viscosities (GEN.955)."""
    if body_type == "g":
        return types.SimpleNamespace(mass=mass_kg, radius=radius_km, body_type=body_type, **fields)
    return rocky(mass_kg, radius_km, **fields)


MERCURY = _body(3.30e23, 2439.7)
VENUS = _body(4.867e24, 6051.8)
EARTH = _body(EARTH_KG, 6371.0)
MOON = _body(7.346e22, 1737.4)
HOT_JUPITER = _body(1.898e27, 1.2 * 71492.0, "g")
JUPITER = _body(1.898e27, 71492.0, "g")


def _evolve(body, host_kg, au, eccentricity, age_gyr=4.5, birth=BIRTH):
    found = body_spin.evolve(body, host_kg, au * pc.AU_TO_M, birth, age_gyr * GYR_S, eccentricity)
    return found.rate, found.ratio


def _hours(rate):
    return 2.0 * math.pi / rate / 3600.0


def test_the_accretion_spin_follows_the_published_law():
    """[R2431] equation 4 with the last impacts' mass ratio (GEN.979): a planet of a mass the paper's runs reached, after a giant impact, has
    its rms spin at `impact_spin.mean_rms_breakup_fraction` of breakup (0.71 of it at the paper's mass ratio of 0.28), and an isotropic axis."""
    draw.set_run_seed(731)
    mass = 1.0e25
    radius_km = (3.0 * mass / (4.0 * math.pi * 3000.0)) ** (1.0 / 3.0) / 1000.0  # the paper's 3000 kg/m3
    samples = [spin.accretion_spin(mass, radius_km) for _ in range(6000)]
    critical = spin.critical_spin_rad_s(mass, radius_km)
    rates = [rate / critical for rate, _ in samples]
    assert all(rate <= spin.BREAKUP_FRACTION + 1e-12 for rate in rates)
    expected = impact_spin.mean_rms_breakup_fraction(mass / pc.EARTH_MASS_TO_KG)
    assert statistics.fmean(rate ** 2 for rate in rates) ** 0.5 == pytest.approx(expected, rel=0.25)   # the cap trims the top
    cosines = [math.cos(math.radians(tilt)) for _, tilt in samples]
    assert statistics.fmean(cosines) == pytest.approx(0.0, abs=0.05)   # isotropic: cos(tilt) is uniform in -1..1 (the rare small-branch draws lean prograde)


def test_the_planets_of_the_generator_turn_in_hours_to_days():
    """GEN.979: the paper's planet turned in 2.7 hours; Vesta (5.3 h), Ceres (9.1 h) and the Earth after the Moon-forming impact (about 5 h) say slower
    and Mars (24.6 h) slower still: the median free day of 0.1 to 5 Earth masses is 5 to 10 hours."""
    draw.set_run_seed(5)
    hours = []
    for _ in range(4000):
        mass_earth = 10.0 ** (-1.0 + 1.7 * draw.random())
        hours.append(_hours(spin.accretion_spin(mass_earth * pc.EARTH_MASS_TO_KG, 6371.0 * mass_earth ** 0.27)[0]))
    assert 4.5 < statistics.median(hours) < 10.0


def test_a_planet_the_size_of_mars_can_turn_as_slowly_as_mars_does():
    """Mars (0.107 Earth masses, a day of 24.6 hours, 0.07 of breakup) is drawn about as often as not; a planet of three Earths seldom does (under one in five)."""
    draw.set_run_seed(9)
    mars_mass, mars_radius_km = 0.107 * 5.972e24, 3389.5
    earth3_mass, earth3_radius_km = 3.0 * 5.972e24, 9400.0
    mars = [_hours(spin.accretion_spin(mars_mass, mars_radius_km)[0]) for _ in range(3000)]
    big = [_hours(spin.accretion_spin(earth3_mass, earth3_radius_km)[0]) for _ in range(3000)]
    assert 0.25 < sum(h >= 24.0 for h in mars) / len(mars) < 0.8
    assert sum(h >= 24.0 for h in big) / len(big) < 0.2 < sum(h >= 24.0 for h in mars) / len(mars)


def test_the_rheology_is_each_bodys_own():
    """`interior_heat.tidal_response`: the Andrade `k2` and `Q` of the body's mantle at the orbital frequency (GEN.955); the fixed pairs of the giant
    and the ice; a stiffer small body, the Moon's `k2` of 0.024 and a `Q` of tens."""
    n = 2.0e-7
    k2, q = interior_heat.tidal_response(EARTH, n)
    assert k2 == pytest.approx(tidal_heating.rigid_k2(EARTH_KG, 6.371e6), rel=0.3) and 1.0 < q < 100.0
    assert interior_heat.tidal_response(MERCURY, n)[0] < interior_heat.tidal_response(EARTH, n)[0]  # a small body is stiffer
    moon_k2, moon_q = interior_heat.tidal_response(MOON, n)
    assert moon_k2 == pytest.approx(0.0242, rel=0.2) and 10.0 < moon_q < 100.0
    assert interior_heat.tidal_response(JUPITER, n) == (tidal_heating.GIANT_K2, tidal_heating.GIANT_TIDAL_Q)
    icy = _body(1.0e22, 1000.0, water_mass_fraction=0.5)
    k2, q = interior_heat.tidal_response(icy, n)      # the ice shell's (GEN.960): the homogeneous ice sphere's k2 and a finite Q
    assert k2 == pytest.approx(icy_tide.love_k2(1.0e22, 1.0e6, n)) and 1.0 < q < math.inf


def test_a_hotter_mantle_dissipates_more_and_a_fast_tide_less():
    cool = rocky(EARTH_KG, 6371.0, eta0_pa_s=1.0e13, activation_kj_mol=300.0)
    warm = rocky(EARTH_KG, 6371.0, eta0_pa_s=2.6e10, activation_kj_mol=200.0)
    assert interior_heat.tide_temperature_k(cool) != interior_heat.tide_temperature_k(warm)
    slow, fast = 1.0e-8, 1.0e-4
    assert interior_heat.tidal_response(cool, fast)[1] > interior_heat.tidal_response(cool, slow)[1]


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
    young_rate, young_ratio = _evolve(MERCURY, SUN, 0.387, 0.2056, age_gyr=0.1)
    assert young_ratio is None
    assert young_rate < BIRTH and _hours(young_rate) < 87.969 * 24.0 / 1.5
    assert 0.1 < locking_time_s(MERCURY, SUN, 0.387 * pc.AU_TO_M, BIRTH, 0.2056) / GYR_S < 10.0


def test_earth_is_not_locked_by_the_sun():
    """Its locking time against the sun's tide alone is tens of Gyr: it keeps its spin. (Venus's, with the Q of a hot mantle, is a few Gyr, which its
    243-day retrograde spin shows is too short: the atmospheric tide that holds it away is the spin torque sum, GEN.957.)"""
    rate, ratio = _evolve(EARTH, SUN, 1.0, 0.0167)
    assert ratio is None
    assert locking_time_s(EARTH, SUN, pc.AU_TO_M, BIRTH, 0.0167) > 5.0 * GYR_S
    assert 0.5 * BIRTH < rate <= BIRTH  # worn down by a tenth to a half, no more
    assert spin.spin_state(_hours(rate), 1.0)[0] == "free"


def test_a_hot_jupiter_is_locked_within_a_few_million_years():
    rate, ratio = _evolve(HOT_JUPITER, SUN, 0.05, 0.0, age_gyr=0.1, birth=2.0 * math.pi / (10.0 * 3600.0))
    assert ratio == 1.0
    assert _hours(rate) == pytest.approx(0.05 ** 1.5 * pc.SECONDS_PER_YEAR / 3600.0, rel=2e-2)
    time = locking_time_s(HOT_JUPITER, SUN, 0.05 * pc.AU_TO_M, 2.0 * math.pi / (10.0 * 3600.0))
    assert 1.0e5 * pc.SECONDS_PER_YEAR < time < 1.0e8 * pc.SECONDS_PER_YEAR
    # and a Jupiter at 5.2 AU keeps its spin
    assert _evolve(JUPITER, SUN, 5.2, 0.048, birth=2.0 * math.pi / (10.0 * 3600.0))[1] is None


def test_the_moon_is_locked_to_the_earth_in_tens_of_millions_of_years():
    time = locking_time_s(MOON, EARTH_KG, 3.84399e8, 2.0 * math.pi / (10.0 * 3600.0), 0.0549)
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
    fields = {**vars(body), "star": types.SimpleNamespace(age=age_gyr)}
    return types.SimpleNamespace(
        is_moon=is_moon, distance=au, period=period_years, eccentricity=eccentricity, orbital_inclination_deg=1.0,
        orbital_ascending_node_deg=30.0, **fields)


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
    generate_spin(mercury, SUN, birth=(BIRTH, 40.0, spin_torques.ALWAYS_CAPTURED))   # a body whose draws take the 3:2 (GEN.1172: capture there is a chance)
    low, high = tidal_spin.CASSINI_COUPLING_RANGE
    assert low * 1.0 <= mercury.axial_tilt_deg <= high * 1.0    # GEN.1165: a share of its 1 degree inclination, not exactly 0
    assert spin.spin_state(mercury.rotation_period_hours, mercury.period) == ("resonant", 1.5)
    earth = _planet(EARTH, 1.0, 1.0, 0.0167, 4.5)
    generate_spin(earth, SUN)
    assert spin.spin_state(earth.rotation_period_hours, earth.period)[0] == "free"
    assert earth.axial_tilt_deg == pytest.approx(40.0, abs=1e-6)  # the impacts' tilt, kept: the tide has not finished
    assert _hours(spin.critical_spin_rad_s(EARTH.mass, EARTH.radius)) / spin.BREAKUP_FRACTION <= earth.rotation_period_hours
    hot = _planet(HOT_JUPITER, 0.05, 0.05 ** 1.5, 0.0, 4.5)
    generate_spin(hot, SUN)
    assert spin.spin_state(hot.rotation_period_hours, hot.period) == ("locked", 1.0)


def test_every_generated_planet_is_in_a_state_the_tide_allows():
    seen = {"locked": 0, "resonant": 0, "free": 0}
    for seed in range(7100, 7112):    # the seeds of the sheet audit: with planets on, a mix of dwarf and Sun-like stars
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            config = SystemConfig()
            config.PLANETS = True
            config.MAX_PLANETS = True
            system = StarSystem(system_config=config, galactic_center_dist_ly=26000.0)
        for planet in system.planets + system.secondary_planets:
            if not hasattr(planet, "moons") or planet.body_type == "a":
                continue
            state, turns = spin.spin_state(planet.rotation_period_hours, planet.period)
            seen[state] += 1
            if state != "free":
                assert planet.axial_tilt_deg <= min(planet.orbital_inclination_deg, 180.0 - planet.orbital_inclination_deg) + 1e-9
            host = planet.primary_mass_kg
            locking = locking_time_s(planet, host, planet.distance * pc.AU_TO_M, 2.0 * math.pi / (planet.rotation_period_hours * 3600.0),
                                                 planet.eccentricity)
            if state == "free":
                assert locking > 0.0
    assert seen["free"] > 0 and seen["locked"] + seen["resonant"] > 0


def test_the_cassini_tilt_is_a_share_of_the_inclination():
    """GEN.1165: `eps = x i` with `x` log-uniform over `CASSINI_COUPLING_RANGE`; 0 in the plane; a retrograde orbit counts by its angle from the pole."""
    low, high = tidal_spin.CASSINI_COUPLING_RANGE
    assert tidal_spin.cassini_tilt_deg(2.0, 0.0) == pytest.approx(2.0 * low)
    assert tidal_spin.cassini_tilt_deg(2.0, 1.0) == pytest.approx(2.0 * high)
    assert tidal_spin.cassini_tilt_deg(2.0, 0.5) == pytest.approx(2.0 * math.sqrt(low * high))
    assert tidal_spin.cassini_tilt_deg(0.0, 0.7) == 0.0
    assert tidal_spin.cassini_tilt_deg(178.0, 0.5) == pytest.approx(tidal_spin.cassini_tilt_deg(2.0, 0.5))


def test_a_locked_or_resonant_planet_keeps_a_small_tilt_not_exactly_zero():
    """GEN.1165 / GEN.753 row 7: a locked planet's obliquity is a fraction of a degree to 2 degrees, no longer exactly 0."""
    tilts = []
    for seed in range(24):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            system = StarSystem(SystemConfig(), galactic_center_dist_ly=26000.0)
        for planet in system.planets + system.secondary_planets:
            if not hasattr(planet, "moons") or planet.body_type == "a":
                continue
            if spin.spin_state(planet.rotation_period_hours, planet.period)[0] in ("locked", "resonant"):
                tilts.append(planet.axial_tilt_deg)
    assert len(tilts) >= 20
    tilts.sort()
    assert sum(1 for t in tilts if t > 0.0) >= 0.9 * len(tilts)
    assert tilts[len(tilts) // 2] < 1.0 and tilts[int(0.95 * len(tilts))] < 5.0


def test_a_generated_giant_tail_reaches_past_sixty_degrees_for_an_ice_giant_only():
    """GEN.1165 / GEN.753 row 10: over generated systems an ice giant of 5 to 30 Earth masses is sometimes tilted past 60 degrees; a heavier giant never is."""
    from planetgen.physics import giant_tilt
    ice, gas = [], []
    for seed in range(7100, 7160):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            config = SystemConfig()
            config.PLANETS = True
            config.MAX_PLANETS = True
            system = StarSystem(system_config=config, galactic_center_dist_ly=26000.0)
        for planet in system.planets:
            if getattr(planet, "body_type", None) != "g":
                continue
            low, high = giant_tilt.MASS_RANGE_EARTH
            (ice if low <= planet.mass / pc.EARTH_MASS_TO_KG <= high else gas).append(planet.axial_tilt_deg)
    assert len(ice) >= 4 and len(gas) >= 4
    assert all(0.0 <= tilt <= 180.0 for tilt in ice + gas)
    assert max(gas) < 90.0
