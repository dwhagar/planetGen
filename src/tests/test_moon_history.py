# tests/test_moon_history.py

"""
GEN.959: the history of a planet's moons applied to the bodies (`planetgen.physics.moon_history`): outward drift, the fall of a moon into a ring, the
planet's day from the exchange, and the moons' tide in the planet's heat (`tidal_heating.Forcing.moons`, `interior_heat`).
"""

import math
import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants as pc
from planetgen.physics import interior_heat, moon_history, planets, radiogenic, ring_origin, rings, tidal_heating
from planetgen.physics.mathcheck import _bulk_silicate_earth
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy, first_seed_with

EARTH = pc.EARTH_MASS_TO_KG
MOON = 7.346e22
DAY = 86164.0
ROCK = radiogenic.Rock(0.675, 0.325, 0.0, radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION)
STIFF = interior_heat.Interior(1.0e12, 250.0)


def _heat(spin_days, distance_m=3.844e8, retrograde=False, moons=True):
    forcing = tidal_heating.Forcing(1.989e30, pc.AU_TO_M, 0.0, ((MOON, distance_m, retrograde, (-1.0 if retrograde else 1.0) * math.sqrt(pc.G * (EARTH + MOON) / distance_m ** 3)),) if moons else (), 2.0 * math.pi / (spin_days * DAY))
    return interior_heat.heat_flow(_bulk_silicate_earth(), STIFF, EARTH, 6371.0, 4.5, ROCK, forcing)


def _system(seed):
    draw.set_run_seed(seed)
    config = SystemConfig()
    config.MOONS = True     # TEST.657: after GEN.962 and GEN.1082 to 1084 a rocky planet has a moon only by a giant impact (1 in 12) or capture, so 16 plain systems held 2 moons
    return StarSystem(config, galactic_center_dist_ly=26000.0)


def _systems(seeds=range(16)):
    for seed in seeds:
        with deterministic_entropy(seed):
            yield _system(seed)


def _planet_where(wanted, what):
    """TEST.714: the first planet, over the first 64 seeds, that `wanted` accepts: a case a seed draws rarely is searched for, not pinned to the seed that once had it."""
    def find(system):
        return next((planet for planet in system.planets if planet.body_type != "a" and wanted(planet)), None)

    _seed, found = first_seed_with(lambda seed: find(_system(seed)), lambda planet: planet is not None, what=what)
    return found


def test_the_moon_flexes_the_earth_by_a_terawatt_or_two():
    """A rocky Earth under the Moon dissipates of the order of the solid-body share of the measured 3.7 TW (the oceans take most of it), steeply less farther out
    and next to nothing when the planet turns with the moon; a retrograde planet is flexed as much."""
    plain = _heat(1.0, moons=False)
    near = _heat(1.0)
    assert plain.tidal_w == 0.0
    assert 5.0e11 < near.tidal_w < 1.0e13
    assert _heat(1.0, 1.5 * 3.844e8).tidal_w == pytest.approx(near.tidal_w / 1.5 ** 6, rel=0.2)
    assert _heat(27.3).tidal_w < 1.0e-3 * near.tidal_w
    assert _heat(-1.0).tidal_w == pytest.approx(near.tidal_w, rel=0.2)
    assert near.total_w > plain.total_w


def test_a_planet_with_moons_has_a_forcing_even_on_a_circle():
    planet = types.SimpleNamespace(mass=EARTH, primary_mass_kg=1.989e30, eccentricity=0.0, distance=1.0, rotation_period_hours=24.0, axial_tilt_deg=170.0, is_moon=False,
                                   moons=[types.SimpleNamespace(mass=MOON, distance=0.00257, orbital_inclination_deg=5.0),
                                          types.SimpleNamespace(mass=MOON, distance=0.005, orbital_inclination_deg=150.0)])
    forcing = tidal_heating.forcing_of(planet)
    assert forcing.eccentricity == 0.0 and [m[2] for m in forcing.moons] == [False, True]
    assert forcing.moons[0][3] > 0.0 > forcing.moons[1][3]
    assert forcing.moons[0][:2] == (MOON, 0.00257 * pc.AU_TO_M)
    assert forcing.spin_rate_rad_s == pytest.approx(-2.0 * math.pi / 86400.0)
    planet.moons = []
    assert tidal_heating.forcing_of(planet) is None
    planet.is_moon, planet.moons = True, [types.SimpleNamespace(mass=MOON, distance=0.1, orbital_inclination_deg=0.0)]
    assert tidal_heating.forcing_of(planet) is None


def test_the_ring_of_a_fallen_moon_is_disrupted_moon_spread_for_the_time_since_and_capped():
    planet = ring_origin.Planet(EARTH, 6371.0, 1.0, 1.989e30)
    for fallen in (1.0e19, 1.0e22, 1.0e30):
        found = ring_origin.ring_origin(planet, 4.0e9, (ring_origin.FallenMoon(fallen, 3.9e9),))
        assert found.exists and found.origin == "disrupted_moon"
        assert found.mass_kg == pytest.approx(min(ring_origin.spread(ring_origin.RING_RETAINED_FRACTION * fallen, 1.0e8, planet), rings.MAX_MASS_OVER_PLANET * EARTH))
        assert found.age_yr == pytest.approx(1.0e8)


def test_every_planet_with_moons_ends_with_moons_it_can_keep():
    """Over generated systems: a moon that stays is outside its planet's Roche limit for its own density and inside the stable Hill fraction, its period and
    distance agree, a planet that lost moons to the Roche limit has a ring of one of the three origins whenever it has one (a dust ring only when it is the dusty kind), and the planet's day is a real one."""
    seen_moons = seen_fallen = 0
    for number, system in enumerate(_systems(range(64))):
        if number >= 16 and seen_moons > 5 and seen_fallen > 0:
            break      # TEST.714: the first 16 seeds always run; more only while the sample has not yet held a fallen moon and enough moons
        for planet in system.planets:
            if planet.body_type == "a":
                continue
            fallen = getattr(planet, "fallen_moons", ()) or ()
            if fallen:
                seen_fallen += 1
                if planet.ring_kind is not None:
                    # The origin is the largest of the ring's events (`ring_origin.ring_origin`), and the steady dust ring of the small moons near the planet (GEN.1233, origin
                    # `moon_erosion`) is one of them: a fallen moon's ring spreads thin with the time since, and a faint dust ring can outweigh what is left of it (one planet of this sample,
                    # a dusty ring of 1e8 kg after Train 21 changed what the seeds draw). So the origin is one of the three, and a `moon_erosion` ring is a dusty one.
                    assert planet.ring_origin in rings.ORIGINS and planet.ring_age_yr > 0.0
                    assert planet.ring_origin != "moon_erosion" or planet.ring_kind == "dusty"
                    assert planet.ring_mass_kg <= rings.MAX_MASS_OVER_PLANET * planet.mass * (1.0 + 1e-9)
            assert planet.rotation_period_hours > 0.0
            for moon in planet.moons:
                seen_moons += 1
                limit_km = rings.roche_limit_km(planet.radius, moon_history.density_kg_m3(planet), moon_history.density_kg_m3(moon))
                assert moon.distance * pc.AU_TO_M / 1000.0 > limit_km
                assert moon.period == pytest.approx(planets.calculate_orbital_period_years(moon.distance, planet.mass), rel=1e-9)
                assert moon.rotation_period_hours > 0.0
    assert seen_moons > 5 and seen_fallen > 0


def test_the_history_runs_again_to_the_same_state_without_new_draws():
    """It starts from the moons the planet was born with and its birth spin, so applying it twice gives one result, and only the spin axis' precession is drawn."""
    planet = _planet_where(lambda p: getattr(p, "_born_moons", None), "a planet with the moons it was born with")
    # the star's age is adjusted for its planets after they are made, so the first run is of the age it had; two runs on the final age agree
    planets.generate_spin(planet, planet.star.mass, birth=planet._birth_spin)
    before = ([m.distance for m in planet.moons], planet.rotation_period_hours, getattr(planet, "fallen_moons", ()))
    draw.set_run_seed(7)
    planets.generate_spin(planet, planet.star.mass, birth=planet._birth_spin)
    after = ([m.distance for m in planet.moons], planet.rotation_period_hours, getattr(planet, "fallen_moons", ()))
    assert after[0] == pytest.approx(before[0], rel=1e-9) and after[1] == pytest.approx(before[1], rel=1e-9) and after[2] == before[2]


def test_a_moon_moves_with_its_planet_tide_and_keeps_its_orbit_plane():
    """An outward-drifting moon is farther than it was born, with its inclination and eccentricity unchanged and its position at its new distance."""
    planet = _planet_where(lambda p: any(abs(m.distance / m._birth_distance_au - 1.0) >= 0.05 for m in p.moons), "a moon moved by 5 percent")
    moon = next(m for m in planet.moons if abs(m.distance / m._birth_distance_au - 1.0) >= 0.05)
    radius = math.sqrt(moon.position_x ** 2 + moon.position_y ** 2 + moon.position_z ** 2)
    assert moon.distance * (1.0 - moon.eccentricity) * (1.0 - 1e-6) <= radius <= moon.distance * (1.0 + moon.eccentricity) * (1.0 + 1e-6)


def test_a_history_no_integrator_can_finish_leaves_the_moons_as_born(monkeypatch):
    """TM8: when `moon_orbits.evolve` gives up, the moons keep their birth distances, nothing falls and the planet's spin is the one it would have without them."""
    from planetgen.physics import moon_orbits

    def fail(*args, **kwargs):
        raise moon_orbits.IntegrationFailed("stiff")

    planet = _planet_where(lambda p: getattr(p, "_born_moons", None), "a planet with the moons it was born with")
    monkeypatch.setattr(moon_orbits, "evolve", fail)
    planets.generate_spin(planet, planet.star.mass, birth=planet._birth_spin)
    assert planet.moons == planet._born_moons
    assert all(m.distance == pytest.approx(m._birth_distance_au) for m in planet.moons)
    assert planet.fallen_moons == () and planet.rotation_period_hours > 0.0
