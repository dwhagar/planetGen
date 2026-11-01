# tests/test_belt_birth.py

"""
The belts of a system from its finished planets (GEN.1068, GEN.1156): `generation/belt_birth.py` and its wiring in `StarSystem._generate_planets`.
"""

import types

import pytest

from planetgen.generation import belt_birth
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import belt_origin, constants
from planetgen.util import draw

from tests.fuzz_support import deterministic_entropy


def _config(**overrides):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MOONS = False
    cfg.BINARY_SYSTEM = False
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return cfg


def _belts(system):
    return [body for body in system.planets if body.body_type == 'a']


def _planet(distance, mass, body_type="t", formed=None, eccentricity=0.0):
    return types.SimpleNamespace(body_type=body_type, distance=distance, mass=mass, formation_distance_au=formed, eccentricity=eccentricity)


def test_the_chain_skips_belts_and_marks_a_giant_that_moved():
    chain = belt_birth.chain_of([_planet(0.05, 1.9e27, "g", formed=4.0), _planet(1.0, 6e24), types.SimpleNamespace(body_type='a', distance=2.0), _planet(5.2, 1.9e27, "g", formed=5.2)])
    assert [body.distance_au for body in chain] == [0.05, 1.0, 5.2]
    assert [body.giant for body in chain] == [True, False, True]
    assert [body.migrated for body in chain] == [True, False, False]


def test_a_nudged_planet_is_not_a_migrated_giant():
    chain = belt_birth.chain_of([_planet(5.2, 1.9e27, "g", formed=5.0)])
    assert not chain[0].migrated


def test_a_belt_is_inner_inside_a_snow_line_and_a_fifth():
    snow = 3.0
    inner = types.SimpleNamespace(upper_limit=3.5)
    outer = types.SimpleNamespace(upper_limit=40.0)
    assert belt_birth.is_inner(inner, snow) and not belt_birth.is_inner(outer, snow)


def test_the_belts_of_a_system_carry_where_they_came_from():
    with deterministic_entropy(1068):
        for seed in range(1, 9):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=_config())
            belts = _belts(system)
            assert [belt.distance for belt in belts] == sorted(belt.distance for belt in belts)
            for belt in belts:
                assert belt.origin in belt_origin.CAUSES and belt.origin_sentence
                assert belt.primordial_mass_earth >= belt_origin.MIN_BELT_MASS_EARTH
                assert belt.cleared_mass_earth >= 0.0 and belt.scattered_mass_earth >= 0.0
                # what is left after the giants and the collisions is never more than what the giants left
                left_earth = belt.mass_kg / constants.EARTH_MASS_TO_KG
                assert left_earth <= belt.primordial_mass_earth - belt.cleared_mass_earth - belt.scattered_mass_earth + 1e-12
                assert left_earth >= belt_birth.MIN_PRESENT_MASS_EARTH


def test_a_stored_belt_reads_back_with_its_origin():
    with deterministic_entropy(1068):
        for seed in range(1, 25):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=_config())
            for belt in _belts(system):
                again = AsteroidBelt.from_dict(belt.to_dict(), system.system_config)
                assert again.origin == belt.origin and again.origin_sentence == belt.origin_sentence
                assert again.co_causes == list(belt.co_causes) and again.primordial_mass_earth == belt.primordial_mass_earth
                assert again.mass_kg == belt.mass_kg and again.stirred == belt.stirred
                return
    pytest.skip("no belt in the sample")


def test_a_floor_above_every_belt_leaves_none(monkeypatch):
    monkeypatch.setattr(belt_birth, "MIN_PRESENT_MASS_EARTH", 1.0e30)
    with deterministic_entropy(1068):
        for seed in range(1, 6):
            draw.set_run_seed(seed)
            assert not _belts(StarSystem(system_config=_config()))


def test_no_belts_when_the_request_says_none_or_no_planets():
    with deterministic_entropy(1068):
        for seed in range(1, 6):
            draw.set_run_seed(seed)
            assert not _belts(StarSystem(system_config=_config(ASTEROID_BELT=False)))
            assert not StarSystem(system_config=_config(PLANETS=False)).planets


def test_a_slots_list_names_its_belts_and_nothing_else_is_added():
    with deterministic_entropy(1068):
        draw.set_run_seed(3)
        system = StarSystem(system_config=_config(NUM_ORBITS=3, SLOTS=[{"type": "planet"}, {"type": "asteroid_belt"}, {"type": "planet"}]))
        assert [body.body_type == 'a' for body in system.planets][:3] == [False, True, False]
        assert len(_belts(system)) == 1
        assert _belts(system)[0].origin is None


def test_a_count_of_orbits_is_the_whole_count_of_bodies():
    with deterministic_entropy(1068):
        for seed in range(1, 5):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=_config(NUM_ORBITS=2))
            assert system.planet_count + system.belt_count == 2
            assert not _belts(system)
            fewest = StarSystem(system_config=_config(MAX_PLANETS=False))
            assert fewest.planet_count + fewest.belt_count == 1


def test_a_belt_asked_for_is_found_by_the_retries():
    with deterministic_entropy(1068):
        for seed in range(1, 5):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=_config(ASTEROID_BELT=True))
            assert _belts(system) or "an asteroid belt" in system.unmet_requirements


def test_a_belt_past_the_final_ceiling_is_cut_back_not_lost():
    from planetgen.generation import validation
    with deterministic_entropy(1068):
        for seed in range(1, 40):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=_config())
            belt = next((b for b in _belts(system) if b.origin is not None and b.upper_limit > 2.0 * b.lower_limit), None)
            if belt is None:
                continue
            before = (belt.primordial_mass_earth, belt.mass_kg, belt.upper_limit)
            ceiling = (belt.lower_limit * belt.upper_limit) ** 0.5
            bodies = [body for body in system.planets if body.distance <= belt.distance] + [belt]
            assert validation.trim_to_orbit_ceiling(bodies, ceiling) == []
            assert bodies[-1] is belt and belt.upper_limit == pytest.approx(ceiling)
            assert belt.primordial_mass_earth < before[0] and belt.mass_kg <= before[1]
            assert belt.primordial_mass_earth - belt.cleared_mass_earth - belt.scattered_mass_earth >= 0.0
            tight = [belt]
            assert validation.trim_to_orbit_ceiling(tight, belt.lower_limit * 1.01) == [belt] and not tight
            return
    pytest.skip("no wide belt in the sample")


def test_a_belt_asked_for_is_still_placed_where_the_chain_leaves_no_stretch(monkeypatch):
    monkeypatch.setattr(belt_birth, "build", lambda *args, **kwargs: [])
    with deterministic_entropy(1068):
        for seed in range(1, 4):
            draw.set_run_seed(seed)
            system = StarSystem(system_config=_config(ASTEROID_BELT=True))
            belts = _belts(system)
            assert belts or "an asteroid belt" in system.unmet_requirements
            for belt in belts:
                assert belt.origin is None and belt.upper_limit >= belt.lower_limit * 1.1


def test_a_stretch_the_giants_cleared_of_everything_is_no_belt(monkeypatch):
    real = belt_origin.build

    def cleared(*args, **kwargs):
        return [item._replace(mass_earth=0.0) for item in real(*args, **kwargs)]

    monkeypatch.setattr(belt_origin, "build", cleared)
    with deterministic_entropy(1068):
        for seed in range(1, 8):
            draw.set_run_seed(seed)
            assert not _belts(StarSystem(system_config=_config()))
