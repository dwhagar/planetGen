"""
System builder internals (TODO TEST.32).

Direct tests of the `StarSystem` placement helpers that the validate
module (`planetgen/generation/validation.py`, ADM.5) wraps or sits beside:
`generate_slot_object`, `calculate_distance_for_class`, `build_planet`,
`_forced_habitable_distance`, `_trim_to_orbit_ceiling`,
`_reconcile_moved_planet`, `_clear_circumbinary_floor`, and the
`from_dict` newer-schema error. Systems are built planet-less so each
helper is exercised on hand-placed bodies.
"""
import pytest

from planetgen.generation import system as systemData
from planetgen import tuning
from planetgen.physics import constants, daylight, planet_classes, planets as planetPhysics
from planetgen.physics import unit_registry as units
from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.system import StarSystem

from tests.fuzz_support import deterministic_entropy

DRAWS = 40
GIANT_EARTH, EARTH = 300.0, 1.0   # masses (Earth masses) that give a gas giant and an Earth-sized rocky planet


@pytest.fixture(autouse=True)
def _seeded():
    with deterministic_entropy(3201):
        yield


def make_system(star_type="G2V", binary=False, **overrides):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.PLANETS = False
    cfg.MOONS = False
    cfg.BINARY_SYSTEM = binary
    if binary:
        cfg.WIDE_BINARY = False
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return StarSystem(cfg)


@pytest.fixture
def system():
    return make_system()


@pytest.fixture
def close_binary():
    system = make_system(binary=True)
    assert system.binary_type == "close"
    return system


def planet_at(system, distance_au, mass_earth, moon_count=0):
    """A planet of the given mass; its class is derived from what it comes out as (GEN.646)."""
    return Planet(system.system_config, system.star, system.star.habitable_zone, distance_au,
                  mass=mass_earth * constants.EARTH_MASS_TO_KG, moon_count=moon_count)


def belt(system, lower, upper):
    return AsteroidBelt(system.system_config, system.star, lower, lower, upper)


def zone_of(system, distance_au):
    inner, outer = system.star.habitable_zone
    return 'h' if distance_au < inner else 'c' if distance_au > outer else 'e'


# --- generate_slot_object -------------------------------------------------

def test_slot_asteroid_belt_spans_from_the_estimated_distance(system):
    for _ in range(DRAWS):
        obj = system.generate_slot_object(system.star, {"type": "asteroid_belt"}, 3.0)
        assert isinstance(obj, AsteroidBelt) and obj.body_type == 'a'
        assert obj.distance == obj.lower_limit == 3.0
        assert (3.0 * tuning.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN <= obj.upper_limit
                <= 3.0 * tuning.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX)


def test_slot_planet_honours_class_and_moon_count(system):
    cold = system.star.habitable_zone[1] * 4
    obj = system.generate_slot_object(system.star, {"type": "planet", "planet_class": "I", "moons": 2}, cold)
    assert isinstance(obj, Planet)
    assert obj.planet_class == "I" and obj.distance == cold
    assert len(obj.moons) <= 2


def test_slot_type_defaults_to_planet(system):
    obj = system.generate_slot_object(system.star, {}, system.star.habitable_zone[1] * 2)
    assert isinstance(obj, Planet) and obj.body_type in ('t', 'g')


def test_a_slot_for_a_class_that_can_carry_life_is_placed_in_the_ecosphere(system):
    """A class that can carry life asked for in a cold-zone slot is placed in the habitable zone instead of failing."""
    letter = sorted(tuning.HABITABLE_PLANET_CLASSES)[0]
    distance = system.calculate_distance_for_class(letter, system.star.habitable_zone[1] * 5)
    assert zone_of(system, distance) == 'e'


@pytest.mark.parametrize("slot_type", ["comet", "", "Planet", None])
def test_slot_with_an_unknown_type_raises(system, slot_type):
    with pytest.raises(ValueError, match=f"Invalid slot type '{slot_type}'; expected 'planet' or 'asteroid_belt'."):
        system.generate_slot_object(system.star, {"type": slot_type}, 1.0)


# --- calculate_distance_for_class -----------------------------------------

@pytest.mark.parametrize("planet_class", [None, "Z", "not-a-class", "J", "A", "T"])
def test_distance_unchanged_for_no_class_or_a_class_no_zone_belongs_to(system, planet_class):
    """A class is derived from a finished body, so only a class that can carry life has a place (GEN.688)."""
    assert system.calculate_distance_for_class(planet_class, 7.25) == 7.25


def test_distance_unchanged_when_the_life_class_already_fits(system):
    inner, outer = system.star.habitable_zone
    distance = (inner + outer) / 2
    for letter in tuning.HABITABLE_PLANET_CLASSES:
        assert system.calculate_distance_for_class(letter, distance) == distance


def test_a_life_class_is_moved_inside_the_zone_with_margin(system):
    inner, outer = system.star.habitable_zone
    margin = min(tuning.MIN_ASTEROID_BELT_SEPARATION, (outer - inner) / 4)
    letter = sorted(tuning.HABITABLE_PLANET_CLASSES)[0]
    for estimated in (inner * 0.1, outer * 10):
        for _ in range(DRAWS):
            distance = system.calculate_distance_for_class(letter, estimated)
            assert inner + margin <= distance <= outer - margin


def test_ecosphere_draw_avoids_a_belt_in_the_zone(system):
    inner, outer = system.star.habitable_zone
    mid = (inner + outer) / 2
    placed = [belt(system, mid - 0.1, mid + 0.1)]
    for _ in range(DRAWS):
        distance = system.calculate_distance_for_class(sorted(tuning.HABITABLE_PLANET_CLASSES)[0], outer * 5, planets=placed)
        assert inner <= distance <= outer
        assert not (mid - 0.1 < distance < mid + 0.1)


def test_distance_returned_is_a_float(system):
    assert isinstance(system.calculate_distance_for_class(sorted(tuning.HABITABLE_PLANET_CLASSES)[0], 100.0), float)


# --- _forced_habitable_distance -------------------------------------------

def test_forced_habitable_distance_is_inside_the_zone(system):
    inner, outer = system.star.habitable_zone
    for _ in range(DRAWS):
        distance = system._forced_habitable_distance()
        assert inner < distance < outer
        assert zone_of(system, distance) == 'e'


def test_forced_habitable_distance_avoids_a_belt(system):
    inner, outer = system.star.habitable_zone
    placed = [belt(system, inner, inner + (outer - inner) * 0.6)]
    for _ in range(DRAWS):
        distance = system._forced_habitable_distance(planets=placed)
        assert placed[0].upper_limit <= distance < outer


def test_forced_habitable_distance_respects_a_floor_inside_the_zone(close_binary, monkeypatch):
    inner, outer = close_binary.star.habitable_zone
    floor = inner + (outer - inner) * 0.5
    monkeypatch.setattr(close_binary, "_orbit_floor_au", lambda star: floor)
    for _ in range(DRAWS):
        assert floor <= close_binary._forced_habitable_distance() <= outer


def test_forced_habitable_distance_ignores_a_floor_past_the_zone(close_binary, monkeypatch):
    """A floor beyond the whole zone can't be honoured; the draw falls back
    to the full zone rather than an empty range."""
    inner, outer = close_binary.star.habitable_zone
    monkeypatch.setattr(close_binary, "_orbit_floor_au", lambda star: outer * 2)
    for _ in range(DRAWS):
        assert inner <= close_binary._forced_habitable_distance() <= outer


# --- _trim_to_orbit_ceiling -----------------------------------------------

def test_trim_on_an_empty_list_is_a_no_op(system):
    planets = []
    assert system._trim_to_orbit_ceiling(planets, 1.0) is None
    assert planets == []


def test_trim_removes_trailing_bodies_past_the_ceiling(system):
    planets = [planet_at(system, 6.0, GIANT_EARTH), planet_at(system, 30.0, GIANT_EARTH), planet_at(system, 60.0, GIANT_EARTH)]
    keep = planets[0]
    system._trim_to_orbit_ceiling(planets, 10.0)
    assert planets == [keep]


def test_trim_counts_a_planets_hill_sphere_not_just_its_center(system):
    planet = planet_at(system, 9.99, GIANT_EARTH)
    assert planet.distance < 10.0 < planet.distance + planet.min_orbit_distance
    planets = [planet]
    system._trim_to_orbit_ceiling(planets, 10.0)
    assert planets == []


def test_trim_uses_a_belts_upper_limit(system):
    inside, straddling = belt(system, 2.0, 4.0), belt(system, 8.0, 12.0)
    planets = [inside, straddling]
    system._trim_to_orbit_ceiling(planets, 10.0)
    assert planets == [inside]


def test_trim_only_works_from_the_outside_in(system):
    """An interior body past the ceiling stays when the outermost one is
    inside it: trimming stops at the first body that fits."""
    far, near = planet_at(system, 50.0, GIANT_EARTH), planet_at(system, 3.0, GIANT_EARTH)
    planets = [far, near]
    system._trim_to_orbit_ceiling(planets, 10.0)
    assert planets == [far, near]


def test_trim_keeps_everything_under_a_generous_ceiling(system):
    planets = [planet_at(system, 6.0, GIANT_EARTH), belt(system, 20.0, 25.0)]
    before = list(planets)
    system._trim_to_orbit_ceiling(planets, 1e6)
    assert planets == before


# --- _reconcile_moved_planet ----------------------------------------------

def test_reconcile_keeps_a_giant_in_the_cold_zone_after_a_move(system):
    planet = planet_at(system, system.star.habitable_zone[1] * 3, GIANT_EARTH)
    planet.distance *= 2
    system._reconcile_moved_planet(planet)
    assert planet.body_type == "g" and planet.zone == 'c'
    assert planet.planet_class == planet_classes.classify(planet)


def test_reconcile_derives_the_class_again_for_a_planet_pushed_out_of_its_zone(system):
    inner, outer = system.star.habitable_zone
    planet = planet_at(system, (inner + outer) / 2, EARTH)
    mass, radius = planet.mass, planet.radius
    planet.distance = outer * 4
    system._reconcile_moved_planet(planet)
    assert planet.zone == 'c'
    assert planet.distance == outer * 4  # the move itself is kept
    assert (planet.mass, planet.radius) == (mass, radius)  # where it sits does not change what it is made of
    assert planet.planet_class == planet_classes.classify(planet)


def test_reconcile_refreshes_the_period_for_the_new_distance(system):
    inner, outer = system.star.habitable_zone
    planet = planet_at(system, (inner + outer) / 2, EARTH)
    old_period = planet.period
    planet.distance = outer * 4
    system._reconcile_moved_planet(planet)
    assert planet.period > old_period


def test_reconcile_carries_moons_into_the_parents_new_zone(system):
    _, outer = system.star.habitable_zone
    planet = planet_at(system, system.star.habitable_zone[0] * 0.5, GIANT_EARTH, moon_count=3)
    planet.distance = outer * 4
    system._reconcile_moved_planet(planet)
    assert planet.moons
    for moon in planet.moons:
        assert moon.zone == planet.zone == 'c'
        assert moon.planet_class == planet_classes.classify(moon)


def test_reconcile_reruns_spin_dose_and_score_for_each_moon_of_a_moved_planet(system, monkeypatch):
    """GEN.704: a moon keeps its class when its parent only moves within a zone, but its dose, UV and score were drawn
    at the old distance, so `generate_spin` (which ends with them) runs again for each moon."""
    planet = planet_at(system, system.star.habitable_zone[1] * 3, GIANT_EARTH, moon_count=3)
    assert planet.moons
    reran = []
    real = planetPhysics.generate_spin
    monkeypatch.setattr(planetPhysics, "generate_spin", lambda body, mass, **kwargs: (reran.append(body), real(body, mass, **kwargs))[1])
    planet.distance *= 2
    system._reconcile_moved_planet(planet)
    assert all(any(body is moon for body in reran) for moon in planet.moons)


def _assert_solar_day_follows_the_motion(planet):
    """GEN.904: the stored solar day is the one the planet's rotation and its orbit round the star give now."""
    expected = daylight.solar_day_hours(planet.rotation_period_hours, units.converter("yr", "h")(planet.period),
                                        prograde=planet.axial_tilt_deg <= 90.0)
    assert planet.solar_day_hours == (pytest.approx(expected) if expected is not None else None)


@pytest.mark.parametrize("mass_earth, move", [(GIANT_EARTH, 2.0), (EARTH, 4.0)], ids=["giant", "rocky"])
def test_reconcile_leaves_the_daylight_of_a_moved_planet_and_its_moons_in_step_with_its_spin(system, mass_earth, move):
    """GEN.904: `generate_spin` runs again for a moved planet (and a reclassified one is regenerated), so its solar day and
    daylight follow the new rotation and orbit, on the planet and on each moon."""
    inner, outer = system.star.habitable_zone
    for _ in range(12):
        giant = mass_earth == GIANT_EARTH
        distance = outer * 3.0 if giant else (inner + outer) / 2
        planet = planet_at(system, distance, mass_earth, moon_count=2)
        planet.distance = distance * move if giant else outer * move
        system._reconcile_moved_planet(planet)
        if planet.rotation_period_hours:
            assert planet.day_share_mean is not None
            assert planet.polar_circle_deg is not None
            _assert_solar_day_follows_the_motion(planet)
        for moon in planet.moons:
            if moon.rotation_period_hours:
                assert moon.day_share_mean is not None
                assert moon.parent_shadow_share is not None


# --- _clear_circumbinary_floor --------------------------------------------

def test_clear_floor_on_an_empty_list_is_a_no_op(close_binary):
    planets = []
    assert close_binary._clear_circumbinary_floor(planets) is None
    assert planets == []


def test_clear_floor_leaves_a_single_star_alone(system):
    assert system._orbit_floor_au(system.star) == 0.0
    planet = planet_at(system, 0.05, GIANT_EARTH)
    system._clear_circumbinary_floor([planet])
    assert planet.distance == 0.05


def test_clear_floor_pushes_an_inner_planet_out_to_the_limit(close_binary):
    floor = close_binary._orbit_floor_au(close_binary.star)
    assert floor > 0
    planet = planet_at(close_binary, floor / 2, GIANT_EARTH)
    outer_planet = planet_at(close_binary, floor * 50, GIANT_EARTH)
    close_binary._clear_circumbinary_floor([planet, outer_planet])
    assert planet.distance == pytest.approx(floor)
    assert outer_planet.distance == floor * 50  # only the innermost body moves


def test_clear_floor_shifts_a_belt_keeping_its_width(close_binary):
    floor = close_binary._orbit_floor_au(close_binary.star)
    inner_belt = belt(close_binary, floor / 4, floor / 2)
    width = inner_belt.upper_limit - inner_belt.lower_limit
    close_binary._clear_circumbinary_floor([inner_belt])
    assert inner_belt.lower_limit == pytest.approx(floor)
    assert inner_belt.upper_limit - inner_belt.lower_limit == pytest.approx(width)
    assert inner_belt.distance == pytest.approx(floor)


def test_clear_floor_leaves_a_body_already_outside_it(close_binary):
    floor = close_binary._orbit_floor_au(close_binary.star)
    planet = planet_at(close_binary, floor * 3, GIANT_EARTH)
    outer_belt = belt(close_binary, floor * 2, floor * 3)
    close_binary._clear_circumbinary_floor([planet])
    close_binary._clear_circumbinary_floor([outer_belt])
    assert planet.distance == floor * 3
    assert (outer_belt.lower_limit, outer_belt.upper_limit) == (floor * 2, floor * 3)


def test_clear_floor_reconciles_a_moved_planet_to_its_new_zone(close_binary, monkeypatch):
    inner, outer = close_binary.star.habitable_zone
    monkeypatch.setattr(close_binary, "_orbit_floor_au", lambda star: outer * 3)
    planet = planet_at(close_binary, (inner + outer) / 2, EARTH)
    close_binary._clear_circumbinary_floor([planet])
    assert planet.distance == outer * 3
    assert planet.zone == 'c' and planet.planet_class == planet_classes.classify(planet)


# --- from_dict ------------------------------------------------------------

def test_from_dict_rejects_a_newer_schema(system):
    data = system.to_dict()
    data["schema_version"] = systemData.SERIALIZATION_SCHEMA_VERSION + 1
    with pytest.raises(ValueError, match=(
            rf"StarSystem.from_dict: schema_version {systemData.SERIALIZATION_SCHEMA_VERSION + 1} is newer "
            rf"than this code understands \(max {systemData.SERIALIZATION_SCHEMA_VERSION}\)\.")):
        StarSystem.from_dict(data)


def test_from_dict_accepts_the_current_and_a_missing_schema(system):
    data = system.to_dict()
    assert data["schema_version"] == systemData.SERIALIZATION_SCHEMA_VERSION
    assert StarSystem.from_dict(data).star.name == system.star.name
    del data["schema_version"]
    assert StarSystem.from_dict(data).star.name == system.star.name


def test_from_dict_rejects_a_newer_schema_before_reading_anything_else():
    """The version check runs first, so an otherwise empty dict still gets
    the clear error rather than a KeyError."""
    with pytest.raises(ValueError, match="is newer than this code understands"):
        StarSystem.from_dict({"schema_version": 10 ** 6})


# ---------------------------------------------------------------------------
# Life data: no chemistry, no timeline
# ---------------------------------------------------------------------------

def test_a_habitable_world_with_no_viable_chemistry_gets_no_timeline(monkeypatch):
    """A habitable-zone planet whose class and star share no life chemical
    is lifeless, so it gets no evolutionary timeline; it used to get one,
    and the population pass then gave it a species no page showed."""
    from planetgen.generation import life as planetLife

    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.HABITABLE_WORLD = True
    cfg.BINARY_SYSTEM = False
    # Not every habitable-zone world draws a life chemistry; take the first that did.
    planet = None
    for _ in range(50):
        system = StarSystem(system_config=cfg)
        planet = next((p for p in system.planets if getattr(p, "zone", None) == "e" and not p.is_moon
                       and p.has_timeline), None)
        if planet is not None:
            break
    assert planet is not None

    monkeypatch.setattr(planetLife.life_classes, "score_classes", lambda *args, **kwargs: {})
    planetLife.apply_life_data(planet)
    assert planet.life_class is None
    assert not planet.has_timeline and planet.evolutionary_pace is None
