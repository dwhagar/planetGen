"""
System builder internals (TODO TEST.32).

Direct tests of the `StarSystem` placement helpers that the validate
module (`stellarObjects/validation.py`, ADM.5) wraps or sits beside:
`generate_slot_object`, `calculate_distance_for_class`,
`_forced_habitable_distance`, `_trim_to_orbit_ceiling`,
`_reconcile_moved_planet`, `_clear_circumbinary_floor`, and the
`from_dict` newer-schema error. Systems are built planet-less so each
helper is exercised on hand-placed bodies.
"""
import pytest

from stellarObjects import program_constants, systemData
from stellarObjects.asteroidData import AsteroidBelt
from stellarObjects.config import SystemConfig
from stellarObjects.planetData import Planet
from stellarObjects.systemData import StarSystem

from tests.fuzz_support import deterministic_entropy

DRAWS = 40


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


def planet_at(system, distance_au, planet_class, moon_count=0):
    return Planet(system.system_config, system.star, system.star.habitable_zone, distance_au,
                  planet_class=planet_class, moon_count=moon_count)


def belt(system, lower, upper):
    return AsteroidBelt(system.system_config, lower, lower, upper)


def zone_of(system, distance_au):
    inner, outer = system.star.habitable_zone
    return 'h' if distance_au < inner else 'c' if distance_au > outer else 'e'


# --- generate_slot_object -------------------------------------------------

def test_slot_asteroid_belt_spans_from_the_estimated_distance(system):
    for _ in range(DRAWS):
        obj = system.generate_slot_object({"type": "asteroid_belt"}, 3.0)
        assert isinstance(obj, AsteroidBelt) and obj.body_type == 'a'
        assert obj.distance == obj.lower_limit == 3.0
        assert (3.0 * program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MIN <= obj.upper_limit
                <= 3.0 * program_constants.ASTEROID_BELT_MAX_DISTANCE_FACTOR_MAX)


def test_slot_planet_honours_class_and_moon_count(system):
    cold = system.star.habitable_zone[1] * 4
    obj = system.generate_slot_object({"type": "planet", "planet_class": "J", "moons": 2}, cold)
    assert isinstance(obj, Planet)
    assert obj.planet_class == "J" and obj.distance == cold
    assert len(obj.moons) <= 2


def test_slot_type_defaults_to_planet(system):
    obj = system.generate_slot_object({}, system.star.habitable_zone[1] * 2)
    assert isinstance(obj, Planet) and obj.body_type in ('t', 'g')


def test_slot_planet_class_is_moved_into_a_zone_that_supports_it(system):
    """An ecosphere-only class asked for in a cold-zone slot is placed in
    the habitable zone instead of failing."""
    obj = system.generate_slot_object({"type": "planet", "planet_class": "M"}, system.star.habitable_zone[1] * 5)
    assert obj.planet_class == "M" and obj.zone == 'e'


@pytest.mark.parametrize("slot_type", ["comet", "", "Planet", None])
def test_slot_with_an_unknown_type_raises(system, slot_type):
    with pytest.raises(ValueError, match=f"Invalid slot type '{slot_type}'; expected 'planet' or 'asteroid_belt'."):
        system.generate_slot_object({"type": slot_type}, 1.0)


# --- calculate_distance_for_class -----------------------------------------

@pytest.mark.parametrize("planet_class", [None, "Z", "not-a-class"])
def test_distance_unchanged_for_no_or_unknown_class(system, planet_class):
    assert system.calculate_distance_for_class(planet_class, 7.25) == 7.25


@pytest.mark.parametrize("planet_class, zone", [("J", 'h'), ("J", 'e'), ("J", 'c'), ("M", 'e'), ("A", 'h'), ("T", 'c')])
def test_distance_unchanged_when_the_class_already_fits(system, planet_class, zone):
    inner, outer = system.star.habitable_zone
    distance = {'h': inner * 0.5, 'e': (inner + outer) / 2, 'c': outer * 2}[zone]
    assert system.calculate_distance_for_class(planet_class, distance) == distance


def test_ecosphere_class_is_moved_inside_the_zone_with_margin(system):
    inner, outer = system.star.habitable_zone
    margin = min(program_constants.MIN_ASTEROID_BELT_SEPARATION, (outer - inner) / 4)
    for estimated in (inner * 0.1, outer * 10):
        for _ in range(DRAWS):
            distance = system.calculate_distance_for_class("M", estimated)
            assert inner + margin <= distance <= outer - margin


def test_hot_only_class_is_moved_into_the_hot_zone(system):
    inner, _ = system.star.habitable_zone
    for _ in range(DRAWS):
        distance = system.calculate_distance_for_class("A", system.star.habitable_zone[1] * 3)
        assert inner * 0.05 <= distance <= inner * 0.95


def test_cold_only_class_is_moved_into_the_cold_zone(system):
    _, outer = system.star.habitable_zone
    for _ in range(DRAWS):
        distance = system.calculate_distance_for_class("T", system.star.habitable_zone[0] * 0.2)
        assert outer * 1.05 <= distance <= outer * 3.0


@pytest.mark.parametrize("planet_class, estimated_factor", [("A", 4.0), ("T", 0.2)])
def test_moved_distance_avoids_an_already_placed_belt(system, planet_class, estimated_factor):
    inner, outer = system.star.habitable_zone
    if planet_class == "A":
        lo, hi = inner * 0.05, inner * 0.95
    else:
        lo, hi = outer * 1.05, outer * 3.0
    # A belt covering the middle 80% of the target range.
    span = hi - lo
    placed = [belt(system, lo + span * 0.1, hi - span * 0.1)]
    for _ in range(DRAWS):
        distance = system.calculate_distance_for_class(planet_class, inner * estimated_factor, planets=placed)
        assert lo <= distance <= hi
        assert not (placed[0].lower_limit < distance < placed[0].upper_limit)


def test_ecosphere_draw_avoids_a_belt_in_the_zone(system):
    inner, outer = system.star.habitable_zone
    mid = (inner + outer) / 2
    placed = [belt(system, mid - 0.1, mid + 0.1)]
    for _ in range(DRAWS):
        distance = system.calculate_distance_for_class("M", outer * 5, planets=placed)
        assert inner <= distance <= outer
        assert not (mid - 0.1 < distance < mid + 0.1)


def test_distance_returned_is_a_float(system):
    assert isinstance(system.calculate_distance_for_class("M", 100.0), float)


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
    planets = [planet_at(system, 6.0, "J"), planet_at(system, 30.0, "J"), planet_at(system, 60.0, "J")]
    keep = planets[0]
    system._trim_to_orbit_ceiling(planets, 10.0)
    assert planets == [keep]


def test_trim_counts_a_planets_hill_sphere_not_just_its_center(system):
    planet = planet_at(system, 9.99, "J")
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
    far, near = planet_at(system, 50.0, "J"), planet_at(system, 3.0, "J")
    planets = [far, near]
    system._trim_to_orbit_ceiling(planets, 10.0)
    assert planets == [far, near]


def test_trim_keeps_everything_under_a_generous_ceiling(system):
    planets = [planet_at(system, 6.0, "J"), belt(system, 20.0, 25.0)]
    before = list(planets)
    system._trim_to_orbit_ceiling(planets, 1e6)
    assert planets == before


# --- _reconcile_moved_planet ----------------------------------------------

def test_reconcile_keeps_a_class_that_still_fits(system):
    planet = planet_at(system, system.star.habitable_zone[1] * 3, "J")
    planet.distance *= 2
    assert system._reconcile_moved_planet(planet) is False
    assert planet.planet_class == "J" and planet.zone == 'c'


def test_reconcile_reclassifies_a_planet_pushed_out_of_its_zone(system):
    inner, outer = system.star.habitable_zone
    planet = planet_at(system, (inner + outer) / 2, "M")
    planet.distance = outer * 4
    assert system._reconcile_moved_planet(planet) is True
    assert planet.zone == 'c'
    assert planet.planet_class != "M"
    assert program_constants.PLANET_CLASSES[planet.planet_class]['c']
    assert planet.distance == outer * 4  # the move itself is kept


def test_reconcile_refreshes_the_period_for_the_new_distance(system):
    inner, outer = system.star.habitable_zone
    planet = planet_at(system, (inner + outer) / 2, "M")
    old_period = planet.period
    planet.distance = outer * 4
    system._reconcile_moved_planet(planet)
    assert planet.period > old_period


def test_reconcile_carries_moons_into_the_parents_new_zone(system):
    _, outer = system.star.habitable_zone
    planet = planet_at(system, system.star.habitable_zone[0] * 0.5, "J", moon_count=3)
    planet.distance = outer * 4
    system._reconcile_moved_planet(planet)
    for moon in planet.moons:
        assert moon.zone == planet.zone == 'c'
        assert program_constants.PLANET_CLASSES[moon.planet_class]['c']


# --- _clear_circumbinary_floor --------------------------------------------

def test_clear_floor_on_an_empty_list_is_a_no_op(close_binary):
    planets = []
    assert close_binary._clear_circumbinary_floor(planets) is None
    assert planets == []


def test_clear_floor_leaves_a_single_star_alone(system):
    assert system._orbit_floor_au(system.star) == 0.0
    planet = planet_at(system, 0.05, "J")
    system._clear_circumbinary_floor([planet])
    assert planet.distance == 0.05


def test_clear_floor_pushes_an_inner_planet_out_to_the_limit(close_binary):
    floor = close_binary._orbit_floor_au(close_binary.star)
    assert floor > 0
    planet = planet_at(close_binary, floor / 2, "J")
    outer_planet = planet_at(close_binary, floor * 50, "J")
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
    planet = planet_at(close_binary, floor * 3, "J")
    outer_belt = belt(close_binary, floor * 2, floor * 3)
    close_binary._clear_circumbinary_floor([planet])
    close_binary._clear_circumbinary_floor([outer_belt])
    assert planet.distance == floor * 3
    assert (outer_belt.lower_limit, outer_belt.upper_limit) == (floor * 2, floor * 3)


def test_clear_floor_reconciles_a_moved_planet_to_its_new_zone(close_binary, monkeypatch):
    inner, outer = close_binary.star.habitable_zone
    monkeypatch.setattr(close_binary, "_orbit_floor_au", lambda star: outer * 3)
    planet = planet_at(close_binary, (inner + outer) / 2, "M")
    close_binary._clear_circumbinary_floor([planet])
    assert planet.distance == outer * 3
    assert planet.zone == 'c' and planet.planet_class != "M"


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
    from stellarObjects import planetLife

    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.HABITABLE_WORLD = True
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    planet = next(p for p in system.planets if getattr(p, "zone", None) == "e" and not p.is_moon)
    assert planet.evolutionary_data

    monkeypatch.setattr(planetLife, "get_viable_life_chemicals", lambda *args, **kwargs: {})
    planet.evolutionary_data = []
    planetLife.apply_life_data(planet)
    assert planet.life_chemical is None
    assert planet.evolutionary_data == []
