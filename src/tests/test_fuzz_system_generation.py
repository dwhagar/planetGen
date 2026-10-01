# tests/test_fuzz_system_generation.py

"""
Property-based fuzzing of whole-system generation (`StarSystem` and
everything it drives: `starData`, `planetData`/`planetPhysics` (moons
included), `asteroidData`, `cometData`, `doubleStar` (P-type/close
binaries), `wideBinary` (S-type/wide binaries), `compactRemnant`
(black-hole/neutron-star anchored systems), `evolution` and `planetLife`),
plus the pure physics helpers underneath it called directly with boundary
inputs. See `tests/fuzz_support.py` for the `ci`/`deep` profiles.

What this adds over the seeded `test_bughunt_*.py` files and
`test_systems.py`:

- configurations are drawn by hypothesis from the whole valid option space
  at once (every spectral letter x subclass x Yerkes class including the
  `D` alias, `NUM_ORBITS` from 0 up to far past
  `ABSOLUTE_MAX_SYSTEM_OBJECTS`, every tri-state flag on/off/None together,
  hand-written `SLOTS` lists, hostile system names) instead of a fixed
  list, and shrink to a minimal failing config;
- generation is made fully reproducible: besides the global `random`
  module, generation draws from `secrets` (`utils.reseed_rng`,
  `planetData`/`planetLife`/`planetPhysics` class/flavor picks), so
  `_deterministic_entropy` routes `secrets.randbits`/`randbelow`/`choice`
  through a `random.Random` seeded from the hypothesis-drawn seed -- a
  failing example replays exactly;
- the serialized tree is walked for every numeric leaf (not a hand-picked
  attribute list) and must survive strict JSON (`allow_nan=False`) and a
  JSON -> `from_dict` -> `to_dict` round trip unchanged;
- orbital sanity is checked with *independent* physics (adjacent Hill
  spheres disjoint, belts disjoint from every other body, moons inside
  their parent's Hill sphere and outside its body, circumbinary bodies
  outside the binary's own orbit) rather than by mirroring
  `validate_system`'s own formula.

Real bugs found here are kept as `xfail(strict=True)` tests named after
the invariant they break, each with a minimal reproduction in its reason.
"""

import json
import math
import re

import pytest
from hypothesis import HealthCheck, assume, example, given, note, settings
from hypothesis import strategies as st

from stellarObjects import evolution, keplerMotion as km, physical_constants as pc
from stellarObjects import planetLife, planetPhysics as pp, program_constants as prog
from stellarObjects import cometData, starData, utils
from stellarObjects.compactRemnant import BlackHole, NeutronStar
from stellarObjects.config import SystemConfig
from stellarObjects.doubleStar import BinaryStarProxy
from stellarObjects.starData import Star
from stellarObjects.systemData import StarSystem

from tests.fuzz_support import deterministic_entropy as _deterministic_entropy, hostile_text, scaled

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------


def _generate(cfg, seed, **kwargs):
    with _deterministic_entropy(seed):
        return StarSystem(system_config=cfg, **kwargs)


# ---------------------------------------------------------------------------
# Configuration space
# ---------------------------------------------------------------------------

SPECTRAL = "OBAFGKM"
YERKES = ["0", "IA+", "IA", "IAB", "IB", "II", "III", "IV", "V", "VI", "VII", "D"]
ALL_STAR_TYPES = [f"{s}{n}{y}" for s in SPECTRAL for n in range(10) for y in YERKES]
TRISTATE = [
    "HABITABLE_WORLD", "ASTEROID_BELT", "COMETS", "LARGE_STAR", "MOONS",
    "MAX_PLANETS", "INTELLIGENT_LIFE", "BINARY_SYSTEM", "WIDE_BINARY", "PLANETS",
]
PLANET_CLASSES = sorted(prog.PLANET_CLASSES)
SEEDS = st.integers(min_value=0, max_value=2**64 - 1)
NUM_ORBITS = st.one_of(
    st.none(),
    st.sampled_from([0, 1, 2, prog.ABSOLUTE_MAX_SYSTEM_OBJECTS, prog.ABSOLUTE_MAX_SYSTEM_OBJECTS + 1, 10**6]),
    st.integers(min_value=0, max_value=40),
)
SLOT_SPEC = st.one_of(
    st.none(),
    st.just({"type": "asteroid_belt"}),
    st.fixed_dictionaries(
        {"type": st.just("planet")},
        optional={
            "planet_class": st.sampled_from(PLANET_CLASSES),
            "moons": st.integers(min_value=-3, max_value=10**6),
        },
    ),
)


def _normalize(cfg):
    """`generate.build_system_config`'s own cross-option normalization,
    applied to a directly-built config, plus `validate_system_args`'s
    rejections as `assume()`s -- so every config drawn here is one the CLI
    could actually produce."""
    if cfg.INTELLIGENT_LIFE is not None:
        assume(cfg.HABITABLE_WORLD is not False)
        cfg.HABITABLE_WORLD = True
    if cfg.PLANETS is False:
        assume(not (cfg.MOONS or cfg.MAX_PLANETS or cfg.HABITABLE_WORLD))
        assume(cfg.NUM_ORBITS is None and not cfg.SLOTS)
    if cfg.HABITABLE_WORLD is True and cfg.ASTEROID_BELT is True:
        assume(cfg.LARGE_STAR is not False)
        cfg.LARGE_STAR = True
    if cfg.HABITABLE_WORLD is False:
        # An explicit habitable-class slot under -habitable_world is a
        # contradiction the generator rejects with a clean ValueError
        # (covered by test_contradictory_slot_and_flag_is_a_clean_value_error).
        assume(not any(spec and spec.get("planet_class") in prog.HABITABLE_PLANET_CLASSES
                       for spec in (cfg.SLOTS or [])))
    return cfg


@st.composite
def system_configs(draw):
    cfg = SystemConfig()
    for attr in TRISTATE:
        setattr(cfg, attr, draw(st.sampled_from([None, True, False])))
    cfg.STAR_TYPE = draw(st.one_of(st.none(), st.sampled_from(ALL_STAR_TYPES)))
    cfg.NUM_ORBITS = draw(NUM_ORBITS)
    cfg.AGE = draw(st.sampled_from([None, "young", "old"]))
    cfg.MARKDOWN = draw(st.booleans())
    cfg.NAME = draw(st.one_of(st.none(), hostile_text.filter(lambda s: s.strip() != "")))
    if draw(st.integers(0, 4)) == 0:
        cfg.SLOTS = draw(st.lists(SLOT_SPEC, max_size=8))
    return _normalize(cfg)


# ---------------------------------------------------------------------------
# Invariant checks
# ---------------------------------------------------------------------------

CMB_K = 2.7  # nothing in the universe is colder than the CMB (2.725 K)

# Serialized keys that must be strictly positive wherever they hold a number.
POSITIVE_KEYS = {
    "mass", "radius", "period", "temperature", "luminosity", "volume", "hill_radius",
    "min_orbit_distance", "system_perimeter", "age", "lifespan", "gravity", "scale_height",
    "atm_density", "atm_molar_density", "nucleus_diameter_km", "perihelion_distance_au",
    "distance_au", "orbital_period_years", "_binary_separation_au", "_effective_mass",
    "separation_au", "binary_mutual_orbital_period_years", "galactic_orbital_period_gy",
    "galactic_orbital_speed_kms", "primary_mass_solar", "rotation_period_hours",
    "orbital_speed_kms", "min_update_interval_years", "upper_limit", "lower_limit", "a_crit_au",
}
NON_NEGATIVE_KEYS = {"atmospheric_pressure", "heliosphere_radius", "eccentricity", "flavor_text_count"}
TEMPERATURE_KEYS = {"temperature", "surface_temperature"}


def _walk(node, path="$"):
    if isinstance(node, dict):
        for key, value in node.items():
            yield from _walk(value, f"{path}.{key}")
    elif isinstance(node, (list, tuple)):
        for index, value in enumerate(node):
            yield from _walk(value, f"{path}[{index}]")
    else:
        yield path, node


def _is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def assert_serialized_numbers_sane(data, *, zero_ok_paths=()):
    """Every numeric leaf finite; sign/range rules by key name. A path in
    `zero_ok_paths` may be exactly 0.0 (a black hole's own, deliberately
    unmodeled, temperature/luminosity) but still never negative."""
    problems = []
    for path, value in _walk(data):
        if not _is_number(value):
            continue
        key = path.rsplit(".", 1)[-1].split("[", 1)[0]
        if not math.isfinite(value):
            problems.append(f"{path} = {value!r} (not finite)")
            continue
        if value == 0 and path in zero_ok_paths:
            continue
        if key in POSITIVE_KEYS and not value > 0:
            problems.append(f"{path} = {value!r} (must be > 0)")
        elif key in NON_NEGATIVE_KEYS and value < 0:
            problems.append(f"{path} = {value!r} (must be >= 0)")
        if key in TEMPERATURE_KEYS and value < CMB_K:
            problems.append(f"{path} = {value!r} K (colder than the CMB)")
    assert not problems, "\n".join(problems[:25])
    # Strict JSON: no NaN/Infinity tokens a standards-conforming parser would reject.
    json.dumps(data, allow_nan=False)


def _edge_out(obj):
    return obj.upper_limit if obj.body_type == "a" else obj.distance


def _hill_au(planet):
    return planet.hill_radius / pc.AU_TO_KM


def assert_orbits_sane(planets):
    """Independent physical checks on one star's ordered body list."""
    distances = [obj.distance for obj in planets]
    assert distances == sorted(distances), f"bodies not in orbital order: {distances}"
    assert len(set(distances)) == len(distances), f"two bodies share an orbit: {distances}"
    for obj in planets:
        if obj.body_type == "a":
            assert 0 < obj.lower_limit <= obj.distance <= obj.upper_limit, (
                f"belt limits out of order: {obj.lower_limit} / {obj.distance} / {obj.upper_limit}")
    for prev, cur in zip(planets, planets[1:]):
        if prev.body_type == "a":
            # The floor every body after a belt keeps; a planet's own Hill
            # sphere clearance is checked by
            # test_planet_hill_sphere_clears_the_belt_inside_it.
            inner_reach = prev.upper_limit + prog.MIN_ASTEROID_BELT_SEPARATION
            outer_reach = cur.distance if cur.body_type != "a" else cur.lower_limit
            assert inner_reach <= outer_reach * (1 + 1e-12), (
                f"body at {cur.distance} AU closer than MIN_ASTEROID_BELT_SEPARATION to the belt "
                f"ending at {prev.upper_limit} AU")
            continue
        inner_reach = prev.distance + _hill_au(prev)
        if cur.body_type == "a":
            outer_reach = cur.lower_limit
        else:
            outer_reach = cur.distance - _hill_au(cur)
        assert inner_reach <= outer_reach * (1 + 1e-12), (
            f"{prev.body_type}@{prev.distance} AU and {cur.body_type}@{cur.distance} AU overlap: "
            f"inner body reaches {inner_reach} AU, outer body starts at {outer_reach} AU")
    belts = [obj for obj in planets if obj.body_type == "a"]
    for belt in belts:
        for other in planets:
            if other is belt:
                continue
            if other.body_type == "a":
                assert not (belt.lower_limit < other.upper_limit and other.lower_limit < belt.upper_limit), (
                    f"belts [{belt.lower_limit}, {belt.upper_limit}] and "
                    f"[{other.lower_limit}, {other.upper_limit}] overlap")
            else:
                assert not (belt.lower_limit < other.distance < belt.upper_limit), (
                    f"planet at {other.distance} AU inside belt [{belt.lower_limit}, {belt.upper_limit}]")


def _all_planet_lists(system):
    return [system.planets, system.secondary_planets]


def assert_flags_honored(system, cfg):
    bodies = system.planets + system.secondary_planets
    real = [b for b in bodies if b.body_type != "a"]
    if cfg.PLANETS is False:
        assert not bodies, "PLANETS=False but bodies were generated"
    explicit_moons = any(spec and spec.get("moons") for spec in (cfg.SLOTS or []))
    if cfg.MOONS is False and not explicit_moons:
        assert all(not p.moons for p in real), "MOONS=False but a planet has moons"
    if cfg.COMETS is False:
        assert not system.comets and not system.secondary_comets
    if cfg.COMETS is True:
        assert len(system.comets) >= prog.SYSTEM_COMET_COUNT_RANGE[0]
    explicit_belt = any(spec and spec.get("type") == "asteroid_belt" for spec in (cfg.SLOTS or []))
    if cfg.ASTEROID_BELT is False and not explicit_belt:
        assert not any(b.body_type == "a" for b in bodies), "ASTEROID_BELT=False but a belt exists"
    if cfg.HABITABLE_WORLD is False:
        assert system.hab_count == 0, "HABITABLE_WORLD=False but a habitable world exists"
    if cfg.BINARY_SYSTEM is False:
        assert system.binary_type is None and len(system.stars) == 1
    if cfg.BINARY_SYSTEM is True:
        expected = {True: "wide", False: "close"}.get(cfg.WIDE_BINARY)
        assert system.binary_type in ({"wide", "close"} if expected is None else {expected})
        assert len(system.stars) == 2
    if cfg.NUM_ORBITS is not None and cfg.PLANETS is not False:
        assert len(system.planets) <= max(cfg.NUM_ORBITS, 0) + 2
    if cfg.NAME:
        # A wide pair's stars add their own word after the system name.
        assert system.name == cfg.NAME


def assert_binary_geometry(system):
    if system.binary_type == "wide":
        wb = system.wide_binary
        assert 0 < system.primary_star.a_crit_au < wb.separation_au
        assert 0 < system.secondary_star.a_crit_au < wb.separation_au
        for star, planets in ((system.primary_star, system.planets),
                              (system.secondary_star, system.secondary_planets)):
            for obj in planets:
                assert _edge_out(obj) <= star.a_crit_au * (1 + 1e-9), (
                    f"body at {_edge_out(obj)} AU past its star's stability limit {star.a_crit_au} AU")
    if system.binary_type == "close":
        assert isinstance(system.star, BinaryStarProxy)
        assert system.star.mass == pytest.approx(system.primary_star.mass + system.secondary_star.mass)
        assert system.star._secondary.mass <= system.star._primary.mass


def assert_life_caps(system, cfg):
    for obj in system.planets + system.secondary_planets:
        if obj.body_type == "a" or not obj.evolutionary_data:
            continue
        stage = evolution.life_stage_from_paragraphs(obj.evolutionary_data)
        if stage is None:
            continue
        cap = prog.PLANET_CLASS_MAX_LIFE_STAGE.get(obj.planet_class)
        if cap is not None:
            assert evolution.MILESTONE_KEYS.index(stage) <= evolution.MILESTONE_KEYS.index(cap), (
                f"class {obj.planet_class} reports {stage}, above its cap {cap}")
        if cfg.INTELLIGENT_LIFE is False:
            assert stage != "technological_civilization"


def assert_round_trip(system):
    before = system.to_dict()
    via_json = json.loads(json.dumps(before))
    rebuilt = StarSystem.from_dict(via_json)
    after = rebuilt.to_dict()
    assert json.dumps(after, sort_keys=True) == json.dumps(before, sort_keys=True)
    # Second hop must be a fixed point too (no drift that only shows up later).
    assert StarSystem.from_dict(after).to_dict() == after
    assert (rebuilt.planet_count, rebuilt.belt_count, rebuilt.moon_count, rebuilt.hab_count,
            rebuilt.m_count, rebuilt.comet_count) == (
        system.planet_count, system.belt_count, system.moon_count, system.hab_count,
        system.m_count, system.comet_count)


def check_system(system, cfg=None, *, round_trip=True):
    data = system.to_dict()
    assert_serialized_numbers_sane(data)
    for planets in _all_planet_lists(system):
        assert_orbits_sane(planets)
    assert_binary_geometry(system)
    if cfg is not None:
        assert_flags_honored(system, cfg)
        assert_life_caps(system, cfg)
    for star in system.stars:
        if star.lifespan != float("inf"):
            assert star.age <= star.lifespan * 1.001
    if round_trip:
        assert_round_trip(system)


# ---------------------------------------------------------------------------
# Whole-system properties
# ---------------------------------------------------------------------------

@settings(max_examples=scaled(80))
@given(seed=SEEDS, cfg=system_configs())
def test_random_valid_configs_generate_sane_systems(seed, cfg):
    note(f"config: {cfg.to_dict()}")
    system = _generate(cfg, seed)
    check_system(system, cfg)
    # Rendering must not crash either, in the format the config asks for.
    assert str(system)


@pytest.mark.parametrize("yerkes", YERKES)
def test_every_star_type_generates_a_sane_system(yerkes):
    """All 70 spectral-letter x subclass combinations for one Yerkes class,
    each with every orbit-producing flag forced on (the most bodies, the
    most chances to break an invariant)."""
    for index, star_type in enumerate(t for t in ALL_STAR_TYPES if t[2:] == yerkes):
        cfg = SystemConfig()
        cfg.STAR_TYPE = star_type
        cfg.BINARY_SYSTEM = False
        cfg.MOONS = True
        cfg.COMETS = True
        cfg.MAX_PLANETS = True
        system = _generate(cfg, seed=YERKES.index(yerkes) * 1000 + index)
        try:
            check_system(system, cfg)
        except AssertionError as exc:
            raise AssertionError(f"{star_type}: {exc}") from exc


@settings(max_examples=scaled(30))
@given(
    seed=SEEDS,
    star_type=st.sampled_from(ALL_STAR_TYPES),
    num_orbits=st.sampled_from([0, 1, prog.ABSOLUTE_MAX_SYSTEM_OBJECTS, 10**6]),
    binary=st.sampled_from([(False, None), (True, False), (True, True)]),
)
def test_extreme_orbit_counts(seed, star_type, num_orbits, binary):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.NUM_ORBITS = num_orbits
    cfg.BINARY_SYSTEM, cfg.WIDE_BINARY = binary
    system = _generate(cfg, seed)
    check_system(system, cfg)
    if num_orbits == 0:
        assert system.planets == []
    assert len(system.planets) <= min(num_orbits, 10**6)


@settings(max_examples=scaled(25))
@given(seed=SEEDS, star_type=st.one_of(st.none(), st.sampled_from(ALL_STAR_TYPES)),
       wide=st.sampled_from([None, True, False]))
def test_every_forced_flag_on_together(seed, star_type, wide):
    cfg = SystemConfig()
    for attr in TRISTATE:
        setattr(cfg, attr, True)
    cfg.WIDE_BINARY = wide
    cfg.STAR_TYPE = star_type
    system = _generate(cfg, seed)
    check_system(system, cfg)


@settings(max_examples=scaled(25))
@given(seed=SEEDS, star_type=st.one_of(st.none(), st.sampled_from(ALL_STAR_TYPES)))
def test_every_forced_flag_off_together(seed, star_type):
    cfg = SystemConfig()
    for attr in TRISTATE:
        setattr(cfg, attr, False)
    cfg.STAR_TYPE = star_type
    system = _generate(cfg, seed)
    check_system(system, cfg)
    assert system.planets == [] and system.comets == [] and system.binary_type is None


@settings(max_examples=scaled(30), suppress_health_check=[HealthCheck.too_slow])
@given(seed=SEEDS, remnant_cls=st.sampled_from([BlackHole, NeutronStar]),
       num_orbits=st.one_of(st.none(), st.integers(0, 30)), moons=st.sampled_from([None, True, False]),
       galactic_dist=st.one_of(st.none(), st.floats(min_value=1.0, max_value=60000.0)))
def test_compact_remnant_anchored_systems(seed, remnant_cls, num_orbits, moons, galactic_dist):
    """`StarSystem(compact_remnant=...)` (`generate.py phenomenon
    --anchor-system`). JSON round trip is documented as unsupported for
    these (no `star_kind` discriminator), so only generation invariants
    are checked; a remnant's own zero-luminosity surface temperature is
    exempt from the CMB floor only if the remnant itself reports it."""
    cfg = SystemConfig()
    cfg.NUM_ORBITS = num_orbits
    cfg.MOONS = moons
    with _deterministic_entropy(seed):
        remnant = remnant_cls(cfg, galactic_center_dist_ly=galactic_dist)
        system = StarSystem(system_config=cfg, compact_remnant=remnant)
    assert system.star is remnant and system.binary_type is None
    zero_ok = ("$.temperature", "$.luminosity") if remnant_cls is BlackHole else ()
    assert_serialized_numbers_sane(remnant.to_dict(), zero_ok_paths=zero_ok)
    assert_serialized_numbers_sane({"planets": [p.to_dict() for p in system.planets]})
    assert_orbits_sane(system.planets)
    # The remnant alone round-trips through its own class.
    assert remnant_cls.from_dict(json.loads(json.dumps(remnant.to_dict())), cfg).to_dict() == remnant.to_dict()


@settings(max_examples=scaled(40))
@given(seed=SEEDS, slots=st.lists(SLOT_SPEC, min_size=1, max_size=10),
       star_type=st.sampled_from(ALL_STAR_TYPES))
def test_explicit_slot_lists(seed, slots, star_type):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.BINARY_SYSTEM = False
    cfg.SLOTS = slots
    cfg.NUM_ORBITS = len(slots)
    system = _generate(cfg, seed)
    check_system(system, cfg)
    for spec, obj in zip(slots, system.planets):
        if spec is not None and spec.get("type") == "asteroid_belt":
            assert obj.body_type == "a"
        if spec is not None and spec.get("type") == "planet":
            assert obj.body_type != "a"
            if spec.get("moons") is not None and spec["moons"] <= 0:
                assert obj.moons == []


def test_invalid_slot_type_and_unknown_planet_class_raise_value_error():
    for bad in ({"type": "comet"}, {"type": "planet", "planet_class": "Z"},
                {"type": "planet", "planet_class": "m"}):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        cfg.SLOTS = [bad]
        cfg.NUM_ORBITS = 1
        with pytest.raises(ValueError):
            _generate(cfg, seed=1)


def test_contradictory_slot_and_flag_is_a_clean_value_error():
    cfg = SystemConfig()
    cfg.HABITABLE_WORLD = False
    cfg.BINARY_SYSTEM = False
    cfg.SLOTS = [{"type": "planet", "planet_class": "M"}]
    cfg.NUM_ORBITS = 1
    with pytest.raises(ValueError):
        _generate(cfg, seed=0)


@settings(max_examples=scaled(40))
@given(bad=st.one_of(
    hostile_text.filter(lambda s: not s[:3].upper().startswith(tuple(f"{a}{d}" for a in SPECTRAL for d in "0123456789"))),
    st.sampled_from(["", "g", "G", "G2", "GV", "X2V", "2GV", " G2V", "G-2V", "G2I", "G2X"]),
))
def test_malformed_star_type_is_a_clean_value_error(bad):
    cfg = SystemConfig()
    cfg.STAR_TYPE = bad
    cfg.BINARY_SYSTEM = False
    if not bad:
        # Empty string means "not forced" (STAR_TYPE is falsy) -- must just work.
        _generate(cfg, seed=0)
        return
    with pytest.raises(ValueError):
        _generate(cfg, seed=0)


@settings(max_examples=scaled(30))
@given(star_type=st.builds(
    lambda head, tail: head + tail,
    st.sampled_from([f"{a}{d}{y}" for a in SPECTRAL for d in "0123456789" for y in YERKES]),
    hostile_text.filter(bool),
))
@example(star_type="G2Vjunk")
@example(star_type="G10V")
@example(star_type="G2V;DROP TABLE sectors")
def test_star_type_with_trailing_garbage_is_rejected(star_type):
    """Regression: generate_star used re.match (a prefix match), so
    'G2Vjunk' was accepted as G2V and 'G10V' parsed as G1 + Yerkes '0' (a
    hypergiant!) with the trailing 'V' dropped. A valid type followed by
    any extra text must be rejected -- unless the whole string happens to
    be another valid type (e.g. 'G2V' + 'I' = 'G2VI')."""
    assume(not re.fullmatch(r"[OBAFGKM][0-9](IA\+|IAB|VII|III|IA|IB|II|IV|VI|0|V|D)", star_type.upper()))
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.BINARY_SYSTEM = False
    with pytest.raises(ValueError):
        _generate(cfg, seed=0)


@settings(max_examples=scaled(20))
@given(seed=SEEDS)
def test_generation_is_reproducible_under_deterministic_entropy(seed):
    """Guards this file's own premise: with `secrets` routed through a
    seeded RNG, the same seed must give the byte-identical system -- if a
    new entropy source sneaks in, every shrunk failure above would stop
    reproducing."""
    cfg_a, cfg_b = SystemConfig(), SystemConfig()
    for cfg in (cfg_a, cfg_b):
        cfg.MOONS = True
        cfg.COMETS = True
    a = _generate(cfg_a, seed).to_dict()
    b = _generate(cfg_b, seed).to_dict()
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


# ---------------------------------------------------------------------------
# Real bugs (strict xfail): moons and circumbinary orbits
# ---------------------------------------------------------------------------

def _moon_pairs(system):
    for planets in _all_planet_lists(system):
        for planet in planets:
            if planet.body_type == "a":
                continue
            for moon in planet.moons:
                yield planet, moon


def _moon_heavy_system(seed):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.BINARY_SYSTEM = False
    cfg.MOONS = True
    cfg.MAX_PLANETS = True
    return _generate(cfg, seed)


def test_moons_orbit_inside_their_parents_hill_sphere():
    bad = []
    for seed in range(5):
        for planet, moon in _moon_pairs(_moon_heavy_system(seed)):
            if moon.distance * pc.AU_TO_KM > planet.hill_radius:
                bad.append((seed, planet.planet_class, moon.distance * pc.AU_TO_KM, planet.hill_radius))
    assert not bad, f"{len(bad)} moons outside their parent's Hill sphere, e.g. {bad[:3]}"


def test_moons_orbit_outside_their_parents_body():
    bad = []
    for seed in range(5):
        for planet, moon in _moon_pairs(_moon_heavy_system(seed)):
            if moon.distance * pc.AU_TO_KM <= planet.radius + moon.radius:
                bad.append((seed, planet.planet_class, planet.radius, moon.distance * pc.AU_TO_KM))
    assert not bad, f"{len(bad)} moons inside their parent's body, e.g. {bad[:3]}"


def test_circumbinary_bodies_orbit_outside_the_binary():
    bad = []
    for seed in range(60):
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = False
        system = _generate(cfg, seed)
        if system.planets:
            inner = system.planets[0]
            edge = inner.lower_limit if inner.body_type == "a" else inner.distance
            if edge <= system.star.binary_separation_au:
                bad.append((seed, edge, system.star.binary_separation_au))
    assert not bad, f"bodies inside the binary orbit (seed, edge AU, separation AU): {bad[:5]}"


def test_planet_hill_sphere_clears_the_belt_inside_it():
    bad = []
    for seed in range(40):
        system = _generate(SystemConfig(), seed)
        for planets in _all_planet_lists(system):
            for prev, cur in zip(planets, planets[1:]):
                if prev.body_type == "a" and cur.body_type != "a" and \
                        cur.distance - _hill_au(cur) < prev.upper_limit:
                    bad.append((seed, prev.upper_limit, cur.distance, _hill_au(cur)))
    assert not bad, f"(seed, belt outer edge AU, planet AU, planet r_H AU): {bad[:5]}"


def test_binary_secondary_is_never_heavier_than_primary():
    bad = []
    for seed in range(200):
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = seed % 2 == 0
        system = _generate(cfg, seed)
        if system.secondary_star.mass > system.primary_star.mass:
            bad.append((seed, system.binary_type, system.primary_star.type, system.secondary_star.type))
    assert not bad, f"secondary heavier than primary: {bad[:5]}"


# ---------------------------------------------------------------------------
# evolution / planetLife driven directly
# ---------------------------------------------------------------------------

@settings(max_examples=scaled(60))
@given(seed=SEEDS, star_type=st.sampled_from(ALL_STAR_TYPES),
       planet_class=st.one_of(st.none(), st.sampled_from(PLANET_CLASSES), hostile_text),
       intelligent=st.sampled_from([None, True, False]),
       age=st.sampled_from([None, "young", "old"]))
def test_evolutionary_timeline_respects_caps_for_any_star(seed, star_type, planet_class, intelligent, age):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.INTELLIGENT_LIFE = intelligent
    cfg.AGE = age
    with _deterministic_entropy(seed):
        star = Star(cfg)
        paragraphs = evolution.get_evolutionary_timeline(star, planet_class)
    assert isinstance(paragraphs, list) and all(isinstance(p, str) and p for p in paragraphs)
    text = " ".join(paragraphs)
    assert "nan" not in text.lower().split() and "inf" not in text.lower().split()
    stage = evolution.life_stage_from_paragraphs(paragraphs)
    if stage is not None:
        cap = prog.PLANET_CLASS_MAX_LIFE_STAGE.get(planet_class)
        if cap is not None:
            assert evolution.MILESTONE_KEYS.index(stage) <= evolution.MILESTONE_KEYS.index(cap)
        if intelligent is False:
            assert stage != "technological_civilization"
        if intelligent is True and cap is None:
            assert stage == "technological_civilization"


@settings(max_examples=scaled(30))
@given(seed=SEEDS)
def test_life_pass_is_stable_on_reload(seed):
    """`apply_life_data`/`decide_flavor_text` run once at generation; the
    reloaded planet must carry exactly the same life data (nothing re-rolled
    by `from_dict`) and re-running the life pass on a copy must not crash."""
    cfg = SystemConfig()
    cfg.HABITABLE_WORLD = True
    cfg.BINARY_SYSTEM = False
    system = _generate(cfg, seed)
    rebuilt = StarSystem.from_dict(system.to_dict())
    for a, b in zip(system.planets, rebuilt.planets):
        if a.body_type == "a":
            continue
        assert (a.life_chemical, a.evolutionary_speed, a.evolutionary_data, a.flavor_text) == \
               (b.life_chemical, b.evolutionary_speed, b.evolutionary_data, b.flavor_text)
        with _deterministic_entropy(seed):
            planetLife.apply_life_data(b)
            planetLife.decide_flavor_text(b)


# ---------------------------------------------------------------------------
# Lower-level physics helpers at their boundaries
# ---------------------------------------------------------------------------
#
# Contract (the one `calculate_orbital_period_years` already documents):
# a helper either returns finite numbers or raises ValueError/TypeError --
# never ZeroDivisionError/OverflowError, and never a silent NaN/inf.

M_SUN = pc.SOLAR_MASS_TO_KG
L_SUN = pc.SOLAR_LUMINOSITY

# name -> (callable, nominal float args)
HELPERS = {
    "calculate_orbital_period_years": (pp.calculate_orbital_period_years, (1.0, M_SUN)),
    "calculate_hill_sphere": (utils.calculate_hill_sphere, (1.5e11, 6e24, M_SUN)),
    "holman_wiegert_critical_semimajor_axis": (utils.holman_wiegert_critical_semimajor_axis, (100.0, 0.5, 0.3)),
    "mutual_hill_radius_au": (utils.mutual_hill_radius_au, (6e24, 6e24, 1.0, 1.5, M_SUN)),
    "mutual_hill_radius_m": (utils.mutual_hill_radius_m, (1.5e11, 2e11, 6e24, 6e24, M_SUN)),
    "snow_line_au": (utils.snow_line_au, (L_SUN,)),
    "disk_surface_density_scale": (utils.disk_surface_density_scale, (M_SUN,)),
    "mmsn_surface_density_gcm2": (utils.mmsn_surface_density_gcm2, (1.0, 2.7, 1.0)),
    "isolation_mass_kg": (utils.isolation_mass_kg, (1.0, 7.0, M_SUN)),
    "calculate_habitable_zone": (utils.calculate_habitable_zone, (L_SUN,)),
    "calculate_galactic_orbit": (utils.calculate_galactic_orbit, (26000.0,)),
    "circular_orbital_speed_kms": (utils.circular_orbital_speed_kms, (1.0, 1.0)),
    "orbital_position_au": (utils.orbital_position_au, (1.0, 10.0, 20.0, 30.0)),
    "minimum_update_interval_years": (utils.minimum_update_interval_years, (1.0,)),
    "gravitational_parameter_au3_yr2": (km.gravitational_parameter_au3_yr2, (1.0,)),
    "mean_motion_per_year": (km.mean_motion_per_year, (1.0, 1.0)),
    "solve_eccentric_anomaly": (km.solve_eccentric_anomaly, (1.0, 0.5)),
    "true_anomaly_and_distance_elliptical": (km.true_anomaly_and_distance_elliptical, (1.0, 0.5, 1.0)),
    "solve_barker_equation": (km.solve_barker_equation, (1.0,)),
    "parabolic_mean_anomaly": (km.parabolic_mean_anomaly, (1.0, 1.0, 1.0)),
    "true_anomaly_and_distance_parabolic": (km.true_anomaly_and_distance_parabolic, (1.0, 1.0)),
    "vis_viva_speed_kms": (km.vis_viva_speed_kms, (1.0, 1.0, 1.0)),
    "_activity_chance": (cometData._activity_chance, (1.0,)),
    "_atmosphere_retention_factor": (pp._atmosphere_retention_factor, (1.0,)),
    "_sample_evolved_star_mass_sol": (starData._sample_evolved_star_mass_sol, (1.0, 5.0)),
    "_calculate_heliosphere_radius_static": (
        Star._calculate_heliosphere_radius_static, (M_SUN, L_SUN, 696000.0, "G2V", "V")),
    "_calculate_system_perimeter_static": (BinaryStarProxy._calculate_system_perimeter_static, (M_SUN, 26000.0)),
    "sample_bounded_bell": (utils.sample_bounded_bell, (1.0, 2.0, 0.5)),
}

# Finite but hostile: zero, signed zero, subnormal, tiny, huge, negative.
FINITE_BOUNDARY = [0.0, -0.0, 5e-324, 1e-300, 1e-12, 1e12, 1e150, 1e300, -1.0, -1e300]
NON_FINITE = [math.nan, math.inf, -math.inf]


def _leaves(value):
    if isinstance(value, (tuple, list)):
        for item in value:
            yield from _leaves(item)
    else:
        yield value


def _outcome(fn, args):
    """None when the call honored the contract, else a short description."""
    try:
        with _deterministic_entropy(0):
            result = fn(*args)
    except (ValueError, TypeError):
        return None
    except Exception as exc:  # noqa: BLE001 - classifying is the point
        return f"{type(exc).__name__}: {exc}"
    for leaf in _leaves(result):
        if _is_number(leaf) and not math.isfinite(leaf):
            return f"silent {leaf!r}"
        if isinstance(leaf, str) and {"nan", "inf", "-inf"} & set(leaf.lower().split()):
            return f"non-finite text {leaf!r}"
    return None


def _grid_failures(name, values):
    fn, nominal = HELPERS[name]
    failures = []
    for index, nominal_value in enumerate(nominal):
        if not isinstance(nominal_value, float):
            continue
        for value in values:
            args = list(nominal)
            args[index] = value
            problem = _outcome(fn, args)
            if problem:
                failures.append(f"arg{index}={value!r}: {problem}")
    return failures


# Regression repros: each of these once broke the contract on a FINITE
# input (now fixed by `utils.finite_domain`); kept as exact examples.
FINITE_REPROS = [
    ("calculate_orbital_period_years", (1.0, 5e-324)),    # kg->Msun underflow: ZeroDivisionError
    ("calculate_orbital_period_years", (1e300, M_SUN)),   # OverflowError
    ("calculate_hill_sphere", (1.5e11, 6e24, 0.0)),       # ZeroDivisionError
    ("calculate_hill_sphere", (1.5e11, 6e24, 5e-324)),    # silent inf
    ("mutual_hill_radius_au", (6e24, 6e24, 1.0, 1.5, 0.0)),
    ("mutual_hill_radius_m", (1.5e11, 2e11, 6e24, 6e24, 0.0)),
    ("disk_surface_density_scale", (1e300,)),             # OverflowError
    ("mmsn_surface_density_gcm2", (0.0, 2.7, 1.0)),
    ("mmsn_surface_density_gcm2", (5e-324, 2.7, 1.0)),
    ("isolation_mass_kg", (1.0, 7.0, 0.0)),
    ("isolation_mass_kg", (1e300, 7.0, M_SUN)),           # silent inf
    ("calculate_galactic_orbit", (5e-324,)),
    ("calculate_galactic_orbit", (1e300,)),
    ("circular_orbital_speed_kms", (1.0, 0.0)),
    ("circular_orbital_speed_kms", (1.0, 5e-324)),
    ("mean_motion_per_year", (0.0, 1.0)),
    ("mean_motion_per_year", (1e300, 1.0)),
    ("solve_barker_equation", (1e300,)),                  # silent nan
    ("parabolic_mean_anomaly", (1.0, 0.0, 1.0)),
    ("parabolic_mean_anomaly", (1.0, 1e300, 1.0)),
    ("true_anomaly_and_distance_parabolic", (1e300, 1.0)),
    ("vis_viva_speed_kms", (0.0, 1.0, 1.0)),
    ("vis_viva_speed_kms", (5e-324, 1.0, 1.0)),
    ("_calculate_heliosphere_radius_static", (0.0, L_SUN, 696000.0, "G2V", "V")),
    ("_calculate_heliosphere_radius_static", (M_SUN, L_SUN, 1e300, "G2V", "V")),
    ("_calculate_system_perimeter_static", (M_SUN, 1e300)),
]


@pytest.mark.parametrize("name", list(HELPERS))
def test_helper_finite_boundary_inputs_are_clean(name):
    failures = _grid_failures(name, FINITE_BOUNDARY)
    assert not failures, f"{name}: " + "; ".join(failures[:8])


@pytest.mark.parametrize("name,args", FINITE_REPROS, ids=[f"{n}{a[:2]}" for n, a in FINITE_REPROS])
def test_helper_finite_repro_is_clean(name, args):
    assert _outcome(HELPERS[name][0], args) is None


@pytest.mark.parametrize("name", list(HELPERS))
def test_helpers_reject_non_finite_inputs(name):
    """Every helper here validates its domain, so NaN/inf input must be a
    ValueError too (regression: `calculate_orbital_period_years(nan, M)`
    returned nan because `nan <= 0` is False) -- or, for an argument the
    helper clamps (holman_wiegert's mu/e, sample_bounded_bell's mode) or
    documents as allowed (vis_viva's a=inf), a finite result."""
    failures = _grid_failures(name, NON_FINITE)
    assert not failures, f"{name}: " + "; ".join(failures[:8])


@given(fn=st.sampled_from(sorted(HELPERS)), index=st.integers(0, 4), bad=st.sampled_from(NON_FINITE))
@example(fn="calculate_orbital_period_years", index=0, bad=math.nan)
@example(fn="calculate_orbital_period_years", index=1, bad=math.nan)
def test_non_finite_argument_never_escapes_as_a_number(fn, index, bad):
    func, nominal = HELPERS[fn]
    assume(index < len(nominal) and isinstance(nominal[index], float))
    args = list(nominal)
    args[index] = bad
    assert _outcome(func, args) is None


# Physically valid domains: every helper must give finite, sensible output.
POS = st.floats(min_value=1e-6, max_value=1e9, allow_nan=False, allow_infinity=False)
MASS = st.floats(min_value=1e20, max_value=1e33)
LUM = st.floats(min_value=1e20, max_value=1e34)
ANGLE = st.floats(min_value=-720.0, max_value=720.0)
ECC = st.floats(min_value=0.0, max_value=0.999)


def _finite(*values):
    return all(math.isfinite(v) for v in _leaves(values) if _is_number(v))


@given(d=POS, m=MASS)
def test_valid_domain_orbital_period(d, m):
    period = pp.calculate_orbital_period_years(d, m)
    assert math.isfinite(period) and period > 0
    # Kepler's third law, inverted, recovers the distance.
    assert (period ** 2 * (m / M_SUN)) ** (1 / 3) == pytest.approx(d, rel=1e-9)


@given(d=POS, m=MASS, big=MASS)
def test_valid_domain_hill_spheres(d, m, big):
    assume(m < big)
    r = utils.calculate_hill_sphere(d, m, big)
    assert math.isfinite(r) and 0 < r < d
    mutual = utils.mutual_hill_radius_au(m, m, d, d, big)
    assert math.isfinite(mutual) and mutual > 0
    assert mutual == pytest.approx(utils.mutual_hill_radius_m(d, d, m, m, big), rel=1e-12)


@given(a=POS, mu=st.floats(0.0, 1.0), e=st.floats(0.0, 1.0))
def test_valid_domain_holman_wiegert(a, mu, e):
    a_crit = utils.holman_wiegert_critical_semimajor_axis(a, mu, e)
    assert math.isfinite(a_crit) and 0 < a_crit < a


@given(lum=LUM, mass=MASS, d=POS)
def test_valid_domain_disk_helpers(lum, mass, d):
    inner, outer = utils.calculate_habitable_zone(lum)
    assert 0 < inner < outer and _finite(inner, outer)
    snow = utils.snow_line_au(lum)
    sigma = utils.mmsn_surface_density_gcm2(d, snow, utils.disk_surface_density_scale(mass))
    iso = utils.isolation_mass_kg(d, sigma, mass)
    assert _finite(snow, sigma, iso) and snow > 0 and sigma > 0 and iso > 0


@given(d=POS, period=POS, inc=ANGLE, node=ANGLE, phase=ANGLE)
def test_valid_domain_orbital_position_and_speed(d, period, inc, node, phase):
    x, y, z = utils.orbital_position_au(d, inc, node, phase)
    assert math.sqrt(x * x + y * y + z * z) == pytest.approx(d, rel=1e-9)
    speed = utils.circular_orbital_speed_kms(d, period)
    interval = utils.minimum_update_interval_years(period)
    assert _finite(speed, interval) and speed > 0 and interval > 0


@given(mean_anomaly=st.floats(-1e6, 1e6), e=ECC, a=POS, m=st.floats(1e-3, 300.0))
def test_valid_domain_kepler_elliptical(mean_anomaly, e, a, m):
    ecc_anomaly = km.solve_eccentric_anomaly(mean_anomaly, e)
    nu, r = km.true_anomaly_and_distance_elliptical(mean_anomaly, e, a)
    assert _finite(ecc_anomaly, nu, r)
    assert a * (1 - e) * (1 - 1e-9) <= r <= a * (1 + e) * (1 + 1e-9)
    speed = km.vis_viva_speed_kms(r, a, m)
    assert math.isfinite(speed) and speed > 0
    assert math.isfinite(km.mean_motion_per_year(a, m))


@given(t=st.floats(-1e5, 1e5), q=st.floats(1e-3, 1e4), m=st.floats(1e-3, 300.0))
def test_valid_domain_kepler_parabolic(t, q, m):
    big_m = km.parabolic_mean_anomaly(t, q, m)
    nu, r = km.true_anomaly_and_distance_parabolic(big_m, q)
    assert _finite(big_m, nu, r)
    assert -math.pi < nu < math.pi and r >= q * (1 - 1e-9)


@given(dist=st.floats(1e-3, 1e6), mass=MASS)
def test_valid_domain_galactic_and_perimeter(dist, mass):
    speed, period = utils.calculate_galactic_orbit(dist)
    perimeter = BinaryStarProxy._calculate_system_perimeter_static(mass, dist)
    assert _finite(speed, period, perimeter) and speed > 0 and period > 0 and perimeter > 0


@given(lo=st.floats(-1e6, 1e6), width=st.floats(1e-6, 1e6), mode=st.floats(0.0, 1.0),
       seed=SEEDS)
def test_valid_domain_bounded_bell(lo, width, mode, seed):
    hi = lo + width
    assume(hi > lo)
    with _deterministic_entropy(seed):
        value = utils.sample_bounded_bell(lo, hi, mode)
    assert math.isfinite(value) and lo <= value <= hi


@given(years=st.floats(1 / (365.25 * 24 * 60), 1e15), age=st.floats(1e-6, 1e4))
def test_valid_domain_text_formatters(years, age):
    for text in (utils.years_to_time_string(years), utils.format_age_string(age)):
        assert isinstance(text, str) and text
        assert not {"nan", "inf", "-inf"} & set(text.lower().split())
