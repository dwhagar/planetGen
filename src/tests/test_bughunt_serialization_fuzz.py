# tests/test_bughunt_serialization_fuzz.py

"""
Tier 1 bug-hunt coverage: `to_dict`/`from_dict` round-trip fuzzing across
every class using the shared allowlist-based serialization pattern
(`serialization.fields_to_dict`/`fields_from_dict` -- see that module's
own docstring). Existing tests (`test_serialization.py`,
`test_db_persistence.py`) already cover the common single/binary-system
case; this file drives many more randomized system/phenomenon shapes
(seeded, via `bughunt_support.run_seeded`) through the same round trip and
compares `to_dict()` before/after -- a lossy round trip for any field (a
type coercion, a dropped nested value, a stale cached value) shows up as a
dict-equality mismatch here regardless of which specific class or seed
triggers it.

Also fails on any instance attribute of a star, planet, moon or belt NOT
in that class's own `SERIALIZABLE_FIELDS` and not on the short known list
of runtime-only attributes (back-references like `Planet.star`, children
serialized on their own like `Planet.moons`) -- the "field drift" risk
`serialization.py`'s own docstring names: an attribute silently isn't
persisted until someone adds it to the allowlist.
"""

import random

import pytest

from tests.bughunt_support import FUZZ_SEEDS

from planetgen.generation.belt import AsteroidBelt
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.comet import Comet
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.config import SystemConfig
from planetgen.generation.binary import BinaryStarProxy
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.planet import Planet
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.star import Star
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.system import StarSystem
from planetgen.generation.wide_binary import WideBinaryPair

# A modest subset of the full 150-seed fuzz list -- a full StarSystem
# generation (with binaries) is far more expensive per-iteration than the
# pure-function fuzz tests, so this trades seed-count for still covering
# every star-type/binary-shape combination below at least once.
_SEEDS = FUZZ_SEEDS[:20]

STANDALONE_PHENOMENON_CLASSES = [Nebula, SupernovaRemnant, RoguePlanet, InterstellarComet, AsteroidField]


def _make_config(**overrides):
    cfg = SystemConfig()
    for attr, value in overrides.items():
        setattr(cfg, attr, value)
    return cfg


def _roundtrip_ok(obj, cls, *from_dict_args):
    before = obj.to_dict()
    reconstructed = cls.from_dict(before, *from_dict_args)
    after = reconstructed.to_dict()
    return before == after, before, after


def _assert_roundtrip(obj, cls, label, *from_dict_args):
    ok, before, after = _roundtrip_ok(obj, cls, *from_dict_args)
    if not ok:
        diffs = {k: (before.get(k), after.get(k)) for k in before if before.get(k) != after.get(k)}
        pytest.fail(f"{label} ({cls.__name__}) round-trip mismatch: {diffs}")


# --- SystemConfig ---------------------------------------------------------

def test_system_config_roundtrip_fuzz():
    def check(seed):
        cfg = SystemConfig()
        cfg.HABITABLE_WORLD = random.choice([True, False, None])
        cfg.ASTEROID_BELT = random.choice([True, False, None])
        cfg.NUM_ORBITS = random.choice([None, 0, 1, 12])
        cfg.STAR_TYPE = random.choice([None, "G2V", "M5V", "O3I"])
        _assert_roundtrip(cfg, SystemConfig, f"SystemConfig seed={seed}")

    from tests.bughunt_support import run_seeded
    run_seeded(check, seeds=_SEEDS)


# --- Standalone phenomena --------------------------------------------------

@pytest.mark.parametrize("cls", STANDALONE_PHENOMENON_CLASSES)
def test_standalone_phenomenon_roundtrip_fuzz(cls):
    def check(seed):
        cfg = _make_config()
        obj = cls(cfg)
        _assert_roundtrip(obj, cls, f"{cls.__name__} seed={seed}", cfg)

    from tests.bughunt_support import run_seeded
    run_seeded(check, seeds=_SEEDS)


@pytest.mark.parametrize("cls", [BlackHole, NeutronStar])
def test_compact_remnant_roundtrip_fuzz(cls):
    def check(seed):
        cfg = _make_config()
        obj = cls(cfg)
        _assert_roundtrip(obj, cls, f"{cls.__name__} seed={seed}", cfg)

    from tests.bughunt_support import run_seeded
    run_seeded(check, seeds=_SEEDS)


# --- Full StarSystem generation: Star/Planet/AsteroidBelt/Comet/binaries --

_STAR_TYPES = ["M5V", "K2V", "G2V", "F5V", "A1V", "B3V", "O5V"]


# Attributes a generated body may hold that `SERIALIZABLE_FIELDS` leaves
# out on purpose: back-references and the config (rebuilt by `from_dict`),
# children `to_dict` writes itself (moons, a belt's composition), and the
# galactic distance the star's serialized perimeter and orbit come from, and
# the `SpatialPosition3D` a body holds (GEN.74), which the stored position
# columns rebuild on load. A star's name is serialized through its `name`
# property, which draws it on first read from `_name_seed` (PERF.43).
_RUNTIME_ONLY = {
    Star: {"system_config", "galactic_center_dist_ly", "spatial", "_name", "_name_seed"},
    Planet: {"system_config", "star", "moons", "spatial", "_staged_au", "_staged_velocity_kms",
             "primary_mass_kg", "_star_distance_au"},
    AsteroidBelt: {"system_config", "composition"},
}


def test_starsystem_children_have_no_unlisted_attributes():
    """Every attribute generation sets is either serialized or known to be
    runtime-only, so a new field can't silently drop out of saves (TEST.4:
    was a Tier 2 report that listed the same expected extras every run)."""
    from tests.bughunt_support import run_seeded

    def check(seed):
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = False
        cfg.STAR_TYPE = random.choice(_STAR_TYPES)
        system = StarSystem(system_config=cfg)
        bodies = list(system.stars)
        for obj in system.planets:
            bodies.append(obj)
            bodies.extend(getattr(obj, "moons", None) or [])
        for body in bodies:
            cls = type(body)
            if cls not in _RUNTIME_ONLY:
                continue
            extra = set(vars(body)) - set(cls.SERIALIZABLE_FIELDS) - _RUNTIME_ONLY[cls]
            assert not extra, f"{cls.__name__} seed={seed}: attrs not in SERIALIZABLE_FIELDS: {sorted(extra)}"

    run_seeded(check, seeds=_SEEDS)


def test_starsystem_star_and_planet_roundtrip_fuzz():
    def check(seed):
        star_type = random.choice(_STAR_TYPES)
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = False
        cfg.STAR_TYPE = star_type
        system = StarSystem(system_config=cfg)

        for star in system.stars:
            _assert_roundtrip(star, Star, f"Star seed={seed} type={star_type}", cfg)

        for obj in system.planets:
            if isinstance(obj, AsteroidBelt):
                _assert_roundtrip(obj, AsteroidBelt, f"AsteroidBelt seed={seed}", cfg)
            elif isinstance(obj, Planet):
                _assert_roundtrip(obj, Planet, f"Planet seed={seed}", system.star, cfg)
                for moon in obj.moons:
                    _assert_roundtrip(moon, Planet, f"Moon seed={seed}", system.star, cfg)

        for comet in getattr(system, "comets", []) or []:
            _assert_roundtrip(comet, Comet, f"Comet seed={seed}", cfg)

    from tests.bughunt_support import run_seeded
    run_seeded(check, seeds=_SEEDS)


def test_starsystem_close_binary_proxy_roundtrip_fuzz():
    def check(seed):
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = False
        star_type = random.choice(_STAR_TYPES)
        cfg.STAR_TYPE = star_type
        system = StarSystem(system_config=cfg)
        assert isinstance(system.star, BinaryStarProxy)
        _assert_roundtrip(system.star, BinaryStarProxy, f"BinaryStarProxy seed={seed}", cfg)

    from tests.bughunt_support import run_seeded
    run_seeded(check, seeds=_SEEDS)


def test_starsystem_wide_binary_pair_roundtrip_fuzz():
    def check(seed):
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = True
        star_type = random.choice(_STAR_TYPES)
        cfg.STAR_TYPE = star_type
        system = StarSystem(system_config=cfg)
        assert system.wide_binary is not None
        _assert_roundtrip(
            system.wide_binary, WideBinaryPair, f"WideBinaryPair seed={seed}",
            cfg, system.primary_star, system.secondary_star,
        )

    from tests.bughunt_support import run_seeded
    run_seeded(check, seeds=_SEEDS)
