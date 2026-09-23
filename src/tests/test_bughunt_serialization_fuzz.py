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

Also reports (Tier 2 -- informational, not a hard failure) any instance
attribute NOT present in that class's own `SERIALIZABLE_FIELDS` -- a
currently-undetectable "field drift" risk (`serialization.py`'s own
docstring: a class attribute silently isn't persisted until someone adds
it to the allowlist). This is intentionally soft: several classes hold
legitimate non-serialized attributes (back-references like `Planet.star`,
recursively-serialized children like `Planet.moons`) that would otherwise
show up as constant false positives.
"""

import random

import pytest

from tests.bughunt_support import FUZZ_SEEDS, tier2

from stellarObjects.asteroidData import AsteroidBelt
from stellarObjects.asteroidFieldData import AsteroidField
from stellarObjects.cometData import Comet
from stellarObjects.compactRemnant import BlackHole, NeutronStar
from stellarObjects.config import SystemConfig
from stellarObjects.doubleStar import BinaryStarProxy
from stellarObjects.nebulaData import Nebula
from stellarObjects.planetData import Planet
from stellarObjects.roguePlanetData import InterstellarComet, RoguePlanet
from stellarObjects.starData import Star
from stellarObjects.supernovaRemnantData import SupernovaRemnant
from stellarObjects.systemData import StarSystem
from stellarObjects.wideBinary import WideBinaryPair

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


@tier2
def test_starsystem_children_roundtrip_fuzz_soft_field_drift_report():
    """Tier 2 companion to the Tier 1 test below -- same generation loop,
    but only responsible for the soft "field not in SERIALIZABLE_FIELDS"
    report, kept separate so a Tier 1 failure there never masks it."""
    from tests.bughunt_support import Tier2Report, run_seeded

    report = Tier2Report()

    def check(seed):
        star_type = random.choice(_STAR_TYPES)
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = False
        cfg.STAR_TYPE = star_type
        system = StarSystem(system_config=cfg)
        for star in system.stars:
            extra = set(vars(star)) - set(Star.SERIALIZABLE_FIELDS)
            if extra:
                report.add(f"Star seed={seed}: attrs not in SERIALIZABLE_FIELDS: {sorted(extra)}")
        for obj in system.planets:
            if isinstance(obj, Planet):
                extra = set(vars(obj)) - set(Planet.SERIALIZABLE_FIELDS)
                if extra:
                    report.add(f"Planet seed={seed}: attrs not in SERIALIZABLE_FIELDS: {sorted(extra)}")

    run_seeded(check, seeds=_SEEDS[:10])
    report.flush("starsystem-field-drift")


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
