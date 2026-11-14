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

# A small subset of the full 150-seed fuzz list -- a full StarSystem
# generation (with binaries) is far more expensive per-iteration than the
# pure-function fuzz tests, so this trades seed-count for time (TEST.734: 5
# seeds, was 20, about 400 s off a full run; the shapes each test names are
# built by the test itself, not drawn by the seed, so five seeds still reach
# each one).
_SEEDS = FUZZ_SEEDS[:5]

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

@pytest.mark.parametrize("cls", [*STANDALONE_PHENOMENON_CLASSES, BlackHole, NeutronStar])
def test_phenomenon_roundtrip_fuzz(cls):
    """The standalone phenomena and the compact remnants share one body (TEST.734: two tests had the same one)."""
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
# property, which draws it on first read from `_name_seed` (PERF.43). A planet's `phi_bio` and `phi_cpx` are the
# likelihoods its life draw and life scores are made from; the stored life scores say what they gave (GEN.649). A moon's
# `moon_disk_system` is the giant's disk while its moons are drawn; the giant keeps the disk's numbers in columns (DB.70), as a moon does its
# `moon_origin`. A planet's `formation_distance_au` is where it formed, which for all but a pulsar's disk planet is the orbit distance it is
# stored with; what the draws read it for is kept as `born_beyond_snow_line` (TEST.427). A planet's `_orbit_birth` is the eccentricity
# and tidal state it was drawn with, which `settle_orbit` takes again once the star's age is final (GEN.985); it is not stored. A giant's `_giant_draw` is the heavy-element draws that
# `giant_draw.settle` keeps (GEN.961), also runtime state. A planet's `_surface_w` is the surface water its air solve drew, kept so the final-age
# refresh (GEN.839) solves again with the same draw. A moon's `_belt_exposure_msv_yr` is the dose its parent's belts give at its orbit before the air
# cuts it (GEN.983); the stored `dose_belt_msv_yr` is that after the air.
# A planet's `_air_replay` is the state of the draws of its air, kept so the end-of-making pass can make the air again with the late accretion it read (GEN.1493); runtime state, spent to None and never stored. Its `_late_accreted_percent` is the late accretion the pass set (GEN.1493), also not stored.
# A star's `initial_mass_sol`, `phase_end_age_gy`, `_model_age_lifespan` and `_model_fe_h` are the population model's
# own state, which a forced type shares (GEN.941); they are runtime state that the stored columns rebuild. `_model_rotation_hours` is the
# spin period a scattered star was placed with (GEN.846), which the built star takes from the `bright_stars` row it was built from.
_RUNTIME_ONLY = {
    Star: {"system_config", "galactic_center_dist_ly", "spatial", "_name", "_name_seed", "initial_mass_sol",
           "phase_end_age_gy", "_model_age_lifespan", "_model_fe_h", "_model_rotation_hours",
           "state", "decay_age_yr",
           "_companions"},    # GEN.1441: the other radiating objects a star's bodies see besides the host, attached by the system from its tree (`radiation_field.attach`), not stored; derived with no draw at __init__ from the class and the phase (GEN.785, GEN.786); the store writes the columns
    Planet: {"system_config", "star", "moons", "spatial", "_staged_au", "_staged_velocity_kms",
             "primary_mass_kg", "_star_distance_au", "phi_bio", "phi_cpx",
             "moon_disk_system", "formation_distance_au", "_orbit_birth", "_giant_draw", "_age_drawn",
             "_surface_w", "_belt_exposure_msv_yr", "_birth_spin", "_born_moons", "fallen_moon_mass_kg", "fallen_moons", "_birth_distance_au",
             "_formation_fractions", "irregulars", "_deferred_climate",
             "_bare_core", "_captured_envelope", "_late_core",
             "_air_choice", "_air_context", "belt_member", "_air_haze", "_co2_sink", "_tidal_radius_km", "_resonance_e", "_kozai", "_retention_cap_pa", "_venus_spin", "oddity_quantities", "_frame_axis", "_air_light_au",
             "_air_replay", "_late_accreted_percent", "_air_stale", "_air_eccentricity", "_tide_parent", "_tide_siblings", "_envelope_helium_y"},    # GEN.1047 to 1051, GEN.1340: the air branch and its context the oddity rows read, the resonance eccentricity, the retention cap the life stage keeps its oxygen inside, the Venus-type spin and the oddity quantities, the helium share of a hydrogen-first loss (GEN.1527, rebuilt by the settle pass); the stored columns (`oddities`, the eccentricity, the spin, the air) hold their result;    # the core the envelope was captured onto, the envelope it kept and the late-core request, which only the draw reads (GEN.1323); the stored columns hold their result; `irregulars` is written by `to_dict` as its own entry (GEN.1083); the solid fractions drawn before the envelope took its share (GEN.1226); the stored fractions are the ones after; `_deferred_climate` is the climate step held until the moons are known (GEN.1282), None once it ran    # kept so the moon history can run again (GEN.959); the ring it builds is stored
    AsteroidBelt: {"system_config", "composition", "makeup", "_members", "cleared_spans", "scatter_giant_au", "captured_mass_kg", "oddity_quantities", "swept_kg", "implanted_kg", "scatter_giant_kg", "scatter_giant_radius_km", "_secular_planets"},    # where the mass an outer disk lost went and the giant that scattered it (GEN.1308, not stored); the giant's distance the members of an outer disk are drawn with (GEN.1092); the stored members hold what it gave; the mass the planets caught (GEN.1085), already out of `mass_kg`
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
