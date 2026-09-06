"""
Full-system generation regression tests.

The true cross product of every SystemConfig tri-state flag x AGE x star
type x NUM_ORBITS x SLOTS is enormous, so this sweeps a representative
design instead: a spread of star types (covering every Yerkes class) each
crossed with every tri-state flag forced True, forced False, and left at
its default (None), plus dedicated tests for AGE and BINARY_SYSTEM. Every
generated system is checked against invariants systemData.py itself claims
to guarantee (no orbital overlap, requested features actually present,
counts matching reality, star age never exceeding lifespan even after
adjust_age_for_planets runs against the system's actual planets).

Run with: pytest tests/test_systems.py
"""
import math

import pytest

from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem
from stellarObjects.doubleStar import BinaryStarProxy
from stellarObjects import program_constants as prog_c

# One representative star type per Yerkes class, spanning several spectral
# letters, so the system-generation sweep exercises every evolutionary track
# (main sequence, giant, subgiant, bright giant, supergiant, hypergiant,
# subdwarf, white dwarf) without re-running the full 770-combo star matrix.
STAR_TYPES = [
    "G2V", "M5V", "K3V", "A0V", "B2V", "O5V",
    "K1III", "F0IB", "M2VII", "K5IV", "A0II", "M3IA+", "M50", "K8VI",
]

TRISTATE_ATTRS = [
    "HABITABLE_WORLD", "ASTEROID_BELT", "LARGE_STAR", "MOONS",
    "MAX_PLANETS", "INTELLIGENT_LIFE", "BINARY_SYSTEM", "PLANETS",
]

TRIALS = 2


def make_config(star_type, **overrides):
    """
    Builds a SystemConfig the same way systemGen.main() would: applying the
    same normalization it does (INTELLIGENT_LIFE implies HABITABLE_WORLD;
    HABITABLE_WORLD + ASTEROID_BELT together imply LARGE_STAR) so configs
    built directly here stay consistent with what the CLI ever actually
    produces.
    """
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    for attr, value in overrides.items():
        setattr(cfg, attr, value)
    if cfg.INTELLIGENT_LIFE is not None:
        cfg.HABITABLE_WORLD = True
    if cfg.HABITABLE_WORLD is True and cfg.ASTEROID_BELT is True:
        cfg.LARGE_STAR = True
    return cfg


def all_stars(system):
    return system.stars


def assert_ages_never_exceed_lifespan(system):
    for star in all_stars(system):
        if star.lifespan == float('inf'):
            continue
        assert star.age <= star.lifespan * 1.001, (
            f"{star.name} ({star.type}): age {star.age} Gy exceeds lifespan {star.lifespan} Gy "
            f"after adjust_age_for_planets"
        )


def assert_no_orbital_overlap(system):
    """
    Mirrors StarSystem.validate_system's own separation case matrix, as a
    check that its correction pass actually converges rather than as
    independent new physics -- the minimums it enforces are its own.
    """
    objects = system.planets
    for i in range(1, len(objects)):
        prev, cur = objects[i - 1], objects[i]
        gap = cur.distance - prev.upper_limit if prev.type == 'a' else cur.distance - prev.distance

        if prev.type == 'a':
            # validate_system requires MIN_ASTEROID_BELT_SEPARATION after a
            # belt regardless of what follows it (belt or planet).
            min_gap = prog_c.MIN_ASTEROID_BELT_SEPARATION
        elif cur.type == 'a':
            min_gap = prev.min_orbit_distance
        else:
            min_gap = max(cur.min_orbit_distance, prev.min_orbit_distance)

        assert gap >= min_gap - 1e-9, (
            f"objects {i - 1} ({prev.type}) and {i} ({cur.type}) overlap: gap={gap}, required>={min_gap}"
        )


def assert_counts_are_consistent(system):
    planet_count, belt_count, moon_count = system.count_objects()
    assert planet_count == system.planet_count
    assert belt_count == system.belt_count
    assert moon_count == system.moon_count
    assert planet_count + belt_count == len(system.planets)
    hab_count, m_count = system.count_habitable()
    assert hab_count == system.hab_count
    assert m_count == system.m_count
    assert hab_count >= m_count >= 0


@pytest.mark.parametrize("star_type", STAR_TYPES)
def test_baseline_random_system_generates_without_error(star_type):
    for _ in range(TRIALS):
        system = StarSystem(system_config=make_config(star_type))
        assert str(system).strip()
        assert_ages_never_exceed_lifespan(system)
        assert_no_orbital_overlap(system)
        assert_counts_are_consistent(system)


@pytest.mark.parametrize("attr", TRISTATE_ATTRS)
@pytest.mark.parametrize("value", [True, False])
@pytest.mark.parametrize("star_type", STAR_TYPES)
def test_each_tristate_flag_forced(star_type, value, attr):
    for _ in range(TRIALS):
        system = StarSystem(system_config=make_config(star_type, **{attr: value}))
        assert str(system).strip()
        assert_ages_never_exceed_lifespan(system)
        assert_no_orbital_overlap(system)
        assert_counts_are_consistent(system)

        if attr == "HABITABLE_WORLD" and value is True:
            assert system.hab_count >= 1, f"{star_type}: HABITABLE_WORLD=True produced no habitable world"
        if attr == "ASTEROID_BELT" and value is True:
            assert system.belt_count >= 1, f"{star_type}: ASTEROID_BELT=True produced no belt"
        if attr == "PLANETS" and value is False:
            assert len(system.planets) == 0
        if attr == "MOONS" and value is False:
            assert system.moon_count == 0
        if attr == "BINARY_SYSTEM" and value is True:
            assert isinstance(system.star, BinaryStarProxy)
            assert len(system.stars) == 2
            assert math.isclose(system.star.mass, sum(s.mass for s in system.stars))
            assert all(s.mass > 0 for s in system.stars)


def test_asteroid_belt_forcing_is_reliable():
    """
    Regression guard for a fixed bug: StarSystem's orbit-placement loop used
    to pick a single random belt_index up front and only place a belt there
    if that slot's estimated_distance also happened to fall outside the
    habitable zone, with no retry -- unlike HABITABLE_WORLD, which actively
    steers/retries placement through the whole loop. Measured failure rates
    before the fix (100 trials each): G2V 2%, M5V 27%, K3V 2%, A0V 9%, B2V 0%,
    K1III 3%, M2VII 4%, K5IV 4%, K8VI 15%.

    Fixed in systemData.py's shared orbit loop by giving ASTEROID_BELT its own
    guaranteed last-resort fallback slot (mirroring HABITABLE_WORLD's), with a
    reserved separate slot when both are forced simultaneously so neither's
    fallback can collide with or overwrite the other's placement. This test
    uses M5V (highest observed pre-fix rate) over 60 trials so it would
    reproduce the bug reliably (< 1e-6 chance of spuriously passing) if it
    ever regressed.
    """
    failures = 0
    trials = 60
    for _ in range(trials):
        system = StarSystem(system_config=make_config("M5V", ASTEROID_BELT=True))
        if system.belt_count == 0:
            failures += 1
    assert failures == 0, f"ASTEROID_BELT=True failed to produce a belt in {failures}/{trials} trials"


@pytest.mark.parametrize("star_type", STAR_TYPES)
def test_asteroid_belt_and_habitable_world_forced_together_both_succeed(star_type):
    """
    Regression guard for a second bug found while fixing the above: the
    pre-existing HABITABLE_WORLD "retroactive replace" logic could overwrite
    the immediately preceding slot with a new M-class planet -- including a
    belt just placed to satisfy ASTEROID_BELT's own guarantee -- since it had
    no awareness of that guarantee. Fixed by protecting a belt from being
    retroactively overwritten whenever ASTEROID_BELT is required.
    """
    failures_belt = 0
    failures_hab = 0
    trials = 15
    for _ in range(trials):
        system = StarSystem(system_config=make_config(star_type, ASTEROID_BELT=True, HABITABLE_WORLD=True))
        if system.belt_count == 0:
            failures_belt += 1
        if system.hab_count == 0:
            failures_hab += 1
    assert failures_belt == 0, f"{star_type}: ASTEROID_BELT+HABITABLE_WORLD failed to produce a belt in {failures_belt}/{trials} trials"
    assert failures_hab == 0, f"{star_type}: ASTEROID_BELT+HABITABLE_WORLD failed to produce a habitable world in {failures_hab}/{trials} trials"


@pytest.mark.parametrize("age", ["young", "old", None])
@pytest.mark.parametrize("star_type", ["G2V", "K1III", "M2VII", "B2V"])
def test_age_flag_propagates_and_respects_lifespan(star_type, age):
    for _ in range(TRIALS):
        system = StarSystem(system_config=make_config(star_type, AGE=age))
        assert_ages_never_exceed_lifespan(system)


def test_intelligent_life_implies_habitable_world_ends_true():
    for value in [True, False]:
        cfg = make_config("G2V", INTELLIGENT_LIFE=value)
        assert cfg.HABITABLE_WORLD is True
        system = StarSystem(system_config=cfg)
        assert system.hab_count >= 1


def test_binary_system_star_properties_are_sane():
    for star_type in ["G2V", "M2VII", "O5V"]:
        for _ in range(TRIALS):
            system = StarSystem(system_config=make_config(star_type, BINARY_SYSTEM=True))
            proxy = system.star
            assert isinstance(proxy, BinaryStarProxy)
            primary, secondary = proxy.stars
            assert secondary.mass <= primary.mass
            # Not a tighter bound than "positive": doubleStar.py draws the
            # secondary's mass as a 0.1-0.8 fraction of the primary's, but for
            # a primary near its Yerkes class's own minimum mass (e.g. a
            # white dwarf near the 0.5 Msun floor), generate_star's own
            # min-mass clamp on mass_override can correctly push the
            # secondary above that naive fraction.
            assert 0 < secondary.mass
            assert proxy.system_perimeter > 0 and math.isfinite(proxy.system_perimeter)
            assert proxy.heliosphere_radius > 0 and math.isfinite(proxy.heliosphere_radius)
