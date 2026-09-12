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
from stellarObjects import physical_constants, program_constants as prog_c
from stellarObjects.utils import (circular_orbital_speed_kms, minimum_update_interval_years,
                                   mutual_hill_radius_m, orbital_position_au)

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

    Defaults `BINARY_SYSTEM` to False: since `StarSystem._should_generate_binary`
    rolls real, spectral-class-dependent chance whenever it's left at its
    own default (None), every test in this module that isn't specifically
    about binary generation would otherwise sometimes get a companion star
    unpredictably -- most incidentally harmless, but a wide (S-type) pair's
    own `a_crit_au` stability ceiling can make a forced HABITABLE_WORLD/
    ASTEROID_BELT guarantee geometrically impossible for an especially
    luminous host (e.g. an O-type supergiant's habitable zone can sit
    beyond any sampled companion's stability limit) -- an explicit
    `BINARY_SYSTEM=...` override (including the dedicated tri-state sweep
    below) still takes precedence over this default.
    """
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.BINARY_SYSTEM = False
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


def _assert_no_overlap_within(objects):
    for i in range(1, len(objects)):
        prev, cur = objects[i - 1], objects[i]
        gap = cur.distance - prev.upper_limit if prev.body_type == 'a' else cur.distance - prev.distance

        if prev.body_type == 'a':
            # validate_system requires MIN_ASTEROID_BELT_SEPARATION after a
            # belt regardless of what follows it (belt or planet).
            min_gap = prog_c.MIN_ASTEROID_BELT_SEPARATION
        elif cur.body_type == 'a':
            min_gap = prev.min_orbit_distance
        else:
            # `cur.star.mass` (not `system.star.mass`): for an S-type (wide)
            # binary's secondary list, `system.star` is still the primary --
            # see `StarSystem._mutual_min_distance_au`'s own docstring for
            # why the owning planet's own `.star` is used instead.
            mutual_radius_m = mutual_hill_radius_m(
                cur.distance * physical_constants.AU_TO_M, prev.distance * physical_constants.AU_TO_M,
                cur.mass, prev.mass, cur.star.mass,
            )
            min_gap = (mutual_radius_m / physical_constants.AU_TO_M) * prog_c.MUTUAL_HILL_RADII_SEPARATION

        # A relative tolerance, not a fixed 1e-9 AU: `min_gap` here and the
        # value `validate_system` actually enforced are the same quantity
        # computed via two different (mathematically but not
        # floating-point-identically equivalent) arithmetic paths, and a
        # system that pushes a body out to hundreds of thousands of AU
        # (an M50 star's own `MAX_PLANETS` case, e.g.) needs a tolerance
        # that scales with that magnitude, not a fixed absolute slop sized
        # for AU-scale gaps.
        tolerance = max(1e-9, abs(min_gap) * 1e-9)
        assert gap >= min_gap - tolerance, (
            f"objects {i - 1} ({prev.body_type}) and {i} ({cur.body_type}) overlap: gap={gap}, required>={min_gap}"
        )


def assert_no_orbital_overlap(system):
    """
    Mirrors StarSystem.validate_system's own separation case matrix, as a
    check that its correction pass actually converges rather than as
    independent new physics -- the minimums it enforces are its own.

    For an S-type (wide) binary, `system.planets`/`system.secondary_planets`
    are two independent orbital sequences around two different stars --
    each is checked for internal overlap on its own; there is no ordering
    relationship between an object in one list and an object in the other
    for this check to apply to (their cross-star clearance is a separate
    concern, see `StarSystem._validate_cross_star_clearance`).
    """
    _assert_no_overlap_within(system.planets)
    _assert_no_overlap_within(system.secondary_planets)


def _all_planets_and_moons(system):
    bodies = []
    for obj in system.planets + system.secondary_planets:
        if obj.body_type == "a":
            continue
        bodies.append(obj)
        bodies.extend(obj.moons)
    return bodies


def assert_positions_and_speeds_are_consistent(system):
    """
    Every planet's/moon's `position_x/y/z` must sit exactly on the sphere
    of radius `distance` (the point `orbital_position_au` derives is, by
    construction, at that fixed radius from the orbital anchor regardless
    of inclination/node/phase), `orbital_speed_kms` must match
    `circular_orbital_speed_kms(distance, period)`, and
    `min_update_interval_years` must match
    `minimum_update_interval_years(period)` -- a guard against the same
    "recomputed value silently drifts from a corrected distance"
    staleness bug `StarSystem.validate_system`'s own docstring warns
    `period` used to have, now extended to position/speed/update-guard
    (see `planetPhysics.update_orbital_position`).
    """
    for body in _all_planets_and_moons(system):
        radius = math.sqrt(body.position_x ** 2 + body.position_y ** 2 + body.position_z ** 2)
        assert radius == pytest.approx(body.distance, rel=1e-9), (
            f"{body.name}: position radius {radius} AU != distance {body.distance} AU"
        )
        expected_position = orbital_position_au(
            body.distance, body.orbital_inclination_deg,
            body.orbital_ascending_node_deg, body.orbital_phase_deg,
        )
        assert (body.position_x, body.position_y, body.position_z) == pytest.approx(expected_position)

        expected_speed = circular_orbital_speed_kms(body.distance, body.period)
        assert body.orbital_speed_kms == pytest.approx(expected_speed, rel=1e-9), (
            f"{body.name}: orbital_speed_kms {body.orbital_speed_kms} != expected {expected_speed}"
        )

        expected_interval = minimum_update_interval_years(body.period)
        assert body.min_update_interval_years == pytest.approx(expected_interval, rel=1e-9), (
            f"{body.name}: min_update_interval_years {body.min_update_interval_years} "
            f"!= expected {expected_interval}"
        )


def assert_counts_are_consistent(system):
    planet_count, belt_count, moon_count = system.count_objects()
    assert planet_count == system.planet_count
    assert belt_count == system.belt_count
    assert moon_count == system.moon_count
    # count_objects() defaults to both stars' combined lists (see its own
    # docstring) -- for a single star or P-type binary,
    # system.secondary_planets is always [], so this is unchanged from
    # `len(system.planets)` alone in those cases.
    assert planet_count + belt_count == len(system.planets) + len(system.secondary_planets)
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
        assert_positions_and_speeds_are_consistent(system)


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
        assert_positions_and_speeds_are_consistent(system)

        if attr == "HABITABLE_WORLD" and value is True:
            assert system.hab_count >= 1, f"{star_type}: HABITABLE_WORLD=True produced no habitable world"
        if attr == "ASTEROID_BELT" and value is True:
            assert system.belt_count >= 1, f"{star_type}: ASTEROID_BELT=True produced no belt"
        if attr == "PLANETS" and value is False:
            assert len(system.planets) == 0
            assert len(system.secondary_planets) == 0
        if attr == "MOONS" and value is False:
            assert system.moon_count == 0
        if attr == "BINARY_SYSTEM" and value is True:
            # BINARY_SYSTEM=True with WIDE_BINARY left at its default (None)
            # picks either binary configuration at random (see
            # SystemConfig.WIDE_BINARY) -- assertions below cover whichever
            # one this trial happened to roll; test_binary_system_star_properties_are_sane
            # and test_wide_binary.py's own tests force each configuration
            # explicitly for configuration-specific checks.
            assert len(system.stars) == 2
            assert all(s.mass > 0 for s in system.stars)
            if isinstance(system.star, BinaryStarProxy):
                assert system.binary_type == "close"
                assert math.isclose(system.star.mass, sum(s.mass for s in system.stars))
            else:
                assert system.binary_type == "wide"
                assert system.star is system.primary_star
                assert system.primary_star.a_crit_au is not None
                assert system.secondary_star.a_crit_au is not None


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


def test_render_is_idempotent_and_does_not_double_roll_flavor(monkeypatch):
    """
    Regression test for the Phase 0 bug (see TODO.md): both
    Planet._generate_life_and_flavor_paragraphs and StarSystem.__str__ used
    to roll flavor text and mutate system_config.system_flavor_count/
    recent_flavor_texts at render time (inside to_paragraph_list()/__str__())
    rather than at generation time, so calling either twice on the same
    object double-rolled and double-mutated shared counters. Flavor text is
    now decided once during StarSystem.__init__ (system-level in
    self.system_flavor_text, planet/moon-level via
    planetLife.decide_flavor_text), and __str__/to_paragraph_list() are pure
    reads.
    """
    monkeypatch.setattr(prog_c, "FLAVOR_CHANCE_SYSTEM", 1.0)
    monkeypatch.setattr(prog_c, "FLAVOR_CHANCE_PLANET", 1.0)

    system = StarSystem(system_config=make_config("G2V", PLANETS=True, MAX_PLANETS=True, HABITABLE_WORLD=True))

    assert system.system_flavor_text is not None
    assert system.system_config.system_flavor_count > 0

    flavor_count_after_generation = system.system_config.system_flavor_count
    recent_flavor_after_generation = list(system.system_config.recent_flavor_texts)

    first_render = str(system)
    second_render = str(system)

    assert first_render == second_render
    assert system.system_config.system_flavor_count == flavor_count_after_generation
    assert system.system_config.recent_flavor_texts == recent_flavor_after_generation


def test_binary_system_star_properties_are_sane():
    # WIDE_BINARY=False forces the P-type (close) configuration this test
    # specifically checks -- see test_wide_binary.py for the S-type
    # equivalent (WideBinaryPair/a_crit-based) checks.
    for star_type in ["G2V", "M2VII", "O5V"]:
        for _ in range(TRIALS):
            system = StarSystem(system_config=make_config(star_type, BINARY_SYSTEM=True, WIDE_BINARY=False))
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
            assert proxy.galactic_orbital_speed_kms > 0 and math.isfinite(proxy.galactic_orbital_speed_kms)
            assert proxy.galactic_orbital_period_gy > 0 and math.isfinite(proxy.galactic_orbital_period_gy)
            # Galactic orbit doesn't depend on mass -- the proxy's combined
            # value should match either constituent star's own value.
            assert proxy.galactic_orbital_speed_kms == pytest.approx(primary.galactic_orbital_speed_kms)
            assert proxy.galactic_orbital_period_gy == pytest.approx(primary.galactic_orbital_period_gy)
