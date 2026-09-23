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

from stellarObjects.asteroidData import AsteroidBelt
from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem
from stellarObjects.doubleStar import BinaryStarProxy
from stellarObjects import physical_constants, program_constants as prog_c
from stellarObjects.utils import (circular_orbital_speed_kms, minimum_update_interval_years,
                                   mutual_hill_radius_au, mutual_hill_radius_m, orbital_position_au)

# One representative star type per Yerkes class, spanning several spectral
# letters, so the system-generation sweep exercises every evolutionary track
# (main sequence, giant, subgiant, bright giant, supergiant, hypergiant,
# subdwarf, white dwarf) without re-running the full 770-combo star matrix.
STAR_TYPES = [
    "G2V", "M5V", "K3V", "A0V", "B2V", "O5V",
    "K1III", "F0IB", "M2VII", "K5IV", "A0II", "M3IA+", "M50", "K8VI",
]

TRISTATE_ATTRS = [
    "HABITABLE_WORLD", "ASTEROID_BELT", "COMETS", "LARGE_STAR", "MOONS",
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


def assert_no_belt_overlap_all_pairs(system):
    """
    Independent, all-pairs safety net -- deliberately NOT a re-derivation
    of `validate_system`'s own spacing formula (that's what
    `assert_no_orbital_overlap`/`_assert_no_overlap_within` already is, an
    adjacent-list-pairs mirror of it). This instead asserts the one
    invariant that must hold regardless of *how* objects got spaced: no
    asteroid belt's own `[lower_limit, upper_limit]` span may overlap
    another belt's span, or contain a planet's `distance`, for ANY pair in
    a list -- not just list-adjacent ones.

    This is the regression guard for a fixed bug where a zone-forced
    placement (a forced `HABITABLE_WORLD`, or an explicit `SLOTS` class
    request) could land a planet or belt at a distance the sequential
    generation loop's own `estimated_distance` progression never implied,
    breaking the "list order == distance order" assumption
    `validate_system`'s adjacent-pairs sweep depends on -- so two belts
    could overlap, or a planet's orbit could sit inside a belt, without
    either being caught, because the object that would have caught it
    wasn't actually its list neighbor. Fixed via belt-aware distance
    selection (`StarSystem._distance_avoiding_belts`) and cross-star belt
    clearance (`StarSystem._validate_cross_star_clearance`).

    Deliberately does NOT compare across `system.planets` and
    `system.secondary_planets` for a wide (S-type) binary: each list's
    `distance` is measured from its OWN star, so the two lists don't share
    a coordinate frame at all -- a primary-list belt at "3 AU" and a
    secondary-list planet at "3 AU" are three AU from two different stars,
    typically separated by tens of AU of binary orbital separation, not
    anywhere near each other in real space. That cross-star relationship
    has its own physically-meaningful check
    (`StarSystem._validate_cross_star_clearance`'s worst-case-gap-vs-
    binary-separation reasoning), which is exercised separately below.
    """
    for planet_list in (system.planets, system.secondary_planets):
        belts = [obj for obj in planet_list if obj.body_type == 'a']

        for i, belt in enumerate(belts):
            for other in belts[i + 1:]:
                overlap = belt.lower_limit < other.upper_limit and other.lower_limit < belt.upper_limit
                assert not overlap, (
                    f"asteroid belts overlap: [{belt.lower_limit}, {belt.upper_limit}] "
                    f"and [{other.lower_limit}, {other.upper_limit}]"
                )
            for planet in planet_list:
                if planet.body_type == 'a':
                    continue
                inside = belt.lower_limit < planet.distance < belt.upper_limit
                assert not inside, (
                    f"{planet.name}'s orbit at {planet.distance} AU sits inside asteroid belt "
                    f"[{belt.lower_limit}, {belt.upper_limit}]"
                )


def assert_cross_star_clearance_holds(system):
    """
    For a wide (S-type) binary, mirrors `StarSystem._validate_cross_star_
    clearance`'s own worst-case-gap-vs-threshold formula (the physically
    correct way to compare the two stars' outermost objects -- each list's
    `distance` is measured from its own star, so a raw interval-overlap
    check across the two lists, the way `assert_no_belt_overlap_all_pairs`
    checks within one list, would be meaningless -- see that function's
    own docstring) to confirm the correction actually converged, the same
    "mirror the formula as a convergence check" role
    `assert_no_orbital_overlap` plays for `validate_system` itself.

    Regression guard for `_validate_cross_star_clearance` unconditionally
    skipping either star's trailing `AsteroidBelt` when finding its
    "outermost planet" -- this now finds the true outermost object (belt
    included) on each side and applies the same belt-aware threshold
    (`MIN_ASTEROID_BELT_SEPARATION` when either side is a belt, the
    Gladman mutual-Hill-radius criterion when both are real planets) the
    fixed method itself uses.
    """
    if system.binary_type != "wide" or not system.planets or not system.secondary_planets:
        return

    def edge_au(obj):
        return obj.upper_limit if obj.body_type == 'a' else obj.distance

    outer_p = max(system.planets, key=edge_au)
    outer_s = max(system.secondary_planets, key=edge_au)

    a_bin = system.wide_binary.separation_au
    worst_case_gap_au = a_bin - edge_au(outer_p) - edge_au(outer_s)

    if outer_p.body_type == 'a' or outer_s.body_type == 'a':
        threshold_au = prog_c.MIN_ASTEROID_BELT_SEPARATION
    else:
        central_mass_kg = system.primary_star.mass + system.secondary_star.mass
        r_h_mutual_au = mutual_hill_radius_au(
            outer_p.mass, outer_s.mass, outer_p.distance, outer_s.distance, central_mass_kg,
        )
        threshold_au = physical_constants.GLADMAN_MUTUAL_HILL_STABILITY_FACTOR * r_h_mutual_au

    tolerance = max(1e-9, abs(threshold_au) * 1e-9)
    assert worst_case_gap_au >= threshold_au - tolerance, (
        f"wide binary's two disks encroach: outermost primary object edge={edge_au(outer_p)} AU, "
        f"outermost secondary object edge={edge_au(outer_s)} AU, a_bin={a_bin} AU, "
        f"worst-case gap={worst_case_gap_au}, required>={threshold_au}"
    )


def test_validate_system_belt_overlap_correction_lands_exactly_past_the_belt():
    """
    Regression test for a bug where `validate_system`'s
    `additional_correction` term double-counted `last_planet.distance` on
    top of the offset that already cancels the negative gap, roughly
    doubling the corrected distance for any overlap correction involving
    an asteroid belt instead of nudging the overlapping body just past
    it. `_assert_no_overlap_within` (this file's own general invariant
    check) can't catch this on its own -- it only asserts the resulting
    gap is AT LEAST the required minimum, which a large overshoot still
    (trivially) satisfies -- so this asserts the exact corrected value.
    """
    system = StarSystem(system_config=make_config("G2V", PLANETS=False))
    belt = AsteroidBelt(system.system_config, distance=2.0, lower_limit=1.8, upper_limit=2.2)
    overlapping_belt = AsteroidBelt(system.system_config, distance=2.0, lower_limit=1.8, upper_limit=2.2)

    system.validate_system([belt, overlapping_belt])

    expected_distance = belt.upper_limit + prog_c.MIN_ASTEROID_BELT_SEPARATION
    assert overlapping_belt.distance == pytest.approx(expected_distance, rel=1e-9)


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
    # comets is its own list, separate from planets/secondary_planets (see
    # StarSystem._generate_comets' docstring) -- count_comets() defaults
    # to both stars' combined lists the same way count_objects() does.
    comet_count = system.count_comets()
    assert comet_count == system.comet_count
    assert comet_count == len(system.comets) + len(system.secondary_comets)


@pytest.mark.parametrize("star_type", STAR_TYPES)
def test_baseline_random_system_generates_without_error(star_type):
    for _ in range(TRIALS):
        system = StarSystem(system_config=make_config(star_type))
        assert str(system).strip()
        assert_ages_never_exceed_lifespan(system)
        assert_no_orbital_overlap(system)
        assert_no_belt_overlap_all_pairs(system)
        assert_cross_star_clearance_holds(system)
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
        assert_no_belt_overlap_all_pairs(system)
        assert_cross_star_clearance_holds(system)
        assert_counts_are_consistent(system)
        assert_positions_and_speeds_are_consistent(system)

        if attr == "HABITABLE_WORLD" and value is True:
            assert system.hab_count >= 1, f"{star_type}: HABITABLE_WORLD=True produced no habitable world"
        if attr == "ASTEROID_BELT" and value is True:
            assert system.belt_count >= 1, f"{star_type}: ASTEROID_BELT=True produced no belt"
        if attr == "COMETS" and value is True:
            assert system.comet_count >= 1, f"{star_type}: COMETS=True produced no comet"
        if attr == "COMETS" and value is False:
            assert system.comet_count == 0, f"{star_type}: COMETS=False still produced a comet"
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
        assert_no_belt_overlap_all_pairs(system)
    assert failures_belt == 0, f"{star_type}: ASTEROID_BELT+HABITABLE_WORLD failed to produce a belt in {failures_belt}/{trials} trials"
    assert failures_hab == 0, f"{star_type}: ASTEROID_BELT+HABITABLE_WORLD failed to produce a habitable world in {failures_hab}/{trials} trials"


def test_wide_binary_asteroid_belts_forced_do_not_overlap_across_stars():
    """
    Regression guard for `_validate_cross_star_clearance` unconditionally
    skipping any star's trailing `AsteroidBelt` when finding its
    "outermost planet" -- meaning a wide (S-type) binary's two disks could
    have belts (or a belt and the other star's outermost planet) overlap
    in real space with no check catching it at all. Forces ASTEROID_BELT
    on both stars' generation (primary via `apply_guarantees`, secondary
    always independently random -- so this also forces `MAX_PLANETS` to
    push both disks as far out as possible, maximizing how often the two
    disks would actually reach each other and exercise the cross-star
    check) and asserts the fixed all-pairs invariant across both stars'
    combined object lists.
    """
    trials = 15
    for _ in range(trials):
        system = StarSystem(system_config=make_config(
            "M2VII", BINARY_SYSTEM=True, WIDE_BINARY=True, ASTEROID_BELT=True, MAX_PLANETS=True,
        ))
        assert system.binary_type == "wide"
        assert_no_belt_overlap_all_pairs(system)
        assert_cross_star_clearance_holds(system)


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


def test_binary_system_markdown_renders_each_stars_table():
    """
    Regression guard for a fixed bug: each binary star's own `###` header
    was joined to its property table by a single '\\n', not a blank line,
    so html/lib/mdconvert.py's blank-line block splitter lumped the header
    and table into one block -- which is neither a valid single-line
    heading nor a valid table -- and rendered as one escaped, literal
    paragraph of '#'/'|' text instead of a real <h3> + <table>. Checks both
    binary configurations (close/P-type and wide/S-type), since the bug was
    duplicated in both code paths (systemData.py's close- and wide-binary
    branches each had their own copy of the same one-newline join).
    """
    import html as html_module
    import os
    import sys

    _src_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, os.path.join(_src_dir, "html", "lib"))
    from mdconvert import markdown_to_html  # same path setup as test_mdconvert.py

    for wide_binary in (False, True):
        system = StarSystem(system_config=make_config("G2V", BINARY_SYSTEM=True, WIDE_BINARY=wide_binary))
        # _db.py's persist_star_system renders the web-facing copy with
        # MARKDOWN=True (ATX headers) -- match that exactly, since MARKDOWN
        # defaults to False (wikitext '===' headers) which mdconvert.py
        # doesn't parse as headings at all.
        system.system_config.MARKDOWN = True
        html = markdown_to_html(str(system))

        star_names = [s.name for s in system.stars]
        assert len(star_names) == 2
        table_count = html.count("<table>")
        assert table_count >= 2, (
            f"expected at least one rendered <table> per star (got {table_count} "
            f"for {star_names}); a binary star's header+table likely collapsed "
            f"into one unrendered paragraph again"
        )
        for name in star_names:
            # mdconvert.py HTML-escapes every heading's text (see its module
            # docstring) -- a name containing '/&/</> (e.g. from an
            # apostrophe-bearing generated name) is expected to come out
            # escaped, not literal.
            escaped_name = html_module.escape(name)
            assert f">{escaped_name}<" in html, f"{name}'s header did not render as a heading"
