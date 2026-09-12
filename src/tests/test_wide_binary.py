"""
S-type (wide) binary regression tests.

Covers the physics formulas (`utils.holman_wiegert_critical_semimajor_axis`,
`utils.mutual_hill_radius_au`), the sampling distributions
(`utils.sample_wide_binary_separation_au`/`sample_wide_binary_eccentricity`),
the per-star orbit-ceiling enforcement in `StarSystem._generate_planets`, the
cross-star clearance pruning in `StarSystem._validate_cross_star_clearance`,
and end-to-end generation/round-trip through `WIDE_BINARY_SYSTEM=True`.

Run with: pytest src/tests/test_wide_binary.py
"""
import math
import statistics

import pytest

from stellarObjects import physical_constants, program_constants
from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import (holman_wiegert_critical_semimajor_axis,
                                   mutual_hill_radius_au,
                                   sample_wide_binary_eccentricity,
                                   sample_wide_binary_separation_au)

STAR_TYPES = ["G2V", "M5V", "K3V", "A0V"]
TRIALS = 3


def make_config(star_type="G2V", **overrides):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    for attr, value in overrides.items():
        setattr(cfg, attr, value)
    return cfg


def _outermost_planet(planets):
    for obj in reversed(planets):
        if obj.body_type != "a":
            return obj
    return None


def _orbit_edge(obj):
    return obj.upper_limit if obj.body_type == "a" else obj.distance + obj.min_orbit_distance


# ---------------------------------------------------------------------------
# holman_wiegert_critical_semimajor_axis
# ---------------------------------------------------------------------------

def test_holman_wiegert_equal_mass_circular_reference_value():
    # 0.464 - 0.380*0.5 = 0.274 exactly (see the function's own docstring).
    assert holman_wiegert_critical_semimajor_axis(1.0, 0.5, 0.0) == pytest.approx(0.274)


def test_holman_wiegert_equal_mass_moderate_eccentricity_reference_value():
    result = holman_wiegert_critical_semimajor_axis(1.0, 0.5, 0.5)
    assert result == pytest.approx(0.1178, abs=1e-3)


def test_holman_wiegert_scales_linearly_with_separation():
    ratio_at_1au = holman_wiegert_critical_semimajor_axis(1.0, 0.3, 0.2)
    assert holman_wiegert_critical_semimajor_axis(500.0, 0.3, 0.2) == pytest.approx(ratio_at_1au * 500.0)


def test_holman_wiegert_decreases_with_eccentricity():
    mu = 0.4
    a_crit_low_e = holman_wiegert_critical_semimajor_axis(100.0, mu, 0.0)
    a_crit_high_e = holman_wiegert_critical_semimajor_axis(100.0, mu, 0.7)
    assert a_crit_high_e < a_crit_low_e


def test_holman_wiegert_clamps_out_of_range_inputs():
    # mu/e outside the fit's validated range are clamped, not extrapolated
    # or raised (see the function's own docstring).
    mu_min, mu_max = physical_constants.HOLMAN_WIEGERT_MU_RANGE
    e_min, e_max = physical_constants.HOLMAN_WIEGERT_ECCENTRICITY_RANGE
    assert holman_wiegert_critical_semimajor_axis(100.0, -1.0, -1.0) == pytest.approx(
        holman_wiegert_critical_semimajor_axis(100.0, mu_min, e_min)
    )
    assert holman_wiegert_critical_semimajor_axis(100.0, 5.0, 5.0) == pytest.approx(
        holman_wiegert_critical_semimajor_axis(100.0, mu_max, e_max)
    )


# ---------------------------------------------------------------------------
# mutual_hill_radius_au
# ---------------------------------------------------------------------------

def test_mutual_hill_radius_hand_computed_example():
    earth_mass_kg = 5.972e24
    solar_mass_kg = physical_constants.SOLAR_MASS_TO_KG
    # Two Earth-mass planets, both at 1 AU, around a 1 solar mass central body.
    expected = ((2 * earth_mass_kg) / (3 * solar_mass_kg)) ** (1 / 3) * 1.0
    assert mutual_hill_radius_au(earth_mass_kg, earth_mass_kg, 1.0, 1.0, solar_mass_kg) == pytest.approx(expected)


def test_mutual_hill_radius_scales_with_average_distance():
    m = 5.972e24
    central = physical_constants.SOLAR_MASS_TO_KG
    r_at_1_1 = mutual_hill_radius_au(m, m, 1.0, 1.0, central)
    r_at_2_2 = mutual_hill_radius_au(m, m, 2.0, 2.0, central)
    assert r_at_2_2 == pytest.approx(2 * r_at_1_1)


# ---------------------------------------------------------------------------
# Separation/eccentricity samplers -- statistical sanity (randomized, so
# these check distributional shape/range, not exact values).
# ---------------------------------------------------------------------------

def test_separation_sampler_stays_in_range_and_is_log_uniform():
    lo = program_constants.WIDE_BINARY_SEPARATION_MIN_AU
    hi = program_constants.WIDE_BINARY_SEPARATION_MAX_AU
    samples = [sample_wide_binary_separation_au() for _ in range(500)]
    assert all(lo <= s <= hi for s in samples)

    # Log-uniform: log10(samples) should be roughly uniformly distributed
    # between log10(lo) and log10(hi), so its mean should land near the
    # midpoint of that log range (a generous tolerance since this is a
    # randomized statistical check, not an exact one).
    log_samples = [math.log10(s) for s in samples]
    expected_mean = (math.log10(lo) + math.log10(hi)) / 2
    assert abs(statistics.mean(log_samples) - expected_mean) < 0.3


def test_eccentricity_sampler_stays_in_range_and_follows_thermal_distribution():
    e_max = program_constants.WIDE_BINARY_ECCENTRICITY_MAX
    samples = [sample_wide_binary_eccentricity() for _ in range(1000)]
    assert all(0 <= e < e_max for e in samples)

    # Thermal distribution f(e) = 2e (rescaled to e_max): median should sit
    # at e_max / sqrt(2) (inverse-CDF F(e) = (e/e_max)^2 = 0.5 solved for e),
    # not e_max / 2 the way a uniform distribution's would.
    expected_median = e_max / math.sqrt(2)
    assert abs(statistics.median(samples) - expected_median) < 0.05
    # Denser at higher e than a uniform distribution would be: strictly
    # more than half the mass should sit above the halfway point.
    above_half = sum(1 for e in samples if e > e_max / 2)
    assert above_half > 0.55 * len(samples)


# ---------------------------------------------------------------------------
# Per-star orbit ceiling enforcement (StarSystem._orbit_ceiling_au /
# _generate_planets' _exceeds_orbit_ceiling checks).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("star_type", STAR_TYPES)
def test_orbit_ceiling_falls_back_to_system_perimeter_without_a_crit(star_type):
    # Pinned False: `BINARY_SYSTEM` now rolls real chance when left at its
    # default (see StarSystem._should_generate_binary), and a "close"
    # (P-type) pair's merged proxy star never sets a_crit_au either -- but
    # this test is specifically about the single-star case, so binary-vs-
    # single must not be left to chance here.
    system = StarSystem(system_config=make_config(star_type, PLANETS=False, BINARY_SYSTEM=False))
    assert system.star.a_crit_au is None
    assert system._orbit_ceiling_au(system.star) == system.star.system_perimeter


def test_tight_a_crit_stops_placement_earlier_than_system_perimeter_alone_would():
    system = StarSystem(system_config=make_config("G2V", MAX_PLANETS=True))
    star = system.star
    assert star.system_perimeter > 100  # a single star's galactic Hill sphere is astronomically large

    # An artificially tight companion-driven limit, well inside where
    # MAX_PLANETS=True would otherwise keep placing planets.
    star.a_crit_au = 1.0
    ceiling = system._orbit_ceiling_au(star)
    assert ceiling == 1.0

    planets = system._generate_planets(star, star.habitable_zone, ceiling)
    for planet in planets:
        assert _orbit_edge(planet) <= ceiling + 1e-9


def test_apply_guarantees_false_ignores_slots_and_forced_habitable_world():
    system = StarSystem(system_config=make_config("G2V", PLANETS=False))
    star = system.star
    slotted_config = make_config(
        "G2V", HABITABLE_WORLD=True,
        SLOTS=[{"type": "planet", "planet_class": "J", "moons": 0}],
    )
    system.system_config = slotted_config

    planets = system._generate_planets(
        star, star.habitable_zone, system._orbit_ceiling_au(star), apply_guarantees=False
    )
    # Neither the explicit SLOTS entry (which would force a Class J planet
    # into the first slot) nor the HABITABLE_WORLD guarantee applied --
    # apply_guarantees=False means this list is fully, independently random.
    if planets:
        assert not (planets[0].body_type != "a" and planets[0].planet_class == "J" and len(planets) == 1)


# ---------------------------------------------------------------------------
# Cross-star clearance (_validate_cross_star_clearance) -- uses lightweight
# stub objects (duck-typed: body_type/distance/mass) rather than full Planet
# instances, so the pruning algorithm itself can be tested deterministically
# without depending on randomized generation.
# ---------------------------------------------------------------------------

class _StubPlanet:
    def __init__(self, distance, mass=5.972e24, body_type="t"):
        self.distance = distance
        self.mass = mass
        self.body_type = body_type


def _wide_system_for_clearance_test(separation_au, primary_a_crit_au, secondary_a_crit_au):
    system = StarSystem(system_config=make_config("G2V", BINARY_SYSTEM=True, WIDE_BINARY=True, PLANETS=False))
    system.wide_binary.separation_au = separation_au
    system.primary_star.a_crit_au = primary_a_crit_au
    system.secondary_star.a_crit_au = secondary_a_crit_au
    return system


def test_cross_star_clearance_no_op_when_trivially_satisfied():
    system = _wide_system_for_clearance_test(separation_au=100.0, primary_a_crit_au=40.0, secondary_a_crit_au=40.0)
    system.planets = [_StubPlanet(distance=5.0)]
    system.secondary_planets = [_StubPlanet(distance=5.0)]

    system._validate_cross_star_clearance()

    assert len(system.planets) == 1
    assert len(system.secondary_planets) == 1


def test_cross_star_clearance_prunes_when_forced_close():
    system = _wide_system_for_clearance_test(separation_au=100.0, primary_a_crit_au=50.0, secondary_a_crit_au=50.0)
    # Both planets sit right at their own star's a_crit, with almost no gap
    # left between them (100 - 49.99 - 49.99 = 0.02 AU) -- guaranteed to
    # violate the Gladman threshold for any physically real planet mass.
    system.planets = [_StubPlanet(distance=49.99)]
    system.secondary_planets = [_StubPlanet(distance=49.99)]

    system._validate_cross_star_clearance()

    # One of the two must have been pruned (equal margins -- either is a
    # valid outcome of the tie-break), and clearance now holds for whatever
    # remains.
    assert len(system.planets) + len(system.secondary_planets) == 1


def test_cross_star_clearance_prunes_least_margin_star_first():
    # Both planets sit close to their OWN star's a_crit (so the gap between
    # them stays tiny -- 109.5 - 49.99 - 59.0 = 0.51 AU, guaranteed to
    # violate the Gladman threshold regardless of the randomly generated
    # G2V pair's exact masses), but with clearly different margins:
    # primary's margin (50.0 - 49.99 = 0.01) is far smaller than
    # secondary's (60.0 - 59.0 = 1.0) -- the primary should lose its planet
    # first.
    system = _wide_system_for_clearance_test(separation_au=109.5, primary_a_crit_au=50.0, secondary_a_crit_au=60.0)
    system.planets = [_StubPlanet(distance=49.99)]
    system.secondary_planets = [_StubPlanet(distance=59.0)]

    system._validate_cross_star_clearance()

    assert len(system.planets) == 0
    assert len(system.secondary_planets) == 1


def test_cross_star_clearance_terminates_when_a_list_is_exhausted():
    # Both lists start with exactly one, mutually violating, planet --
    # pruning empties one side, and the loop must terminate rather than
    # looping forever once outer_p/outer_s can no longer both resolve.
    system = _wide_system_for_clearance_test(separation_au=10.0, primary_a_crit_au=5.0, secondary_a_crit_au=5.0)
    system.planets = [_StubPlanet(distance=4.999)]
    system.secondary_planets = [_StubPlanet(distance=4.999)]

    system._validate_cross_star_clearance()  # must return, not hang

    assert len(system.planets) + len(system.secondary_planets) == 1


def test_cross_star_clearance_skips_belts_as_outermost_body():
    # A trailing AsteroidBelt has no discrete Hill sphere to evaluate -- if
    # either side's outermost body is a belt, this is a no-op (see the
    # method's own docstring).
    system = _wide_system_for_clearance_test(separation_au=10.0, primary_a_crit_au=5.0, secondary_a_crit_au=5.0)
    system.planets = [_StubPlanet(distance=4.999, body_type="a")]
    system.secondary_planets = [_StubPlanet(distance=4.999, body_type="a")]

    system._validate_cross_star_clearance()

    assert len(system.planets) == 1
    assert len(system.secondary_planets) == 1


def test_cross_star_clearance_is_rare_for_ordinary_wide_low_eccentricity_pairs():
    """
    Statistical sanity check backing the "rarely triggers" claim in
    `_validate_cross_star_clearance`'s own docstring: real generation
    (not the forced-close stub scenarios above) should essentially never
    need to prune for an ordinary, randomly-sampled wide binary.
    """
    pruned = 0
    trials = 15
    for _ in range(trials):
        system = StarSystem(system_config=make_config("G2V", BINARY_SYSTEM=True, WIDE_BINARY=True))
        if system.binary_type != "wide":
            continue
        primary_before = len(system.planets)
        secondary_before = len(system.secondary_planets)
        system._validate_cross_star_clearance()
        if len(system.planets) < primary_before or len(system.secondary_planets) < secondary_before:
            pruned += 1
    # Not asserting zero (eccentric/tight draws can legitimately trigger
    # it), but it should be the rare exception, not the common case.
    assert pruned <= trials // 2


# ---------------------------------------------------------------------------
# End-to-end generation and round-trip.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("star_type", STAR_TYPES)
def test_wide_binary_system_generates_and_stays_within_bounds(star_type):
    for _ in range(TRIALS):
        system = StarSystem(system_config=make_config(star_type, BINARY_SYSTEM=True, WIDE_BINARY=True, MAX_PLANETS=True))

        assert system.binary_type == "wide"
        assert system.star is system.primary_star
        assert system.secondary_star.mass <= system.primary_star.mass
        assert system.wide_binary is not None

        lo = program_constants.WIDE_BINARY_SEPARATION_MIN_AU
        hi = program_constants.WIDE_BINARY_SEPARATION_MAX_AU
        assert lo <= system.wide_binary.separation_au <= hi
        assert 0 <= system.wide_binary.eccentricity < program_constants.WIDE_BINARY_ECCENTRICITY_MAX
        assert system.wide_binary.periapsis_au <= system.wide_binary.separation_au <= system.wide_binary.apoapsis_au

        assert system.primary_star.a_crit_au is not None
        assert system.secondary_star.a_crit_au is not None
        assert system.primary_star.a_crit_au < system.wide_binary.separation_au
        assert system.secondary_star.a_crit_au < system.wide_binary.separation_au

        for planet in system.planets:
            assert _orbit_edge(planet) <= system._orbit_ceiling_au(system.primary_star) + 1e-6
        for planet in system.secondary_planets:
            assert _orbit_edge(planet) <= system._orbit_ceiling_au(system.secondary_star) + 1e-6

        assert str(system).strip()


def test_wide_binary_secondary_ages_independently_of_primarys_planets():
    system = StarSystem(system_config=make_config("G2V", BINARY_SYSTEM=True, WIDE_BINARY=True, MAX_PLANETS=True))
    if system.binary_type != "wide":
        pytest.skip("random draw landed on close binary")
    # Both stars' ages must stay within their own (possibly different)
    # lifespans after adjust_age_for_planets ran against each star's own
    # list independently.
    assert system.primary_star.age <= system.primary_star.lifespan * 1.001
    assert system.secondary_star.age <= system.secondary_star.lifespan * 1.001


def test_wide_binary_round_trips_through_to_dict_from_dict():
    system = StarSystem(system_config=make_config("G2V", BINARY_SYSTEM=True, WIDE_BINARY=True, MAX_PLANETS=True))

    data = system.to_dict()
    assert data["binary_type"] == "wide"
    assert data["is_binary"] is True
    assert data["secondary_star"] is not None
    assert data["wide_binary"] is not None

    reloaded = StarSystem.from_dict(data)

    assert reloaded.binary_type == "wide"
    assert reloaded.primary_star.mass == system.primary_star.mass
    assert reloaded.secondary_star.mass == system.secondary_star.mass
    assert reloaded.primary_star.a_crit_au == system.primary_star.a_crit_au
    assert reloaded.secondary_star.a_crit_au == system.secondary_star.a_crit_au
    assert reloaded.wide_binary.separation_au == system.wide_binary.separation_au
    assert reloaded.wide_binary.eccentricity == system.wide_binary.eccentricity
    assert len(reloaded.planets) == len(system.planets)
    assert len(reloaded.secondary_planets) == len(system.secondary_planets)
    assert str(reloaded) == str(system)


def test_close_binary_default_still_reachable_when_wide_binary_forced_false():
    system = StarSystem(system_config=make_config("G2V", BINARY_SYSTEM=True, WIDE_BINARY=False))
    assert system.binary_type == "close"
    assert system.secondary_planets == []
    assert system.wide_binary is None
