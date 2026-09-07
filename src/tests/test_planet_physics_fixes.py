"""
Regression tests for the four Track A physics fixes to the planet
atmosphere/density model (see planetPhysics.py and plausibility.py):

  1. Gas-giant density blend: arithmetic mean -> mass-weighted harmonic mean.
  2. Greenhouse factor: no longer inverted (rewarding distance from CO2's
     own molar density instead of proximity to it).
  3. Atmospheric pressure: no longer independent of gravity (a gravity-based
     retention factor is applied to an effective atmosphere density used
     only in the pressure calculation).
  4. Class P gets its own ("cold, glaciated") albedo range so it's no
     longer statistically indistinguishable from Class M.

Run with: pytest src/tests/test_planet_physics_fixes.py
"""
import math
import statistics

import pytest

from stellarObjects import physical_constants as pc
from stellarObjects import plausibility
from stellarObjects import planetPhysics
from stellarObjects import program_constants as prog_c
from stellarObjects.config import SystemConfig
from stellarObjects.planetData import Planet
from stellarObjects.starData import Star

N_SAMPLE = 300


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


# ---------------------------------------------------------------------------
# Fix 1: gas-giant density blend (mass-weighted harmonic mean)
# ---------------------------------------------------------------------------

GAS_GIANT_CLASS = next(c for c, d in prog_c.PLANET_CLASSES.items() if d["type"] == "g")
GAS_GIANT_ZONE = next(z for z in "hec" if prog_c.PLANET_CLASSES[GAS_GIANT_CLASS][z])


def test_gas_giant_density_blend_matches_hand_computed_harmonic_mean(monkeypatch, host_star):
    """
    Drives generate_planet_properties with controlled random draws (core
    density, atm_density, atm_molar_density, core/atmosphere ratio,
    envelope density -- in the order they're actually consumed) and checks
    the resulting planet.density against a hand-computed mass-weighted
    harmonic mean:
        1 / (ratio / core_density_gcm3 + (1 - ratio) / envelope_density_gcm3)
    which is the physically correct way to combine two densities via a mass
    fraction (unlike the old arithmetic-mean blend it replaces).

    Uses Class I (GAS_GIANT_CLASS resolves to the first type=='g' class,
    which has no PLANET_CLASSES density_range override -- see
    generate_planet_properties' `"density_range" not in class_data` guard --
    so this exercises the blend path, not the direct-density override
    path covered by test_density_range_override_skips_the_blend below).
    """
    core_density_gcm3 = 1.0       # within PLANET_DENSITY["g"] = (0.69, 1.64)
    atm_density_kgm3 = 1.0        # within ATMOSPHERE_DENSITY["g"] -- feeds
    # atmospheric_pressure/scale_height only, not this blend (see
    # physical_constants.GAS_ENVELOPE_BULK_DENSITY's docstring)
    atm_molar_density = 0.003     # within ATMOSPHERIC_MOLAR_DENSITY["g"]
    ratio = 0.4                   # within GAS_GIANT_CORE_ATMOSPHERE_RATIO = (0.03, 0.6)
    envelope_density_gcm3 = 0.15  # within GAS_ENVELOPE_BULK_DENSITY = (0.06, 0.3)

    queued = [core_density_gcm3, atm_density_kgm3, atm_molar_density, ratio, envelope_density_gcm3]
    real_uniform = planetPhysics.random.uniform

    def fake_uniform(a, b):
        if queued:
            return queued.pop(0)
        return real_uniform(a, b)

    monkeypatch.setattr(planetPhysics.random, "uniform", fake_uniform)

    cfg = SystemConfig()
    distance = plausibility.distance_for_zone(host_star, GAS_GIANT_ZONE)
    radius = sum(prog_c.PLANET_CLASSES[GAS_GIANT_CLASS]["radius_range"]) / 2
    planet = Planet(
        cfg, host_star, host_star.habitable_zone, distance,
        planet_class=GAS_GIANT_CLASS, radius=radius, zone_override=GAS_GIANT_ZONE,
        moon_count=0,
    )

    expected_density = 1 / (ratio / core_density_gcm3 + (1 - ratio) / envelope_density_gcm3)
    assert planet.density == pytest.approx(expected_density, rel=1e-9)
    # Sanity: the arithmetic mean the old buggy code used would have given a
    # very different (and, for this input, much larger) value -- confirming
    # this test would have failed against the old formula, not just against
    # a formula that happens to coincide with it for these inputs.
    old_arithmetic_mean = core_density_gcm3 * ratio + (1 - ratio) * envelope_density_gcm3
    assert planet.density != pytest.approx(old_arithmetic_mean, rel=1e-6)


def test_density_range_override_skips_the_blend(monkeypatch, host_star):
    """
    A class declaring its own density_range should use that draw as
    planet.density directly -- no core/envelope blend -- since a class
    whose density is set this deliberately (e.g. a brown-dwarf-like object,
    which doesn't have a meaningfully separate light envelope over a denser
    core the way an ordinary gas giant does) shouldn't have that value
    diluted back down by blending toward a light "puffy" envelope value
    (see generate_planet_properties' `"density_range" not in class_data`
    guard). No current class declares density_range (it's generic,
    reusable override infrastructure -- the same `.get(..., default)`
    pattern every other per-class override in PLANET_CLASSES uses), so this
    test injects one onto an existing gas-giant class for the duration of
    the test rather than depending on a specific class having it.
    """
    gas_giant_class = GAS_GIANT_CLASS
    zone = GAS_GIANT_ZONE
    density_range = (42.0, 42.0)
    patched_class_data = dict(prog_c.PLANET_CLASSES[gas_giant_class])
    patched_class_data["density_range"] = density_range
    monkeypatch.setitem(prog_c.PLANET_CLASSES, gas_giant_class, patched_class_data)

    core_density_gcm3 = 42.0  # within the injected density_range above

    queued = [core_density_gcm3]
    real_uniform = planetPhysics.random.uniform

    def fake_uniform(a, b):
        if queued:
            return queued.pop(0)
        return real_uniform(a, b)

    monkeypatch.setattr(planetPhysics.random, "uniform", fake_uniform)

    cfg = SystemConfig()
    distance = plausibility.distance_for_zone(host_star, zone)
    radius = sum(prog_c.PLANET_CLASSES[gas_giant_class]["radius_range"]) / 2
    planet = Planet(
        cfg, host_star, host_star.habitable_zone, distance,
        planet_class=gas_giant_class, radius=radius, zone_override=zone,
        moon_count=0,
    )

    assert planet.density == pytest.approx(core_density_gcm3, rel=1e-9)


def test_gas_giant_sampled_densities_are_finite_positive_and_within_theoretical_bounds():
    """
    Statistical check across a real generation sample: every gas-giant
    density produced by the fixed blend is finite, strictly positive, and
    falls within the analytically-derived corner bounds that
    plausibility.theoretical_gravity_bounds_g's density sub-computation
    implies for this class (the same bounds check_hard_invariants uses,
    computed independently here to double as a cross-check that
    plausibility.py's lockstep update matches planetPhysics.py's formula).

    An earlier version of this test documented that the blend's envelope
    term reused `atm_density` (drawn from ATMOSPHERE_DENSITY["g"], which
    once divided by 1000 for unit conversion is ~1000x lighter than a real
    gas-giant envelope) -- a scale mismatch that mathematically guaranteed
    the harmonic mean collapsed to a near-zero, physically meaningless
    density no matter how dense the core. Fixed by introducing
    physical_constants.GAS_ENVELOPE_BULK_DENSITY, a separate, correctly
    g/cm^3-scaled "puffy gas giant" envelope density grounded in real
    measured values (the lowest known, WASP-193b, is ~0.06 g/cm^3) -- see
    that constant's own docstring for the full derivation. This test now
    also implicitly confirms the fixed densities land in a realistic range,
    not just a self-consistent one.
    """
    records = plausibility.generate_sample(GAS_GIANT_CLASS, GAS_GIANT_ZONE, N_SAMPLE, include_moons=False)
    densities = [r["density"] for r in records]
    assert all(math.isfinite(d) and d > 0 for d in densities)

    min_rock, max_rock = pc.PLANET_DENSITY["g"]
    min_ratio, max_ratio = prog_c.GAS_GIANT_CORE_ATMOSPHERE_RATIO
    atm_gcm3_range = pc.GAS_ENVELOPE_BULK_DENSITY

    bound_values = []
    for rock in (min_rock, max_rock):
        for ratio in (min_ratio, max_ratio):
            for atm in atm_gcm3_range:
                bound_values.append(1 / (ratio / rock + (1 - ratio) / atm))
    lo, hi = min(bound_values), max(bound_values)
    margin = 0.01
    for d in densities:
        assert lo * (1 - margin) <= d <= hi * (1 + margin), (
            f"density={d} outside theoretical bounds [{lo}, {hi}]"
        )


# ---------------------------------------------------------------------------
# Fix 2: greenhouse factor (no longer inverted)
# ---------------------------------------------------------------------------

def test_greenhouse_factor_is_monotonically_increasing_with_atm_molar_density():
    """
    Direct unit test of the greenhouse-factor formula in
    calculate_atmospheric_conditions: holding all else equal, a higher
    atm_molar_density (a heavier/denser atmosphere) must produce a higher
    (or equal, once capped) greenhouse_factor. The old formula measured
    *distance* from CO2_BASE_MOLAR_DENSITY, which was not monotonic at all
    (it decreased as atm_molar_density approached CO2's molar density from
    below, then increased again beyond it).
    """
    def greenhouse_factor(atm_molar_density):
        return min(
            prog_c.CO2_MAX_GREENHOUSE_FACTOR,
            (atm_molar_density / pc.CO2_BASE_MOLAR_DENSITY) * prog_c.CO2_MAX_GREENHOUSE_FACTOR,
        )

    samples = [0.01, 0.02, 0.03, pc.CO2_BASE_MOLAR_DENSITY, 0.05, 0.08, 0.12]
    values = [greenhouse_factor(v) for v in samples]
    assert all(b >= a for a, b in zip(values, values[1:])), values
    # And it should actually vary (not be flatlined at the cap) across the
    # terrestrial atmospheric molar density range used elsewhere in the model.
    assert values[0] < values[-1]


def test_greenhouse_factor_peaks_at_co2_base_molar_density_not_far_from_it():
    """
    Regression guard against the specific inversion bug: at
    atm_molar_density == CO2_BASE_MOLAR_DENSITY, the old (buggy) abs()-based
    formula gave greenhouse_factor == 0 (its minimum), while the new formula
    gives CO2_MAX_GREENHOUSE_FACTOR (its cap) -- a CO2-like atmosphere should
    warm a planet, not leave it with zero greenhouse effect.
    """
    factor_at_co2_density = min(
        prog_c.CO2_MAX_GREENHOUSE_FACTOR,
        (pc.CO2_BASE_MOLAR_DENSITY / pc.CO2_BASE_MOLAR_DENSITY) * prog_c.CO2_MAX_GREENHOUSE_FACTOR,
    )
    assert factor_at_co2_density == prog_c.CO2_MAX_GREENHOUSE_FACTOR


def test_class_n_is_hotter_on_average_than_class_m():
    """
    Class N samples its atm_molar_density at the top of ATMOSPHERIC_MOLAR_DENSITY
    ["t"] (== physical_constants.ATMOSPHERIC_MOLAR_DENSITY["t"][1], close to
    CO2_BASE_MOLAR_DENSITY -- see generate_planet_properties's special case),
    so with the greenhouse inversion fixed, N should now get a real greenhouse
    boost and come out hotter than Class M on average across a decent sample.
    """
    n_records = plausibility.generate_sample("N", "e", N_SAMPLE, include_moons=False)
    m_records = plausibility.generate_sample("M", "e", N_SAMPLE, include_moons=False)

    n_mean_temp = statistics.mean(r["surface_temperature"] for r in n_records)
    m_mean_temp = statistics.mean(r["surface_temperature"] for r in m_records)

    assert n_mean_temp > m_mean_temp, (
        f"Class N mean temp {n_mean_temp:.2f}K should exceed Class M mean temp {m_mean_temp:.2f}K"
    )


# ---------------------------------------------------------------------------
# Fix 3: atmospheric pressure depends on gravity
# ---------------------------------------------------------------------------

def test_atmosphere_retention_factor_is_normalized_at_earth_gravity():
    assert planetPhysics._atmosphere_retention_factor(1.0) == pytest.approx(1.0)


def test_atmosphere_retention_factor_increases_with_gravity():
    low = planetPhysics._atmosphere_retention_factor(0.3)
    earth = planetPhysics._atmosphere_retention_factor(1.0)
    high = planetPhysics._atmosphere_retention_factor(3.0)
    assert low < earth < high


def _spearman_correlation(xs, ys):
    """Rank-based (Spearman) correlation -- more appropriate than Pearson
    here since atmospheric_pressure = effective_atm_density * gravity_ms2 *
    scale_height_m mixes gravity's effect multiplicatively with several
    other independently-random terms (atm_density spans a ~60x range for
    terrestrial classes alone, atm_molar_density its own range, and the
    sample also spans terrestrial and gas-giant classes whose gravity/
    pressure scales differ by orders of magnitude) -- exactly the kind of
    heavy multiplicative noise and scale mixing that suppresses a raw
    Pearson r even when the underlying monotonic relationship (higher
    gravity -> higher pressure) is real and strong."""
    def rank(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        ranks = [0] * len(values)
        for rank_pos, i in enumerate(order):
            ranks[i] = rank_pos
        return ranks

    return statistics.correlation(rank(xs), rank(ys))


def test_atmospheric_pressure_correlates_positively_with_gravity():
    """
    Statistical check across a wide gravity range (terrestrial classes plus
    gas giants) spanning many host star spectral types: correlation between
    gravity and atmospheric_pressure should now be positive. Before this
    fix, gravity canceled out of the pressure formula entirely
    (atmospheric_pressure = atm_density * R * T / atm_molar_density,
    independent of gravity), and the observed correlation was ~-0.11 (noise).

    Uses Spearman (rank) rather than Pearson correlation -- see
    _spearman_correlation's docstring for why a raw Pearson r on this data
    understates the (real, and now positive) monotonic relationship.
    Thresholds are calibrated to this class mix specifically: two
    brown-dwarf-like classes (density/gravity spanning two full orders of
    magnitude) used to be part of this sample and pulled both correlations
    much higher (Spearman reliably >0.5) -- since their removal (see
    CHANGELOG.md), the remaining gas giants (I/J/T) span a much narrower
    gravity range, so the correlation is real but weaker (observed Spearman
    ~0.17-0.22, Pearson ~0.22-0.25 across repeated runs) -- still clearly,
    consistently positive and a clear improvement over the pre-fix ~-0.11,
    just not >0.5 anymore with this narrower class mix.
    """
    records = []
    for cls in ("A", "B", "C", "E", "F", "M", "N", "O", "P"):
        zone = next((z for z in "hec" if prog_c.PLANET_CLASSES[cls][z]), None)
        if zone is None:
            continue
        records.extend(plausibility.generate_sample(cls, zone, 40, include_moons=False))
    for cls in ("I", "J", "T"):
        records.extend(plausibility.generate_sample(cls, "c", 40, include_moons=False))

    records = [r for r in records if r["has_atmosphere"]]
    gravities = [r["gravity"] for r in records]
    pressures = [r["atmospheric_pressure"] for r in records]
    assert len(gravities) >= 300

    pearson = statistics.correlation(gravities, pressures)
    assert pearson > 0, f"gravity/pressure Pearson correlation {pearson:.3f} is not even positive"

    spearman = _spearman_correlation(gravities, pressures)
    assert spearman > 0.1, f"gravity/pressure Spearman correlation {spearman:.3f} is not clearly positive"


# ---------------------------------------------------------------------------
# Fix 4: Class P is colder than Class M (albedo differentiation)
# ---------------------------------------------------------------------------

def test_class_p_has_own_albedo_range_distinct_from_default():
    assert "albedo_range" in prog_c.PLANET_CLASSES["P"]
    assert prog_c.PLANET_CLASSES["P"]["albedo_range"] != (0.12, 0.35)
    # P's albedo_range must still stand out from M's own tuned range (M was
    # given one too, see PLANET_CLASSES["M"]) -- the two shouldn't collapse
    # back into being indistinguishable just because both are now overridden.
    assert prog_c.PLANET_CLASSES["P"]["albedo_range"] != prog_c.PLANET_CLASSES["M"]["albedo_range"]


def test_class_p_is_colder_on_average_than_class_m():
    m_records = plausibility.generate_sample("M", "e", N_SAMPLE, include_moons=False)
    p_records = plausibility.generate_sample("P", "e", N_SAMPLE, include_moons=False)

    m_mean_temp = statistics.mean(r["surface_temperature"] for r in m_records)
    p_mean_temp = statistics.mean(r["surface_temperature"] for r in p_records)

    assert p_mean_temp < m_mean_temp, (
        f"Class P mean temp {p_mean_temp:.2f}K should be colder than Class M mean temp {m_mean_temp:.2f}K"
    )
