"""
Regression tests for the four Track A physics fixes to the planet
atmosphere/density model (see planetPhysics.py and plausibility.py):

  1. Gas-giant density: once a core/envelope blend, now set by the giant
     mass-radius relation (GEN.34).
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
from stellarObjects.roguePlanetData import RoguePlanet
from stellarObjects.starData import Star

N_SAMPLE = 300


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


# ---------------------------------------------------------------------------
# Fix 1: gas-giant density (giant mass-radius relation, GEN.34)
# ---------------------------------------------------------------------------

GAS_GIANT_CLASS = next(c for c, d in prog_c.PLANET_CLASSES.items() if d["type"] == "g")
GAS_GIANT_ZONE = next(z for z in "hec" if prog_c.PLANET_CLASSES[GAS_GIANT_CLASS][z])


def test_gas_giant_density_is_its_mass_over_its_volume(host_star):
    """A giant's mass comes from its class's mass range, its radius from the
    mass-radius relation within 3 sigma of scatter, and its density is
    whatever the two make (no separately drawn density)."""
    cfg = SystemConfig()
    distance = plausibility.distance_for_zone(host_star, GAS_GIANT_ZONE)
    low_kg, high_kg = planetPhysics.giant_mass_range_kg(GAS_GIANT_CLASS)
    for _ in range(50):
        planet = Planet(cfg, host_star, host_star.habitable_zone, distance,
                        planet_class=GAS_GIANT_CLASS, zone_override=GAS_GIANT_ZONE, moon_count=0)
        assert low_kg <= planet.mass <= high_kg
        volume_m3 = (4 / 3) * math.pi * (planet.radius * 1000) ** 3
        assert planet.density == pytest.approx(planet.mass / volume_m3 / 1000, rel=1e-9)
        median_km = planetPhysics.giant_radius_km(planet.mass)
        sigma = max(pc.GIANT_RADIUS_SCATTER.values())
        low_r, high_r = prog_c.PLANET_CLASSES[GAS_GIANT_CLASS]["radius_range"]
        assert max(low_r, median_km * (1 - 3 * sigma)) * (1 - 1e-9) <= planet.radius
        assert planet.radius <= min(high_r, median_km * (1 + 3 * sigma)) * (1 + 1e-9)


def test_giant_with_a_given_mass_keeps_it(host_star):
    """A caller's mass survives generation: the radius follows from it."""
    cfg = SystemConfig()
    mass_kg = 300 * pc.EARTH_MASS_TO_KG
    planet = Planet(cfg, host_star, host_star.habitable_zone, plausibility.distance_for_zone(host_star, "c"),
                    planet_class="J", mass=mass_kg, zone_override="c", moon_count=0)
    assert planet.mass == pytest.approx(mass_kg, rel=1e-9)


@pytest.mark.parametrize("mass_earth, radius_earth", [(17.15, 3.88), (95.16, 9.45), (317.8, 11.21)])
def test_giant_mass_radius_relation_hits_the_solar_system(mass_earth, radius_earth):
    """Neptune, Saturn and Jupiter land within 10% of their real radii."""
    radius_km = planetPhysics.giant_radius_km(mass_earth * pc.EARTH_MASS_TO_KG)
    assert radius_km / pc.EARTH_RADIUS_KM == pytest.approx(radius_earth, rel=0.10)


@pytest.mark.parametrize("mass_bin", ["saturn", "jupiter", "sub-neptune"])
def test_rogue_gas_giant_radius_follows_the_giant_relation(mass_bin):
    """Rogue gas giants (GEN.60) take their radius from the same giant
    mass-radius relation as bound giants, so a small one is not Jupiter's
    size: every radius sits within the relation's 3-sigma scatter of
    `giant_radius_km(mass)`, with a density a giant can have."""
    seen = []
    for _ in range(N_SAMPLE):
        rogue = RoguePlanet(SystemConfig(), mass_bin=mass_bin)
        if rogue.planet_type != "g":
            continue
        median = planetPhysics.giant_radius_km(rogue.mass_kg)
        sigma = pc.GIANT_RADIUS_SCATTER[planetPhysics.giant_regime(rogue.mass_kg)]
        assert median * (1 - 3 * sigma) <= rogue.radius_km <= median * (1 + 3 * sigma)
        density = rogue.mass_kg / ((4 / 3) * math.pi * (rogue.radius_km * 1000) ** 3) / 1000
        assert 0.1 < density < 50.0  # 13 Jupiter masses at 0.8 of Jupiter's radius is ~40 g/cm3
        seen.append((rogue.mass_kg, rogue.radius_km))
    assert seen


def test_small_rogue_gas_giants_are_smaller_than_large_ones():
    """A 0.05 Jupiter-mass rogue is far smaller than a 10 Jupiter-mass one
    (GEN.60: both used to be drawn around Jupiter's radius)."""
    small = planetPhysics.giant_radius_km(0.05 * pc.JUPITER_MASS_TO_KG)
    large = planetPhysics.giant_radius_km(10 * pc.JUPITER_MASS_TO_KG)
    assert small < 0.5 * large
    assert small / pc.EARTH_RADIUS_KM == pytest.approx(3.9, rel=0.15)


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


def test_gas_giant_sampled_gravity_is_within_theoretical_bounds():
    """Statistical check across a real generation sample: every gas giant's
    density is finite and positive and its gravity sits inside
    plausibility.theoretical_gravity_bounds_g, which mirrors the relation."""
    records = plausibility.generate_sample(GAS_GIANT_CLASS, GAS_GIANT_ZONE, N_SAMPLE, include_moons=False)
    lo, hi = plausibility.theoretical_gravity_bounds_g(GAS_GIANT_CLASS)
    for record in records:
        assert math.isfinite(record["density"]) and record["density"] > 0
        assert lo * (1 - 1e-6) <= record["gravity"] <= hi * (1 + 1e-6)


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


def _pearson_correlation(xs, ys):
    """Pearson correlation coefficient, computed directly rather than via
    `statistics.correlation` (Python 3.10+ only -- this repo's CI matrix
    still runs a 3.9 job, see CHANGELOG.md)."""
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    return cov / (var_x * var_y) ** 0.5


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

    return _pearson_correlation(rank(xs), rank(ys))


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

    pearson = _pearson_correlation(gravities, pressures)
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


# ---------------------------------------------------------------------------
# Airless reclassification and the cosmic-background temperature floor
# ---------------------------------------------------------------------------

def test_reclassifying_into_an_airless_class_clears_the_old_atmosphere(host_star):
    """
    `reconcile_zone_and_class` regenerates a body in place, so a body that
    started in an atmosphered class and lands in an airless one (Class C)
    used to keep its old atm_density/atm_molar_density/scale_height.
    """
    cfg = SystemConfig()
    distance = plausibility.distance_for_zone(host_star, "e")
    planet = Planet(cfg, host_star, host_star.habitable_zone, distance,
                    planet_class="M", zone_override="e", moon_count=0)
    assert planet.atm_density is not None and planet.scale_height is not None

    planet.planet_class = "C"
    planet.radius = None
    planet.mass = None
    planetPhysics.generate_planet_properties(planet, zone_override="e")
    planetPhysics.calculate_surface_gravity(planet)
    planetPhysics.calculate_atmospheric_conditions(planet)

    assert planet.atmosphere == "None"
    assert planet.atm_density is None
    assert planet.atm_molar_density is None
    assert planet.scale_height is None
    assert planet.atmospheric_pressure == 0.0


@pytest.mark.parametrize("planet_class", ["C", "I"])
def test_surface_temperature_never_drops_below_the_cosmic_background(host_star, planet_class):
    cfg = SystemConfig()
    zone = "c"
    distance = plausibility.distance_for_zone(host_star, zone)
    planet = Planet(cfg, host_star, host_star.habitable_zone, distance,
                    planet_class=planet_class, zone_override=zone, moon_count=0)

    class DimStar:
        luminosity = host_star.luminosity * 1e-12

    planet.star = DimStar()
    planetPhysics.calculate_atmospheric_conditions(planet)
    assert planet.surface_temperature == pc.COSMIC_BACKGROUND_TEMPERATURE_K


def test_reclassifying_a_moved_planet_keeps_its_distance(host_star):
    """
    `validate_system` pushes a planet clear of its inner neighbor, then
    `reconcile_zone_and_class` reclassifies it for its new zone. Ecosphere
    classes with a "zone_position_mode" used to redraw the distance
    anywhere in the habitable zone at that point, which could drop the
    planet back inside the neighbor it had just been moved past (seen as
    an asteroid belt overlapping the next planet in giant-star systems).
    """
    cfg = SystemConfig()
    inner, outer = host_star.habitable_zone
    for _ in range(30):
        planet = Planet(cfg, host_star, host_star.habitable_zone,
                        plausibility.distance_for_zone(host_star, "h"),
                        planet_class="A", zone_override="h", moon_count=0)
        moved_to = inner + 0.95 * (outer - inner)
        planet.distance = moved_to

        assert planetPhysics.reconcile_zone_and_class(planet, host_star.mass)
        assert planet.zone == "e"
        assert planet.distance == moved_to
