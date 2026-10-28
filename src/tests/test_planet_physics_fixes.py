"""
Regression tests for the Track A physics fixes to the planet model (see planets.py and plausibility.py):

  1. Gas-giant density: once a core/envelope blend, now set by the giant
     mass-radius relation (GEN.34).
  2. Class P gets its own ("cold, glaciated") albedo range so it's no
     longer statistically indistinguishable from Class M.

The greenhouse factor and the pressure of an air are not tested here: since GEN.411 and GEN.422 they are the solution of the fixed point
of `surface_air` (tests/test_surface_air.py), not class ranges.

Run with: pytest src/tests/test_planet_physics_fixes.py
"""
import copy
import math
import statistics

import pytest

from planetgen.physics import constants as pc
from planetgen.generation import plausibility
from planetgen.physics import planets
from planetgen import tuning as prog_c
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.star import Star

N_SAMPLE = 300

# TEST.410: a rocky body is a solve of its air, about 2 s on a 4-core machine (5 s on a loaded one) since GEN.411, so the sample that compares
# two classes' temperatures is 40 a class; measured over 60 each the medians were 312 K (M) and 245 K (P), and 300 K and 238 K over the first 40.
CLIMATE_SAMPLE = 40


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
    """A giant's mass comes from its class's mass range, its radius stays in its class's radius range (the interior
    tables of GEN.460 set it, and a hot young giant swells well past the median mass-radius relation, so no tight band
    around that median holds), and its density is whatever the two make (no separately drawn density)."""
    cfg = SystemConfig()
    distance = plausibility.distance_for_zone(host_star, GAS_GIANT_ZONE)
    low_kg, high_kg = planets.giant_mass_range_kg(GAS_GIANT_CLASS)
    for _ in range(50):
        planet = Planet(cfg, host_star, host_star.habitable_zone, distance,
                        planet_class=GAS_GIANT_CLASS, zone_override=GAS_GIANT_ZONE, moon_count=0)
        assert low_kg <= planet.mass <= high_kg
        volume_m3 = (4 / 3) * math.pi * (planet.radius * 1000) ** 3
        assert planet.density == pytest.approx(planet.mass / volume_m3 / 1000, rel=1e-9)
        low_r, high_r = prog_c.PLANET_CLASSES[GAS_GIANT_CLASS]["radius_range"]
        assert low_r <= planet.radius <= high_r


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
    radius_km = planets.giant_radius_km(mass_earth * pc.EARTH_MASS_TO_KG)
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
        median = planets.giant_radius_km(rogue.mass_kg)
        sigma = pc.GIANT_RADIUS_SCATTER[planets.giant_regime(rogue.mass_kg)]
        assert median * (1 - 3 * sigma) <= rogue.radius_km <= median * (1 + 3 * sigma)
        density = rogue.mass_kg / ((4 / 3) * math.pi * (rogue.radius_km * 1000) ** 3) / 1000
        assert 0.1 < density < 50.0  # 13 Jupiter masses at 0.8 of Jupiter's radius is ~40 g/cm3
        seen.append((rogue.mass_kg, rogue.radius_km))
    assert seen


def test_small_rogue_gas_giants_are_smaller_than_large_ones():
    """A 0.05 Jupiter-mass rogue is far smaller than a 10 Jupiter-mass one
    (GEN.60: both used to be drawn around Jupiter's radius)."""
    small = planets.giant_radius_km(0.05 * pc.JUPITER_MASS_TO_KG)
    large = planets.giant_radius_km(10 * pc.JUPITER_MASS_TO_KG)
    assert small < 0.5 * large
    assert small / pc.EARTH_RADIUS_KM == pytest.approx(3.9, rel=0.15)


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
# Class P is colder than Class M (albedo differentiation)
# ---------------------------------------------------------------------------

def test_class_p_has_own_albedo_range_distinct_from_default():
    assert "albedo_range" in prog_c.PLANET_CLASSES["P"]
    assert prog_c.PLANET_CLASSES["P"]["albedo_range"] != (0.12, 0.35)
    # P's albedo_range must still stand out from M's own tuned range (M was
    # given one too, see PLANET_CLASSES["M"]) -- the two shouldn't collapse
    # back into being indistinguishable just because both are now overridden.
    assert prog_c.PLANET_CLASSES["P"]["albedo_range"] != prog_c.PLANET_CLASSES["M"]["albedo_range"]


def test_class_p_is_colder_than_class_m():
    """Medians, not means: since GEN.411 a world on the edge of its liquid range can land on a hot, steam-thick branch of the fixed point
    (the air is bistable there), and a few of those at 900 K say nothing about the typical body of the class."""
    m_records = plausibility.generate_sample("M", "e", CLIMATE_SAMPLE, include_moons=False)
    p_records = plausibility.generate_sample("P", "e", CLIMATE_SAMPLE, include_moons=False)

    m_median_temp = statistics.median(r["surface_temperature"] for r in m_records)
    p_median_temp = statistics.median(r["surface_temperature"] for r in p_records)

    assert p_median_temp < m_median_temp, (
        f"Class P median temp {p_median_temp:.2f}K should be colder than Class M median temp {m_median_temp:.2f}K"
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
    planets.generate_planet_properties(planet, zone_override="e")
    planets.calculate_surface_gravity(planet)
    planets.calculate_atmospheric_conditions(planet)

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

    dim_star = copy.copy(host_star)
    dim_star.luminosity = host_star.luminosity * 1e-12
    planet.star = dim_star
    planets.calculate_atmospheric_conditions(planet)
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

        assert planets.reconcile_zone_and_class(planet, host_star.mass)
        assert planet.zone == "e"
        assert planet.distance == moved_to
