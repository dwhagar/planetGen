"""
Regression tests for the Track A physics fixes to the planet model (see planets.py and plausibility.py):

  1. Gas-giant density: once a core/envelope blend, now set by the giant
     mass-radius relation (GEN.34).
  2. A body redrawn in place keeps no trace of its old air or class.

The greenhouse factor and the pressure of an air are not tested here: since GEN.411 and GEN.422 they are the solution of the fixed point
of `surface_air` (tests/test_surface_air.py), not class ranges.

Run with: pytest src/tests/test_planet_physics_fixes.py
"""
import copy
import math

import pytest

from planetgen.physics import constants as pc
from planetgen.generation import plausibility
from planetgen.physics import planets
from planetgen import tuning as prog_c
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.star import Star
from tests import class_support

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

GIANT_ZONE = "c"


@pytest.mark.parametrize("mass_earth", [25.0, 100.0, 318.0, 1500.0])
def test_giant_density_is_its_mass_over_its_volume(host_star, mass_earth):
    """A giant of a given mass has a radius the interior tables of GEN.460 set (a hot young giant swells well past the median
    mass-radius relation, so no tight band around that median holds), and its density is whatever the two make (no separately
    drawn density)."""
    cfg = SystemConfig()
    distance = plausibility.distance_for_zone(host_star, GIANT_ZONE)
    mass_kg = mass_earth * pc.EARTH_MASS_TO_KG
    for seed in range(3):
        planet = class_support.draw_planet(host_star, distance, seed, cfg, mass=mass_kg, zone_override=GIANT_ZONE, moon_count=0)
        assert planet.body_type == "g"
        assert planet.mass == pytest.approx(mass_kg, rel=1e-9)
        volume_m3 = (4 / 3) * math.pi * (planet.radius * 1000) ** 3
        assert planet.density == pytest.approx(planet.mass / volume_m3 / 1000, rel=1e-9)
        assert planet.radius >= prog_c.GIANT_MIN_RADIUS_KM


def test_giant_with_a_given_mass_keeps_it(host_star):
    """A caller's mass survives generation: the radius follows from it."""
    cfg = SystemConfig()
    mass_kg = 300 * pc.EARTH_MASS_TO_KG
    planet = Planet(cfg, host_star, host_star.habitable_zone, plausibility.distance_for_zone(host_star, "c"),
                    mass=mass_kg, zone_override="c", moon_count=0)
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


def test_sampled_gravity_is_within_theoretical_bounds():
    """Statistical check across a real generation sample: every body's density is finite and positive and its gravity sits inside
    plausibility.theoretical_gravity_bounds_g for its type, which mirrors the relation."""
    records = plausibility.generate_sample("c", 40, include_moons=False)
    for record in records:
        lo, hi = plausibility.theoretical_gravity_bounds_g(record["body_type"])
        assert math.isfinite(record["density"]) and record["density"] > 0
        assert lo * (1 - 1e-6) <= record["gravity"] <= hi * (1 + 1e-6)


# ---------------------------------------------------------------------------
# A body drawn again in place and the cosmic-background temperature floor
# ---------------------------------------------------------------------------

def test_drawing_a_body_again_clears_the_old_atmosphere_and_class(host_star):
    """
    `reconcile_zone` regenerates a body in place, so a body that held an air must not keep its old atm_density/atm_molar_density/
    scale_height or its class while the new body is built.
    """
    distance = plausibility.distance_for_zone(host_star, "e")
    planet = class_support.planets_where(host_star, distance, lambda p: p.atm_density is not None, 1, moon_count=0)[0]
    assert planet.scale_height is not None and planet.planet_class is not None

    planet.radius = None
    planet.mass = None
    planets.generate_planet_properties(planet, zone_override="e")

    assert planet.atmosphere is None
    assert planet.atm_density is None
    assert planet.atm_molar_density is None
    assert planet.scale_height is None
    assert planet.planet_class is None

    planets.calculate_surface_gravity(planet)
    planets.calculate_atmospheric_conditions(planet)
    planets.apply_class(planet)
    assert planet.planet_class in prog_c.PLANET_CLASSES
    assert (planet.atmosphere == "None") == (planet.atmospheric_pressure == 0.0)


@pytest.mark.parametrize("mass_earth", [0.05, 20.0])
def test_surface_temperature_never_drops_below_the_cosmic_background(host_star, mass_earth):
    """A rock with no air round a starless-dim star sits at the background; an ice giant keeps its own heat (GEN.816), so it sits
    above the background at the temperature that heat radiates, tens of kelvin."""
    cfg = SystemConfig()
    zone = "c"
    distance = plausibility.distance_for_zone(host_star, zone)
    planet = class_support.draw_planet(host_star, distance, 0, cfg, mass=mass_earth * pc.EARTH_MASS_TO_KG,
                                       zone_override=zone, moon_count=0)

    dim_star = copy.copy(host_star)
    dim_star.luminosity = host_star.luminosity * 1e-12
    planet.star = dim_star
    planets.calculate_atmospheric_conditions(planet)
    if planet.body_type == "g":
        assert 30.0 < planet.surface_temperature < 100.0
    else:
        assert planet.surface_temperature == pc.COSMIC_BACKGROUND_TEMPERATURE_K


def test_redrawing_a_moved_planet_keeps_its_distance(host_star):
    """
    `validate_system` pushes a planet clear of its inner neighbor, then `reconcile_zone` redraws it for its new zone. The distance must
    stay where the planet was moved to, or the planet could drop back inside the neighbor it had just been moved past (seen as an
    asteroid belt overlapping the next planet in giant-star systems).
    """
    cfg = SystemConfig()
    inner, outer = host_star.habitable_zone
    for seed in range(8):
        planet = class_support.draw_planet(host_star, plausibility.distance_for_zone(host_star, "h"), seed, cfg,
                                           zone_override="h", moon_count=0)
        moved_to = inner + 0.95 * (outer - inner)
        planet.distance = moved_to

        assert planets.reconcile_zone(planet, host_star.mass)
        assert planet.zone == "e"
        assert planet.distance == moved_to
