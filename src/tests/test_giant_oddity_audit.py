# tests/test_giant_oddity_audit.py

"""
TEST.493, the giant oddity audit of docs/design/systems/weird-worlds.md (rows G1 to G11): the super-Jupiter mass tail (G9), the
hot Neptune desert and its survivors (G7), and that the rows which are derived outputs are not draws, with the ice-giant heat
spread left to GEN.1033. Over samples of the pure mass and architecture laws, drawn from named streams so the sample is the same every run.
"""

import math

import pytest

from planetgen import tuning
from planetgen.physics import constants, mass_boundaries, oddities, planet_architecture, planet_mass
from planetgen.util import draw

EARTH_KG = constants.EARTH_MASS_TO_KG
JUPITER_EARTH = 317.83
SUN_KG = 1.989e30
BROWN_DWARF_EDGE_KG = mass_boundaries.kg("deuterium_burning_limit")
DERIVED_ONLY = {"G2", "G3", "G4", "G5", "G6", "G8", "G10", "G11"}      # section 8: outputs of the interior, spin, ring and cloud models, not draws


def _runaway(star_solar, count):
    return [planet_mass.runaway_mass_kg(star_solar * SUN_KG, draw.Stream(f"g9:{star_solar}:{i}")) for i in range(count)]


@pytest.mark.parametrize("star_solar", [1.0, 1.5])
def test_the_runaway_mass_tail_reaches_ten_jupiter_masses_and_stops_at_the_brown_dwarf_edge(star_solar):
    masses = _runaway(star_solar, 20000)
    jupiters = [m / EARTH_KG / JUPITER_EARTH for m in masses]
    assert max(masses) <= BROWN_DWARF_EDGE_KG
    heavy = sum(1 for j in jupiters if j >= 10.0) / len(jupiters)
    slope = tuning.GIANT_MASS_FUNCTION_SLOPE
    low = tuning.GIANT_RUNAWAY_MIN_MASS_EARTH / JUPITER_EARTH
    ceiling_kg = min(BROWN_DWARF_EDGE_KG, tuning.PLANET_MAX_STAR_MASS_RATIO * star_solar * SUN_KG)
    high = ceiling_kg / EARTH_KG / JUPITER_EARTH
    expected = (10.0 ** -slope - high ** -slope) / (low ** -slope - high ** -slope)
    assert heavy == pytest.approx(expected, abs=0.004)
    assert heavy > 0.002                      # a minority tail (HAT-P-2 b, Kepler-39 b and the like) that exists
    # the ceiling is 0.01 of the star's mass until that passes the 13 Jupiter-mass edge, from about 1.24 solar masses up
    assert max(jupiters) > 0.9 * high


def test_a_light_star_cannot_make_a_super_jupiter():
    ceiling = max(tuning.GIANT_RUNAWAY_MIN_MASS_EARTH * EARTH_KG, tuning.PLANET_MAX_STAR_MASS_RATIO * 0.15 * SUN_KG)      # the floor wins where no giant fits
    assert max(_runaway(0.15, 2000)) <= ceiling * (1.0 + 1e-9)


def test_the_architecture_of_ordinary_stars_makes_super_jupiters_but_none_above_the_brown_dwarf_edge():
    heaviest, heavy = 0.0, 0
    for seed in range(900):
        fe = draw.Stream(f"g9fe:{seed}").gauss(0.0, 0.25)
        plan = planet_architecture.build(SUN_KG, 2.7, 10.0 ** fe, fe, 30.0, draw.Stream(f"g9:{seed}"))
        for planet in plan.planets:
            heaviest = max(heaviest, planet.mass_kg)
            heavy += planet.mass_kg / EARTH_KG / JUPITER_EARTH >= 8.0
    assert heavy >= 1
    assert heaviest <= BROWN_DWARF_EDGE_KG


def test_no_giant_of_a_heavy_metal_rich_star_is_above_the_brown_dwarf_edge():
    for star_solar in (1.5, 2.0):
        for i in range(3000):
            stream = draw.Stream(f"g9over:{star_solar}:{i}")
            out = planet_mass.planet_mass(5.0 + i % 20, star_solar * SUN_KG, 2.7 * star_solar, 0.3, 2.0, stream)
            assert out.mass_kg <= BROWN_DWARF_EDGE_KG, (star_solar, i, out.mass_kg / EARTH_KG / JUPITER_EARTH)


def _period_days(distance_au, star_solar=1.0):
    return math.sqrt(distance_au ** 3 / star_solar) * 365.25


def _chains(count):
    found = []
    for seed in range(count):
        fe = draw.Stream(f"g7fe:{seed}").gauss(0.0, 0.2)
        plan = planet_architecture.build(SUN_KG, 2.7, 10.0 ** fe, fe, 30.0, draw.Stream(f"g7:{seed}"))
        found.extend((_period_days(p.distance_au), p.mass_kg / EARTH_KG, p.runaway) for p in plan.planets)
    return found


def test_few_neptunes_sit_inside_three_days_and_a_few_survive_there():
    planets = _chains(1500)
    neptunes = [(period, mass) for period, mass, _ in planets if 10.0 <= mass < 100.0]
    inside = [mass for period, mass in neptunes if period < 3.0]
    beside = [mass for period, mass in neptunes if 3.0 <= period < 10.0]
    assert len(neptunes) > 400
    assert len(inside) >= 1                      # the survivors: LTT 9779 b (about 29 Earth masses at 0.8 d) is the real case
    assert len(inside) * 3 < len(beside)         # the desert: a third or fewer of what the next octave holds (Mazeh et al. 2016)
    assert max(inside) > 25.0 or len(inside) >= 3


def test_hot_jupiters_inside_three_days_are_rare_and_stay_jupiters():
    planets = _chains(1500)
    jupiters = [(period, mass) for period, mass, _ in planets if mass >= 100.0]
    hot = [mass for period, mass in jupiters if period < 3.0]
    assert jupiters and len(hot) / len(jupiters) < 0.05          # about 1 percent of Sun-like stars host a hot Jupiter, a few percent of giants


def test_the_oddity_table_has_no_row_for_a_derived_giant_output():
    assert not DERIVED_ONLY & {row.id for row in oddities.ROWS}
    for row in oddities.ROWS:
        if row.id.startswith("G"):
            assert row.kind == "giant" and row.id in {"G1", "G7", "G9"}      # section 8: only G1 and the survivor outliers need a rate row


def test_the_ice_giant_heat_spread_is_left_to_gen_1033():
    giant_heat = [row.id for row in oddities.ROWS
                  if row.kind == "giant" and getattr(row.apply, "input", None) == "internal_heat_flux"]
    assert giant_heat == [], "G2 (the Uranus heat anomaly) is GEN.1033's; an oddity row on a giant's internal heat would duplicate it"
