"""
The volatile inventory of the rocky outgassing against the totals of the Earth, Venus and Mars (GEN.813,
docs/design/composition/volatile-inventory-against-planet-totals.md): the drawn carbon and nitrogen of a mantle, in ppm by weight, must
contain what those bodies are known to hold.
"""

import math

import numpy
import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import mixture_rocky as md
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw

C_G_MOL = 12.011
N_G_MOL = 14.007
SI_G_MOL = 28.086
SI_WEIGHT_FRACTION_BSE = 0.21           # the bulk silicate Earth (tier 3: recalled, 45 weight percent SiO2)

# Bergin et al. 2015 (arXiv:1507.04756) Table S1: atomic ratios to silicon of the bulk silicate Earth
BSE_C_PER_SI = 1.1e-3
BSE_N_PER_SI = 2.5e-5
BSE_C_OVER_N = 49.0

BSE_C_PPM = BSE_C_PER_SI * SI_WEIGHT_FRACTION_BSE / SI_G_MOL * C_G_MOL * 1e6
BSE_N_PPM = BSE_N_PER_SI * SI_WEIGHT_FRACTION_BSE / SI_G_MOL * N_G_MOL * 1e6

# Venus: 92 bar, 96.5 percent CO2 and 3.5 percent N2 by volume, g 8.87 m/s2, radius 6051.8 km (NASA planetary fact sheet); mantle 3.3e24 kg
VENUS_AIR_KG = 92e5 * 4.0 * math.pi * 6.0518e6 ** 2 / 8.87
VENUS_N2_MASS_SHARE = 0.035 * 28.013 / (0.035 * 28.013 + 0.965 * 44.01)
VENUS_MANTLE_KG = 3.3e24
VENUS_DEGASSED_SHARE = 0.43     # the share of its mantle carbon that the 4.6e20 kg of CO2 of its air needs (surface-pressure-supply-law.md 2.2)
VENUS_C_PPM = VENUS_AIR_KG * (1.0 - VENUS_N2_MASS_SHARE) * C_G_MOL / 44.01 / VENUS_MANTLE_KG * 1e6
VENUS_N_PPM = VENUS_AIR_KG * VENUS_N2_MASS_SHARE / VENUS_MANTLE_KG * 1e6

# Mars: the carbon dioxide of the source magma, 5 ppm at IW-1 to 1000 ppm at IW+1 (Hirschmann and Withers 2008; Grott et al. 2011,
# as the model of arXiv:2302.04241 adopts them), and the nitrogen of the source magma 1.9 ppm (same paper)
MARS_C_PPM = (5.0 * C_G_MOL / 44.01, 1000.0 * C_G_MOL / 44.01)
MARS_N_PPM = 1.9

# Nitrogen per carbon atom: the bulk silicate Earth 1/49, the CI chondrites 1/17, the air of Venus 1/13.5 and of Mars 1/18 (Table S1);
# Kargel and Lewis 1993 through Schaefer, Lodders and Fegley 2012: 0.013
N_PER_C = {"bulk silicate Earth": 1.0 / BSE_C_OVER_N, "CI chondrites": 1.0 / 17.0, "Venus air": 1.0 / 13.5, "Mars air": 1.0 / 18.0,
           "Kargel and Lewis": 0.013}

EARTH = dict(mass_earth=1.0, ts_k=288.0, pressure_bar=1.0, age_gy=4.5, long_lived_geology=False, ocean_present=False,
             snow_ratio=0.5, zone="e", host_type="G", is_moon=False, insolation=1.0, luminosity_solar=1.0, water_lost=False,
             runaway_ratio=0.5, teq_k=255.0, xuv_index=1.0, v_esc_kms=11.2, radius_km=6371.0, host_mass_solar=1.0)
SNOW_RATIO = {"Venus": 0.72 / 2.69, "Earth": 1.0 / 2.69, "Mars": 1.52 / 2.69}      # a over the snow line of 2.69 AU at one solar luminosity


def _ppm(snow_ratio, count=240):
    carbon, nitrogen = [], []
    for seed in range(count):
        subject = Subject("terrestrial", {**EARTH, "snow_ratio": snow_ratio})
        with draw.bound(seed):
            recipe = RECIPES.draw_all(subject, categories=md.ROCKY_CATEGORIES)
            inventory = md.rocky_inventory(subject, recipe, {})
        if recipe.groups["b"].name == "b4":
            continue        # the water-rich class is its own item (GEN.814)
        carbon.append(inventory["C"] * C_G_MOL / 100.0 * 1e6)
        nitrogen.append(inventory["N"] * N_G_MOL / 100.0 * 1e6)
    return numpy.array(carbon), numpy.array(nitrogen)


def _inside(value, sample, low=5.0, high=95.0):
    lo, hi = numpy.percentile(sample, [low, high])
    return lo <= value <= hi


def test_the_dry_class_carbon_is_the_bulk_silicate_earths_not_sixteen_times_it():
    low, high, shape = md.INVENTORY["co2_dry_wt_percent"]
    assert shape == "log"
    middle_ppm_c = math.sqrt(low * high) / 44.01 * C_G_MOL * 1e4
    # GEN.1005 replaced the bracket tuned to the bulk silicate Earth (log middle 115 ppm C) with one derived from Blanchard et al. 2022 (80 to 294 ppm C of
    # the body). The test followed the tuning: it asked the middle to sit within 20 percent of the 99 ppm of this table. The derived bracket's middle is the
    # model's own answer for the mantle's carbon, which that paper gives as 140 +/- 40 ppm; the Earth's 99 ppm lies inside the bracket (next line) and
    # inside the paper's spread of 100 to 180, which is what the derivation can be held to without tuning it back to one planet.
    assert 100.0 <= middle_ppm_c <= 180.0, middle_ppm_c
    assert low / 44.01 * C_G_MOL * 1e4 < BSE_C_PPM < high / 44.01 * C_G_MOL * 1e4
    assert 0.6 / 44.01 * C_G_MOL * 1e4 > 10.0 * BSE_C_PPM, "the delivered 0.6 weight percent of Ortenzi's table is what the bracket replaces"


def test_the_nitrogen_per_carbon_bracket_holds_every_measured_ratio_with_the_earths_in_the_log_middle():
    low, high, _ = md.INVENTORY["nitrogen_per_carbon"]
    for name, ratio in N_PER_C.items():
        assert low < ratio < high, name
    assert math.sqrt(low * high) == pytest.approx(N_PER_C["bulk silicate Earth"], rel=0.2)


@pytest.mark.parametrize("body", ["Venus", "Earth", "Mars"])
def test_an_earth_like_mantle_draws_the_bulk_silicate_earth_carbon_and_nitrogen_inside_its_spread(body):
    carbon, nitrogen = _ppm(SNOW_RATIO[body])
    assert _inside(BSE_C_PPM, carbon)
    assert _inside(BSE_N_PPM, nitrogen)
    assert numpy.median(carbon) == pytest.approx(BSE_C_PPM, rel=0.5)


def test_the_air_of_venus_is_inside_what_the_inventory_can_supply():
    carbon, nitrogen = _ppm(SNOW_RATIO["Venus"])
    assert VENUS_C_PPM == pytest.approx(38.0, rel=0.1) and VENUS_N_PPM == pytest.approx(3.3, rel=0.1)
    # the air is what the degassed share of the mantle inventory has brought up (no ocean, no carbonate sink, no nitrogen sink)
    assert _inside(VENUS_C_PPM / VENUS_DEGASSED_SHARE, carbon) and _inside(VENUS_N_PPM / VENUS_DEGASSED_SHARE, nitrogen)


def test_mars_is_inside_the_carbon_range_of_its_source_magma_and_the_nitrogen_of_the_silicate_earth():
    carbon, nitrogen = _ppm(SNOW_RATIO["Mars"])
    assert MARS_C_PPM[0] < numpy.median(carbon) < MARS_C_PPM[1]
    assert numpy.percentile(carbon, 95) < 2.0 * MARS_C_PPM[1]
    assert _inside(MARS_N_PPM, nitrogen)


def test_the_earth_class_draw_no_longer_puts_a_thick_nitrogen_air_in_a_tenth_of_the_worlds():
    """Finding 2 of reports/air-fixed-point-findings-2026-10-24.md: with 0.6 weight percent of carbon, 37 percent of the draws held ten
    times the bulk silicate Earth's nitrogen (one seed: 87 bar of N2 on a 1.8 Earth-mass planet); now none do."""
    _, nitrogen = _ppm(SNOW_RATIO["Earth"], count=600)
    assert numpy.mean(nitrogen > 10.0 * BSE_N_PPM) < 0.01
