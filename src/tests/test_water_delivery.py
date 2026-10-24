# tests/test_water_delivery.py

"""GEN.510: the water-rich share of row b4 against mass is the logistic fitted to the Bern population synthesis."""

import json
import math
import os

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import elements
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw

with open(os.path.join(elements.DATA_DIRECTORY, "ngpps_water_fraction.json"), encoding="utf-8") as _handle:
    NGPPS = json.load(_handle)


def _subject(mass, snow_ratio=1.5, plateau=1.0):
    return Subject("V", {"mass_earth": mass, "snow_ratio": snow_ratio, "const:b:c_b4": plateau, "core_radius_earth": mass ** 0.27})


@pytest.mark.parametrize("bin_", NGPPS["bins"], ids=lambda b: f"{b['mass_low']}-{b['mass_high']}")
def test_the_water_rich_probability_follows_the_population_synthesis_bins_within_two_sigma(bin_):
    mass = math.sqrt(bin_["mass_low"] * bin_["mass_high"])
    assert RECIPES.probability_of("b4", _subject(mass)) == pytest.approx(bin_["fraction"], abs=2 * NGPPS["sigma"])


def test_the_probability_is_the_logistic_of_the_fit_and_scales_with_the_plateau():
    centre, slope = NGPPS["logistic"]["centre_log10_earth_masses"], NGPPS["logistic"]["slope_per_dex"]
    for mass in (0.8, 3.0, 6.1, 20.0):
        logistic = 1 / (1 + math.exp(-slope * (math.log10(mass) - centre)))
        assert RECIPES.probability_of("b4", _subject(mass)) == pytest.approx(logistic, rel=1e-9)
        assert RECIPES.probability_of("b4", _subject(mass, plateau=0.5)) == pytest.approx(0.5 * logistic, rel=1e-9)
    assert RECIPES.probability_of("b4", _subject(6.1)) == pytest.approx(0.5, abs=0.01), "half of the planets at the centre, 6.1 Earth masses"
    assert RECIPES.probability_of("b4", _subject(5.0, plateau=0.0)) == 0.0


def test_a_planet_born_inside_the_snow_line_is_never_water_rich_and_a_heavy_one_cannot_reach_probability_one():
    assert "b4" not in dict(RECIPES.probabilities("b", _subject(10.0, snow_ratio=0.5)))
    assert RECIPES.probability_of("b4", _subject(1.0e4)) < 1.0


def test_the_plateau_bracket_covers_the_published_shares_and_the_row_is_a_labelled_model():
    low, high = RECIPES.document["groups"]["b"]["constants"]["c_b4"][:2]
    assert low == 0.0 and high == 1.0 and high >= 0.87
    entry = RECIPES.entry("b4")
    assert entry.state == "partly verified" and "MODEL" in entry.limit and "not an observation" in entry.limit


def test_a_drawn_plateau_stays_inside_its_bracket():
    with draw.bound(7):
        values = [RECIPES.draw_group("b", Subject("V", {"mass_earth": 8.0, "snow_ratio": 1.5, "core_radius_earth": 1.8})).constants[("b", "c_b4")]
                  for _ in range(20)]
    assert all(0.0 <= v <= 1.0 for v in values) and len(set(values)) > 1
