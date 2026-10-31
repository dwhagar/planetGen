# tests/test_water_delivery.py

"""GEN.802: the delivered-water rows b1 to b4 follow one law in the formation distance and the mass, fitted to the Bern population synthesis."""

import json
import math
import os
import random

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import elements
from planetgen.physics.atmosphere_recipes import RECIPES

with open(os.path.join(elements.DATA_DIRECTORY, "ngpps_formation_water.json"), encoding="utf-8") as _handle:
    NGPPS = json.load(_handle)

ROWS = ("b1", "b2", "b3", "b4")


def _subject(mass, snow_ratio, plateau=1.0):
    return Subject("terrestrial", {"mass_earth": mass, "snow_ratio": snow_ratio, "const:b:c_b4": plateau, "core_radius_earth": mass ** 0.27})


def _shares(mass, snow_ratio, plateau=1.0):
    table = dict(RECIPES.probabilities("b", _subject(mass, snow_ratio, plateau)))
    return [table.get(row, 0.0) for row in ROWS]


@pytest.mark.parametrize("bin_", [b for b in NGPPS["bins"] if b["planets"] >= 100], ids=lambda b: f"{b['x_low']}-{b['x_high']}")
def test_the_class_shares_follow_the_population_against_the_formation_distance(bin_):
    """The mean of the rows over the bin (log-uniform in x, 0.5 to 30 Earth masses log-uniform) is the observed share within 0.12."""
    low, high = max(bin_["x_low"], 0.03), min(bin_["x_high"], 4.0)
    xs = [math.exp(math.log(low) + (i + 0.5) / 10 * (math.log(high) - math.log(low))) for i in range(10)]
    masses = [10 ** (math.log10(0.5) + (j + 0.5) / 8 * (math.log10(30.0) - math.log10(0.5))) for j in range(8)]
    mean = [sum(_shares(m, x)[k] for x in xs for m in masses) / (len(xs) * len(masses)) for k in range(4)]
    assert mean == pytest.approx(bin_["observed"], abs=0.12)


def test_the_shares_are_a_probability_table_that_is_ordered_in_the_water():
    for x in (0.03, 0.2, 0.6, 1.0, 2.5, 4.0, 9.0):
        for mass in (0.2, 1.0, 5.0, 30.0):
            shares = _shares(mass, x)
            assert sum(shares) == pytest.approx(1.0, abs=2e-5) and min(shares) >= 0.0
    for mass in (1.0, 8.0):
        water_rich = [_shares(mass, x)[3] for x in (0.05, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 3.0)]
        dry = [_shares(mass, x)[0] for x in (0.05, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 3.0)]
        assert water_rich == sorted(water_rich) and dry == sorted(dry, reverse=True)


def test_a_planet_inside_the_snow_line_can_be_water_rich_as_the_population_has_it():
    """The gate x >= 1 of the earlier row b4 contradicted the population: a quarter of the planets at x = 0.55 and a half at 0.75 are water-rich."""
    assert _shares(3.0, 0.55)[3] == pytest.approx(0.25, abs=0.08)
    assert _shares(3.0, 0.75)[3] == pytest.approx(0.55, abs=0.1)
    assert _shares(3.0, 0.75)[0] < 0.1


def test_the_mass_law_of_gen_510_is_not_applied_on_top_of_the_distance_law():
    """The earlier mass logistic (centre 6.1 Earth masses, slope 6.28 per dex) took the water-rich share from 2 percent at 1 Earth mass to 91
    percent at 12; here, at a fixed distance, the mass moves the water-rich share by a few percent (0.13 per dex in the logit)."""
    for x in (0.5, 1.0, 1.5):
        assert abs(_shares(12.0, x)[3] - _shares(1.0, x)[3]) < 0.08, x
    assert _shares(1.0, 1.5)[3] > 0.95                      # not 0.02 as the mass law alone gave


def test_the_mass_matters_for_the_dry_share_at_a_fixed_distance_as_in_the_population():
    """Small planets are less often dry at a given x (the population's -0.84 per dex in the logit of 'not dry')."""
    assert _shares(0.6, 0.3)[0] < _shares(20.0, 0.3)[0]


def test_a_sample_of_a_thousand_planets_follows_the_population_in_distance_and_mass():
    rng = random.Random(802)
    tally = {}
    for _ in range(1000):
        x = math.exp(rng.uniform(math.log(0.05), math.log(3.0)))
        mass = 10 ** rng.uniform(math.log10(0.5), math.log10(30.0))
        shares = _shares(mass, x)
        pick = rng.choices(range(4), weights=shares)[0]
        band = 0 if x < 0.3 else (1 if x < 0.7 else (2 if x < 1.2 else 3))
        row = tally.setdefault(band, [0, 0, 0, 0, 0])
        row[pick] += 1
        row[4] += 1
    observed = {0: (0.05, 0.3), 1: (0.3, 0.7), 2: (0.7, 1.2), 3: (1.2, 3.0)}
    for band, (low, high) in observed.items():
        total = tally[band][4]
        assert total > 60
        covered = [b for b in NGPPS["bins"] if b["x_low"] >= low - 1e-9 and b["x_high"] <= high + 1e-9]
        weight = sum(b["planets"] for b in covered)
        expect = [sum(b["observed"][k] * b["planets"] for b in covered) / weight for k in range(4)]
        got = [tally[band][k] / total for k in range(4)]
        assert got == pytest.approx(expect, abs=0.16), (band, got, expect)


def test_the_plateau_bracket_covers_the_published_shares_and_the_rows_are_labelled_models():
    low, high = RECIPES.document["groups"]["b"]["constants"]["c_b4"][:2]
    assert low == 0.0 and high == 1.0 and high >= 0.87
    for row in ROWS:
        entry = RECIPES.entry(row)
        # b1's dry carbon is derived from the building blocks and the core (GEN.1005), so it is partly verified like the others; its limit stays the model's
        assert entry.state == "partly verified" and "MODEL" in entry.limit and "not an observation" in entry.limit, row
        assert entry.status == "DERIVED", row


def test_the_plateau_moves_the_water_rich_planets_it_removes_into_the_high_class():
    full, half, none = _shares(5.0, 0.75), _shares(5.0, 0.75, 0.5), _shares(5.0, 0.75, 0.0)
    assert half[3] == pytest.approx(0.5 * full[3]) and none[3] == 0.0
    assert half[2] == pytest.approx(full[2] + 0.5 * full[3]) and none[2] == pytest.approx(full[2] + full[3])
    assert half[0] == pytest.approx(full[0]) and half[1] == pytest.approx(full[1])


def test_a_drawn_plateau_stays_inside_its_bracket():
    from planetgen.util import draw
    with draw.bound(7):
        values = [RECIPES.draw_group("b", Subject("terrestrial", {"mass_earth": 8.0, "snow_ratio": 1.5, "core_radius_earth": 1.8})).constants[("b", "c_b4")]
                  for _ in range(20)]
    assert all(0.0 <= v <= 1.0 for v in values) and len(set(values)) > 1


def test_the_data_file_carries_its_fit_and_the_coefficients_are_nested():
    model = NGPPS["model"]
    assert NGPPS["label"] == "COMPUTED" and NGPPS["fit"]["planets"] == 10506
    assert model["x_range"] == [0.03, 4.0] and model["mass_range_earth"] == [0.5, 30.0] and model["knot"] == 0.6
    assert len(model["coefficients"]) == 3 and all(len(c) == 4 for c in model["coefficients"])
    assert sum(b["planets"] for b in NGPPS["bins"]) == NGPPS["fit"]["planets"]
    for bin_ in NGPPS["bins"]:
        assert sum(bin_["observed"]) == pytest.approx(1.0, abs=0.01) and sum(bin_["model"]) == pytest.approx(1.0, abs=0.01)
