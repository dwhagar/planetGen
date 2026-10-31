"""
The gaps of the atmosphere recipe draw that GEN.674 closed: a comet beyond the coma's perihelion has no gas, the weights of a group whose
fallback is "none applies" are probabilities of occurrence. (The carbon and sulfur an ocean and a long-lived geology take out of the air
are the sinks of GEN.726 and GEN.757, tested in test_surface_pressure.py and test_mixture_sinks.py.)
"""

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import mixture_draw as md
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw


DRY_AIR = {"mass_earth": 0.002, "host_type": "G", "age_gy": 4.5, "zone": "c", "redox_diw": 0.0, "water_lost": False, "runaway_ratio": 0.5, "is_moon": False, "radius_km": 6000.0, "insolation": 1.0}
"""dict: What the gates of the evolved-air rows (table j) read besides the quantities a test sets, for a body that passes none of them."""

SOLAR_HOST = {"n_over_c": 1.0, "ar_over_c": 1.0, "s_over_c": 1.0}
"""dict: The host ratios a comet's species rows read (GEN.1005, k7), at the proto-solar values."""


def test_a_comet_beyond_the_coma_perihelion_has_no_gas_and_a_near_one_has_a_coma():
    far = Subject("comet", {"perihelion_au": 5.0, "family": "x", "age_gy": 4.5, "ts_k": 120.0, **SOLAR_HOST})
    near = Subject("comet", {"perihelion_au": 1.0, "family": "x", "age_gy": 4.5, "ts_k": 200.0, **SOLAR_HOST})
    for seed in range(4):
        with draw.bound(seed):
            result = md.draw_mixture(far)
        assert result.family == "airless" and result.mole_fractions == {} and result.pressure_bar == 0.0
        with draw.bound(seed):
            assert md.draw_mixture(near).family == "coma"


def test_a_row_alone_in_a_none_group_is_drawn_as_often_as_its_weight_says():
    # j5 (an abiotic-oxygen world) is the only row of group j open to this subject; its weight is the constant c_j5 in 0.05 to 0.5
    subject = Subject("terrestrial", {**DRY_AIR, "host_type": "M", "age_gy": 3.0, "mass_earth": 1.0, "zone": "e"})
    table = dict(RECIPES.probabilities("j", subject))
    assert set(table) == {"j5", "none"}
    assert table["j5"] + table["none"] == pytest.approx(1.0)
    assert 0.05 <= table["j5"] <= 0.5
    with_row = sum(1 for seed in range(300) if _drawn_row("j", subject, seed) == "j5")
    assert 0.05 * 300 * 0.4 < with_row < 0.5 * 300 * 1.3        # a minority, as the weight says, and not every eligible world


def _drawn_row(group, subject, seed):
    with draw.bound(seed):
        row = RECIPES.draw_group(group, subject).row
    return row.id if row else "none"


def test_a_frost_row_far_from_its_window_is_not_normalised_to_one():
    warm = Subject("terrestrial", {**DRY_AIR, "ts_k": 181.0, "drawn.d": "d5"})
    table = dict(RECIPES.probabilities("j", warm))
    assert table.get("j4", 0.0) < 1e-6 and table["none"] > 1.0 - 1e-6
    cold = Subject("terrestrial", {**DRY_AIR, "ts_k": 45.0, "drawn.d": "d5"})
    assert dict(RECIPES.probabilities("j", cold))["j4"] > 0.99


def test_weights_that_add_up_to_more_than_one_are_still_shares_in_a_none_group():
    coma = Subject("comet", {"perihelion_au": 1.0, "family": "x", **SOLAR_HOST})
    table = dict(RECIPES.probabilities("k", coma))
    assert sum(table.values()) == pytest.approx(1.0) and "none" not in table
