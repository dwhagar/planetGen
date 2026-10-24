"""
The gaps of the atmosphere recipe draw that GEN.674 closed: a comet beyond the coma's perihelion has no gas, the weights of a group whose
fallback is "none applies" are probabilities of occurrence, and a body with an ocean and a long-lived geology keeps little of its
carbon and sulfur in the air.
"""

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import mixture_draw as md
from planetgen.physics import mixture_rocky
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw
from tests.test_atmosphere_mixture_rocky import EARTH


def test_a_comet_beyond_the_coma_perihelion_has_no_gas_and_a_near_one_has_a_coma():
    far = Subject("comet", {"perihelion_au": 5.0, "family": "x", "age_gy": 4.5, "ts_k": 120.0})
    near = Subject("comet", {"perihelion_au": 1.0, "family": "x", "age_gy": 4.5, "ts_k": 200.0})
    for seed in range(4):
        with draw.bound(seed):
            result = md.draw_mixture(far)
        assert result.family == "airless" and result.mole_fractions == {} and result.pressure_bar == 0.0
        with draw.bound(seed):
            assert md.draw_mixture(near).family == "coma"


def test_a_row_alone_in_a_none_group_is_drawn_as_often_as_its_weight_says():
    # j5 (an abiotic-oxygen world) is the only row of group j open to this subject; its weight is the constant c_j5 in 0.05 to 0.5
    subject = Subject("Z", {"host_type": "M", "age_gy": 3.0, "mass_earth": 1.0, "zone": "e"})
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
    warm = Subject("C", {"ts_k": 181.0, "drawn.d": "d5"})
    table = dict(RECIPES.probabilities("j", warm))
    assert table.get("j4", 0.0) < 1e-6 and table["none"] > 1.0 - 1e-6
    cold = Subject("C", {"ts_k": 45.0, "drawn.d": "d5"})
    assert dict(RECIPES.probabilities("j", cold))["j4"] > 0.99


def test_weights_that_add_up_to_more_than_one_are_still_shares_in_a_none_group():
    coma = Subject("comet", {"perihelion_au": 1.0, "family": "x"})
    table = dict(RECIPES.probabilities("k", coma))
    assert sum(table.values()) == pytest.approx(1.0) and "none" not in table


def _inventory(subject, seed):
    constants = {}
    with draw.bound(seed):
        recipe = RECIPES.draw_all(subject, categories=mixture_rocky.ROCKY_CATEGORIES)
        return mixture_rocky.rocky_inventory(subject, recipe, constants), constants


def test_an_ocean_and_a_long_lived_geology_leave_little_carbon_and_sulfur_and_all_the_nitrogen_and_argon():
    wet = Subject("M", EARTH)
    dry = Subject("M", {**EARTH, "ocean_present": False})
    for seed in range(6):
        sunk, constants = _inventory(wet, seed)
        kept, _ = _inventory(dry, seed)
        sink = constants["weathering_sink"]
        low, high, _shape = mixture_rocky.INVENTORY["weathering_sink"]
        assert low <= sink <= high
        assert sunk["C"] == pytest.approx(kept["C"] * sink) and sunk["S"] == pytest.approx(kept["S"] * sink)
        assert sunk["N"] == kept["N"] and sunk["Ar"] == kept["Ar"] and sunk["H"] == kept["H"]


def test_an_earth_like_air_keeps_little_carbon_dioxide_and_sulfur_next_to_a_dry_world():
    def medians(subject):
        rows = []
        for seed in range(15):
            with draw.bound(seed):
                rows.append(md.draw_mixture(subject).mole_fractions)
        carbon = sorted(f.get("CO2", 0.0) + f.get("CO", 0.0) for f in rows)[7]
        sulfur = sorted(f.get("SO2", 0.0) + f.get("S2", 0.0) + f.get("OCS", 0.0) + f.get("H2S", 0.0) for f in rows)[7]
        return carbon, sulfur
    wet_carbon, wet_sulfur = medians(Subject("M", EARTH))
    dry_carbon, dry_sulfur = medians(Subject("M", {**EARTH, "ocean_present": False}))
    assert wet_carbon < 0.05 and wet_sulfur < 0.05
    assert wet_carbon < 0.1 * dry_carbon and wet_sulfur < 0.1 * dry_sulfur
