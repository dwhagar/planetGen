# tests/test_acid_freezing.py

"""GEN.1644: the freezing curve of the sulfuric acid hydrates (`physics/acid_freezing.py`), read from Figure 1 of Zeleznik 1991 [R4877]."""

import numpy as np
import pytest

from planetgen.physics import acid_freezing as af
from planetgen.physics import sulfuric_vapor

PHASE_MOLE_FRACTION = {"ice": 0.0, "acid": 1.0, "1H": 1 / 2, "2H": 1 / 3, "3H": 1 / 4, "4H": 1 / 5, "6.5H": 1 / 7.5}


def test_the_congruent_melting_points_of_table_6_are_on_the_curve():
    for phase in ("ice", "acid", "1H", "4H"):
        assert af.freezing_point_k(PHASE_MOLE_FRACTION[phase]) == pytest.approx(af.ANCHORS[phase], abs=0.05), phase


def test_the_incongruent_reference_points_are_within_a_kelvin_of_the_curve_at_their_composition():
    # a dihydrate, trihydrate and 6.5 hydrate melt incongruently: the curve runs near their Table 6 temperature at the composition of the solid
    for phase in ("2H", "3H", "6.5H"):
        assert abs(af.freezing_point_k(PHASE_MOLE_FRACTION[phase]) - af.ANCHORS[phase]) < 3.0, phase


def test_the_stable_eutectics_are_the_low_points_between_the_peaks():
    for name, (x, temperature) in af.EUTECTICS.items():
        assert af.freezing_point_k(x) == pytest.approx(temperature, abs=0.2), name
        assert af.freezing_point_k(max(x - 0.03, 0.0)) > temperature - 20 or name == "ice+6.5H"
    assert af.LOWEST_K == pytest.approx(af.EUTECTICS["ice+6.5H"][1], abs=0.01)
    assert af.HIGHEST_K == pytest.approx(af.ANCHORS["acid"])


def test_the_ice_branch_falls_with_acid_and_the_monohydrate_peak_is_the_highest_of_the_hydrates():
    ice = [af.freezing_point_k(x) for x in np.linspace(0.0, 0.1, 11)]
    assert all(a > b for a, b in zip(ice, ice[1:]))
    hydrates = np.linspace(0.13, 0.73, 601)
    peak = max(hydrates, key=af.freezing_point_k)
    assert peak == pytest.approx(0.5, abs=0.02) and af.freezing_point_k(peak) == pytest.approx(af.ANCHORS["1H"], abs=0.1)


def test_the_phase_that_freezes_out_follows_the_diagram():
    assert af.solid_phase(0.0) == "ice" and af.solid_phase(0.05) == "ice"
    assert af.solid_phase(0.2) == "4H" and af.solid_phase(0.5) == "1H" and af.solid_phase(1.0) == "acid"
    assert af.solid_phase(0.9) == "acid" and af.solid_phase(0.6) == "1H"


def test_known_strengths_freeze_where_the_literature_puts_them():
    # 38 weight percent is the ice eutectic (-62 C), 57.6 the tetrahydrate (-28 C), 75 percent the Venus cloud acid (about -30 C)
    assert af.freezing_point_of_weight_percent(38.0) == pytest.approx(211.1, abs=1.0)
    assert af.freezing_point_of_weight_percent(57.6) == pytest.approx(244.9, abs=0.6)
    assert 235.0 < af.freezing_point_of_weight_percent(75.0) < 250.0
    assert af.freezing_point_of_weight_percent(100.0) == pytest.approx(283.45, abs=0.05)
    assert af.freezing_point_of_weight_percent(0.0) == pytest.approx(273.15, abs=0.1)


def test_is_frozen_is_the_comparison_with_the_curve():
    x = sulfuric_vapor.mole_fraction(75.0)
    boundary = af.freezing_point_k(x)
    assert af.is_frozen(boundary - 1.0, x) and af.is_frozen(boundary, x) and not af.is_frozen(boundary + 1.0, x)
    assert not af.is_frozen(300.0, 0.3) and af.is_frozen(200.0, 0.3)


def test_the_curve_stays_inside_the_published_range_everywhere():
    temperatures = [af.freezing_point_k(x) for x in np.linspace(0.0, 1.0, 1001)]
    assert min(temperatures) >= af.LOWEST_K - 1e-9 and max(temperatures) <= af.HIGHEST_K + 1e-9
    assert 200.0 < af.LOWEST_K < 215.0


def test_a_composition_outside_the_range_is_refused():
    for bad in (-0.01, 1.01):
        with pytest.raises(af.AcidFreezingError):
            af.freezing_point_k(bad)
    with pytest.raises(af.AcidFreezingError):
        af.freezing_point_of_weight_percent(101.0)
