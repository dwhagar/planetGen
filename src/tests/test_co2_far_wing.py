# tests/test_co2_far_wing.py

"""
GEN.1019: the far wing of carbon dioxide (`planetgen.physics.co2_far_wing`): the Perrin and Hartmann factor, the Lorentz tail, and the
pure-function checks of the wing a few band strengths make. The Venus check against the k-tables is in test_radiative_convective_real.py.
"""

import json

import numpy as np
import pytest

from planetgen.physics import co2_far_wing as fw

EDGES = np.array([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])
"""ndarray: Ten bands of 100 cm-1."""


def one_band(index, strength=1.0e-18):
    strengths = np.zeros((1, len(EDGES) - 1))
    strengths[0, index] = strength
    return strengths


def test_chi_is_one_inside_sigma1_and_falls_beyond():
    assert fw.chi(0.0, 300.0, "co2_co2") == 1.0 and fw.chi(3.0, 300.0, "co2_co2") == 1.0
    distance = np.array([4.0, 10.0, 30.0, 60.0, 120.0, 300.0])
    values = fw.chi(distance, 300.0, "co2_co2")
    assert np.all(np.diff(values) < 0.0) and np.all((values > 0.0) & (values <= 1.0))


def test_chi_matches_the_hand_calculation_at_300_k():
    """b1 = 0.0888 - 0.160 exp(-0.0041 T), b2 = 0.0526 exp(-0.00152 T), b3 = 0.0232 (the HELIOS-K file of [R3490])."""
    t = 300.0
    b1 = 0.0888 - 0.160 * np.exp(-0.0041 * t)
    b2 = 0.0526 * np.exp(-0.00152 * t)
    expected = np.exp(-(b1 * 27.0 + b2 * 90.0 + 0.0232 * 80.0))
    assert fw.chi(200.0, t, "co2_co2") == pytest.approx(expected, rel=1e-9)


def test_chi_never_exceeds_one_when_cold_and_holds_the_measured_range():
    cold = fw.chi(np.array([10.0, 50.0]), 100.0, "co2_n2")
    assert np.all(cold <= 1.0)
    assert fw.chi(80.0, 100.0, "co2_co2") == fw.chi(80.0, fw.number("temperature_low_k"), "co2_co2")
    assert fw.chi(80.0, 2000.0, "co2_co2") == fw.chi(80.0, fw.number("temperature_high_k"), "co2_co2")


def test_nitrogen_drops_the_line_faster_than_carbon_dioxide_does():
    """The CO2-N2 sigma2 is 10 cm-1 against 30, so the nitrogen-broadened wing is the lower one beyond the first break."""
    assert fw.chi(50.0, 300.0, "co2_n2") < fw.chi(50.0, 300.0, "co2_co2")


def test_the_half_width_grows_with_pressure_and_falls_with_temperature():
    assert fw.half_width_cm1(296.0, 101325.0, 1.0) == pytest.approx(fw.number("gamma_self_cm1_atm"))
    assert fw.half_width_cm1(296.0, 101325.0, 0.0) == pytest.approx(fw.number("gamma_foreign_cm1_atm"))
    assert fw.half_width_cm1(296.0, 2.0e5, 1.0) > fw.half_width_cm1(296.0, 1.0e5, 1.0)
    assert fw.half_width_cm1(600.0, 1.0e5, 1.0) < fw.half_width_cm1(300.0, 1.0e5, 1.0)


def test_the_wing_is_zero_inside_the_cutoff_and_at_no_pressure():
    strengths = one_band(3)
    result = fw.cross_section_cm2(strengths, EDGES, [300.0], [1.0e6], [1.0])[0]
    assert result[3] > 0.0                                   # a 100 cm-1 band reaches itself beyond 25 cm-1
    narrow = np.array([0.0, 10.0, 20.0, 30.0])
    assert fw.cross_section_cm2(np.array([[1.0e-18, 0.0, 0.0]]), narrow, [300.0], [1.0e6], [1.0])[0, 0] == 0.0
    assert not fw.cross_section_cm2(strengths, EDGES, [300.0], [0.0], [1.0]).any()
    assert not fw.cross_section_cm2(strengths, EDGES, [300.0], [1.0e6], [0.0]).any()


def test_the_wing_is_linear_in_the_strength_and_in_the_pressure():
    base = fw.cross_section_cm2(one_band(2), EDGES, [300.0], [1.0e6], [1.0])
    assert fw.cross_section_cm2(one_band(2, 3.0e-18), EDGES, [300.0], [1.0e6], [1.0]) == pytest.approx(3.0 * base, rel=1e-12)
    double = fw.cross_section_cm2(one_band(2), EDGES, [300.0], [2.0e6], [1.0])
    assert double == pytest.approx(2.0 * base, rel=0.01)    # the 1/(gamma^2 + d^2) shape is flat while gamma is far under the distance


def test_the_wing_falls_with_the_distance_from_the_band():
    result = fw.cross_section_cm2(one_band(0), EDGES, [300.0], [5.0e6], [1.0])[0]
    assert np.all(np.diff(result[1:6]) < 0.0)


def test_the_wing_ends_at_the_reach(tmp_path, monkeypatch):
    """The widest range of the two truncations in question (120 cm-1 and the whole grid): more reach, more wing, and none past it."""
    totals = []
    for reach in (120.0, 500.0, 1.0e5):
        monkeypatch.setitem(fw.NUMBERS["reach_cm1"], "value", reach)
        totals.append(fw.cross_section_cm2(one_band(0), EDGES, [300.0], [5.0e6], [1.0])[0])
    assert np.all(totals[0] <= totals[1] + 1e-30) and np.all(totals[1] <= totals[2] + 1e-30)
    assert totals[0][4:].sum() < 1e-30 and totals[1][7:].sum() < 1e-30 and totals[2][9] > totals[1][9] + 1e-34


def test_a_nitrogen_air_has_a_smaller_wing_than_a_pure_carbon_dioxide_one():
    pure = fw.cross_section_cm2(one_band(2), EDGES, [300.0], [1.0e6], [1.0])
    diluted = fw.cross_section_cm2(one_band(2), EDGES, [300.0], [1.0e6], [0.1])
    assert diluted.sum() < 0.5 * pure.sum()


def test_band_strengths_are_the_mean_times_the_width():
    assert fw.band_strengths(np.array([1.0e-20, 2.0e-20]), np.array([0.0, 10.0, 40.0])) == pytest.approx([1.0e-19, 6.0e-19])


def test_every_number_of_the_rule_file_has_a_tier_and_a_tweak_and_a_bad_file_is_refused(tmp_path, monkeypatch):
    for name, entry in fw.NUMBERS.items():
        assert entry["tier"].startswith("T") and entry["tweak"].startswith("L") and entry["basis"], name
    with open(fw.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    del document["numbers"]["reach_cm1"]["tier"]
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    monkeypatch.setattr(fw, "DATA_FILE", str(path))
    with pytest.raises(fw.FarWingError):
        fw._load()
