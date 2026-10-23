# tests/test_oxygen_fugacity.py

"""GEN.498: one oxygen-fugacity scale, against the figures of the G3 research (docs/design/body-composition-research.md section 15)."""

import pytest

from planetgen.physics import oxygen_fugacity as fo2

# T K: (Frost IW, Hirschmann IW) log10 fO2 at 1 bar.
IW_1BAR = {1500: (-11.624, -11.621), 1800: (-8.570, -8.571), 2000: (-7.043, -7.054), 2200: (-5.793, -5.817)}

# (T K, P GPa): QFM minus IW with Frost's FMQ and the reference IW.
QFM_MINUS_IW = {
    (1500, 1e-4): 3.625, (1800, 1e-4): 3.364, (2000, 1e-4): 3.240, (2200, 1e-4): 3.145,
    (1500, 1.0): 3.970, (1800, 1.0): 3.651, (2000, 1.0): 3.499, (2200, 1.0): 3.379,
    (1500, 3.0): 4.671, (1800, 3.0): 4.235, (2000, 3.0): 4.024, (2200, 3.0): 3.856,
}

LOG_K = {1500: (5.311, 5.724), 1800: (3.687, 4.269), 2173: (2.303, 3.018), 2200: (2.221, 2.943)}
CROSSOVER = {1500: (1.00, 0.17), 1800: (1.20, 0.03), 2173: (1.37, -0.06), 2200: (1.38, -0.07)}

# dIW: (CO2/CO, H2O/H2) at 2173 K.
RATIOS_2173 = {-4.7: (9.3e-4, 4.8e-3), -3.0: (6.6e-3, 0.034), -2.0: (0.021, 0.108), 0.0: (0.208, 1.08), 0.5: (0.369, 1.92), 2.0: (2.08, 10.8), 3.0: (6.6, 34.0)}


@pytest.mark.parametrize("temperature", sorted(IW_1BAR))
def test_the_reference_and_frost_iron_wustite_agree_with_the_research_at_one_bar(temperature):
    frost, hirschmann = IW_1BAR[temperature]
    assert fo2.frost_log_fo2("IW", temperature) == pytest.approx(frost, abs=0.002)
    assert fo2.iw_log_fo2(temperature) == pytest.approx(hirschmann, abs=0.002)
    assert abs(fo2.frost_log_fo2("IW", temperature) - fo2.iw_log_fo2(temperature)) < 0.03


@pytest.mark.parametrize("temperature", [1500.0, 2000.0, 3000.0])
def test_the_two_iron_phases_join_to_within_the_fit_at_the_boundary(temperature):
    boundary = -18.64 + 0.04359 * temperature - 5.069e-6 * temperature ** 2
    assert abs(fo2.iw_log_fo2(temperature, boundary + 1e-6) - fo2.iw_log_fo2(temperature, boundary - 1e-6)) < 0.03


@pytest.mark.parametrize("key", sorted(QFM_MINUS_IW))
def test_qfm_stands_three_to_five_log_units_above_iw(key):
    temperature, pressure = key
    assert fo2.buffer_offset("FMQ", temperature, pressure) == pytest.approx(QFM_MINUS_IW[key], abs=0.005)


def test_qfm_minus_iw_is_3_1_to_3_7_at_one_bar_from_1500_to_2200_k():
    values = [fo2.buffer_offset("FMQ", t) for t in range(1500, 2201, 100)]
    assert all(3.1 <= v <= 3.7 for v in values)
    assert values == sorted(values, reverse=True)


@pytest.mark.parametrize("temperature", sorted(LOG_K))
def test_the_equilibrium_constants_match_the_research(temperature):
    k1, k2 = LOG_K[temperature]
    assert fo2.log_k("co2_co", temperature) == pytest.approx(k1, abs=0.002)
    assert fo2.log_k("h2o_h2", temperature) == pytest.approx(k2, abs=0.002)


@pytest.mark.parametrize("temperature", sorted(CROSSOVER))
def test_the_gases_are_equal_at_the_research_crossovers(temperature):
    co, h = CROSSOVER[temperature]
    assert fo2.delta_iw_where_equal("co2_co", temperature) == pytest.approx(co, abs=0.01)
    assert fo2.delta_iw_where_equal("h2o_h2", temperature) == pytest.approx(h, abs=0.01)
    assert fo2.co2_to_co(fo2.delta_iw_where_equal("co2_co", temperature), temperature) == pytest.approx(1.0)


@pytest.mark.parametrize("delta", sorted(RATIOS_2173))
def test_the_gas_ratios_at_2173_k_match_the_research_within_one_percent(delta):
    co, h = RATIOS_2173[delta]
    assert fo2.co2_to_co(delta, 2173.0) == pytest.approx(co, rel=0.012)
    assert fo2.h2o_to_h2(delta, 2173.0) == pytest.approx(h, rel=0.012)


def test_sossi_gas_ratios_at_his_iw_plus_half_are_the_reference_at_zero():
    sossi_value = 0.5
    reference = fo2.convert_delta_iw(sossi_value, "sossi2020", "hirschmann2021")
    assert reference == pytest.approx(0.0)
    # Sossi's own ratios at 2173 K: H2O/H2 1.05, CO2/CO 0.21.
    assert fo2.h2o_to_h2(reference, 2173.0) == pytest.approx(1.05, rel=0.04)
    assert fo2.co2_to_co(reference, 2173.0) == pytest.approx(0.21, rel=0.02)
    # Treating his +0.5 as a reference value would have moved the ratios by 3.2 times.
    assert fo2.h2o_to_h2(0.5, 2173.0) / fo2.h2o_to_h2(0.0, 2173.0) == pytest.approx(10 ** 0.25, rel=1e-6)


@pytest.mark.parametrize("scale", fo2.SCALES)
@pytest.mark.parametrize("other", fo2.SCALES)
def test_converting_there_and_back_is_the_identity(scale, other):
    there = fo2.convert_delta_iw(1.7, scale, other, 1800.0, 2.0)
    assert fo2.convert_delta_iw(there, other, scale, 1800.0, 2.0) == pytest.approx(1.7, abs=1e-12)


def test_the_frost_scale_differs_from_the_reference_by_the_buffers_difference():
    expected = fo2.frost_log_fo2("IW", 2000.0) - fo2.iw_log_fo2(2000.0)
    assert fo2.convert_delta_iw(0.0, "frost1991", "hirschmann2021", 2000.0, 1e-4) == pytest.approx(expected)
    assert abs(expected) < 0.03


def test_delta_iw_and_log_fo2_are_inverses():
    assert fo2.delta_iw(fo2.log_fo2(-1.4, 1700.0, 3.0), 1700.0, 3.0) == pytest.approx(-1.4)


def test_the_earth_mantle_value_of_the_generator_sits_between_the_buffers():
    # REDOX_EARTH_DELTA_IW (3.5) is QFM, which is 3.1 to 3.7 above IW at magma-ocean temperatures.
    from planetgen.physics import atmosphere

    assert 3.1 <= atmosphere.REDOX_EARTH_DELTA_IW <= 3.7


@pytest.mark.parametrize(
    "fraction, expected",
    [(0.014, -1.26), (0.02, -0.64), (0.037, 0.45), (0.05, 1.00), (0.10, 2.28), (0.20, 3.68), (0.44, 5.66)],
)
def test_melt_iron_to_delta_iw_matches_the_research_table(fraction, expected):
    assert fo2.delta_iw_from_fe3_fraction(fraction) == pytest.approx(expected, abs=0.01)
    assert fo2.fe3_fraction(fo2.delta_iw_from_fe3_fraction(fraction)) == pytest.approx(fraction)


def test_a_fraction_outside_zero_to_one_raises():
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError):
            fo2.delta_iw_from_fe3_fraction(bad)


def test_an_unknown_scale_or_temperature_raises():
    with pytest.raises(ValueError):
        fo2.convert_delta_iw(0.0, "mystery", "hirschmann2021")
    with pytest.raises(ValueError):
        fo2.iw_log_fo2(-1.0)
    with pytest.raises(ValueError):
        fo2.log_k("co2_co", 50.0)


def test_every_constant_set_has_a_source_row_and_the_unchecked_ones_say_stand_in():
    assert all(len(row) == 4 and all(row) for row in fo2.SOURCES)
    assert {row[2] for row in fo2.SOURCES} <= {"VERIFIED", "DERIVED", "CROSS-CHECKED", "STAND-IN"}
    assert [row[0] for row in fo2.SOURCES if row[2] == "STAND-IN"] == ["Frost 1991 FMQ and NNO constants"]


def test_nno_stands_only_a_few_tenths_above_fmq_at_one_gigapascal():
    # Anenburg 2026: "at 1 GPa, NNO is only 0.3 log units above FMQ"; at 1 bar the gap is about 0.7 to 0.8.
    for temperature in (1000.0, 1273.0, 1473.0, 1873.0):
        gap_1gpa = fo2.frost_log_fo2("NNO", temperature, 1e4) - fo2.frost_log_fo2("FMQ", temperature, 1e4)
        gap_1bar = fo2.frost_log_fo2("NNO", temperature) - fo2.frost_log_fo2("FMQ", temperature)
        assert 0.1 <= gap_1gpa <= 0.4
        assert 0.7 <= gap_1bar <= 0.8


def test_the_fit_statistics_are_the_papers():
    assert fo2.IW_FIT_RMS == 0.00647 and fo2.IW_FIT_MAX_MISMATCH == 0.0283
