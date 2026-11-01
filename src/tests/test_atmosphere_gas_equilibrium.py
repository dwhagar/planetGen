# tests/test_atmosphere_gas_equilibrium.py

"""The ideal-gas equilibrium solver (GEN.451, GEN.454, GEN.499) against the speciation of the research (part 05 section 3.4,
computed there with the same NASA Glenn data) and against the laws it must obey."""

import math

import pytest

from planetgen.physics import gas_equilibrium as ge, oxygen_fugacity as fo2, species_thermo as st

BUDGET = {"H": 1.0, "C": 0.1, "N": 0.02, "S": 0.02}
"""dict: Part 05 section 3.4's illustration: per H atom C 0.1, N 0.02, S 0.02."""


def outgassed(temperature, delta_iw, pressure):
    return ge.solve(BUDGET, temperature, pressure, fixed_o2_log_fugacity=fo2.log_fo2(delta_iw, temperature))


# Part 05 section 3.4, "1200 K, 1 bar, dIW = 0" and "dIW = +4"; the research's buffer is not the Hirschmann one, so 12 percent.
@pytest.mark.parametrize("delta_iw, expected", [
    (0.0, {"H2": 0.48, "H2O": 0.30, "CO": 0.11, "CO2": 0.051, "H2S": 0.033}),
    (4.0, {"H2O": 0.78, "CO2": 0.16, "SO2": 0.021, "H2S": 0.011}),
])
def test_the_outgassed_gas_at_1200_k_is_the_research_speciation(delta_iw, expected):
    result = outgassed(1200.0, delta_iw, 1.0)
    for name, value in expected.items():
        assert result.mole_fractions[name] == pytest.approx(value, rel=0.12), name


def test_the_reduced_gas_has_methane_only_as_a_trace_at_1200_k_and_one_bar():
    fractions = outgassed(1200.0, 0.0, 1.0).mole_fractions
    assert 5e-6 < fractions["CH4"] < 5e-5
    assert fractions["SO2"] < 1e-6


def test_methane_takes_the_carbon_in_cool_dense_reduced_gas():
    # Part 05 3.4: 1000 K, 100 bar, dIW <= +1: CH4 about 0.24, CO 1e-6 to 1e-3 of... the carbon is methane.
    fractions = outgassed(1000.0, 1.0, 100.0).mole_fractions
    assert fractions["CH4"] > 0.15
    assert fractions["CH4"] > 10 * fractions["CO2"] and fractions["CH4"] > 10 * fractions["CO"]


def test_ammonia_appears_at_800_k_and_100_bar_near_the_research_value():
    # Part 05 3.4: NH3 about 1.5e-2 to 2e-2 (mole fraction) at 800 K and 100 bar.
    fractions = outgassed(800.0, 0.0, 100.0).mole_fractions
    assert 1.2e-2 < fractions["NH3"] < 2.5e-2


def test_the_element_balance_and_the_mole_fractions_close():
    result = outgassed(1500.0, 2.0, 10.0)
    assert sum(result.mole_fractions.values()) == pytest.approx(1.0, abs=1e-12)
    atoms = {symbol: 0.0 for symbol in BUDGET}
    for name, moles in result.moles.items():
        for symbol, count in st.species(name).elements.items():
            if symbol in atoms:
                atoms[symbol] += count * moles
    for symbol, value in BUDGET.items():
        assert atoms[symbol] == pytest.approx(value, rel=1e-8), symbol


def test_the_oxygen_fugacity_of_the_gas_is_the_one_held():
    for delta_iw in (-2.0, 1.0, 5.0):
        result = outgassed(1300.0, delta_iw, 3.0)
        held = fo2.log_fo2(delta_iw, 1300.0)
        assert math.log10(result.mole_fractions["O2"] * 3.0) == pytest.approx(held, abs=1e-6)


def test_the_total_pressure_does_not_move_the_ratio_of_reactions_that_keep_the_mole_number():
    low = outgassed(1300.0, 1.0, 0.01).mole_fractions
    high = outgassed(1300.0, 1.0, 100.0).mole_fractions
    assert low["CO2"] / low["CO"] == pytest.approx(high["CO2"] / high["CO"], rel=1e-6)
    assert low["H2O"] / low["H2"] == pytest.approx(high["H2O"] / high["H2"], rel=1e-6)
    # but methane, which lowers the mole number, grows with pressure
    assert high["CH4"] > 100 * low["CH4"]


def test_the_water_gas_shift_follows_its_equilibrium_constant():
    for temperature in (900.0, 1500.0):
        fractions = outgassed(temperature, 1.0, 5.0).mole_fractions
        k = st.log10_k({"CO2": 1, "H2": 1, "CO": -1, "H2O": -1}, temperature)
        assert math.log10(fractions["CO2"] * fractions["H2"] / (fractions["CO"] * fractions["H2O"])) == pytest.approx(k, abs=1e-6)


def test_solar_gas_follows_zahnle_and_marleys_methane_constant():
    # Z and M 2014: K = pCH4 pH2O / (pCO pH2^3) = 5.24e-14 exp(27285 / T) atm^-2, here in bar and a free oxygen budget.
    abundances = {"H": 1.0, "He": 0.0851 * 0.5, "C": 2.9e-4 / 2, "O": 4.9e-4 / 2, "N": 6.8e-5 / 2}
    for temperature, pressure in ((800.0, 1.0), (1000.0, 10.0)):
        result = ge.solve(abundances, temperature, pressure)
        x = result.mole_fractions
        ratio = x["CH4"] * x["H2O"] / (x["CO"] * x["H2"] ** 3) * pressure ** -2
        expected = 5.24e-14 * math.exp(27285.0 / temperature) * 1.01325 ** 2
        assert ratio == pytest.approx(expected, rel=0.35), temperature


def test_solar_gas_holds_its_carbon_as_methane_when_cool_and_as_carbon_monoxide_when_hot():
    abundances = {"H": 1.0, "He": 0.0425, "C": 1.45e-4, "O": 2.45e-4, "N": 3.4e-5}
    cool = ge.solve(abundances, 600.0, 1.0).mole_fractions
    hot = ge.solve(abundances, 1800.0, 1.0).mole_fractions
    assert cool["CH4"] > 100 * cool["CO"]
    assert hot["CO"] > 100 * hot["CH4"]
    assert cool["NH3"] > cool["N2"] and hot["N2"] > hot["NH3"]


def test_a_solution_is_deterministic_and_independent_of_the_scale_of_the_abundances():
    a = ge.solve({"H": 1.0, "C": 0.1, "O": 0.3}, 1000.0, 2.0)
    b = ge.solve({"H": 10.0, "C": 1.0, "O": 3.0}, 1000.0, 2.0)
    assert a.mole_fractions == pytest.approx(b.mole_fractions, rel=1e-7)
    assert ge.solve({"H": 1.0, "C": 0.1, "O": 0.3}, 1000.0, 2.0).mole_fractions == a.mole_fractions


def test_a_species_whose_fit_stops_below_the_temperature_is_named_as_skipped():
    result = ge.solve({"H": 1.0, "O": 0.45}, 10000.0, 1.0)
    assert "H2O" in result.skipped and "H2O" not in result.mole_fractions
    assert result.mole_fractions["H"] > 0.5
    cold = ge.solve({"H": 1.0, "O": 0.45}, 250.0, 1.0)
    assert "H2O" not in cold.skipped and "O3" in cold.skipped      # O3 and H2O2 are fitted from 300 K
    assert cold.mole_fractions["H2O"] == pytest.approx(0.9, rel=1e-6) and cold.mole_fractions["H2"] == pytest.approx(0.1, rel=1e-6)


def test_refused_inputs():
    with pytest.raises(ValueError):
        ge.solve({"H": 1.0}, -5.0, 1.0)
    with pytest.raises(ValueError):
        ge.solve({"H": 1.0}, 1000.0, 0.0)
    with pytest.raises(ValueError):
        ge.solve({}, 1000.0, 1.0)
    with pytest.raises(ValueError):
        ge.solve({"H": 1.0, "O": 1.0}, 1000.0, 1.0, fixed_o2_log_fugacity=-12.0)
    with pytest.raises(ge.EquilibriumError):
        ge.solve({"Zz": 1.0}, 1000.0, 1.0)


@pytest.mark.parametrize("abundances", [
    {"H": 2.0000000259769655, "O": 1.0, "Ar": 1.0608501544615678e-12},
    {"H": 2.0000000024410274, "O": 1.0, "Ar": 1.01370169311403e-12, "C": 3.021494751116587e-11},
    {"H": 1.9980520712128589, "O": 1.0, "Ar": 1.2121667733025634e-12, "N": 2.043517676950264e-06},
])
def test_water_with_a_trace_of_argon_at_726_k_converges(abundances):
    """GEN.1003: the generator reached 'no convergence at 725.97 K and 1.0 bar' on a mixture of hydrogen and oxygen at the proportion of water with
    argon at 1e-12 of it. The Newton step jittered at the 1e-3 level for ever (the total row of the scaled system was rounding noise), so the solve
    is retried with that row scaled by the total."""
    result = ge.solve(abundances, 725.97, 1.0)
    assert result.mole_fractions["H2O"] > 0.999
    assert result.mole_fractions["Ar"] == pytest.approx(abundances["Ar"] / sum(abundances.values()), rel=1e-2)
    assert result.iterations < 150
    balance = sum(result.moles[name] * st.SPECIES[name].elements.get("H", 0) for name in result.moles) / abundances["H"]
    assert balance == pytest.approx(1.0, abs=1e-8)


def test_a_nitrogen_gas_with_1e_15_of_everything_else_leaves_the_rarest_elements_out_when_it_will_not_settle_with_them():
    """GEN.819: the late nitrogen leaves an air of N2 with the carbon, hydrogen, oxygen, sulfur and argon at 1e-9 to 1e-16 of it; the Newton step
    of that gas stalls with all five in ("the element balance is off"), and settles when the two rarest (under `TRACE_ELEMENT_FRACTION`) are out."""
    abundances = {"H": 3.735440154784931e-16, "Ar": 3.379812146652516e-21, "O": 6.93043683615438e-15, "C": 6.583546815986791e-15,
                  "N": 1.1914940104724934e-05, "S": 6.894484200178456e-20}
    with pytest.raises(ge.EquilibriumError):
        ge._solve(abundances, 434.2659729417367, 1.0000000000000002)
    result = ge.solve(abundances, 434.2659729417367, 1.0000000000000002)
    assert result.mole_fractions["N2"] > 0.999999
    assert "Ar" not in result.moles and "H2S" not in result.moles and "SO2" not in result.moles
