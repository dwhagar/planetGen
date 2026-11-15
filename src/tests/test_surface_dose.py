"""
Surface dose kernels (GEN.444, `planetgen.physics.surface_dose`): the Band spectrum of a stellar particle event
and the shielding of a mixed column. The Solar System bodies
are test targets, never inputs of the code
(docs/design/composition/radiation-model.md sections 11 and 13).
"""

import pytest

from planetgen.physics import gases, proton_range, radiation, solar_particles
from planetgen.physics import surface_dose as sd
from planetgen.physics.elements import DataError

def test_band_spectrum_is_continuous_at_its_break_and_falls_with_rigidity():
    for event in solar_particles.SPECTRA:
        below, above = event.flux(event.r1_gv * (1 - 1e-9)), event.flux(event.r1_gv * (1 + 1e-9))
        assert below == pytest.approx(above, rel=1e-3)
        fluxes = [event.flux(r) for r in (0.1, 0.3, 1.0, 3.0, 10.0)]
        assert fluxes == sorted(fluxes, reverse=True)


def test_band_spectrum_refuses_nonsense():
    with pytest.raises(DataError):
        sd.BandSpectrum(1e8, 2.0, 1.0, 0.2)
    with pytest.raises(DataError):
        sd.BandSpectrum(1e8, 1.0, 5.0, 0.0)


def test_the_interaction_length_of_a_gas_follows_from_its_elements():
    """The rule `1 / lambda = sum w / lambda_i` over the elements gives the Particle Data Group's own molecules (GEN.444, [R4926])."""
    def length(gas):
        return sd.interaction_length_g_cm2(sd.element_mass_fractions({gas: 1.0}))
    for gas, printed in (("CO2", 88.9), ("H2O", 83.3), ("CH4", 73.8), ("NH3", 79.5), ("N2", 89.7), ("O2", 90.2), ("H2", 52.0)):
        assert length(gas) == pytest.approx(printed, rel=0.002)
    assert sd.AIR_INTERACTION_LENGTH_G_CM2 == pytest.approx(90.1, rel=0.001)    # the printed value for dry air


def test_a_mixed_column_is_as_thick_as_the_sum_of_its_gases_in_air_units():
    assert sd.effective_air_column_g_cm2({"air": 100.0}) == pytest.approx(100.0)
    assert sd.effective_air_column_g_cm2({"CO2": 88.9}) == pytest.approx(90.1, rel=0.003)
    assert sd.effective_air_column_g_cm2({"H2": 52.0}) == pytest.approx(90.1, rel=0.002)
    assert sd.effective_air_column_g_cm2({"N2": 10.0, "O2": 5.0}) == pytest.approx(15.0, rel=0.01)
    assert sd.effective_air_column_g_cm2({"N2": 0.0}) == 0.0


def test_every_gas_of_the_table_has_an_interaction_length():
    for gas in gases.GASES:
        length = sd.interaction_length_g_cm2(sd.element_mass_fractions({gas.formula: 1.0}))
        assert 50.0 < length < 210.0, gas.formula


def test_krypton_and_bromine_lie_between_iron_and_xenon():
    fe, xe = sd.ELEMENT_INTERACTION_LENGTH_G_CM2["Fe"], sd.ELEMENT_INTERACTION_LENGTH_G_CM2["Xe"]
    for symbol in ("Kr", "Br"):
        assert fe < sd.ELEMENT_INTERACTION_LENGTH_G_CM2[symbol] < xe


def test_an_element_without_an_interaction_length_is_an_error():
    with pytest.raises(DataError):
        sd.effective_air_column_g_cm2({"ZnS": 1.0})


def test_a_stopping_mixture_uses_the_rows_of_the_range_table_and_the_nearest_for_the_rest():
    mix = sd.stopping_mixture({"CO2": 1.0, "SO2": 0.1, "HDO": 0.01}, proton_range.MATERIALS)
    assert set(mix) <= set(proton_range.MATERIALS)
    assert mix["Si"] > 0.0 and mix["H"] > 0.0           # sulfur takes silicon's row, deuterium hydrogen's
    assert sum(mix.values()) == pytest.approx(1.11)    # nothing is lost: the mass of the columns
    air_like = sd.stopping_mixture({"N2": 78.0, "O2": 21.0, "Ar": 1.0}, proton_range.MATERIALS)
    assert radiation.proton_cutoff_gv(100.0, air_like) == pytest.approx(radiation.proton_cutoff_gv(100.0), rel=0.005)


def test_the_air_of_a_body_is_split_into_columns_by_mass():
    columns = radiation.gas_columns_g_cm2({"N2": 78000.0, "O2": 21000.0}, 100.0)
    assert sum(columns.values()) == pytest.approx(100.0)
    assert columns["N2"] > columns["O2"]
    assert radiation.gas_columns_g_cm2({}, 100.0) == {}
    assert radiation.gas_columns_g_cm2({"N2": 1.0}, 0.0) == {}
