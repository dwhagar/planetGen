"""
Surface dose kernels (GEN.444, `planetgen.physics.surface_dose`): the Band spectrum of a stellar particle event
and the shielding of a mixed column. The Solar System bodies
are test targets, never inputs of the code
(docs/design/composition/radiation-model.md sections 11 and 13).
"""

import pytest

from planetgen.physics import solar_particles
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


def test_a_mixed_column_is_as_thick_as_the_sum_of_its_gases_in_air_units():
    assert sd.effective_air_column_g_cm2({"air": 100.0}) == pytest.approx(100.0)
    assert sd.effective_air_column_g_cm2({"CO2": 88.9}) == pytest.approx(90.1)
    assert sd.effective_air_column_g_cm2({"H2": 52.0}) == pytest.approx(90.1)
    assert sd.effective_air_column_g_cm2({"N2": 10.0, "O2": 5.0}) == pytest.approx(15.0)


def test_a_gas_without_an_interaction_length_is_an_error():
    with pytest.raises(DataError):
        sd.effective_air_column_g_cm2({"SF6": 1.0})
