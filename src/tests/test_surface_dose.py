"""
The Band spectrum of a stellar particle event (GEN.444, `planetgen.physics.surface_dose`); the share past a cutoff, the shielding of a mixed column and the weights it once held were retired (GEN.1559). The Solar System bodies
are test targets, never inputs of the code
(docs/design/composition/radiation-model.md sections 11 and 13).
"""

import pytest

from planetgen.physics import proton_range
from planetgen.physics import surface_dose as sd
from planetgen.physics.elements import DataError

def test_band_spectrum_is_continuous_at_its_break_and_falls_with_rigidity():
    event = sd.EVENT_2006_12_13
    below, above = event.flux(event.r1_gv * (1 - 1e-9)), event.flux(event.r1_gv * (1 + 1e-9))
    assert below == pytest.approx(above, rel=1e-3)
    rigidities = (0.1, 0.3, 1.0, 3.0, 10.0)
    fluxes = [event.flux(r) for r in rigidities]
    assert fluxes == sorted(fluxes, reverse=True)


def test_band_ratio_of_30_to_200_mev_is_the_design_value():
    event = sd.EVENT_2006_12_13
    assert event.flux(proton_range.rigidity_gv(30.0)) / event.flux(proton_range.rigidity_gv(200.0)) == pytest.approx(28.0, rel=0.01)


def test_band_spectrum_refuses_nonsense():
    with pytest.raises(DataError):
        sd.BandSpectrum(1e8, 2.0, 1.0, 0.2)
    with pytest.raises(DataError):
        sd.BandSpectrum(1e8, 1.0, 5.0, 0.0)


def test_the_share_past_a_column_and_a_field_comes_from_the_measured_events():
    """The retired kernel of one event (GEN.1559) is `solar_particles.dose_share` of the 47: a column and a field each cut it."""
    from planetgen.physics import solar_particles as sp
    mars = proton_range.rigidity_gv(proton_range.energy_mev(21.0, "CO2"))
    earth = proton_range.rigidity_gv(proton_range.energy_mev(1033.23, "air"))
    assert sp.dose_share(0.0, 0.0) == pytest.approx(1.0)
    assert sp.dose_share(earth, 0.0) < sp.dose_share(mars, 0.0) < 1.0
    assert sp.dose_share(mars, 14.34) < sp.dose_share(mars, 0.0)
