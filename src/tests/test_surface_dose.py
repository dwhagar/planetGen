"""
Surface dose kernels (GEN.444 slice 1, `planetgen.physics.surface_dose`): the Band spectrum of a stellar particle
event and the share past a cutoff, the shielding of a mixed column and the radiation weights. The Solar System bodies
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
    assert event.flux(sd.rigidity_gv(30.0)) / event.flux(sd.rigidity_gv(200.0)) == pytest.approx(28.0, rel=0.01)


def test_band_spectrum_refuses_nonsense():
    with pytest.raises(DataError):
        sd.BandSpectrum(1e8, 2.0, 1.0, 0.2)
    with pytest.raises(DataError):
        sd.BandSpectrum(1e8, 1.0, 5.0, 0.0)


def test_no_cutoff_lets_the_whole_event_through():
    assert sd.share_past_cutoffs(sd.EVENT_2006_12_13, 0.0, 0.0) == pytest.approx(1.0)


def test_a_column_and_a_field_each_cut_the_share():
    event = sd.EVENT_2006_12_13
    mars = sd.rigidity_gv(proton_range.energy_mev(21.0, "CO2"))
    earth = sd.rigidity_gv(proton_range.energy_mev(1033.23, "air"))
    assert sd.share_past_cutoffs(event, mars) == pytest.approx(0.057, rel=0.02)
    assert sd.share_past_cutoffs(event, earth) < sd.share_past_cutoffs(event, mars)
    assert sd.share_past_cutoffs(event, 0.0, 14.34) < 1.0
    assert sd.share_past_cutoffs(event, mars, 14.34) < sd.share_past_cutoffs(event, mars)


def test_earth_air_and_field_stop_nearly_every_particle_of_an_event():
    earth = sd.rigidity_gv(proton_range.energy_mev(1033.23, "air"))
    assert sd.share_past_cutoffs(sd.EVENT_2006_12_13, earth, 14.34) < 1e-5


def test_the_rigidity_of_a_proton():
    assert sd.rigidity_gv(0.0) == 0.0
    assert sd.rigidity_gv(1000.0) == pytest.approx(1.696, rel=0.001)


def test_a_mixed_column_is_as_thick_as_the_sum_of_its_gases_in_air_units():
    assert sd.effective_air_column_g_cm2({"air": 100.0}) == pytest.approx(100.0)
    assert sd.effective_air_column_g_cm2({"CO2": 88.9}) == pytest.approx(90.1)
    assert sd.effective_air_column_g_cm2({"H2": 52.0}) == pytest.approx(90.1)
    assert sd.effective_air_column_g_cm2({"N2": 10.0, "O2": 5.0}) == pytest.approx(15.0)


def test_a_gas_without_an_interaction_length_is_an_error():
    with pytest.raises(DataError):
        sd.effective_air_column_g_cm2({"SF6": 1.0})


def test_neutron_weight_peaks_near_one_mev_and_never_drops_below_the_floor():
    assert sd.neutron_weight(1.0) == pytest.approx(20.7, rel=0.005)
    peak = max(sd.neutron_weight(e) for e in (0.3, 0.5, 1.0, 1.5, 2.0))
    assert peak == pytest.approx(20.7, abs=0.5)
    assert min(sd.neutron_weight(e) for e in (1e-9, 1e-6, 1e3, 1e6)) >= 2.5
    assert sd.neutron_weight(1e-9) == pytest.approx(2.5, abs=0.05)


def test_neutron_weight_is_continuous_at_its_edges():
    assert sd.neutron_weight(1.0 - 1e-9) == pytest.approx(sd.neutron_weight(1.0 + 1e-9), rel=1e-3)
    assert sd.neutron_weight(50.0 - 1e-9) == pytest.approx(sd.neutron_weight(50.0 + 1e-9), rel=0.001)


def test_neutron_weight_refuses_a_nonpositive_energy():
    with pytest.raises(ValueError):
        sd.neutron_weight(0.0)
