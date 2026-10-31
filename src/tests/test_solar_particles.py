"""
Solar energetic particle events (GEN.747, `planetgen.physics.solar_particles`): the 47 measured spectra, the yearly dose they bring and the
share of it a column of air or a magnetic cutoff lets through, against the NOAA list of 255 events.
"""

import math

import pytest

from planetgen.physics import proton_range, solar_particles as sp
from planetgen.physics.surface_dose import BandSpectrum


def test_the_event_set_is_the_47_fits_of_the_two_tables():
    assert len(sp.EVENTS) == len(sp.SPECTRA) == 47
    assert sum(1 for event in sp.EVENTS if event[0] == "GLE") == 21 and sum(1 for event in sp.EVENTS if event[0] == "subGLE") == 26
    assert sum(1 for event in sp.EVENTS if "esp" in event[1]) == 4
    for kind, _, _, j0, gamma1, gamma2, r0 in sp.EVENTS:
        assert j0 > 0 and r0 > 0 and gamma2 > gamma1
    first = sp.EVENTS[0]
    assert (first[2], first[3], first[4], first[5], first[6]) == ("1997-Nov-06", 8.15e8, 0.28, 5.38, 0.116)  # GLE 55, as typed
    assert isinstance(sp.SPECTRA[0], BandSpectrum)


def test_the_set_is_the_events_of_100_pfu_and_above_in_the_noaa_list():
    """NOAA's list has 2.6 events a year at 100 pfu or more (109 in 42 years); the 43 events of the set (the four `esp` rows are second fits of
    four of them) over its 18 years are 2.4 a year, so the set misses almost none of the events that carry the fluence."""
    events = sum(1 for event in sp.EVENTS if "esp" not in event[1])
    assert events == 43
    assert events / sp.RECORD_YEARS == pytest.approx(sp.events_per_year_above(100.0), rel=0.15)


@pytest.mark.parametrize("peak, count", sp.NOAA_CUMULATIVE_COUNTS)
def test_the_event_rate_reproduces_the_noaa_counts(peak, count):
    assert sp.events_per_year_above(peak) * sp.NOAA_YEARS == pytest.approx(count, rel=1e-9)


def test_the_noaa_rate_falls_with_size_and_the_tail_is_steeper_than_the_record():
    peaks = [10.0, 20.0, 50.0, 100.0, 500.0, 2000.0, 1e4, 3e4, 1.2e5]
    rates = [sp.events_per_year_above(p) for p in peaks]
    assert rates == sorted(rates, reverse=True)
    assert 255.0 / 42.0 == pytest.approx(6.07, abs=0.01)
    # Usoskin and Kovaltsov 2012: four events of 1e10 cm-2 in 600 years of annual 10Be data, 0.0077 a year (0.0032 to 0.0150, 90 percent)
    assert sp.events_per_year_above(1.2e5) == pytest.approx(0.0077, rel=1e-9)
    assert 0.0032 < sp.events_per_year_above(1.2e5) < 0.0150
    assert sp.TAIL_SLOPE == pytest.approx(1.34, abs=0.01)
    # no event above 2e10 in 600 years: the one-sided 90 percent bound is 0.0027 a year, and the power law of the tail keeps under it
    assert sp.events_per_year_above(2.4e5) < 0.0027
    assert sp.events_per_year_above(1.2e5) / sp.events_per_year_above(2.4e5) == pytest.approx(2.0 ** sp.EXTREME_SLOPE)


def test_the_yearly_fluence_above_30_mev_is_a_few_1e8_to_1e9_and_the_dose_is_hundreds_of_msv_in_free_space():
    """18 years of GLE and sub-GLE events: 1.6e10 protons per cm2 above 30 MeV in all, 9e8 a year (cycle 23's big years and cycle 24's weak ones)."""
    fluence = sp.yearly_fluence(proton_range.rigidity_gv(30.0))
    assert fluence == pytest.approx(9e8, rel=0.2)
    assert 100.0 < sp.FREE_MSV_YR < 600.0
    assert sp.FREE_MSV_YR == pytest.approx(265.0, rel=0.05)


def test_the_largest_events_carry_the_dose():
    """Nov 2000, Jul 2000, Nov 2001 and Oct 2003 each add a tenth of the dose or more."""
    def dose_of(index):
        spectrum = sp.SPECTRA[index]
        total = 0.0
        flux = [(t, spectrum.flux(proton_range.rigidity_gv(t))) for t in [sp.proton_range.energy_mev(5.0, "H2O liquid") * 1.01 ** i for i in range(2000)]]
        for (t0, f0), (t1, f1) in zip(flux, flux[1:]):
            total += (f0 - f1) * sp.PROTON_WEIGHT * sp.GY_PER_MEV_PER_G * 1e3 * sp._stopping_mev_cm2_g(math.sqrt(t0 * t1))
        return total / sp.RECORD_YEARS
    doses = sorted(((dose_of(i), sp.EVENTS[i][2]) for i in range(len(sp.EVENTS))), reverse=True)
    assert doses[0][1] == "2000-Nov-09" and doses[0][0] > 0.1 * sp.FREE_MSV_YR
    assert sum(d for d, _ in doses) == pytest.approx(sp.FREE_MSV_YR, rel=0.03)
    assert sum(d for d, _ in doses[:8]) > 0.5 * sp.FREE_MSV_YR


def test_a_proton_needs_78_mev_to_reach_the_organs_and_the_dose_is_the_stopping_power_weighted_fluence():
    assert proton_range.energy_mev(sp.BFO_DEPTH_G_CM2, "H2O liquid") == pytest.approx(78.4, abs=1.0)  # PSTAR CSDA range 5 g/cm2
    assert sp.ORGAN_RIGIDITY_GV == pytest.approx(proton_range.rigidity_gv(78.4), rel=0.01)
    assert proton_range.stopping_power_mev_cm2_g(100.0, "H2O liquid") == pytest.approx(7.289, rel=0.01)  # PSTAR, 100 MeV in water


def test_the_dose_share_is_one_without_air_and_field_and_falls_with_column_and_cutoff():
    assert sp.dose_share(0.0) == 1.0 and sp.dose_share(0.0, 0.0) == 1.0 and sp.dose_share(sp.ORGAN_RIGIDITY_GV) == 1.0
    columns = [0.0, 5.0, 10.0, 22.6, 60.0, 100.0, 300.0, 1033.0]
    shares = [sp.dose_share(proton_range.rigidity_gv(proton_range.energy_mev(x))) for x in columns]
    assert shares == sorted(shares, reverse=True) and shares[-1] < 1e-4
    mars = sp.dose_share(proton_range.rigidity_gv(proton_range.energy_mev(22.6)))
    assert 0.03 < mars < 0.15  # only the protons above 172 MeV reach the ground
    earth = proton_range.rigidity_gv(proton_range.energy_mev(1033.0))
    assert sp.dose_share(0.0, 14.34) < sp.dose_share(0.0, 1.0) < 1.0  # a stronger field cuts more
    assert sp.dose_share(earth, 14.34) <= sp.dose_share(earth, 0.0)


def test_the_polar_caps_let_a_weak_field_through_and_the_latitude_average_is_between_the_pole_and_the_equator():
    low = sp.dose_share(0.0, 1.0)
    pole = sp.dose_share(0.0, 0.0)
    equator = sp.dose_above_msv_yr(1.0) / sp.FREE_MSV_YR
    assert equator < low < pole


def test_the_dose_above_a_rigidity_is_continuous_and_falls_to_nothing_above_the_grid():
    rigidities = [0.3 * 1.07 ** i for i in range(80)]
    doses = [sp.dose_above_msv_yr(r) for r in rigidities]
    assert doses == sorted(doses, reverse=True) and doses[0] == sp.FREE_MSV_YR
    for a, b in zip(doses, doses[1:]):
        assert b >= a * 0.5 or b < 1e-6 * sp.FREE_MSV_YR
    assert sp.dose_above_msv_yr(1e4) == 0.0
