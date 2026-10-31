# tests/test_oxygen_fugacity.py

"""GEN.498: one oxygen-fugacity scale, against the figures of the G3 research (docs/design/composition/body-composition-research.md section 15)."""

import pytest

from planetgen.physics import oxygen_fugacity as fo2

# T K: Hirschmann (2021) iron-wustite log10 fO2 at 1 bar.
IW_1BAR = {1500: -11.621, 1800: -8.571, 2000: -7.054, 2200: -5.817}


@pytest.mark.parametrize("temperature", sorted(IW_1BAR))
def test_the_reference_iron_wustite_agrees_with_the_research_at_one_bar(temperature):
    assert fo2.iw_log_fo2(temperature) == pytest.approx(IW_1BAR[temperature], abs=0.002)


@pytest.mark.parametrize("temperature", [1500.0, 2000.0, 3000.0])
def test_the_two_iron_phases_join_to_within_the_fit_at_the_boundary(temperature):
    boundary = -18.64 + 0.04359 * temperature - 5.069e-6 * temperature ** 2
    assert abs(fo2.iw_log_fo2(temperature, boundary + 1e-6) - fo2.iw_log_fo2(temperature, boundary - 1e-6)) < 0.03


def test_delta_iw_and_log_fo2_are_inverses():
    assert fo2.delta_iw(fo2.log_fo2(-1.4, 1700.0, 3.0), 1700.0, 3.0) == pytest.approx(-1.4)


@pytest.mark.parametrize(
    "delta, fraction",
    [(-1.26, 0.014), (-0.64, 0.02), (0.45, 0.037), (1.00, 0.05), (2.28, 0.10), (3.68, 0.20), (5.66, 0.44)],
)
def test_melt_iron_matches_the_research_table(delta, fraction):
    assert fo2.fe3_fraction(delta) == pytest.approx(fraction, abs=0.0025, rel=0.03)


def test_a_bad_temperature_raises():
    with pytest.raises(ValueError):
        fo2.iw_log_fo2(-1.0)


def test_every_constant_set_has_a_source_row():
    assert all(len(row) == 4 and all(row) for row in fo2.SOURCES)
    assert {row[2] for row in fo2.SOURCES} <= {"VERIFIED", "DERIVED", "CROSS-CHECKED", "STAND-IN"}
