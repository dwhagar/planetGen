# tests/test_brown_dwarf.py

"""GEN.252: a brown dwarf as a star: its mass, how bright and hot it is at an age (Burrows et al. 2001's two cooling
laws), the L, T or Y type that temperature makes (Kirkpatrick et al. 2021, Table 13), and where those reach a page
(a free-floating brown dwarf, a brown dwarf member of a multi-star tree).

The 20 pc census's densities (Kirkpatrick et al. 2021, Table 15) are test targets that show the cooling laws give a
plausible spread of temperatures. None of them is an input of the generator."""

import math

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.physics import brown_dwarf, constants
from planetgen.util import draw
from planetgen.web import scatter_members_view
from tests import render_support


def test_the_masses_run_from_the_deuterium_limit_to_the_hydrogen_limit():
    low, high = brown_dwarf.mass_range_solar()
    assert low == pytest.approx(13 * constants.JUPITER_MASS_TO_KG / constants.SOLAR_MASS_TO_KG)
    assert high == pytest.approx(80 * constants.JUPITER_MASS_TO_KG / constants.SOLAR_MASS_TO_KG)
    with draw.bound(252):
        masses = [brown_dwarf.sample_mass_solar() for _ in range(4000)]
    assert low <= min(masses) and max(masses) <= high


def test_the_mass_function_is_the_censuss_power_law():
    """dN/dM ~ M^-0.6: the share below a mass M is (M^0.4 - low^0.4) / (high^0.4 - low^0.4)."""
    low, high = brown_dwarf.mass_range_solar()
    with draw.bound(252):
        masses = [brown_dwarf.sample_mass_solar() for _ in range(20000)]
    power = 1.0 - tuning.BROWN_DWARF_MASS_FUNCTION_SLOPE
    for mass in (0.02, 0.04, 0.06):
        expected = (mass ** power - low ** power) / (high ** power - low ** power)
        assert sum(m < mass for m in masses) / len(masses) == pytest.approx(expected, abs=0.015)


def test_the_temperature_is_burrows_second_law_and_falls_with_age():
    assert brown_dwarf.old_age_temperature_k(0.05, 1.0) == pytest.approx(1550.0)
    assert brown_dwarf.old_age_temperature_k(0.05, 8.0) == pytest.approx(1550.0 * 8.0 ** -0.32)
    assert brown_dwarf.old_age_temperature_k(0.1, 1.0) == pytest.approx(1550.0 * 2.0 ** 0.83)
    assert brown_dwarf.effective_temperature_k(0.05, 8.0) < brown_dwarf.effective_temperature_k(0.05, 2.0)
    assert brown_dwarf.effective_temperature_k(0.03, 5.0) < brown_dwarf.effective_temperature_k(0.07, 5.0)


def test_a_young_brown_dwarf_is_no_hotter_than_a_late_m_dwarf():
    assert brown_dwarf.old_age_temperature_k(0.075, 0.05) > tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K
    assert brown_dwarf.effective_temperature_k(0.075, 0.05) == tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K
    assert brown_dwarf.spectral_type(brown_dwarf.effective_temperature_k(0.075, 0.05))[0] == "M"


def test_the_radius_is_about_jupiters_and_falls_with_mass():
    """The laws give 0.088 solar radii at 0.05 solar masses (Burrows et al. 2001's figure 3: older radii within about
    30 percent of Jupiter's over a broad mass range)."""
    assert brown_dwarf.radius_sol(0.05) == pytest.approx(0.088, abs=0.002)
    radii = [brown_dwarf.radius_sol(m) for m in (0.0125, 0.03, 0.05, 0.07, 0.0764)]
    assert radii == sorted(radii, reverse=True)
    jupiter = constants.JUPITER_RADIUS_KM * constants.KM_TO_M_FACTOR / constants.SOLAR_RADIUS_M
    assert all(0.7 * jupiter < r < 1.5 * jupiter for r in radii)


def test_the_luminosity_agrees_with_stefan_boltzmann_where_the_temperature_is_not_capped():
    """The radius is the 1 Gy one, so the two laws' `t^-0.02` in `R^2` leaves a few percent at 8 Gy."""
    for mass, age in ((0.03, 3.0), (0.05, 1.0), (0.07, 8.0)):
        found = brown_dwarf.state(mass, age)
        assert found["luminosity_sol"] == pytest.approx(brown_dwarf.luminosity_sol(mass, age), rel=0.05)


def test_the_types_start_where_kirkpatrick_et_al_put_them():
    assert brown_dwarf.type_temperature_k(0) == pytest.approx(2237.5)
    assert brown_dwarf.type_temperature_k(10) == pytest.approx(1254.8, abs=0.5)
    assert brown_dwarf.type_temperature_k(20) == pytest.approx(460.6, abs=0.5)
    # Y0 is inside the 371 to 510 K they measured.
    assert 371.0 < brown_dwarf.type_temperature_k(20) < 510.0
    # The fits meet where one range ends and the next begins.
    for number in (8.75, 14.75):
        below, above = brown_dwarf.type_temperature_k(number - 1e-9), brown_dwarf.type_temperature_k(number + 1e-9)
        assert below == pytest.approx(above, abs=0.5)


@pytest.mark.parametrize("temperature, expected", [
    (3000.0, "M"), (2300.0, "M"), (2200.0, "L0"), (1600.0, "L5"), (1300.0, "L8"), (1200.0, "T2"), (900.0, "T6"),
    (700.0, "T7"), (500.0, "T9"), (455.0, "Y0"), (300.0, "Y1"), (200.0, "Y2"),
])
def test_the_type_is_read_from_the_temperature(temperature, expected):
    assert brown_dwarf.spectral_type(temperature).startswith(expected)


def test_a_cooler_brown_dwarf_never_has_an_earlier_type():
    order = "MLTY"
    previous = None
    for temperature in range(3000, 150, -25):
        letter = brown_dwarf.spectral_type(float(temperature))[0]
        number = brown_dwarf.spectral_number(float(temperature))
        key = (order.index(letter), -1.0 if number is None else number)
        assert previous is None or key >= previous
        previous = key


def test_the_temperatures_of_a_census_spread_as_the_20_pc_census_does():
    """Table 15: the density per 150 K bin (1e-3 per pc3) at 600 to 1200 K is 2.80, 1.99, 1.72, 1.11 (complete there).
    The cooling laws, a constant birth rate over 10 Gy and the mass function give those shares within a factor 2."""
    low, high = brown_dwarf.mass_range_solar()
    census = {600: 2.80, 750: 1.99, 900: 1.72, 1050: 1.11, 1200: 1.95, 1350: 0.94, 1500: 0.81}
    counts = dict.fromkeys(census, 0)
    with draw.bound(252):
        for _ in range(60000):
            temperature = brown_dwarf.effective_temperature_k(brown_dwarf.sample_mass_solar(), draw.uniform(0.0, 10.0))
            for edge in counts:
                if edge <= temperature < edge + 150:
                    counts[edge] += 1
    shares = {edge: n / sum(counts.values()) for edge, n in counts.items()}
    observed = {edge: n / sum(census.values()) for edge, n in census.items()}
    for edge in (600, 750, 900, 1050):
        assert 0.5 <= shares[edge] / observed[edge] <= 2.0, (edge, shares[edge], observed[edge])


@pytest.mark.xfail(strict=True, reason="GEN.252: Burrows et al. 2001's old-age power law is too cool for the heaviest "
                   "brown dwarfs, so L dwarfs hotter than 1500 K come out at 0.2 to 0.5 of the census's share "
                   "(docs/design/brown-dwarf-state.md, 'What it does not do')")
def test_the_hot_l_dwarfs_reach_the_censuss_share():
    census = {1500: 0.81, 1650: 0.78, 1800: 0.50, 1950: 0.72, 1200: 1.95, 1350: 0.94}
    counts = dict.fromkeys(census, 0)
    with draw.bound(252):
        for _ in range(60000):
            temperature = brown_dwarf.effective_temperature_k(brown_dwarf.sample_mass_solar(), draw.uniform(0.0, 10.0))
            for edge in counts:
                if edge <= temperature < edge + 150:
                    counts[edge] += 1
    for edge in (1500, 1650, 1800, 1950):
        share = counts[edge] / sum(counts.values())
        assert 0.5 <= share / (census[edge] / sum(census.values())) <= 2.0


def test_the_free_floating_brown_dwarf_has_the_laws_mass_radius_and_a_type():
    with draw.bound(252):
        dwarfs = [RoguePlanet(SystemConfig(), mass_bin="brown-dwarf") for _ in range(30)]
    for dwarf in dwarfs:
        solar = dwarf.mass_kg / constants.SOLAR_MASS_TO_KG
        assert dwarf.radius_km == pytest.approx(brown_dwarf.radius_sol(solar) * constants.SOLAR_RADIUS_M
                                                / constants.KM_TO_M_FACTOR)
        assert dwarf.spectral_type == brown_dwarf.spectral_type(dwarf.effective_temperature_k)
        section = render_support.phenomenon_section(dwarf)
        assert f"{dwarf.spectral_type} brown dwarf" in render_support.prose_text(section)
        assert f"type {dwarf.spectral_type}" in render_support.facts_text(section)
    assert RoguePlanet(SystemConfig(), mass_bin="terrestrial").spectral_type is None


def test_a_brown_dwarf_member_has_a_type_and_the_radius_the_laws_give():
    with draw.bound(252):
        leaf = multistar._star_leaf(0.05, 4.0e9, draw)
    found = brown_dwarf.state(0.05, 4.0)
    assert leaf.kind == "brown-dwarf" and leaf.state == found["spectral_type"]
    assert leaf.radius_au == pytest.approx(found["radius_sol"] * constants.SOLAR_RADIUS_M / constants.AU_M)
    assert leaf.state[0] in "LTY"


def test_the_members_page_names_a_brown_dwarfs_type():
    tree = {"members": [{"member_no": 0, "kind": "main-sequence", "state": "V", "mass_solar": 0.6},
                        {"member_no": 1, "kind": "brown-dwarf", "state": "T5", "mass_solar": 0.05},
                        {"member_no": 2, "kind": "brown-dwarf", "state": None, "mass_solar": 0.05}],
            "nodes": []}
    labels = [row[1] for row in scatter_members_view.members_view(tree)["members"]]
    assert labels == ["Main-sequence star", "Brown dwarf (T5)", "Brown dwarf"]
