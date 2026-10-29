# tests/test_datasets_readers.py

"""The three datasets of GEN.548 read through their functions, on fixtures the tests build offline: a small MIST
archive, a white-dwarf sequence archive and an equation-of-state and outer-crust set. No test reaches the network."""

import numpy as np
import pytest

from planetgen.datasets import catalogue, fetch, isochrones, neutron_star_matter, white_dwarf_models
from planetgen.physics import constants
from tests import datasets_support as support


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    support.reset_caches(tmp_path / "cache", monkeypatch)
    return tmp_path / "cache"


# ---- MIST ----

@pytest.fixture
def mist(monkeypatch):
    return support.install_isochrones(monkeypatch)


def test_the_derived_file_keeps_only_main_sequence_rows(mist, cache):
    assert not (cache / "isochrones" / isochrones.ARCHIVE).exists()
    with np.load(cache / "isochrones" / isochrones.DERIVED) as data:
        assert list(data["fe_h"]) == [-1.0, 0.0] and data["rows_0"].shape == (120, 6)
        assert data["rows_0"][:, 1].max() <= 10 ** ((10 - 8.0) / 2.5) + 1e-3


def test_a_main_sequence_point_interpolates_in_mass_age_and_metallicity(mist):
    log_l, log_teff, log_r, present_mass = isochrones.main_sequence_point(0.5, 1.0, 0.0)
    assert log_l == pytest.approx(4 * np.log10(0.5), abs=0.02)
    assert log_teff == pytest.approx(3.7 + 0.5 * np.log10(0.5), abs=0.005)
    assert log_r == pytest.approx(0.8 * np.log10(0.5), abs=0.01)
    assert isochrones.main_sequence_point(0.5, 1.0, -0.5)[0] == pytest.approx(4 * np.log10(0.5) - 0.05, abs=0.02)
    between = isochrones.main_sequence_point(0.5, 0.5, 0.0)
    assert between is not None and between[0] == pytest.approx(4 * np.log10(0.5), abs=0.02)


def test_a_star_that_left_the_main_sequence_or_is_outside_the_ages_gets_none(mist):
    assert isochrones.main_sequence_point(2.0, 5.0, 0.0) is None     # turn-off at 5 Gy is 1.2 Msun
    assert isochrones.main_sequence_point(0.5, 1e-5, 0.0) is None    # younger than the grid
    assert isochrones.main_sequence_point(0.5, 200.0, 0.0) is None   # older than the grid
    assert isochrones.main_sequence_point(0.1, 1.0, 0.0) is None     # lighter than the grid's lowest mass


def test_the_arrival_on_the_main_sequence_is_read_from_the_contracting_rows(mist):
    """GEN.947: the age at which the heaviest star still in phase -1 has fallen to the mass, in the log of the mass."""
    assert isochrones.pre_main_sequence_arrival_yr(0.9, 0.0) == pytest.approx(1.0e8, rel=1e-6)
    assert isochrones.pre_main_sequence_arrival_yr(0.35, 0.0) == pytest.approx(1.0e9, rel=1e-6)
    halfway = (np.log10(0.5) - np.log10(0.9)) / (np.log10(0.35) - np.log10(0.9))
    assert isochrones.pre_main_sequence_arrival_yr(0.5, 0.0) == pytest.approx(10 ** (8.0 + halfway), rel=1e-6)
    assert isochrones.pre_main_sequence_arrival_yr(2.0, 0.0) is None     # heavier than any contracting star of the rows
    assert isochrones.pre_main_sequence_arrival_yr(0.05, 0.0) is None    # lighter than the oldest isochrone's


def test_the_lifetime_follows_the_turn_off_and_metallicity_is_clamped(mist):
    assert isochrones.main_sequence_lifetime_gy(10 ** (1 / 2.5), 0.0) == pytest.approx(1.0, rel=0.02)  # a grid age
    between = isochrones.main_sequence_lifetime_gy(2.0, 0.0)
    assert 1.0 < between < 10.0
    assert isochrones.main_sequence_lifetime_gy(0.5, 0.0) is None
    assert isochrones.main_sequence_point(0.5, 1.0, 3.0) == isochrones.main_sequence_point(0.5, 1.0, 0.0)
    assert isochrones.metallicities() == (-1.0, 0.0)
    assert isochrones.abundances(-0.5) == pytest.approx((0.27, 0.014))


def test_the_parser_refuses_a_file_that_is_not_a_mist_isochrone():
    with pytest.raises(ValueError):
        isochrones.parse_isochrone_file([b"# nothing\n", b"1 2 3\n"])


def test_the_isochrone_readers_say_so_when_the_archive_was_never_downloaded(cache):
    with pytest.raises(fetch.DataMissing, match=r"missing \(SC1,"):
        isochrones.main_sequence_point(1.0, 1.0, 0.0)


# ---- white dwarfs ----

@pytest.fixture
def dwarfs(monkeypatch):
    return support.install_white_dwarfs(monkeypatch)


def test_the_sequence_parser_reads_sixteen_numbers_per_model(dwarfs):
    first = white_dwarf_models.sequence(0.6)
    assert len(first.teff_k) == 6 and first.teff_k[0] == pytest.approx(100000.0)
    assert first.radius_cm[0] == pytest.approx(7e8)
    assert first.log_mass_fractions["hydrogen"][0] == pytest.approx(-4.0)
    assert first.log_mass_fractions["oxygen"][0] == pytest.approx(-0.3054)


def test_the_radius_interpolates_in_temperature_and_mass(dwarfs):
    one_solar_radius_cm = constants.SOLAR_RADIUS_M * 100.0
    exact = white_dwarf_models.radius_solar(0.6, 100000.0)
    assert exact == pytest.approx(7e8 / one_solar_radius_cm)
    middle = white_dwarf_models.radius_solar(0.55, 100000.0)
    assert white_dwarf_models.radius_solar(0.6, 100000.0) < middle < white_dwarf_models.radius_solar(0.5, 100000.0)
    assert white_dwarf_models.radius_solar(0.6, 1e9) == pytest.approx(exact)  # hotter than the sequence: its hot end
    assert white_dwarf_models.radius_solar(0.6, 1.0) < exact


def test_masses_outside_the_grid_get_none_and_the_masses_are_listed(dwarfs):
    assert white_dwarf_models.masses_sol() == (0.5, 0.6, 0.7)
    assert white_dwarf_models.radius_solar(0.4, 1e4) is None and white_dwarf_models.radius_solar(0.8, 1e4) is None


def test_the_layers_are_the_cases_and_the_core_is_split_evenly(dwarfs):
    thick = white_dwarf_models.layer_mass_fractions(0.6)
    thin = white_dwarf_models.layer_mass_fractions(0.6, "thin")
    assert thick["hydrogen"] == pytest.approx(1e-4) and thin["hydrogen"] == pytest.approx(1e-10)
    assert thick["helium"] == pytest.approx(1e-2) and thick["carbon"] == thick["oxygen"]
    assert sum(thick.values()) == pytest.approx(1.0, abs=0.02)


def test_a_sequence_file_with_a_broken_model_is_refused():
    with pytest.raises(ValueError):
        white_dwarf_models.parse_sequence("\n".join(["x"] * 5 + ["1.0 2.0 3.0"]), 0.6, "thick")


def test_the_white_dwarf_readers_say_so_when_nothing_was_downloaded(cache):
    with pytest.raises(fetch.DataMissing, match=r"missing \(SC2,"):
        white_dwarf_models.radius_solar(0.6, 1e4)


# ---- neutron stars ----

@pytest.fixture
def stars(monkeypatch):
    return support.install_neutron_stars(monkeypatch)


def test_all_equations_of_state_are_read_and_ascend_in_density(stars):
    tables = neutron_star_matter.equations_of_state()
    assert [table.index for table in tables] == list(range(1, 22))
    assert all(np.all(np.diff(table.baryon_density_fm3) > 0) for table in tables)
    assert tables[0].pressure_mev_fm3[-1] < tables[-1].pressure_mev_fm3[-1]


def test_an_equation_of_state_with_falling_density_is_refused():
    with pytest.raises(ValueError):
        neutron_star_matter.parse_equation_of_state(1, "# h\n2 1 1\n1 1 1\n")


def test_only_the_unmagnetised_compositions_are_read(stars):
    assert neutron_star_matter.crust_models() == ("modelA", "modelB")
    layers = neutron_star_matter.outer_crust("modelA")
    assert [(layer.z, layer.a) for layer in layers] == [(26, 56), (28, 62), (28, 64), (36, 86)]


def test_the_crust_nucleus_follows_the_density_and_ends_at_drip(stars):
    grams_per_cm3 = lambda n_fm3: n_fm3 / 1e-39 * 1.66053906660e-24
    assert neutron_star_matter.crust_nucleus("modelA", grams_per_cm3(1e-9)) == (26, 56)
    assert neutron_star_matter.crust_nucleus("modelA", grams_per_cm3(1e-7)) == (28, 62)
    assert neutron_star_matter.crust_nucleus("modelA", grams_per_cm3(5e-7)) == (28, 64)
    assert neutron_star_matter.crust_nucleus("modelB", grams_per_cm3(5e-5)) == (28, 66)
    assert neutron_star_matter.crust_nucleus("modelA", grams_per_cm3(1e-3)) is None
    assert neutron_star_matter.baryon_density_fm3(grams_per_cm3(1e-9)) == pytest.approx(1e-9)
    assert neutron_star_matter.mass_density_g_cm3(1e-9) == pytest.approx(grams_per_cm3(1e-9))
    assert neutron_star_matter.mass_density_g_cm3(neutron_star_matter.baryon_density_fm3(4e11)) == pytest.approx(4e11)


def test_the_neutron_star_readers_say_so_when_nothing_was_downloaded(cache):
    with pytest.raises(fetch.DataMissing, match=r"missing \(SC3,"):
        neutron_star_matter.crust_models()


# ---- the catalogue ----

def test_the_catalogue_lists_the_datasets_with_distinct_fills():
    assert [dataset.key for dataset in catalogue.ALL] == [
        "isochrones", "white_dwarf_models", "neutron_star_matter", "ozone_uv", "ir_ktables", "ir_cia", "harris_clusters", "noaa_proton_events",
        "baraffe_isochrones", "stellar_spectra", "mamajek_dwarf_sequence", "hhe_eos"]
    assert [dataset.fill for dataset in catalogue.ALL] == ["SC1", "SC2", "SC3", "AD9", "AD10", "AD11", "SC4", "SC5", "SC6", "SL1", "PS12", "GS1"]
    assert len(neutron_star_matter.DATASET.pieces) == 1 + neutron_star_matter.EOS_COUNT
    assert all(piece.url.startswith("https://") for dataset in catalogue.ALL for piece in dataset.pieces)
