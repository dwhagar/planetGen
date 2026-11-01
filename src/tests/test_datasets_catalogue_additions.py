# tests/test_datasets_catalogue_additions.py

"""The three datasets that setup now downloads for the sourcing chain (GEN.297, GEN.747, GEN.936 and GEN.948): the Harris globular
cluster catalog, the NOAA solar proton events and the [R0254] isochrones. The parsers and the readers run on small fixtures
built offline (`tests/datasets_support.py`); the real files are checked only on a machine that has fetched them."""

import datetime
import math

import pytest

from planetgen.datasets import baraffe_isochrones, catalogue, fetch, harris_clusters, noaa_proton_events
from tests import datasets_support as support


@pytest.fixture(autouse=True)
def cache(request, tmp_path, monkeypatch):
    """An empty cache folder, except for the tests on the real files, which read the machine's own."""
    if request.node.get_closest_marker("real_datasets"):
        return None
    support.reset_caches(tmp_path / "cache", monkeypatch)
    return tmp_path / "cache"


# ---- the catalogue ----

def test_the_three_datasets_are_in_the_catalogue_with_a_registry_row_and_a_fallback():
    for module, fill in ((harris_clusters, "SC4"), (noaa_proton_events, "SC5"), (baraffe_isochrones, "SC6")):
        assert module.DATASET in catalogue.ALL
        assert module.DATASET.fill == fill and module.DATASET.fallback and module.DATASET.replaces
        assert all(piece.url.startswith("https://") for piece in module.DATASET.pieces)


def test_a_missing_file_is_refused_by_name(monkeypatch):
    for read in (harris_clusters.clusters, noaa_proton_events.events, baraffe_isochrones.isochrones):
        with pytest.raises(fetch.DataMissing):
            read()


# ---- Harris ----

@pytest.fixture
def harris(monkeypatch):
    support.install_harris(monkeypatch)
    return harris_clusters.clusters()


def test_a_cluster_joins_its_three_parts(harris):
    assert [cluster.name for cluster in harris] == ["NGC 104", "NGC 362", "Whiting 1"]
    tuc = harris[0]
    assert tuc.alias == "47 Tuc"
    assert tuc.ra_deg == pytest.approx((0 + 24 / 60 + 5.67 / 3600) * 15)
    assert tuc.dec_deg == pytest.approx(-(72 + 4 / 60 + 52.6 / 3600))
    assert (tuc.l_deg, tuc.b_deg, tuc.r_sun_kpc, tuc.r_gc_kpc) == (305.89, -44.89, 4.5, 7.4)
    assert (tuc.fe_h, tuc.fe_h_weight, tuc.m_v_total, tuc.spectral_type) == (-0.72, 10.0, -9.42, "G4")
    assert (tuc.v_r_km_s, tuc.sigma_v_km_s, tuc.concentration, tuc.core_radius_arcmin, tuc.half_light_radius_arcmin) == (
        -18.0, 11.0, 2.07, 0.36, 3.17)
    assert (tuc.log_core_relaxation_yr, tuc.log_median_relaxation_yr, tuc.core_collapsed) == (7.84, 9.55, False)


def test_a_collapsed_core_is_flagged_and_a_blank_cell_is_none(harris):
    assert harris[1].core_collapsed and harris[1].concentration == 1.76
    bare = harris[2]
    assert bare.alias == "" and bare.spectral_type == "" and not bare.core_collapsed
    assert bare.dec_deg == pytest.approx(-(3 + 15 / 60 + 10 / 3600))
    assert bare.fe_h == -0.7 and bare.ellipticity is None and bare.sigma_v_km_s is None and bare.mu_v is None
    assert bare.log_median_relaxation_yr == 7.93


def test_parts_that_list_different_clusters_are_refused():
    broken = support.HARRIS_TEXT.replace(" Whiting 1   -0.70", " Whiting 2   -0.70")
    with pytest.raises(ValueError, match="same clusters"):
        harris_clusters.parse(broken)


# ---- NOAA proton events ----

@pytest.fixture
def noaa(monkeypatch):
    support.install_noaa(monkeypatch)
    return noaa_proton_events.events()


def test_an_event_has_its_year_from_the_heading_and_its_peak_flux(noaa):
    first = noaa[0]
    assert first.start == datetime.datetime(1976, 4, 30, 21, 20)
    assert first.maximum == datetime.datetime(1976, 5, 1, 17, 0)
    assert first.peak_flux_pfu == 12 and first.flare == "X2/2B"
    assert len(noaa) == 5


def test_a_maximum_in_january_after_a_december_start_is_the_next_year(noaa):
    assert noaa[1].start == datetime.datetime(1989, 12, 31, 23, 10)
    assert noaa[1].maximum == datetime.datetime(1990, 1, 1, 4, 0)


def test_the_odd_spellings_of_the_newest_rows_are_read(noaa):
    assert noaa[2].maximum == datetime.datetime(2017, 8, 10, 4, 35)
    assert noaa[3].maximum == datetime.datetime(2017, 6, 19, 4, 50)
    assert noaa[4].maximum == datetime.datetime(2017, 9, 8, 0, 35) and noaa[4].peak_flux_pfu == 844


def test_a_listing_with_a_row_before_any_year_is_refused():
    with pytest.raises(ValueError, match="before any year"):
        noaa_proton_events.parse("<tr><td>Apr 30/2120<td>May 01/1700<td> 12")


# ---- Baraffe ----

@pytest.fixture
def baraffe(monkeypatch):
    support.install_baraffe(monkeypatch)
    return baraffe_isochrones.isochrones()


def test_each_age_block_is_an_isochrone_ascending_in_mass(baraffe):
    assert [iso.age_gyr for iso in baraffe] == [0.001, 1.0]
    assert list(baraffe[0].mass_solar) == [0.01, 0.1, 1.0]
    assert list(baraffe[1].radius_solar) == [0.1, 0.13, 0.9]


def test_a_grid_mass_at_a_grid_age_returns_the_row(baraffe):
    state = baraffe_isochrones.interpolate(0.1, 0.001)
    assert state == pytest.approx({"teff_k": 2900.0, "log_luminosity_solar": -1.0, "log_g": 3.4, "radius_solar": 1.1})


def test_between_ages_it_is_linear_in_log_age_and_between_masses_in_log_mass(baraffe):
    halfway = baraffe_isochrones.interpolate(0.1, math.sqrt(0.001 * 1.0))
    assert halfway["radius_solar"] == pytest.approx((1.1 + 0.13) / 2)
    assert halfway["log_g"] == pytest.approx((3.4 + 5.1) / 2)
    light = baraffe_isochrones.interpolate(0.03, 0.001)
    assert light["radius_solar"] == pytest.approx(0.4 + (1.1 - 0.4) * math.log(3) / math.log(10))


def test_outside_the_grid_it_answers_none(baraffe):
    assert baraffe_isochrones.interpolate(0.1, 1e-5) is None, "younger than the grid"
    assert baraffe_isochrones.interpolate(0.1, 20.0) is None, "older than the grid"
    assert baraffe_isochrones.interpolate(5.0, 0.001) is None, "heavier than the grid"
    assert baraffe_isochrones.interpolate(0.02, 0.5) is None, "lighter than the older block reaches"
    assert baraffe_isochrones.interpolate(-1.0, 0.5) is None


def test_a_row_of_the_wrong_width_is_refused():
    with pytest.raises(ValueError, match="nine columns"):
        baraffe_isochrones.parse("!  t (Gyr) =   0.0010\n 0.010 2400. -2.40\n")


# ---- the real files, where this machine has them ----

REAL = pytest.mark.skipif(any(fetch.locate(module.DATASET).missing for module in (harris_clusters, noaa_proton_events, baraffe_isochrones)),
                          reason="the Harris, NOAA and Baraffe files are not in this machine's cache folder")


@REAL
@pytest.mark.real_datasets
def test_the_real_harris_catalog_has_157_clusters_and_omega_centauri_is_the_biggest():
    clusters = harris_clusters.clusters()
    assert len(clusters) == 157
    brightest = min((c for c in clusters if c.m_v_total is not None), key=lambda c: c.m_v_total)
    assert brightest.alias == "omega Cen"
    assert 25 <= sum(c.core_collapsed for c in clusters) <= 35


@REAL
@pytest.mark.real_datasets
def test_the_real_proton_events_run_from_1976_to_2017_and_the_biggest_is_over_forty_thousand_pfu():
    events = noaa_proton_events.events()
    assert len(events) == 265 and events[0].start.year == 1976 and events[-1].start.year == 2017
    biggest = max(events, key=lambda e: e.peak_flux_pfu)
    assert biggest.peak_flux_pfu > 40000
    assert all(event.peak_flux_pfu >= 10 and event.maximum >= event.start for event in events)


@REAL
@pytest.mark.real_datasets
def test_the_real_baraffe_grid_gives_a_solar_mass_star_its_pre_main_sequence_radius():
    grid = baraffe_isochrones.isochrones()
    assert len(grid) == 30 and grid[0].age_gyr == 0.0005 and grid[-1].age_gyr == 10.0
    young, old = baraffe_isochrones.interpolate(1.0, 0.001), baraffe_isochrones.interpolate(1.0, 5.0)
    assert young["radius_solar"] > 2.0 and 0.9 < old["radius_solar"] < 1.1
