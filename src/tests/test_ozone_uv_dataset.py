"""
The downloaded ozone cross sections (GEN.443, AD9, `planetgen.datasets.ozone_uv` and `physics.uv_cross_sections`), on a fixture
archive the tests build offline: the derived file, the cell means, the temperature law, the fallback and the refusal.
"""

import numpy as np
import pytest

from planetgen.datasets import fetch, ozone_uv
from planetgen.physics import uv_cross_sections as uvx
from tests import datasets_support as support


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    support.reset_caches(tmp_path / "cache", monkeypatch)
    return tmp_path / "cache"


@pytest.fixture
def installed(monkeypatch):
    return support.install_ozone(monkeypatch)


def _index(edges, nm):
    return int(np.searchsorted(edges, nm, side="right") - 1)


def test_the_archive_is_replaced_by_one_compact_file_with_the_four_tables(installed, cache):
    assert not (cache / "ozone_uv" / ozone_uv.ARCHIVE).exists()
    with np.load(cache / "ozone_uv" / ozone_uv.DERIVED) as data:
        assert sorted(data.files) == sorted(f"{kind}_{key}" for kind in ("nm", "sigma") for key in ozone_uv.MEMBERS)
        assert data["nm_218"][0] == 195.0 and data["nm_218"][-1] == 345.0 and data["nm_vis"][-1] == 830.0


def test_cells_are_means_of_the_table_and_the_gas_bins_fill_what_it_does_not_reach(installed):
    edges = uvx.standard_edges()
    bins = uvx.band_mean_cross_section_cm2("O3", edges)
    sigma = uvx.ozone_cross_section_cm2(edges, 295.0)
    hartley = _index(edges, 255.5)
    assert sigma[hartley] == pytest.approx(support.ozone_fixture_sigma(255.5, 295.0), rel=0.01)
    assert sigma[_index(edges, 150.5)] == bins[_index(edges, 150.5)]         # under 195 nm: the gas table
    assert sigma[_index(edges, 1500.0)] == bins[_index(edges, 1500.0)]       # over 830 nm: the gas table
    chappuis = _index(edges, 605.0)
    assert sigma[chappuis] == pytest.approx(support.ozone_fixture_sigma(605.0, 295.0), rel=0.05)
    assert np.all(np.isfinite(sigma)) and np.all(sigma >= 0.0)


def test_the_temperature_is_linear_between_the_three_tables_and_held_outside_them(installed):
    edges = uvx.standard_edges()
    cell = _index(edges, 300.5)
    cold, warm, hot = (uvx.ozone_cross_section_cm2(edges, t)[cell] for t in (218.0, 243.0, 295.0))
    assert cold < warm < hot
    assert uvx.ozone_cross_section_cm2(edges, 230.5)[cell] == pytest.approx(cold + (warm - cold) * (230.5 - 218.0) / 25.0, rel=1e-6)
    assert uvx.ozone_cross_section_cm2(edges, 270.0)[cell] == pytest.approx(warm + (hot - warm) * (270.0 - 243.0) / 52.0, rel=1e-6)
    assert uvx.ozone_cross_section_cm2(edges, 150.0)[cell] == cold
    assert uvx.ozone_cross_section_cm2(edges, 500.0)[cell] == hot
    assert uvx.ozone_cross_section_cm2(edges, None)[cell] == pytest.approx(uvx.ozone_cross_section_cm2(edges, 220.0)[cell])


def test_without_the_download_a_run_is_refused_unless_it_disregards_missing_data(cache):
    edges = uvx.standard_edges()
    with pytest.raises(fetch.DataMissing):
        uvx.ozone_cross_section_cm2(edges)
    with fetch.disregard_missing_data():
        taken = uvx.ozone_cross_section_cm2(edges)
        assert taken == pytest.approx(uvx.band_mean_cross_section_cm2("O3", edges))
        assert ("ozone_uv", "ozone cross sections from the 20 nm bins of the gas table") in fetch.stand_ins_used()


def test_the_gas_table_bins_land_in_the_right_cells():
    edges = uvx.standard_edges()
    o3 = uvx.band_mean_cross_section_cm2("O3", edges)
    assert o3[_index(edges, 250.5)] == pytest.approx(1.1e-17, rel=0.1)         # the 240 to 260 nm bin
    assert o3[_index(edges, 1500.0)] == 0.0                                     # nothing beyond 1100 nm
    assert np.all(uvx.band_mean_cross_section_cm2("N2", edges) == 0.0)
    assert uvx.cross_section_cm2("SO2", edges)[_index(edges, 285.0)] == pytest.approx(7.5e-19, rel=0.1)
