"""
The downloaded methane absorption table (GEN.949, AD12, `planetgen.datasets.methane_visible` and `physics.uv_cross_sections.methane_optical_depth`),
on a fixture the tests build offline: the units, the depth of the mean transmission, the cells it serves, the fallback and the refusal.
"""

import numpy as np
import pytest

from planetgen.datasets import fetch, methane_visible
from planetgen.physics import uv_cross_sections as uvx
from tests import datasets_support as support


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    support.reset_caches(tmp_path / "cache", monkeypatch)
    return tmp_path / "cache"


@pytest.fixture
def installed(monkeypatch):
    return support.install_methane(monkeypatch)


def _cell(edges, nm):
    return int(np.searchsorted(edges, nm, side="right") - 1)


def test_the_table_is_read_in_cross_section_units(installed):
    nm, sigma = methane_visible.samples(fetch.locate(methane_visible.DATASET))
    band = sigma[(nm > 720.0) & (nm < 734.0)]
    assert band.size and band == pytest.approx(40.0 / (1.0e5 * methane_visible.LOSCHMIDT_CM3))     # per km-amagat to cm2 per molecule
    assert 1.4e-23 < band[0] < 1.6e-23


def test_a_cell_is_the_depth_of_the_mean_transmission_not_the_mean_depth(installed):
    edges = uvx.standard_edges()
    located = fetch.locate(methane_visible.DATASET)
    nm, sigma = methane_visible.samples(located)
    column = 1.0e24
    depth = methane_visible.cell_optical_depth(located, edges, column)
    cell = _cell(edges, 725.0)                                      # the 720 to 730 cell sits fully in the 719 to 735 band
    inside = sigma[(nm >= edges[cell]) & (nm < edges[cell + 1])]
    assert depth[cell] == pytest.approx(column * inside.mean(), rel=1e-9)         # a uniform cell: the same depth either way
    mixed = _cell(edges, 735.5)                                     # the 730 to 740 cell is band for half, shoulder for half
    samples = sigma[(nm >= edges[mixed]) & (nm < edges[mixed + 1])]
    mean_depth = column * samples.mean()
    transmission_depth = -np.log(np.mean(np.exp(-column * samples)))
    assert depth[mixed] == pytest.approx(transmission_depth, rel=1e-9)
    assert depth[mixed] < mean_depth                                # the windows let the light through
    assert column * samples.min() < depth[mixed]                    # but not as clear as the shoulder alone


def test_a_thin_column_gives_the_mean_depth_and_a_zero_column_none(installed):
    edges = uvx.standard_edges()
    located = fetch.locate(methane_visible.DATASET)
    nm, sigma = methane_visible.samples(located)
    mixed = _cell(edges, 735.5)
    samples = sigma[(nm >= edges[mixed]) & (nm < edges[mixed + 1])]
    assert methane_visible.cell_optical_depth(located, edges, 1.0e18)[mixed] == pytest.approx(1.0e18 * samples.mean(), rel=1e-4)
    assert methane_visible.cell_optical_depth(located, edges, 0.0)[mixed] == 0.0
    with pytest.raises(ValueError):
        methane_visible.cell_optical_depth(located, edges, -1.0)


def test_the_table_serves_the_cells_between_400_and_1050_and_the_gas_bins_stand_elsewhere(installed):
    edges = uvx.standard_edges()
    column = 1.0e24
    depth = uvx.methane_optical_depth(column, edges)
    bins = column * uvx.band_mean_cross_section_cm2("CH4", edges)
    for nm in (150.5, 350.5, 1075.0, 1500.0):                       # under 400 nm and over 1050 nm: the gas table
        assert depth[_cell(edges, nm)] == bins[_cell(edges, nm)]
    for nm in (455.0, 725.0, 900.0):
        assert depth[_cell(edges, nm)] != bins[_cell(edges, nm)]
    assert np.all(np.isfinite(depth)) and np.all(depth >= 0.0)


def test_the_run_without_the_download_takes_the_bins_only_when_it_disregards_missing_data():
    edges = uvx.standard_edges()
    with pytest.raises(fetch.DataMissing):
        uvx.methane_optical_depth(1.0e24, edges)
    with fetch.disregard_missing_data():
        depth = uvx.methane_optical_depth(1.0e24, edges)
        assert ("methane_visible", "methane visible bands from the 35 bins of the gas table") in fetch.stand_ins_used()
    assert depth == pytest.approx(1.0e24 * uvx.band_mean_cross_section_cm2("CH4", edges))


def test_the_dataset_names_its_source_replacement_and_fallback():
    dataset = methane_visible.DATASET
    assert dataset.key == "methane_visible" and dataset.fill == "AD12" and dataset.fallback
    assert [piece.name for piece in dataset.pieces] == [methane_visible.TABLE] and dataset.pieces[0].url.startswith("https://")
