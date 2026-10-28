"""
The downloaded infrared opacity data of GEN.715 (`planetgen.datasets.ir_ktables`, `planetgen.datasets.hitran_cia`) on fixture files the
tests build offline: the FITS reader, the reduction of a k-distribution, the CIA parser and band means, and the download step's derived files.
"""

import gzip
import io

import numpy as np
import pytest

from planetgen.datasets import catalogue, fetch, hitran_cia, ir_ktables
from tests import ir_datasets_support as support


@pytest.fixture
def cache(tmp_path, monkeypatch):
    support.install(monkeypatch, tmp_path / "cache")
    return tmp_path / "cache"


# ---- the reduction ----

def test_the_g_grid_is_a_distribution_with_a_finer_top():
    g, weights, edges = ir_ktables.g_grid()
    assert len(g) == sum(ir_ktables.G_POINTS) == 10
    assert weights.sum() == pytest.approx(1.0)
    assert edges[0] == 0.0 and edges[-1] == pytest.approx(1.0) and np.all(np.diff(g) > 0)
    split = np.searchsorted(edges, ir_ktables.G_SPLIT)
    assert edges[split] == pytest.approx(ir_ktables.G_SPLIT)
    assert weights[ir_ktables.G_POINTS[0]:].sum() == pytest.approx(1.0 - ir_ktables.G_SPLIT)


def test_a_constant_opacity_reduces_to_itself():
    _g, _w, edges = ir_ktables.g_grid()
    k = np.full((2, 3, 40), 4.0e-20)
    reduced = ir_ktables.reduce_distribution(k, np.full(40, 1.0 / 40), edges)
    assert reduced.shape == (2, 3, 10)
    assert reduced == pytest.approx(4.0e-20)


def test_reduction_keeps_the_mean_and_sorts_the_values():
    """Half the spectrum at k = 1 and half at k = 3: the lower five g intervals read 1 up to the cut, the rest 3, and the weighted mean is 2."""
    _g, weights, edges = ir_ktables.g_grid()
    k = np.concatenate([np.full(50, 3.0), np.full(50, 1.0)])[None, :]       # unsorted on purpose
    reduced = ir_ktables.reduce_distribution(k, np.full(100, 0.01), edges)[0]
    assert (reduced * weights).sum() == pytest.approx(2.0)
    assert np.all(np.diff(reduced) >= -1e-12)
    assert reduced[0] == pytest.approx(1.0) and reduced[-1] == pytest.approx(3.0)


def test_unequal_weights_move_the_cut():
    _g, weights, edges = ir_ktables.g_grid()
    k = np.array([[1.0, 5.0]])
    reduced = ir_ktables.reduce_distribution(k, np.array([0.8, 0.2]), edges)[0]
    assert (reduced * weights).sum() == pytest.approx(0.8 * 1.0 + 0.2 * 5.0)
    assert reduced[0] == pytest.approx(1.0) and reduced[-1] == pytest.approx(5.0)          # the weak value fills g up to 0.8


# ---- the FITS reader ----

def test_the_reader_returns_the_axes_and_the_data_of_a_table():
    stream = io.BytesIO(gzip.decompress(support.ktable_gz("CO2")))
    data, temperature, pressure, wavelength = ir_ktables.read_ktable(stream)
    assert data.shape == (len(support.PRESSURES_BAR), len(support.TEMPERATURES_K), 20, len(support.wavelengths_cm()))
    assert temperature == pytest.approx(support.TEMPERATURES_K) and pressure == pytest.approx(support.PRESSURES_BAR)
    nodes, _w = np.polynomial.legendre.leggauss(20)
    expected = support.line_k("CO2", pressure, temperature, (nodes + 1.0) / 2.0, wavelength)
    assert data == pytest.approx(expected.astype(np.float32), rel=1e-6)


def test_a_table_that_ends_early_is_refused():
    raw = gzip.decompress(support.ktable_gz("CO2"))
    with pytest.raises(ValueError):
        ir_ktables.read_ktable(io.BytesIO(raw[:20000]))
    with pytest.raises(ValueError):
        ir_ktables.read_ktable(io.BytesIO(b"not a fits file" * 400))


def test_the_reduced_table_has_the_bands_the_temperatures_and_the_line():
    stream = io.BytesIO(gzip.decompress(support.ktable_gz("CO2")))
    k, temperature, _pressure = ir_ktables.reduce_table(*ir_ktables.read_ktable(stream))
    assert k.shape == (129, 3, len(support.PRESSURES_BAR), 10)
    assert list(temperature) == [100.0, 300.0, 1000.0]          # 2500 K is above T_MAX_K
    edges = ir_ktables.BAND_EDGES_CM1
    centre = (edges[1:] + edges[:-1]) / 2.0
    line = int(np.argmin(np.abs(centre - support.LINE_CENTRE_CM1["CO2"])))
    far = int(np.argmin(np.abs(centre - 3000.0)))
    assert k[line, 1, -1, -1] > 1e4 * k[far, 1, -1, -1]
    assert np.all(k[..., 1:] >= k[..., :-1] * (1.0 - 1e-6))          # g sorted by strength
    assert k[line, 1, -1, -1] > k[line, 1, 0, -1]                    # the fixture's pressure law


# ---- collision-induced absorption ----

def test_a_cia_file_is_parsed_into_blocks_and_noise_below_zero_is_dropped():
    text = support.cia_text("N2-N2").decode()
    blocks = hitran_cia.parse(text)
    assert [temperature for temperature, _data in blocks] == [200.0, 300.0]
    assert blocks[0][1].shape == (449, 2)
    lines = text.split("\n")
    lines[3] = f" {10.0 * 3:10.4f}  {-1.0e-48:10.3E}"
    assert all(np.all(data[:, 1] >= 0.0) for _t, data in hitran_cia.parse("\n".join(lines)))
    with pytest.raises(ValueError):
        hitran_cia.parse("")
    with pytest.raises(ValueError):
        hitran_cia.parse("   N2-N2  not a header\n")


def test_band_means_are_means_over_the_covered_part_and_zero_where_nothing_covers():
    nu = np.arange(100.0, 1001.0, 1.0)
    blocks = [(200.0, np.column_stack([nu, np.full_like(nu, 4.0e-46)])), (300.0, np.column_stack([nu, np.full_like(nu, 2.0e-46)]))]
    edges = np.array([0.0, 150.0, 400.0, 1000.0, 2000.0])
    temperatures, k = hitran_cia.band_means(blocks, edges)
    assert list(temperatures) == [200.0, 300.0]
    assert k[0, 1] == pytest.approx(4.0e-46) and k[1, 2] == pytest.approx(2.0e-46)       # fully covered bands
    assert k[0, 0] == 0.0                                                                  # under half covered: not used
    assert k[0, 3] == 0.0                                                                  # beyond the file
    assert k[0, 1] > k[1, 1]


def test_blocks_at_different_temperatures_and_ranges_are_filled_per_band():
    """One band is measured at 200 and 300 K, another at 250 K only: the second is held at that value at every temperature."""
    low = np.arange(10.0, 100.0, 1.0)
    high = np.arange(500.0, 601.0, 1.0)
    blocks = [(200.0, np.column_stack([low, np.full_like(low, 1.0)])), (300.0, np.column_stack([low, np.full_like(low, 3.0)])),
              (250.0, np.column_stack([high, np.full_like(high, 7.0)]))]
    temperatures, k = hitran_cia.band_means(blocks, np.array([0.0, 100.0, 450.0, 650.0]))
    assert list(temperatures) == [200.0, 250.0, 300.0]
    assert k[:, 0] == pytest.approx([1.0, 2.0, 3.0])                                       # linear in temperature between the two blocks that cover it
    assert k[:, 2] == pytest.approx([7.0] * 3)                                       # the mean over the half the file covers


# ---- the download step ----

def test_both_datasets_are_in_the_catalogue_with_registry_ids():
    keys = {dataset.key: dataset for dataset in catalogue.ALL}
    assert keys["ir_ktables"].fill == "AD10" and keys["ir_cia"].fill == "AD11"
    assert all(piece.transient for dataset in (ir_ktables.DATASET, hitran_cia.DATASET) for piece in dataset.pieces)
    assert {piece.url.split("/")[2] for piece in ir_ktables.DATASET.pieces} == {"exomol.com"}
    assert {piece.url.split("/")[2] for piece in hitran_cia.DATASET.pieces} == {"hitran.org"}


def test_the_archives_are_replaced_by_one_compact_file_each(cache):
    assert sorted(path.name for path in (cache / "ir_ktables").iterdir()) == [".derived.json", "ir_ktables.npz"]
    assert sorted(path.name for path in (cache / "ir_cia").iterdir()) == [".derived.json", "ir_cia.npz"]
    with np.load(cache / "ir_ktables" / "ir_ktables.npz") as data:
        assert sorted(data.files) == sorted(["g", "weights", "edges_cm1", "temperature_k", "pressure_bar", "k_H2O", "k_CO2", "k_CH4"])
        assert data["k_H2O"].shape == (129, 3, 5, 10)
    with np.load(cache / "ir_cia" / "ir_cia.npz") as data:
        assert sorted(data.files) == sorted(f"{kind}_{pair}" for kind in ("temperature", "k") for pair in hitran_cia.FILES)
        assert data["k_N2-N2"].shape[1] == 130


def test_the_pair_keys_are_sorted_formulas():
    assert hitran_cia.pair_key("N2", "CH4") == "CH4-N2" and hitran_cia.pair_key("O2", "N2") == "N2-O2"
    assert all(key == hitran_cia.pair_key(*key.split("-")) for key in hitran_cia.FILES)


def test_a_wrong_hash_keeps_the_dataset_missing(tmp_path, monkeypatch):
    files = support.payloads()
    name = ir_ktables.GASES["CO2"][0]
    files[name] = files[name][:-10]
    with pytest.raises(AssertionError):
        support.install(monkeypatch, tmp_path / "cache", {**files, name: b"x" + files[name]})
    assert fetch.locate(ir_ktables.DATASET).missing
