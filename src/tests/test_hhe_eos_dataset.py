"""
The downloaded hydrogen-helium equation of state (GEN.733, GS1, `planetgen.datasets.hhe_eos`), on a fixture archive the tests build
offline: the derived file, the interpolation, the helium shift, the edges of the grid, the fallback and the refusal.
"""

import math

import numpy as np
import pytest

from planetgen.datasets import fetch, hhe_eos
from tests import datasets_support as support


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    support.reset_caches(tmp_path / "cache", monkeypatch)
    return tmp_path / "cache"


@pytest.fixture
def installed(monkeypatch):
    return support.install_hhe_eos(monkeypatch)


def test_the_archive_is_replaced_by_one_compact_file_with_the_mixture_and_the_volume_difference(installed, cache):
    assert not (cache / "hhe_eos" / hhe_eos.ARCHIVE).exists()
    with np.load(cache / "hhe_eos" / hhe_eos.DERIVED) as data:
        assert sorted(data.files) == ["ivl_slope", "log_p", "log_rho", "log_s", "log_t"]
        assert data["log_rho"].shape == (len(support.HHE_LOG_T), len(support.HHE_LOG_P))
        assert data["log_t"][0] == 2.0 and data["log_p"][-1] == pytest.approx(13.0)


def test_the_grid_reads_back_as_the_files_axes_and_flat_rows(installed):
    grid = hhe_eos.grid()
    assert (grid.n_t, grid.n_p) == (len(support.HHE_LOG_T), len(support.HHE_LOG_P))
    assert grid.log_t_range == pytest.approx((2.0, 8.0)) and grid.log_p_range_gpa == pytest.approx((-9.0, 13.0))
    assert len(grid.log_rho) == len(grid.ivl_slope) == grid.n_t * grid.n_p == grid.log_s.size


def test_a_node_is_read_back_exactly_and_between_nodes_the_log_density_is_bilinear(installed):
    grid = hhe_eos.grid()
    mu = 0.725 * support.HHE_HYDROGEN_MU + 0.275 * support.HHE_HELIUM_MU
    log_t, log_p = 3.75, 0.5
    expected = math.log10(support.hhe_fixture_density_g_cm3(log_t, log_p, mu))
    assert grid.sample(log_t, log_p)[0] == pytest.approx(expected, abs=1e-5)
    middle = grid.sample(log_t + 0.125, log_p + 0.125)[0]
    corners = [grid.sample(log_t + a, log_p + b)[0] for a in (0.0, 0.25) for b in (0.0, 0.25)]
    assert middle == pytest.approx(sum(corners) / 4.0, abs=1e-5)


def test_the_table_helium_share_gives_the_mixture_and_the_shift_follows_the_ideal_volume_law(installed):
    temperature_k, pressure_pa = 10 ** 3.75, 10 ** 9.5
    grid = hhe_eos.grid()
    log_rho, slope = grid.sample(3.75, 0.5)
    assert hhe_eos.density_kg_m3(temperature_k, pressure_pa) == pytest.approx(1000.0 * 10.0 ** log_rho, rel=1e-9)
    assert slope < 0.0                                                          # helium is the denser fluid
    pure_hydrogen = hhe_eos.density_kg_m3(temperature_k, pressure_pa, 0.0)
    pure_helium = hhe_eos.density_kg_m3(temperature_k, pressure_pa, 1.0)
    assert pure_hydrogen < hhe_eos.density_kg_m3(temperature_k, pressure_pa) < pure_helium
    volume = 10.0 ** -log_rho * 1e-3                                            # m3/kg of the mixture
    assert 1.0 / pure_helium == pytest.approx(volume * (1.0 + (1.0 - hhe_eos.TABLE_HELIUM_SHARE) * slope), rel=1e-9)
    assert 1.0 / pure_hydrogen == pytest.approx(volume * (1.0 - hhe_eos.TABLE_HELIUM_SHARE * slope), rel=1e-9)


def test_the_slope_is_held_to_its_bounds_where_the_pure_tables_are_unphysical(installed):
    grid = hhe_eos.grid()
    low, high = hhe_eos.IVL_SLOPE_RANGE
    assert all(low - 1e-6 <= value <= high + 1e-6 for value in grid.ivl_slope)


def test_a_point_off_the_grid_takes_the_value_on_its_edge(installed):
    grid = hhe_eos.grid()
    assert grid.sample(1.0, 0.5) == grid.sample(2.0, 0.5)
    assert grid.sample(9.0, 0.5) == grid.sample(8.0, 0.5)
    assert grid.sample(4.0, -12.0) == grid.sample(4.0, -9.0)
    assert grid.sample(4.0, 20.0) == grid.sample(4.0, 13.0)
    assert grid.sample(8.0, 13.0)[0] == grid.log_rho[-1]


def test_the_entropy_is_read_back_and_bilinear(installed):
    grid = hhe_eos.grid()
    assert grid.log_entropy(3.75, 0.5) == pytest.approx(support.hhe_fixture_log_entropy(3.75, 0.5), abs=1e-5)
    assert grid.log_entropy(3.875, 0.625) == pytest.approx(support.hhe_fixture_log_entropy(3.875, 0.625), abs=1e-5)
    assert hhe_eos.entropy_mj_kg_k(10 ** 3.75, 10 ** 9.5) == pytest.approx(10.0 ** support.hhe_fixture_log_entropy(3.75, 0.5), rel=1e-4)


def test_an_adiabat_follows_the_line_of_constant_entropy(installed):
    grid = hhe_eos.grid()
    log_s0 = support.hhe_fixture_log_entropy(3.0, -4.0)           # 1000 K at one bar
    nodes = grid.adiabat_log_t(log_s0)
    assert len(nodes) == grid.n_p
    for log_p in (-4.0, 0.0, 4.0):
        k = round((log_p - grid.log_p0) / grid.d_log_p)
        assert nodes[k] == pytest.approx(3.0 + 0.4 * (log_p + 4.0), abs=1e-3)     # T ~ P^0.4 in the fixture
    k = round((0.0 - grid.log_p0) / grid.d_log_p)
    assert grid.adiabat_log_t(log_s0 + 0.05)[k] > nodes[k]                          # a higher entropy is hotter


def test_the_one_bar_temperature_and_the_entropy_of_an_adiabat_invert_each_other(installed):
    grid = hhe_eos.grid()
    for temperature_k in (150.0, 300.0, 1000.0):
        log_s0 = grid.log_entropy_of_surface_temperature(temperature_k)
        assert grid.surface_temperature_k(log_s0) == pytest.approx(temperature_k, rel=0.02)
    cold = grid.log_entropy_of_surface_temperature(50.0)           # under the table's edge, carried up the molecular gas's adiabat
    assert cold < grid.log_entropy_of_surface_temperature(100.0)
    assert 40.0 < grid.surface_temperature_k(cold) < 60.0


def test_a_table_that_does_not_fill_its_grid_is_refused():
    with pytest.raises(ValueError, match="do not fill"):
        hhe_eos._read_table(b"2.0 -9.0 -3.0 0 0 0 0 0 0 0.3\n2.0 -8.0 -3.0 0 0 0 0 0 0 0.3\n3.0 -9.0 -3.0 0 0 0 0 0 0 0.3\n", "broken")
    with pytest.raises(ValueError, match="columns"):
        hhe_eos._read_table(b"2.0 -9.0 -3.0\n", "narrow")


def test_the_reader_says_so_when_nothing_was_downloaded_and_a_run_that_disregards_it_is_told_the_fallback(cache):
    with pytest.raises(fetch.DataMissing, match=r"missing \(GS1,"):
        hhe_eos.grid()
    located = fetch.locate(hhe_eos.DATASET)
    with pytest.raises(fetch.DataMissing):
        located.use_stand_in("anything")
    with fetch.disregard_missing_data():
        assert located.use_stand_in("giant layer pressures from the polytrope") is True
        assert ("hhe_eos", "giant layer pressures from the polytrope") in fetch.stand_ins_used()
