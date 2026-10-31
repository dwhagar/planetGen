"""
The infrared optical depth of a stack of layers, band by band and g point by g point (GEN.715, `planetgen.physics.line_opacity`), on the
offline fixture tables of `tests/ir_datasets_support.py`: the table lookup and its interpolation, the random overlap of two gases, the
continuum, the collision-induced pairs, the stand-in path of a gas without a table and the refusal when the data is missing. The same module
against measured numbers on the real tables is `test_radiative_convective_real.py`.
"""

import math

import numpy as np
import pytest

from planetgen.datasets import fetch, ir_ktables
from planetgen.physics import constants, gases
from planetgen.physics import line_opacity as ir
from planetgen.physics.air_column import AVOGADRO
from tests import ir_datasets_support as support


@pytest.fixture(autouse=True)
def tables(tmp_path, monkeypatch):
    """The fixture k-tables and pair files in an empty cache folder."""
    support.install(monkeypatch, tmp_path / "cache")
    yield
    fetch._STAND_INS.clear()


def _molecules_cm2(pressure_pa, gravity=9.81, molar_mass=0.0289):
    """The molecules above one square centimetre of a layer with `pressure_pa` of air on it."""
    return pressure_pa / gravity / molar_mass * AVOGADRO * 1.0e-4


def _depth(mix, temperature, pressure_pa, water_fraction=0.0, column_pa=None, gravity=9.81):
    """The depth (band, g) of one layer of `mix` holding `column_pa` (default: the pressure itself) of air."""
    prepared = ir.prepare(mix)
    one = lambda value: np.array([float(value)])
    return ir.layer_depth(prepared, one(temperature), one(pressure_pa), one(water_fraction),
                          one(_molecules_cm2(pressure_pa if column_pa is None else column_pa, gravity)))[0]


def _mean(depth):
    return (depth * ir.G_WEIGHTS).sum(-1)


# ---- the Planck function ----

def test_planck_fraction_matches_a_quadrature():
    xs = np.linspace(1e-9, 60.0, 600001)
    density = xs ** 3 / np.expm1(xs)
    cumulative = np.cumsum((density[1:] + density[:-1]) / 2.0 * np.diff(xs))
    cumulative = np.concatenate([[0.0], cumulative]) * 15.0 / math.pi ** 4
    for x in (0.1, 0.5, 1.0, 1.99, 2.01, 3.0, 5.0, 10.0):
        index = int(np.searchsorted(xs, x))
        assert float(ir.planck_fraction_below(x)) == pytest.approx(float(cumulative[index]), abs=2e-5)
    assert float(ir.planck_fraction_below(0.0)) == 0.0
    assert float(ir.planck_fraction_below(80.0)) == pytest.approx(1.0, abs=1e-12)


def test_the_bands_cover_the_spectrum_a_planet_radiates_in():
    assert ir.BANDS == 130 and ir.EDGES_CM1[0] == 0.0 and ir.EDGES_CM1[-1] >= 12000.0
    assert np.all(np.diff(ir.EDGES_CM1) > 0.0)
    assert ir.G_WEIGHTS.sum() == pytest.approx(1.0) and len(ir.G_POINTS) == ir.GRID_POINTS == 10


@pytest.mark.parametrize("temperature", [60.0, 288.0, 735.0, 2500.0])
def test_band_emission_adds_up_to_stefan_boltzmann(temperature):
    power = ir.band_emission(np.array([temperature]))[0]
    assert power.shape == (ir.BANDS,) and np.all(power >= 0.0)
    lost = 1.0 - ir.planck_fraction_below(ir.SECOND_RADIATION_CM_K * ir.EDGES_CM1[-1] / temperature)
    assert power.sum() == pytest.approx(constants.STEFAN_BOLTZMANN_CONSTANT * temperature ** 4 * (1.0 - lost), rel=2e-3)


def test_band_emission_peaks_where_wien_says_and_takes_any_array_shape():
    power = ir.band_emission(np.array([[288.0, 300.0], [250.0, 3000.0]]))
    assert power.shape == (2, 2, ir.BANDS)
    peak = ir.CENTRES_CM1[int(np.argmax(power[0, 0] / np.diff(ir.EDGES_CM1)))]
    assert peak == pytest.approx(2.821 * 288.0 / ir.SECOND_RADIATION_CM_K, rel=0.25)          # the exitance per wavenumber peaks at h c nu / k T = 2.821
    assert ir.band_emission(np.array([1.0]))[0].sum() == pytest.approx(ir.band_emission(np.array([5.0]))[0].sum())   # the ends hold


# ---- the k-distributions ----

def test_a_distribution_keeps_its_mean_when_it_is_cut_again():
    rng = np.random.default_rng(3)
    k = rng.lognormal(0.0, 2.0, (4, 7, 57))
    weights = rng.random(57)
    weights /= weights.sum()
    cut = ir.reduce_distribution(k, weights, ir.G_EDGES)
    assert cut.shape == (4, 7, 10)
    assert (cut * ir.G_WEIGHTS).sum(-1) == pytest.approx((k * weights).sum(-1))
    assert np.all(np.diff(cut, axis=-1) >= -1e-12 * cut.max())


def test_two_flat_gases_overlap_into_their_sum():
    first, second = np.full((3, ir.GRID_POINTS), 2.0), np.full((3, ir.GRID_POINTS), 5.0)
    assert ir.resort(first, second) == pytest.approx(7.0)
    assert ir.overlap(first, second) == pytest.approx(7.0)


def test_random_overlap_keeps_the_mean_and_is_between_the_larger_gas_and_the_sum_of_the_lines():
    """A weak gas that fills only the bottom of g and a strong one: the sum of the means is the mean of the overlap, and the strongest line
    of the pair is no stronger than the two strongest lines added, nor weaker than the stronger of the two."""
    g = np.linspace(0.0, 1.0, ir.GRID_POINTS)
    first = (10.0 ** (3.0 * g - 2.0))[None, :]
    second = (10.0 ** (4.0 * g - 3.0))[None, :]
    both = ir.overlap(first, second)
    assert _mean(both) == pytest.approx(_mean(first) + _mean(second), rel=1e-6)
    assert np.all(np.diff(both, axis=-1) >= -1e-12)
    assert both[0, -1] <= first[0, -1] + second[0, -1] + 1e-9
    assert np.all(both >= np.maximum(first, second) - 1e-9 * both.max())


def test_a_thin_gas_adds_only_its_mean():
    first = (10.0 ** np.linspace(-1.0, 2.0, ir.GRID_POINTS))[None, :]
    thin = np.full_like(first, 1.0e-3)
    assert ir.overlap(first, thin) == pytest.approx(first + 1.0e-3)
    assert ir.overlap(thin, first) == pytest.approx(first + 1.0e-3)


# ---- the table lookup ----

def _table(gas):
    return ir.prepare({"N2": 1.0e5, gas: 1.0e3}).tables[gas]


@pytest.mark.parametrize("gas", ["H2O", "CO2", "CH4"])
def test_a_lookup_at_a_grid_point_is_the_table(gas):
    prepared = ir.prepare({"N2": 1.0e5, gas: 1.0e3})
    table = prepared.tables[gas]
    temperatures = prepared.table_temperatures
    pressures_bar = np.exp(table.log_pressure_first + table.log_pressure_step * np.arange(table.n_pressures))
    for it in (0, len(temperatures) - 1):
        for ip in (0, 2, table.n_pressures - 1):
            k = ir._lookup(table, temperatures, np.array([temperatures[it]]), np.array([pressures_bar[ip]]))[0]
            assert k.shape == (ir.BANDS, ir.GRID_POINTS)
            assert k[1:] == pytest.approx(10.0 ** table.logk[:, it, ip, :].astype(float), rel=1e-4)
            assert k[0] == pytest.approx(k[1])                  # the 0 to 200 cm-1 band takes the first band's distribution (stand-in IR5)


def test_a_lookup_interpolates_log_k_linearly_in_inverse_temperature_and_in_ln_pressure():
    prepared = ir.prepare({"N2": 1.0e5, "CO2": 1.0e3})
    table, temperatures = prepared.tables["CO2"], prepared.table_temperatures
    inverse = 0.5 * (1.0 / temperatures[0] + 1.0 / temperatures[1])
    p_low = math.exp(table.log_pressure_first)
    p_mid = p_low * math.exp(0.5 * table.log_pressure_step)
    k = ir._lookup(table, temperatures, np.array([1.0 / inverse]), np.array([p_mid]))[0, 1:]
    corners = table.logk[:, :2, :2, :].astype(float)
    assert np.log10(k) == pytest.approx(corners.mean(axis=(1, 2)), abs=1e-4)


def test_a_lookup_off_the_table_holds_the_edge():
    prepared = ir.prepare({"N2": 1.0e5, "H2O": 1.0e3})
    table, temperatures = prepared.tables["H2O"], prepared.table_temperatures
    edge = ir._lookup(table, temperatures, np.array([temperatures[-1]]), np.array([100.0]))
    off = ir._lookup(table, temperatures, np.array([temperatures[-1] + 500.0]), np.array([1000.0]))
    low = ir._lookup(table, temperatures, np.array([20.0]), np.array([1.0e-9]))
    assert off == pytest.approx(edge) and np.all(np.isfinite(low)) and np.all(low > 0.0)


# ---- the depth of a layer ----

def test_a_layer_of_a_line_gas_is_opaque_in_its_band_and_clear_elsewhere():
    depth = _depth({"N2": 1.0e5, "CO2": 1.0e3}, 300.0, 1.0e5)
    mean = _mean(depth)
    line = int(np.argmin(np.abs(ir.CENTRES_CM1 - support.LINE_CENTRE_CM1["CO2"])))
    far = int(np.argmin(np.abs(ir.CENTRES_CM1 - 3000.0)))
    assert depth.shape == (ir.BANDS, ir.GRID_POINTS) and np.all(np.isfinite(depth)) and np.all(depth >= 0.0)
    assert np.all(np.diff(depth, axis=-1) >= -1e-9 * depth.max())
    assert mean[line] > 1.0e3 * mean[far]


def test_more_gas_in_the_layer_is_more_depth_and_is_exactly_linear_for_one_gas_without_pairs():
    mix = {"Ar": 9.0e4, "CO2": 1.0e3}
    one = _depth(mix, 300.0, 1.0e5, column_pa=1.0e4)
    two = _depth(mix, 300.0, 1.0e5, column_pa=2.0e4)
    assert two == pytest.approx(2.0 * one, rel=1e-9)
    assert _mean(_depth({"Ar": 9.0e4, "CO2": 4.0e3}, 300.0, 1.0e5)).sum() > 3.0 * _mean(one).sum() * (1.0e5 / 1.0e4)          # four times the gas in the layer (the pairs more)


def test_a_gas_the_air_does_not_hold_adds_nothing():
    clear = ir.layer_depth(ir.prepare({"Ar": 1.0e5}), np.array([300.0]), np.array([1.0e5]), np.array([0.0]), np.array([1.0e22]))
    assert not ir.prepare({"Ar": 1.0e5}).absorbers
    assert clear.shape == (1, ir.BANDS, ir.GRID_POINTS) and not clear.any()


def test_the_water_fraction_of_a_layer_is_water_and_the_rest_is_the_dry_mixture():
    prepared = ir.prepare({"N2": 7.0e4, "CO2": 1.0e4, "H2O": 1.0e4})
    assert prepared.has_water and prepared.dry["CO2"] == pytest.approx(1.0e4 / 8.0e4)
    assert sum(prepared.dry.values()) == pytest.approx(1.0) and set(prepared.tabled) == {"H2O", "CO2"}
    dry_layer = ir.layer_depth(prepared, np.array([300.0]), np.array([1.0e5]), np.array([0.0]), np.array([1.0e22]))
    wet_layer = ir.layer_depth(prepared, np.array([300.0]), np.array([1.0e5]), np.array([0.2]), np.array([1.0e22]))
    h2o = int(np.argmin(np.abs(ir.CENTRES_CM1 - support.LINE_CENTRE_CM1["H2O"])))
    co2 = int(np.argmin(np.abs(ir.CENTRES_CM1 - support.LINE_CENTRE_CM1["CO2"])))
    assert _mean(wet_layer[0])[h2o] > 100.0 * _mean(dry_layer[0])[h2o]
    assert _mean(wet_layer[0])[co2] < _mean(dry_layer[0])[co2]          # the vapour displaces dry gas


def test_adding_carbon_dioxide_to_a_wet_air_adds_to_its_band_and_takes_nothing_from_the_water_band():
    args = (np.array([300.0]), np.array([1.0e5]), np.array([0.2]), np.array([3.0e22]))
    wet = _mean(ir.layer_depth(ir.prepare({"Ar": 8.0e4, "H2O": 2.0e4}), *args)[0])
    both = _mean(ir.layer_depth(ir.prepare({"Ar": 7.0e4, "CO2": 1.0e4, "H2O": 2.0e4}), *args)[0])
    carbon = int(np.argmin(np.abs(ir.CENTRES_CM1 - support.LINE_CENTRE_CM1["CO2"])))
    assert both[carbon] > 100.0 * wet[carbon]
    assert np.all(both >= wet * (1.0 - 1.0e-9))


# ---- the water vapour continuum ----

def _continuum(fraction, temperature, pressure_pa=1.0e5):
    column = _molecules_cm2(pressure_pa)
    return ir._continuum_depth(np.array([temperature]), np.array([pressure_pa]), np.array([fraction]), np.array([column]))[0]


def test_the_water_continuum_closes_the_window_of_a_humid_air_and_grows_faster_than_the_vapour():
    window = int(np.argmin(np.abs(ir.CENTRES_CM1 - 1000.0)))
    assert _continuum(2.0e-3, 288.0)[window] / _continuum(1.0e-3, 288.0)[window] > 2.5          # the self continuum goes as the vapour squared, the foreign as the vapour
    assert _continuum(4.0e-3, 270.0)[window] > _continuum(4.0e-3, 310.0)[window]          # the same vapour absorbs more when colder
    assert 0.02 < _continuum(1.0e-2, 288.0)[window] < 2.0                                  # a humid window at 10 microns: a tenth of a unit
    assert not _continuum(0.0, 288.0).any()


def test_the_water_continuum_is_confined_to_its_nodes_and_never_negative():
    depth = _continuum(0.03, 300.0)
    outside = (ir.CENTRES_CM1 < 200.0) | (ir.CENTRES_CM1 > 3600.0)
    assert not depth[outside].any() and np.all(depth[~outside] > 0.0) and np.all(np.isfinite(depth))
    assert _continuum(0.03, 300.0)[int(np.argmin(np.abs(ir.CENTRES_CM1 - 400.0)))] > depth[int(np.argmin(np.abs(ir.CENTRES_CM1 - 1000.0)))]


# ---- collision-induced pairs ----

def test_the_pairs_of_the_air_are_the_pairs_of_its_gases_and_missing_pairs_are_not_listed():
    assert {key for key, _first, _second in ir.prepare({"N2": 9.0e4, "H2": 1.0e4}).pairs} == {"N2-N2", "H2-N2", "H2-H2"}
    assert not ir.prepare({"Ar": 1.0e5}).pairs
    assert {key for key, _first, _second in ir.prepare({"N2": 1.0e5}).pairs} == {"N2-N2"}


def test_a_pair_absorbs_with_the_square_of_the_density():
    """Twice the pressure at one temperature and the same air: twice the density, twice the molecules above, four times the depth."""
    one = _depth({"N2": 1.0e5}, 250.0, 1.0e5)
    two = _depth({"N2": 1.0e5}, 250.0, 2.0e5)
    assert two.sum() == pytest.approx(4.0 * one.sum(), rel=1e-9)
    assert np.all(one == one[:, :1])                    # collision-induced absorption has no g structure: the same at every point


def test_the_pair_coefficient_is_interpolated_between_its_temperatures():
    cold, warm, between = (_depth({"N2": 1.0e5}, t, 1.0e5).sum() for t in (200.0, 300.0, 250.0))
    assert cold > between > warm                         # the fixture's k falls with the temperature (and the density with it)
    below = _depth({"N2": 1.0e5}, 120.0, 1.0e5)
    assert np.all(np.isfinite(below)) and below.sum() > 0.0


# ---- the far wings (co2_far_wing) ----

def _wing(partials, temperature, pressure):
    prepared = ir.prepare(partials)
    water = np.array([0.0])
    return ir._wing_depth(prepared, np.array([temperature]), np.array([pressure]), water, np.array([_molecules_cm2(pressure)]))[0]


def test_a_dense_carbon_dioxide_air_has_a_far_wing_and_a_thin_one_has_none():
    thick = _wing({"CO2": 9.0e6, "N2": 3.0e5}, 600.0, 9.3e6)
    thin = _wing({"N2": 0.78e5, "O2": 0.21e5, "CO2": 28.0}, 288.0, 1.0e5)
    assert thick.sum() > 0.0 and not thin.any()


def test_the_wing_grows_faster_than_the_pressure():
    """Its cross section goes with the pressure and the column with it: the depth with the square."""
    one = _wing({"CO2": 3.0e6, "N2": 1.0e5}, 500.0, 3.1e6)
    ten = _wing({"CO2": 3.0e7, "N2": 1.0e6}, 500.0, 3.1e7)
    assert ten.sum() > 50.0 * one.sum()


def test_the_wing_falls_between_the_bands_and_is_not_flat():
    wing = _wing({"CO2": 9.0e6, "N2": 3.0e5}, 600.0, 9.3e6)
    assert wing.max() > 10.0 * wing[wing > 0.0].min()


def test_the_wing_ends_where_the_carbon_dioxide_pressure_is_under_the_floor():
    below = _wing({"CO2": 0.9 * ir.WING_MIN_CO2_PA, "N2": 1.0e5}, 300.0, 1.0e5)
    assert not below.any()


def test_the_water_wing_adds_to_a_carbon_dioxide_air_and_not_to_a_thin_one():
    dry = _wing({"CO2": 9.0e6, "N2": 3.0e5}, 600.0, 9.3e6)
    prepared = ir.prepare({"CO2": 9.0e6, "N2": 3.0e5, "H2O": 3.0e3})
    water = np.array([3.0e3 / 9.303e6])
    wet = ir._wing_depth(prepared, np.array([600.0]), np.array([9.3e6]), water, np.array([_molecules_cm2(9.3e6)]))[0]
    assert wet.sum() > 0.99 * dry.sum()


# ---- gases without a table ----

def test_a_gas_without_a_table_takes_the_stand_in_path_and_still_absorbs():
    prepared = ir.prepare({"N2": 1.0e5, "NH3": 1.0e3})
    assert "NH3" in prepared.standin and prepared.tabled == ()
    depth = _depth({"N2": 1.0e5, "NH3": 1.0e3}, 300.0, 1.0e5)
    assert depth.shape == (ir.BANDS, ir.GRID_POINTS) and np.all(np.isfinite(depth)) and depth.sum() > 0.0
    assert np.all(np.diff(depth, axis=-1) >= -1e-12)


def test_no_gas_of_the_table_gives_a_nan_or_a_negative_depth():
    for row in gases.GASES:
        for temperature in (60.0, 300.0, 1500.0):
            depth = _depth({"N2": 1.0e5, row.formula: 1.0e3}, temperature, 1.0e5)
            assert np.all(np.isfinite(depth)) and np.all(depth >= 0.0), (row.formula, temperature)


def test_a_very_thick_air_does_not_underflow():
    depth = _depth({"CO2": 9.2e9}, 735.0, 9.2e9, gravity=8.87)
    assert np.all(np.isfinite(depth)) and depth.max() > 1.0e5


# ---- missing data ----

def test_a_missing_table_dataset_is_refused_by_name(tmp_path, monkeypatch):
    support.forget(monkeypatch, tmp_path / "empty")
    with pytest.raises(fetch.DataMissing, match="ir_ktables"):
        ir.prepare({"N2": 1.0e5, "CO2": 1.0e3})


def test_a_run_that_disregards_missing_data_takes_the_stand_ins_and_says_so(tmp_path, monkeypatch):
    support.forget(monkeypatch, tmp_path / "empty")
    with fetch.disregard_missing_data():
        prepared = ir.prepare({"N2": 1.0e5, "CO2": 1.0e3})
        assert prepared.tabled == () and "CO2" in prepared.standin and not prepared.pairs
        assert _depth({"N2": 1.0e5, "CO2": 1.0e3}, 300.0, 1.0e5).sum() > 0.0
        assert {key for key, _what in fetch.stand_ins_used()} == {"ir_ktables", "ir_cia"}


def test_tables_made_for_other_bands_are_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(ir_ktables, "BAND_EDGES_CM1", ir_ktables.BAND_EDGES_CM1 * 1.01)
    ir._build_tables.cache_clear()
    with pytest.raises(ValueError, match="other bands"):
        ir.prepare({"N2": 1.0e5, "CO2": 1.0e3})


def test_bad_input_is_refused():
    with pytest.raises(ValueError):
        ir.prepare({})
