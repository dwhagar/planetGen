# tests/test_volatile_depletion.py

"""
GEN.1361: the volatility depletion law (`physics/volatile_depletion.py`) against its own 26-row table.

The table is the bulk Earth over CI chondrites for 26 elements against their 50 percent condensation temperature; the law is a hinge in log space.
The tests refit the hinge from the table (so the constants in the data file cannot drift from it), hold the law to the table's residuals and to
the rows the research set names (sodium, potassium and zinc), and check the distance lever and the refusals.
"""

import json
import math

import numpy as np
import pytest

from planetgen.physics import volatile_depletion as vd
from planetgen.physics.elements import DataError

SOLAR_LUMINOSITY_W = 3.828e26


def _hinge_fit():
    """Least squares of log10 f = min(0, b (t50 - c)) over the table: for each cut on a fine grid the best slope is closed-form."""
    t50 = np.array([row[0] for row in vd.ROWS.values()])
    logged = np.array([row[2] for row in vd.ROWS.values()])
    best = None
    for cut in np.arange(1250.0, 1400.0, 0.1):
        shape = np.minimum(0.0, t50 - cut)
        slope = float(np.dot(shape, logged) / np.dot(shape, shape))
        rms = float(np.sqrt(np.mean((slope * shape - logged) ** 2)))
        if best is None or rms < best[0]:
            best = (rms, slope, float(cut))
    return best


def test_the_table_holds_26_elements_in_t50_order_and_each_log_is_the_log_of_its_ratio():
    assert len(vd.ROWS) == 26
    t50 = [vd.t50_k(symbol) for symbol in vd.ROWS]
    assert t50 == sorted(t50) and t50[0] == 664.0
    for symbol, (_t, ratio, logged, _phase) in vd.ROWS.items():
        assert abs(math.log10(ratio) - logged) < 0.011, symbol
    assert vd.earth_over_ci("Mg") == 1.0
    assert vd.ROWS["Na"][:3] == (958.0, 0.28, -0.55)


def test_the_law_in_the_data_file_is_the_hinge_the_table_gives():
    rms, slope, cut = _hinge_fit()
    assert slope == pytest.approx(vd.SLOPE_PER_K, abs=3e-5)
    assert cut == pytest.approx(vd.CUT_K, abs=1.0)
    assert rms == pytest.approx(0.079, abs=0.002)


def test_the_law_reproduces_the_table_to_its_residual_and_the_named_rows_untuned():
    residuals = [vd.log_retention(vd.t50_k(symbol)) - vd.ROWS[symbol][2] for symbol in vd.ROWS if vd.t50_k(symbol) >= vd.MIN_T50_EARTH_K]
    assert math.sqrt(sum(r * r for r in residuals) / len(residuals)) < 0.085
    assert max(abs(r) for r in residuals) < 0.25
    for symbol in ("Na", "K", "Zn"):                    # the alkalis and zinc the research note names: within a factor of 1.5
        assert vd.retained_fraction(vd.t50_k(symbol)) == pytest.approx(vd.earth_over_ci(symbol), rel=0.5), symbol


def test_the_law_is_the_one_the_research_thread_quotes_to_a_hundredth_of_a_dex():
    for t50 in (664.0, 800.0, 958.0, 1100.0, 1300.0):
        quoted = -2.150 + 0.001623 * t50
        assert abs(vd.log_retention(t50) - quoted) < 0.015


def test_refractory_elements_are_kept_whole_and_the_share_falls_with_a_cooler_t50():
    for symbol in ("Ca", "Al", "Ti", "V", "Sr", "Y", "Zr"):
        assert vd.retained_fraction(vd.t50_k(symbol)) == 1.0, symbol
    shares = list(vd.retained_fractions().values())
    assert shares == sorted(shares) and shares[0] < 0.1 and shares[-1] == 1.0
    assert vd.retained_fraction(vd.CUT_K) == 1.0 and vd.retained_fraction(vd.CUT_K - 100.0) < 1.0


def test_the_cut_follows_the_blackbody_temperature_of_the_formation_distance():
    assert vd.cut_temperature_k(1.0, SOLAR_LUMINOSITY_W) == pytest.approx(vd.CUT_K)
    assert vd.REFERENCE_BLACKBODY_K == pytest.approx(278.3, abs=0.2)
    assert vd.cut_temperature_k(4.0, SOLAR_LUMINOSITY_W) == pytest.approx(0.5 * vd.CUT_K)        # d^-1/2
    assert vd.cut_temperature_k(1.0, 16.0 * SOLAR_LUMINOSITY_W) == pytest.approx(2.0 * vd.CUT_K)  # L^1/4
    mars = vd.cut_temperature_k(1.524, SOLAR_LUMINOSITY_W)
    assert 1050.0 < mars < 1080.0
    # a body formed farther out keeps more of the same element (Mars keeps about 2.6 times the Earth law's sodium; the research note has about twice)
    earth_na = vd.retained_fraction(958.0)
    assert 2.0 < vd.retained_fraction(958.0, mars) / earth_na < 3.2
    # a closer or hotter zone has a higher cut and loses more (Mercury is the known counterexample: its K/U is not lower than the Earth's)
    assert vd.retained_fraction(958.0, vd.cut_temperature_k(0.387, SOLAR_LUMINOSITY_W)) < 0.05 * earth_na
    # beyond the zone where the cut falls under the T50 of the element the element is kept whole
    assert vd.retained_fraction(958.0, vd.cut_temperature_k(5.0, SOLAR_LUMINOSITY_W)) == 1.0
    assert vd.retained_fraction(958.0, vd.cut_temperature_k(0.5, 0.01 * SOLAR_LUMINOSITY_W)) == 1.0


def test_a_colder_cut_keeps_more_than_the_earth_cut_and_a_hotter_cut_less():
    cold, hot = 0.8 * vd.CUT_K, 1.2 * vd.CUT_K
    for t50 in (700.0, 958.0, 1200.0):
        assert vd.retained_fraction(t50, cold) >= vd.retained_fraction(t50) >= vd.retained_fraction(t50, hot)
    assert vd.retained_fraction(958.0, cold) > vd.retained_fraction(958.0, hot)


def test_the_ices_and_noble_gases_and_bad_numbers_are_refused():
    for t50 in (40.0, 180.0, 252.0, 499.0):
        with pytest.raises(ValueError):
            vd.retained_fraction(t50)
    for bad in (0.0, -5.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            vd.retained_fraction(bad)
    with pytest.raises(ValueError):
        vd.retained_fraction(900.0, 0.0)
    with pytest.raises(ValueError):
        vd.cut_temperature_k(0.0, SOLAR_LUMINOSITY_W)
    with pytest.raises(DataError):
        vd.t50_k("Unobtainium")


def test_a_malformed_file_is_refused(tmp_path):
    good = json.load(open(vd.DATA_FILE, encoding="utf-8"))

    def refused(change):
        document = json.loads(json.dumps(good))
        change(document)
        path = tmp_path / "bad.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        with pytest.raises(DataError):
            vd.load(str(path))

    refused(lambda d: d["law"].update(slope_per_k=-1.0))
    refused(lambda d: d["law"].pop("cut_k"))
    refused(lambda d: d["reference"].update(luminosity_w=0))
    refused(lambda d: d["elements"][0].update(log10_depletion=0.5))
    refused(lambda d: d["elements"].__setitem__(1, dict(d["elements"][0])))
    refused(lambda d: d.update(elements=d["elements"][:5]))
    with pytest.raises(DataError):
        vd.load(str(tmp_path / "missing.json"))
