"""
Climate-tuning regression suite.

Companion to `climate_tuning_cli.py` (the human-driven iteration tool used to
arrive at the values below) and `test_physical_plausibility.py` (which
deliberately does NOT assert on class-to-class differentiation, since that
was the whole gap this session's tuning pass closed -- see
docs/TODO.md/CHANGELOG.md for the greenhouse-formula fix and per-class
`albedo_range`/`atm_molar_density_range`/`atm_density_range`/
`greenhouse_multiplier_range` tuning it describes).

Unlike climate_tuning_cli.py's human-reviewed report, this module DOES
assert -- but only generously-toleranced bands and relative orderings
between classes, never tight exact-value checks, since generation is
inherently stochastic (host star type, and every tuned range, are randomized
per sample; `stellarObjects.utils.reseed_rng()` reseeds from `secrets` on
every generation call, so this cannot be made deterministic via a seed).
The point is to catch a real regression (e.g. a future edit that
accidentally collapses two classes back together, or drags Class M's
pressure back below realistic Earth values) without being brittle against
the normal sample-to-sample variance these ranges are designed to produce.

Run with: pytest tests/test_climate_tuning.py
"""
import statistics

import pytest

from stellarObjects import plausibility

N_SAMPLE = 200

TUNED_CLASSES = ("E", "F", "G", "H", "K", "L", "M", "N", "O", "V")


@pytest.fixture(scope="module")
def class_means():
    """
    One shared {class: (mean_surface_temperature_k, mean_atmospheric_pressure_pa)}
    dict for the whole module, built from a single N_SAMPLE-sized batch per
    tuned class (zone 'e', moons excluded -- moons would mix in other
    classes' own tuning and aren't the subject of this suite).
    """
    means = {}
    for cls in TUNED_CLASSES:
        records = plausibility.generate_sample(cls, "e", N_SAMPLE, include_moons=False)
        means[cls] = (
            statistics.mean(r["surface_temperature"] for r in records),
            statistics.mean(r["atmospheric_pressure"] for r in records),
        )
    return means


# --- Real-world-analog band checks (generous -- see module docstring) ---

def test_class_m_lands_within_earth_like_band(class_means):
    temp, pressure = class_means["M"]
    assert 265 <= temp <= 305, f"Class M mean surface_temperature {temp:.1f}K outside Earth-like band"
    assert 60_000 <= pressure <= 140_000, f"Class M mean atmospheric_pressure {pressure:.1f}Pa outside Earth-like band"


def test_class_n_lands_within_venus_like_band(class_means):
    temp, pressure = class_means["N"]
    assert 650 <= temp <= 820, f"Class N mean surface_temperature {temp:.1f}K outside Venus-like band"
    assert 5e6 <= pressure <= 1.5e7, f"Class N mean atmospheric_pressure {pressure:.1f}Pa outside Venus-like band"


def test_class_k_lands_within_mars_like_band(class_means):
    temp, pressure = class_means["K"]
    # Colder/thinner than every other tuned terrestrial class, but not
    # literally at Mars' own absolute values -- K is generated in the same
    # ecosphere zone as every other habitable class, not at Mars' real,
    # farther orbital distance, so its baseline runs warmer than real Mars
    # regardless of tuning (see PLANET_CLASSES["K"]'s docstring note).
    assert 200 <= temp <= 270, f"Class K mean surface_temperature {temp:.1f}K outside Mars-like band"
    assert 100 <= pressure <= 3000, f"Class K mean atmospheric_pressure {pressure:.1f}Pa outside Mars-like band"


# --- Relative ordering / differentiation checks ---

def test_class_n_is_hottest_and_highest_pressure_of_the_tuned_classes(class_means):
    n_temp, n_pressure = class_means["N"]
    for cls in TUNED_CLASSES:
        if cls == "N":
            continue
        temp, pressure = class_means[cls]
        assert n_temp > temp, f"Class N ({n_temp:.1f}K) should be hotter than Class {cls} ({temp:.1f}K)"
        assert n_pressure > pressure, f"Class N ({n_pressure:.1f}Pa) should be higher-pressure than Class {cls} ({pressure:.1f}Pa)"


def test_class_k_is_colder_and_thinner_than_class_l(class_means):
    # "usually has vegetation, K does not" -- L should read as modestly more
    # hospitable (warmer, thicker atmosphere) than K, not identical to it.
    k_temp, k_pressure = class_means["K"]
    l_temp, l_pressure = class_means["L"]
    assert l_temp > k_temp, f"Class L ({l_temp:.1f}K) should be warmer than Class K ({k_temp:.1f}K)"
    assert l_pressure > k_pressure, f"Class L ({l_pressure:.1f}Pa) should be higher-pressure than Class K ({k_pressure:.1f}Pa)"


def test_class_l_is_colder_than_class_m(class_means):
    # L is "marginally habitable" -- it shouldn't reach M's Earth-like norm.
    l_temp, _ = class_means["L"]
    m_temp, _ = class_means["M"]
    assert l_temp < m_temp, f"Class L ({l_temp:.1f}K) should be colder than Class M ({m_temp:.1f}K)"


def test_class_h_is_hotter_and_drier_than_class_o(class_means):
    # H = hot/dry desert extreme, O = warm/wet ocean extreme, both anchored
    # near Class M -- they should diverge from each other, not collapse.
    h_temp, h_pressure = class_means["H"]
    o_temp, o_pressure = class_means["O"]
    assert h_temp > o_temp, f"Class H ({h_temp:.1f}K) should be hotter than Class O ({o_temp:.1f}K)"
    assert h_pressure < o_pressure, f"Class H ({h_pressure:.1f}Pa) should be lower-pressure (drier) than Class O ({o_pressure:.1f}Pa)"


def test_class_e_f_g_form_a_cooling_progression(class_means):
    # "Class E through G are very similar to Class M but ... younger" and
    # "F to G move on a progression of a cooling planet toward an M" --
    # E should be hottest, G coolest, strictly ordered.
    e_temp, _ = class_means["E"]
    f_temp, _ = class_means["F"]
    g_temp, _ = class_means["G"]
    assert e_temp > f_temp > g_temp, (
        f"Expected E > F > G cooling progression, got E={e_temp:.1f}K F={f_temp:.1f}K G={g_temp:.1f}K"
    )


def test_class_g_lands_near_class_m(class_means):
    # G is the last step of the E->F->G progression, converging toward
    # M/O/K/L/N -- it should land close to M, not still running hot like E/F.
    g_temp, _ = class_means["G"]
    m_temp, _ = class_means["M"]
    assert abs(g_temp - m_temp) <= 25, (
        f"Class G ({g_temp:.1f}K) should have converged near Class M ({m_temp:.1f}K) by the end of the E->F->G progression"
    )


def test_class_v_is_hotter_and_thicker_than_m_but_cooler_than_n(class_means):
    # "a thick atmosphere with high surface temperature and pressure" --
    # between M's Earth-like norm and N's full Venus extreme.
    m_temp, m_pressure = class_means["M"]
    v_temp, v_pressure = class_means["V"]
    n_temp, n_pressure = class_means["N"]
    assert m_temp < v_temp < n_temp, (
        f"Expected M < V < N surface_temperature, got M={m_temp:.1f}K V={v_temp:.1f}K N={n_temp:.1f}K"
    )
    assert m_pressure < v_pressure < n_pressure, (
        f"Expected M < V < N atmospheric_pressure, got M={m_pressure:.1f}Pa V={v_pressure:.1f}Pa N={n_pressure:.1f}Pa"
    )
