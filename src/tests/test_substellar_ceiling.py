# tests/test_substellar_ceiling.py

"""GEN.494: the hottest a brown dwarf can be at its age, from the Baraffe et al. 2015 isochrones (`physics/substellar_ceiling.py`).

The grid below is rows of the BHAC15 file itself (the 0.07 and 0.072 solar mass rows at six of its ages), so the numbers are the table's.
"""

import numpy as np
import pytest

from planetgen.datasets.baraffe_isochrones import Isochrone
from planetgen.physics import substellar_ceiling as ceiling

# age (Gyr): (Teff at 0.07, Teff at 0.072) from BHAC15_iso.2mass
ROWS = {0.0005: (2835, 2835), 0.001: (2897, 2896), 0.01: (2946, 2954), 0.1: (2794, 2813), 1.0: (2047, 2109), 10.0: (1626, 1790)}


def grid():
    out = []
    for age, (teff_07, teff_072) in ROWS.items():
        out.append(Isochrone(age_gyr=age, mass_solar=np.array([0.07, 0.072]), teff_k=np.array([teff_07, teff_072], dtype=float),
                             log_luminosity_solar=np.array([-3.0, -3.0]), log_g=np.array([4.0, 4.0]), radius_solar=np.array([0.1, 0.1])))
    return out


G = grid()


def test_the_limit_mass_is_hottest_at_ten_million_years_as_the_paper_says():
    """Baraffe 2015 section 3.3: Teff about 2900 K at 1 Myr and 2950 K at 10 Myr for the 0.07 solar mass limit."""
    assert ceiling.ceiling_k(0.001, G) == pytest.approx(2900, abs=10)
    assert ceiling.ceiling_k(0.01, G) == pytest.approx(2950, abs=10)
    assert max(ceiling.ceiling_k(age, G) for age in ROWS) == ROWS[0.01][0]
    assert ceiling.HOTTEST_K == pytest.approx(ROWS[0.01][0], abs=10)


def test_the_ceiling_falls_with_age_after_the_peak():
    ages = [0.01, 0.1, 1.0, 10.0]
    values = [ceiling.ceiling_k(a, G) for a in ages]
    assert values == sorted(values, reverse=True)
    assert values[-1] < 1700 < 1800 < values[-2]


def test_an_old_brown_dwarf_cannot_be_late_m_like():
    """The late-M-like row (2238 K and up) is a young brown dwarf: at 1 Gyr the limit mass is already near 2050 K."""
    assert ceiling.can_be_this_hot(2300.0, 0.05, G)
    assert not ceiling.can_be_this_hot(2300.0, 1.0, G)
    assert ceiling.can_be_this_hot(1500.0, 10.0, G)


def test_outside_the_table_it_answers_none():
    assert ceiling.ceiling_k(2.0e-4, G) is None
    assert ceiling.ceiling_k(20.0, G) is None
    assert ceiling.can_be_this_hot(2000.0, 20.0, G) is None


def test_the_age_a_temperature_stops_being_possible_is_found_on_the_falling_side():
    age = ceiling.oldest_age_gyr(2400.0, G)
    assert 0.1 < age < 1.0
    assert ceiling.ceiling_k(age, G) == pytest.approx(2400.0, abs=60.0)
    assert ceiling.oldest_age_gyr(3200.0, G) == 0.0
    assert ceiling.oldest_age_gyr(1000.0, G) == 10.0


def test_the_real_table_gives_the_papers_numbers_where_it_is_held():
    from planetgen.datasets import baraffe_isochrones, fetch
    try:
        baraffe_isochrones.isochrones()
    except fetch.DataMissing:
        pytest.skip("the BHAC15 file is not downloaded on this machine")
    assert ceiling.ceiling_k(0.001) == pytest.approx(2900, abs=15)
    assert ceiling.ceiling_k(0.01) == pytest.approx(2950, abs=15)
    assert ceiling.ceiling_k(10.0) < 1700


def test_the_cap_follows_the_age_below_the_limit_mass_and_is_the_peak_at_or_above_it():
    assert ceiling.cap_k(0.05, 0.001, G) == ceiling.hottest_k(0.001, G)
    assert ceiling.cap_k(0.0699, 5.0, G) == ceiling.hottest_k(5.0, G)
    assert ceiling.cap_k(0.07, 5.0, G) == ceiling.HOTTEST_K == ceiling.cap_k(0.0764, 0.0001, G)


def test_the_cap_holds_the_age_to_the_table_and_never_asks_for_a_cooler_answer_than_the_table_has():
    youngest, oldest = G[0].age_gyr, G[-1].age_gyr
    assert ceiling.hottest_k(youngest / 100.0, G) == ceiling.ceiling_k(youngest, G)
    assert ceiling.hottest_k(oldest * 5.0, G) == ceiling.ceiling_k(oldest, G)
