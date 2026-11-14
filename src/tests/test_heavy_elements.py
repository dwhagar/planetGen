# tests/test_heavy_elements.py

"""
GEN.1371: the s- and r-process split of a star's heavy elements and their share between the silicate and the metal of a body (`physics/heavy_elements`). Solar System bodies are test targets only.
"""

import math

import pytest

from planetgen import tuning
from planetgen.physics import heavy_elements as he, stellar_abundances as sa


def _solar():
    return sa.composition(0.0, sa.abundance_offsets(0.0, sa.BRANCH_THIN, (0.0,) * len(sa.RESIDUAL_ELEMENTS)))


def test_a_solar_scale_star_has_the_solar_heavy_elements_to_within_its_own_trend():
    """At [Fe/H] 0 the split gives each element the solar abundance, moved to proto-solar, times a mix of the r and s scales that is close to 1 (the Eu and Y trends sit near zero there)."""
    found = he.birth_heavy(_solar(), sa.BRANCH_THIN)
    for symbol, log_eps in tuning.STELLAR_HEAVY_LOG_EPS.items():
        solar = 10.0 ** (log_eps + tuning.STELLAR_PROTO_METAL_OFFSET_DEX - 12.0)
        assert found["abundances"][symbol] / solar == pytest.approx(1.0, abs=0.2), symbol


def test_an_r_process_enhanced_star_raises_europiumlike_elements_more_than_barium():
    star = _solar()
    base = he.birth_heavy(star, sa.BRANCH_THIN)
    high = he.birth_heavy(star, sa.BRANCH_THIN, eu_normal=2.0)
    gold_gain = high["abundances"]["Au"] / base["abundances"]["Au"]
    barium_gain = high["abundances"]["Ba"] / base["abundances"]["Ba"]
    assert gold_gain > 1.2 and 1.0 < barium_gain < gold_gain          # Au is 94 percent r-process, Ba 15 percent


def test_an_s_process_enhanced_star_raises_barium_more_than_gold():
    star = _solar()
    base = he.birth_heavy(star, sa.BRANCH_THIN)
    high = he.birth_heavy(star, sa.BRANCH_THIN, s_normal=2.0)
    assert high["abundances"]["Ba"] / base["abundances"]["Ba"] > 1.5 > high["abundances"]["Au"] / base["abundances"]["Au"] > 1.0


def test_the_trend_of_the_s_process_is_held_at_the_edge_of_its_fit():
    assert he.s_process_xfe(0.0) == pytest.approx(0.030)
    assert he.s_process_xfe(-3.0) == he.s_process_xfe(tuning.STELLAR_FE_H_FIT_RANGE[sa.BRANCH_THIN][0])


def test_the_siderophile_elements_go_to_the_core_and_the_lithophile_ones_stay_in_the_mantle():
    earth_like = dict(rock_fraction=0.675, metal_fraction=0.325)
    assert he.core_share(he.partition_of("Au"), **earth_like) > 0.99
    assert he.core_share(he.partition_of("Nd"), **earth_like) < 0.01
    assert he.core_share(1.0, **earth_like) == pytest.approx(0.325)          # D 1: shared by mass
    assert he.core_share(1000.0, 0.0, 1.0) == 0.0 and he.core_share(1000.0, 1.0, 0.0) == 0.0


def test_the_mass_of_an_element_is_conserved_between_the_silicate_and_the_metal():
    ratios = he.star_to_mg(0.0, 4.6, "test-key-1371")
    rock, metal = 0.675, 0.325
    ppm = he.bulk_ppm(ratios, rock, metal)
    for symbol in he.ELEMENTS:
        total = (ppm[symbol][0] * rock + ppm[symbol][1] * metal) * 1.0e-6
        assert total == pytest.approx(ratios[symbol] * he.ROCK_MG_MASS_FRACTION * rock, rel=1e-9), symbol


def test_a_body_with_no_metal_keeps_everything_in_its_silicate_and_the_same_star_gives_the_same_numbers():
    ratios = he.star_to_mg(0.1, 3.0, "test-key-1371")
    assert he.star_to_mg(0.1, 3.0, "test-key-1371") == ratios
    ppm = he.bulk_ppm(ratios, 1.0, 0.0)
    assert all(core == 0.0 and mantle > 0.0 for mantle, core in ppm.values())


def test_a_metal_rich_star_scales_the_elements_with_its_offset_from_the_solar_scale():
    offsets = he.offset_dex(0.0, 4.6, "test-key-1371")
    assert all(abs(value) < 1.0 for value in offsets.values())
    assert set(offsets) == set(he.ELEMENTS)
    groups = he.group_ppm(he.bulk_ppm(he.star_to_mg(0.0, 4.6, "test-key-1371"), 0.675, 0.325))
    assert groups["platinum group"][1] > 100.0 * groups["platinum group"][0] and groups["rare earths"][0] > groups["rare earths"][1]
    assert math.isfinite(groups["gold"][1])


class _Conn:
    """Answers the queries `describe/heavy_rows` makes, by the table each names."""

    def __init__(self, planet, star):
        self.planet, self.star = planet, star

    def execute(self, sql, params=()):
        self.last = sql
        return self

    def fetchone(self):
        if "FROM planets" in self.last:
            return self.planet
        if "FROM stars" in self.last:
            return self.star
        return None          # no sector: the default distance from the center


def test_the_sheet_rows_of_a_planet_give_the_core_the_siderophile_metals_and_the_rock_the_rare_earths(monkeypatch):
    from planetgen.describe import heavy_rows
    planet = {"star_system_id": 1, "star_id": 7, "metal_mass_fraction": 0.325, "rock_mass_fraction": 0.675}
    star = {"mass_kg": 1.989e30, "metallicity_fe_h": 0.0, "age_gy": 4.6, "star_system_id": 1}
    found = {row[1]: row for row in heavy_rows.rows(_Conn(planet, star), "planet", 5)}
    assert set(found) == set(heavy_rows.KEYS)
    assert found["heavy_pge_core_ppm"][4] > 100.0 * found["heavy_pge_rock_ppm"][4] > 0.0
    assert found["heavy_gold_core_ppm"][4] > 100.0 * found["heavy_gold_rock_ppm"][4]
    assert found["heavy_ree_rock_ppm"][4] > 1.0
    assert all(heavy_rows.is_derived_key(key) for key in found)
    assert heavy_rows.rows(_Conn(planet, {**star, "metallicity_fe_h": None}), "planet", 5) == []
    assert heavy_rows.rows(_Conn({**planet, "metal_mass_fraction": 0.0}, star), "planet", 5) == []
    assert heavy_rows.rows(_Conn(planet, star), "belt", 5) == []
