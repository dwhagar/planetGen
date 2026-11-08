# tests/test_accessory_minerals.py

"""
GEN.1368: the accessory minerals of the asteroid taxon mixes, scaled by the host star's bulk.
"""

import dataclasses
import math

import pytest

from planetgen.physics import accessory_minerals as am, bulk_draw, compaction, materials, small_body_composition as sbc
from planetgen.util import draw


def _host(**changes):
    logs = am.tuning.STELLAR_SOLAR_LOG_EPS
    solar = {name: am.reference_ratio(element) for name, element in (("cr", "Cr"), ("ti", "Ti"), ("p", "P"), ("al", "Al"), ("ni", "Ni"), ("mn", "Mn"), ("na", "Na"))}
    base = dict(fe=0.8, si=0.9, ca=0.06, c_over_o=0.5, y_envelope=0.25, z=0.014, k=am.reference_ratio("K") if "K" in logs else 0.0, **solar)
    base.update(changes)
    return bulk_draw.HostRatios(**base)


def test_every_accessory_is_a_row_of_the_table_with_a_grain_density():
    for taxon, rows in am.ACCESSORIES.items():
        assert taxon in sbc.TAXA
        for row, share in rows.items():
            assert materials.material(row), row
            assert 0.0 < share < 0.2, (taxon, row)
            assert compaction.grain_density_kg_m3({row: 1.0}) > 0.0, row
    assert set(am.CONTROL) <= {row for rows in am.ACCESSORIES.values() for row in rows} | {"whitlockite", "hydroxylapatite", "alabandite"}


def test_a_mix_adds_up_to_one_and_keeps_the_base_in_proportion():
    base = sbc.TAXON_DATA["S"].fractions
    mixed = am.mix("S", base)
    assert math.fsum(mixed.values()) == pytest.approx(1.0)
    assert mixed["chromite"] > 0.0 and mixed["merrillite"] > 0.0
    ratio = {name: mixed[name] / base[name] for name in base}
    assert max(ratio.values()) == pytest.approx(min(ratio.values()))     # the base minerals give way together


def test_the_ordinary_chondrite_taxon_carries_chromite_and_phosphates_and_the_cv_taxa_the_inclusion_phases():
    s, k, l = am.mix("S", sbc.TAXON_DATA["S"].fractions), am.mix("K", sbc.TAXON_DATA["K"].fractions), am.mix("L", sbc.TAXON_DATA["L"].fractions)
    assert 0.003 < s["chromite"] < 0.01 and s["merrillite"] > 0.0 and s["chlorapatite"] > 0.0
    assert {"melilite", "hibonite", "perovskite", "spinel"} <= set(l) and {"melilite", "hibonite", "spinel"} <= set(k)
    e = am.mix("E", sbc.TAXON_DATA["E"].fractions)
    assert {"oldhamite", "daubreelite", "niningerite"} <= set(e)
    c = am.mix("C", sbc.TAXON_DATA["C"].fractions)
    assert {"cronstedtite", "saponite", "dolomite", "epsomite", "gypsum", "pentlandite", "tochilinite"} <= set(c)
    assert c["carbon_diamond"] > 0.0 and c["moissanite"] > 0.0           # nanodiamond and presolar silicon carbide


def test_a_star_rich_in_an_element_gives_more_of_the_minerals_that_element_sets():
    base = sbc.TAXON_DATA["S"].fractions
    rich, poor = _host(cr=2.0 * am.reference_ratio("Cr"), p=0.5 * am.reference_ratio("P")), _host()
    assert am.mix("S", base, rich)["chromite"] > 1.8 * am.mix("S", base, poor)["chromite"] * 0.99
    assert am.mix("S", base, rich)["merrillite"] < am.mix("S", base, poor)["merrillite"]
    assert am.mix("S", base, poor)["chromite"] == pytest.approx(am.mix("S", base)["chromite"])    # the reference star is the reference
    assert am.scale("Cr", None) == 1.0


def test_a_mineral_no_tracked_element_sets_is_not_scaled():
    host = _host(ni=3.0 * am.reference_ratio("Ni"))
    assert am.accessories("C", host)["dolomite"] == am.ACCESSORIES["C"]["dolomite"]
    assert am.accessories("C", host)["pentlandite"] == pytest.approx(3.0 * am.ACCESSORIES["C"]["pentlandite"])


def test_a_belt_draws_with_its_accessories_and_the_host_changes_them():
    draw.set_run_seed(21)
    plain = sbc.draw_belt(2.7, 3.828e26, 5.2e-5)
    assert math.fsum(plain.components.values()) == pytest.approx(1.0)
    assert any(row in plain.components for row in ("chromite", "merrillite", "melilite", "cronstedtite"))
    draw.set_run_seed(21)
    rich = sbc.draw_belt(2.7, 3.828e26, 5.2e-5, host=_host(cr=3.0 * am.reference_ratio("Cr")))
    assert rich.components.get("chromite", 0.0) > plain.components.get("chromite", 0.0)


def test_the_inclusion_phases_count_as_the_refractory_axis_and_the_clay_phases_as_matrix():
    assert {"melilite", "hibonite", "perovskite"} <= set(sbc.REFRACTORY_ROWS)
    assert {"cronstedtite", "saponite", "tochilinite"} <= set(sbc.MATRIX_ROWS)
