# tests/test_hot_condensates.py

"""
GEN.1375: the cloud temperatures and trace gases of hot atmospheres that follow the host star's elements (`physics/hot_condensates.py`).

At the Sun's composition the fits must give the table of the research part (06 section 2.2, computed from [R3165, R4562]); a host poorer or
richer in one element moves only the clouds that element sets; the trace gases are held to the least of their elements.
"""

import json
import math

import pytest

from planetgen.physics import hot_condensates as hc
from planetgen.physics import stellar_abundances as sa
from planetgen.physics.elements import DataError

PRESSURES = (1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0)

TABLE = {
    "Cr": (1174, 1245, 1325, 1416, 1521, 1642, 1784),
    "MnS": (1096, 1149, 1207, 1271, 1343, 1423, 1514),
    "Na2S": (774, 819, 871, 929, 996, 1072, 1162),
    "ZnS": (665, 694, 725, 760, 798, 841, 888),
    "KCl": (625, 662, 702, 749, 801, 862, 933),
}
"""dict: Part 06 section 2.2, [Fe/H] 0, K at 1e-4 to 100 bar."""


def _host(fe_h=0.0, **offsets):
    return sa.composition(fe_h, offsets)


@pytest.mark.parametrize("species", sorted(TABLE))
def test_the_suns_clouds_are_the_research_table(species):
    host = _host()
    for pressure, expected in zip(PRESSURES, TABLE[species]):
        assert hc.cloud_temperature_k(species, pressure, host) == pytest.approx(expected, abs=1.0)


def test_a_metal_rich_host_with_every_element_following_iron_gives_the_printed_plus_half_row():
    host = _host(0.5)
    expected = (807, 857, 914, 978, 1052, 1138, 1240)                 # Na2S, [Fe/H] +0.5
    for pressure, value in zip(PRESSURES, expected):
        assert hc.cloud_temperature_k("Na2S", pressure, host) == pytest.approx(value, abs=1.0)


def test_a_cloud_follows_its_own_elements_and_no_other():
    sun = _host()
    poor_zinc = _host(0.0, Zn=-0.3)
    assert hc.cloud_temperature_k("ZnS", 1.0, poor_zinc) < hc.cloud_temperature_k("ZnS", 1.0, sun) - 10.0
    for species in ("Cr", "MnS", "Na2S", "KCl"):
        assert hc.cloud_temperature_k(species, 1.0, poor_zinc) == hc.cloud_temperature_k(species, 1.0, sun)
    rich_manganese = _host(0.0, Mn=0.2)
    assert hc.cloud_temperature_k("MnS", 1.0, rich_manganese) > hc.cloud_temperature_k("MnS", 1.0, sun)
    assert hc.cloud_temperature_k("Cr", 1.0, _host(0.0, Cr=0.2)) > hc.cloud_temperature_k("Cr", 1.0, sun)
    # sulfur is a term of MnS, Na2S and ZnS but not of chromium
    sulfur_rich = _host(0.0, S=0.2)
    assert hc.cloud_temperature_k("Cr", 1.0, sulfur_rich) == hc.cloud_temperature_k("Cr", 1.0, sun)
    assert all(hc.cloud_temperature_k(s, 1.0, sulfur_rich) > hc.cloud_temperature_k(s, 1.0, sun) for s in ("MnS", "Na2S", "ZnS"))


def test_potassium_chloride_is_limited_by_potassium_alone():
    sun = _host()
    assert hc.cloud_temperature_k("KCl", 1.0, _host(0.0, K=0.3)) > hc.cloud_temperature_k("KCl", 1.0, sun)
    assert hc.cloud_temperature_k("KCl", 1.0, _host(0.0, Ca=0.3)) == hc.cloud_temperature_k("KCl", 1.0, sun)
    assert hc.KCL["element"] == "K"


def test_a_cloud_is_hotter_deeper_and_in_a_richer_star():
    for species in hc.CLOUD_SPECIES:
        values = [hc.cloud_temperature_k(species, pressure, _host()) for pressure in PRESSURES]
        assert values == sorted(values) and len(set(values)) == len(values)
        assert hc.cloud_temperature_k(species, 1.0, _host(0.3)) > hc.cloud_temperature_k(species, 1.0, _host(-0.3))


def test_the_trace_elements_follow_iron_from_their_solar_abundance():
    sun = _host()
    assert hc.abundance_per_h("Rb", sun) == pytest.approx(10.0 ** (2.32 - 12.0))
    assert hc.abundance_per_h("Cs", sun) == pytest.approx(10.0 ** (1.08 - 12.0))
    assert hc.abundance_per_h("Rb", _host(0.3)) == pytest.approx(2.0 * hc.abundance_per_h("Rb", sun), rel=0.01)
    assert hc.abundance_per_h("Zn", _host(0.0, Zn=-0.1)) == pytest.approx(10.0 ** (4.56 - 12.1))
    assert set(hc.TRACE_ELEMENTS) == {"Rb", "Cs", "Ge", "As", "Se"}


def test_the_halogens_take_the_rocky_inventory_offsets():
    from planetgen.physics import mixture_rocky
    sun = _host()
    assert hc.abundance_per_h("Cl", sun) == pytest.approx(10.0 ** (5.31 - 12.0 + mixture_rocky.STELLAR_CL_OVER_FE_DEX))
    assert hc.abundance_per_h("F", sun) == pytest.approx(10.0 ** (4.4 - 12.0 + mixture_rocky.STELLAR_F_OVER_FE_DEX))


def test_the_trace_gases_are_held_to_their_scarcest_element():
    sun = _host()
    assert set(hc.HOT_TRACE_GASES) == {"GeH4", "AsH3", "H2Se", "MgH", "CaH", "NaH", "LiCl", "LiF", "AlCl", "PN", "HCP"}
    assert hc.gas_ceiling("GeH4", sun) == pytest.approx(hc.abundance_per_h("Ge", sun))
    assert hc.gas_ceiling("LiCl", sun) == pytest.approx(hc.abundance_per_h("Li", sun))              # lithium is scarcer than chlorine
    assert hc.gas_ceiling("HCP", sun) == pytest.approx(hc.abundance_per_h("P", sun))                # phosphorus is scarcer than carbon
    assert hc.gas_ceiling("PN", sun) == pytest.approx(hc.abundance_per_h("P", sun))
    assert hc.gas_ceiling("AlCl", sun) == pytest.approx(hc.abundance_per_h("Cl", sun))              # chlorine is scarcer than aluminium
    assert hc.gas_ceiling("MgH", sun) > hc.gas_ceiling("CaH", sun) > hc.gas_ceiling("AsH3", sun)
    assert hc.gas_ceiling("NaH", _host(0.4)) == pytest.approx(2.51 * hc.gas_ceiling("NaH", sun), rel=0.01)
    for formula in hc.HOT_TRACE_GASES:
        assert 0.0 < hc.gas_ceiling(formula, sun) < 1e-3


def test_refusals():
    sun = _host()
    for pressure in (0.0, -1.0, float("nan")):
        with pytest.raises(ValueError):
            hc.cloud_temperature_k("Cr", pressure, sun)
    with pytest.raises(ValueError):
        hc.cloud_temperature_k("FeS", 1.0, sun)
    with pytest.raises(ValueError):
        hc.gas_ceiling("H2O", sun)
    with pytest.raises(ValueError):
        hc.abundance_per_h("Xe", sun)


def test_a_malformed_file_is_refused(tmp_path):
    good = json.load(open(hc.DATA_FILE, encoding="utf-8"))

    def refused(change):
        document = json.loads(json.dumps(good))
        change(document)
        path = tmp_path / "bad.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        with pytest.raises(DataError):
            hc.load(str(path))

    refused(lambda d: d["fits"][0].update(a=-1.0))
    refused(lambda d: d["fits"][0].update(terms=[]))
    refused(lambda d: d["fits"].append(dict(d["fits"][0])))                                  # a species twice
    refused(lambda d: d["kcl"].pop("gas_per_total"))
    refused(lambda d: d["trace_gases"][0].update(atoms={"Ge": 0}))
    refused(lambda d: d["trace_gases"].append(dict(d["trace_gases"][0])))
    refused(lambda d: d.update(trace_elements=[]))
    with pytest.raises(DataError):
        hc.load(str(tmp_path / "missing.json"))
    assert math.isfinite(hc.cloud_temperature_k("Cr", 1.0, _host()))
