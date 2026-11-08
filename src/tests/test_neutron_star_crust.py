# tests/test_neutron_star_crust.py

"""
GEN.1358, GEN.1427: the crust of a neutron star by density (`physics/neutron_star_crust.py`, `data/neutron_star_crust.json`) and its reader on the page.

The table is [R0722: Tables 2 and 3] (READ): 15 nuclei of the outer crust with the density each holds to, neutron drip at 4.27e11 g/cm3 and
11 cells of the inner crust, and the predicted nuclear pasta layer at its bottom (GEN.1427). The tests check the sequence against the printed one, that every row adds up, the refusals, and that the
neutron star page shows it.
"""

import json

import pytest

from planetgen.describe import describe, model
from planetgen.describe.quantity import inline_text
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import NeutronStar
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from planetgen.physics import neutron_star_crust as crust
from planetgen.physics.elements import DataError


def test_the_outer_crust_runs_from_iron_to_krypton_118_in_the_printed_order():
    names = [nucleus.name for nucleus in crust.OUTER_CRUST]
    assert names == ["iron-56", "nickel-62", "nickel-64", "nickel-66", "krypton-86", "selenium-84", "germanium-82", "zinc-80", "zinc-82",
                     "palladium-128", "ruthenium-126", "molybdenum-124", "zirconium-122", "strontium-120", "krypton-118"]
    assert crust.OUTER_CRUST[0].max_density_g_cm3 == 8.02e6 and crust.OUTER_CRUST[-1].max_density_g_cm3 == 4.27e11
    assert crust.NEUTRON_DRIP_G_CM3 == 4.27e11
    for nucleus in crust.OUTER_CRUST[4:9]:                  # the nuclei with 50 neutrons, zinc-82 with 52
        assert nucleus.neutrons in (50, 52)
    assert all(nucleus.neutrons == 82 for nucleus in crust.OUTER_CRUST[9:])
    assert [n.mass_data for n in crust.OUTER_CRUST][:9] == ["experimental masses (Ruster+2006)"] + ["experimental"] * 8


def test_the_inner_crust_is_zirconium_then_tin_then_zirconium_and_germanium():
    assert [cell.name for cell in crust.INNER_CRUST] == ["zirconium"] * 5 + ["tin"] * 4 + ["zirconium", "germanium"]
    assert crust.INNER_CRUST[0].density_g_cm3 == 4.67e11 and crust.CRUST_BOTTOM_G_CM3 == 1.32e14
    assert crust.INNER_CRUST[0].nucleons_in_cell == 180 and crust.INNER_CRUST[0].protons == 40
    assert crust.inner_runs() == [("zirconium", 4.67e11, 2.66e12), ("tin", 6.24e12, 3.41e13), ("zirconium", 7.94e13, 7.94e13), ("germanium", 1.32e14, 1.32e14)]


def test_every_row_adds_up_and_the_densities_rise():
    for nucleus in crust.OUTER_CRUST:
        assert nucleus.protons + nucleus.neutrons == nucleus.mass_number
    for cell in crust.INNER_CRUST:
        assert cell.protons + cell.neutrons_in_cell == cell.nucleons_in_cell
    densities = [row[2] for row in crust.sequence()]
    assert densities == sorted(densities) and len(crust.sequence()) == 15 + 1 + 11


def test_a_malformed_file_is_refused(tmp_path):
    good = json.load(open(crust.DATA_FILE, encoding="utf-8"))

    def refused(change):
        document = json.loads(json.dumps(good))
        change(document)
        path = tmp_path / "bad.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        with pytest.raises(DataError):
            crust.load(str(path))

    refused(lambda d: d["outer_crust"][0].update(neutrons=31))                       # 26 + 31 is not 56
    refused(lambda d: d["outer_crust"][0].update(symbol="Co"))                       # Z = 26 is iron
    refused(lambda d: d["inner_crust"][1].update(density_g_cm3=1.0))                 # not rising
    refused(lambda d: d.update(neutron_drip_g_cm3=1e12))                             # drip is the last outer layer's density
    refused(lambda d: d.update(outer_crust=d["outer_crust"][:5]))
    with pytest.raises(DataError):
        crust.load(str(tmp_path / "missing.json"))


def test_the_pasta_layer_is_five_shapes_inside_the_crust_core_boundary():
    pasta = crust.PASTA
    assert [item.name for item in pasta.shapes] == ["gnocchi", "spaghetti", "lasagna", "anti-spaghetti", "anti-gnocchi"]
    assert pasta.evidence == "predicted only"
    assert pasta.lower_g_cm3 == pytest.approx(0.2 * 0.165 * crust.G_CM3_PER_FM3)
    assert crust.INNER_CRUST[-1].density_g_cm3 < pasta.upper_g_cm3 <= 1.4e14
    assert pasta.lower_g_cm3 < pasta.upper_g_cm3 and 0.0 < pasta.crust_mass_share_at_most <= 1.0 and pasta.thickness_m > 0.0


def test_a_malformed_pasta_block_is_refused(tmp_path):
    good = json.load(open(crust.DATA_FILE, encoding="utf-8"))
    for change in (lambda p: p.update(shapes=p["shapes"][:2]),                      # fewer than three shapes
                   lambda p: p.update(lower_fraction_of_saturation=9.0),            # starts above the crust-core boundary
                   lambda p: p.update(crust_mass_share_at_most=1.5)):               # more than the whole crust
        document = json.loads(json.dumps(good))
        change(document["pasta"])
        path = tmp_path / "bad.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        with pytest.raises(DataError):
            crust.load(str(path))
    document = json.loads(json.dumps(good))
    del document["pasta"]
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(DataError):
        crust.load(str(path))


def test_the_neutron_star_page_shows_the_crust_by_density():
    draw.set_run_seed(3)
    config = SystemConfig()
    star = NeutronStar(config)
    document = describe(StarSystem(system_config=config, compact_remnant=star))
    model.check(document)
    section = next(s for s in model.sections(document) if s["kind"] == "neutron-star")
    block = next(b for b in section["blocks"] if b["type"] == "facts" and b["caption"] == "Crust")
    rows = {row["key"]: row["value"] for row in block["rows"]}
    assert list(rows) == ["crust_outer", "crust_drip", "crust_inner", "crust_pasta"]
    outer = [item["qty"] for item in rows["crust_outer"] if isinstance(item, dict)]
    assert len(outer) == 15 and all(q["kind"] == "density" and q["unit"] == "g_cm3" for q in outer)
    assert outer[0]["value"] == 8.02e6 and outer[-1]["value"] == 4.27e11
    text = inline_text(rows["crust_outer"])
    assert "iron-56 to" in text and "krypton-118 to" in text
    inner = inline_text(rows["crust_inner"])
    assert "zirconium from" in inner and "tin from" in inner and "germanium at" in inner
    assert [item["qty"]["value"] for item in rows["crust_drip"] if isinstance(item, dict)] == [4.27e11]
    pasta = inline_text(rows["crust_pasta"])
    assert pasta.startswith("predicted, not observed") and "100 m" in pasta
    assert pasta.index("gnocchi") < pasta.index("spaghetti") < pasta.index("lasagna") < pasta.index("anti-spaghetti") < pasta.rindex("anti-gnocchi")
