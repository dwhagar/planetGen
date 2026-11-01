# tests/test_brown_dwarf_types.py

"""
GEN.517: the effective temperature against spectral type of brown dwarfs (`planetgen.physics.brown_dwarf_types`, [R2374: Table 13]): the map and its refusals, the temperatures of the type boundaries, the chance of each class at a temperature, and the
agreement of rows g1 to g5 of the atmosphere recipes with it.
"""

import json

import pytest

from planetgen.generation.possibility_table import Subject
from planetgen.physics import brown_dwarf_types as bd
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.physics.elements import DataError


def written(tmp_path, mutate):
    with open(bd.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    mutate(document)
    path = tmp_path / "map.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return str(path)


def test_the_map_gives_the_published_temperatures_of_the_types():
    for s, teff in ((0, 2238), (5, 1613), (8, 1336), (10, 1255), (14, 1182), (16, 966), (18, 686), (19, 567), (20, 460), (21, 368), (22, 288)):
        assert bd.teff_of_type(s) == pytest.approx(teff, abs=1.0)
    assert bd.scatter_k(3.0) == 134.0 and bd.scatter_k(12.0) == 79.0 and bd.scatter_k(20.0) == 79.0
    assert bd.teff_of_type(8.75) == pytest.approx(1277.7, abs=0.1), "the pieces meet"
    with pytest.raises(DataError, match="outside the map"):
        bd.teff_of_type(23.0)


def test_the_map_is_flat_between_l9_and_t5_so_the_type_is_loose_there():
    slope = bd.teff_of_type(14.0) - bd.teff_of_type(13.0)
    assert slope == pytest.approx(-18.3, abs=0.1) and abs(bd.scatter_k(12.0) / slope) > 4, "about four subtypes"


def test_the_loader_refuses_a_broken_map(tmp_path):
    with pytest.raises(DataError, match="not contiguous"):
        bd.load(written(tmp_path, lambda d: d["pieces"][1].update(s_min=9.0)))
    with pytest.raises(DataError, match="no finite coefficients"):
        bd.load(written(tmp_path, lambda d: d["pieces"][0].update(coefficients=[])))
    with pytest.raises(DataError, match="no scatter"):
        bd.load(written(tmp_path, lambda d: d["pieces"][0].update(rms_k=0)))
    with pytest.raises(DataError, match="outside the map"):
        bd.load(written(tmp_path, lambda d: d["boundaries"][0].update(s=40.0)))
    with pytest.raises(DataError, match="cannot read"):
        bd.load(str(tmp_path / "none.json"))
    with pytest.raises(DataError, match="none of"):
        bd.boundary_teff("x-y")


def test_the_boundaries_are_the_temperatures_of_the_types_they_lie_at():
    assert bd.boundary_teff("m-l") == pytest.approx(2237.5)
    assert bd.boundary_teff("l-t") == pytest.approx(1273.1, abs=0.1)
    assert bd.boundary_teff("t-late-t") == pytest.approx(686.2, abs=0.1)
    assert bd.boundary_teff("t-y") == pytest.approx(460.3, abs=0.1)


def test_the_chance_of_each_class_is_one_in_all_and_half_at_a_boundary():
    for teff in (300.0, 460.0, 686.0, 1000.0, 1270.0, 2238.0, 2700.0):
        probabilities = bd.type_probabilities(teff)
        assert sum(probabilities.values()) == pytest.approx(1.0) and all(p >= 0 for p in probabilities.values())
    at_y = bd.type_probabilities(bd.boundary_teff("t-y"))
    assert at_y["Y"] == pytest.approx(0.5, abs=0.01) and at_y["late T"] == pytest.approx(0.5, abs=0.01)
    assert bd.type_probabilities(1613.0)["L"] > 0.99 and bd.type_probabilities(966.0)["T"] > 0.99 and bd.type_probabilities(300.0)["Y"] > 0.97


@pytest.mark.parametrize("teff", [250.0, 400.0, 460.0, 560.0, 686.0, 800.0, 1000.0, 1273.0, 1500.0, 2000.0, 2238.0, 2500.0])
def test_rows_g1_to_g5_are_the_map_seen_from_the_temperature(teff):
    p = dict(RECIPES.probabilities("g", Subject("brown-dwarf", {"teff_k": teff})))
    expected = bd.type_probabilities(teff)
    for row, name in (("g2", "L"), ("g3", "T"), ("g4", "late T"), ("g5", "Y")):
        assert p.get(row, 0.0) == pytest.approx(expected[name], abs=2e-3), (teff, row)
    if teff < 2700:
        assert p.get("g1", 0.0) == pytest.approx(expected["late M"], abs=2e-3)
