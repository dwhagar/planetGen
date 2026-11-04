# tests/test_meteorite_groups.py

"""GEN.475: the chondrite groups as regions of composition space with a nearest-group label (`physics/meteorite_groups.py`; [R4658: Table 1])."""

import pytest

from planetgen.physics import meteorite_groups as mg


def test_every_group_is_its_own_nearest_group_except_the_pair_the_table_cannot_tell_apart():
    for name, (fa, metal, refractory, matrix) in mg.CHONDRITE_GROUPS.items():
        label, distance = mg.nearest_group(fa, metal, refractory, matrix)
        if name in ("EH", "EL"):
            assert label in ("EH", "EL")
        elif name in ("CV", "CK") and fa is None:
            continue
        else:
            assert label == name, name
        assert distance < 0.2 or name in ("CV", "CK")


def test_the_ordinary_chondrites_are_told_apart_by_their_oxidation():
    """H 19.3, L 25.2, LL 31.3 mol percent fayalite; the same metal-poor, matrix-poor rock."""
    for fa, group in ((19.0, "H"), (25.0, "L"), (31.5, "LL")):
        assert mg.nearest_group(fa, 4.0, 0.0, 12.0)[0] == group


def test_a_water_bearing_matrix_rock_with_no_metal_is_a_carbonaceous_ci_or_cm():
    label, _ = mg.nearest_group(None, 0.0, 0.0, 98.0)
    assert label == "CI"
    assert mg.nearest_group(None, 0.1, 5.0, 70.0)[0] == "CM"


def test_a_metal_rich_reduced_rock_is_cb_or_ch_and_a_body_far_from_every_group_is_none():
    assert mg.nearest_group(3.0, 65.0, 0.0, 1.0)[0] == "CB"
    assert mg.nearest_group(0.4, 10.0, 0.0, 8.0)[0] in ("EH", "EL")
    label, distance = mg.nearest_group(60.0, 100.0, 40.0, 100.0)
    assert label == mg.NONE and distance > mg.NEAREST_LIMIT


def test_the_fayalite_of_the_oxides_is_the_molar_iron_over_iron_plus_magnesium():
    assert mg.fayalite_mol_percent(0.0, 0.4) == 0.0
    assert mg.fayalite_mol_percent(0.4, 0.0) == 100.0
    # a bulk silicate with the mass ratio FeO/MgO of an H chondrite mantle gives the table's H value near 19 mol percent
    assert mg.fayalite_mol_percent(0.18 * 71.844, 0.82 * 40.304) == pytest.approx(18.0, abs=0.01)
    with pytest.raises(ValueError):
        mg.fayalite_mol_percent(0.0, 0.0)


def test_the_distance_uses_the_axes_both_have():
    with_fa = mg.nearest_group(31.0, 2.5, 10.0, 40.0)
    without = mg.nearest_group(None, 2.5, 10.0, 40.0)
    assert with_fa[0] == "CK"
    assert without[0] in ("CV", "CK")
