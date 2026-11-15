# tests/test_asteroid_vocabulary.py

"""GEN.1369: each of the 84 asteroid component names has a host mineral in a primitive body's mix or is tagged vocabulary-only, and every key a mix uses is a row id."""

import pytest

from planetgen import tuning
from planetgen.physics import accessory_minerals, asteroid_vocabulary as av, materials, small_body_composition as sbc

DRAWN = {row for taxon in sbc.TAXON_DATA.values() for row in taxon.fractions} | {row for rows in accessory_minerals.ACCESSORIES.values() for row in rows}


def test_every_one_of_the_84_names_is_hosted_or_vocabulary_only_and_never_both():
    assert len(tuning.ASTEROID_COMPONENTS) == 84
    for name in tuning.ASTEROID_COMPONENTS:
        row = av.canonical(name)
        assert (row in av.HOSTED) != (row in av.VOCABULARY_ONLY), f"{name} must be in exactly one of the host table and the vocabulary-only list"
    assert av.HOSTED and av.VOCABULARY_ONLY


def test_the_tables_name_no_row_the_vocabulary_lacks():
    rows = {av.canonical(name) for name in tuning.ASTEROID_COMPONENTS}
    assert set(av.HOSTED) | set(av.VOCABULARY_ONLY) == rows


def test_every_host_is_a_row_a_taxon_mix_or_accessory_list_draws():
    for name, hosts in av.HOSTED.items():
        assert hosts, name
        missing = [row for row in hosts if row not in DRAWN]
        assert not missing, f"{name} is hosted by rows no taxon mix draws: {missing}"


def test_a_vocabulary_only_name_is_not_drawn_and_carries_a_reason():
    for name, why in av.VOCABULARY_ONLY.items():
        assert name not in DRAWN, f"{name} is in a taxon mix: host it instead of tagging it"
        assert len(why) > 40


def test_every_key_a_mix_uses_is_a_row_id_and_not_an_alias():
    keys = {row for taxon in sbc.TAXON_DATA.values() for row in taxon.fractions} | {row for rows in accessory_minerals.ACCESSORIES.values() for row in rows}
    keys |= set(sbc.VOLATILE_ROWS.values()) | set(sbc.METAL_ROWS) | set(sbc.REFRACTORY_ROWS) | set(sbc.MATRIX_ROWS)
    aliased = {key: materials.material(key).id for key in keys if materials.material(key).id != key}
    assert not aliased, f"written by an alias: {aliased}"


def test_the_platinum_group_elements_sit_in_the_metal_and_the_organics_hold_the_pahs():
    for name in ("iridium", "osmium", "platinum", "palladium", "ruthenium", "rhodium"):
        assert av.HOSTED[av.canonical(name)] == ("kamacite", "taenite")
    assert av.HOSTED["coronene"] == ("complex_organics",)
    assert av.canonical("silicon dioxide") == "quartz" and "quartz" in av.VOCABULARY_ONLY
    assert av.HOSTED[av.canonical("plagioclase feldspars")] == ("albite", "anorthite")


def test_hosted_lists_the_names_a_composition_can_carry_with_the_largest_host():
    found = av.hosted({"kamacite": 0.1, "taenite": 0.3, "complex_organics": 0.05, "forsterite": 0.5})
    assert found["iridium"] == "taenite" and found["coronene"] == "complex_organics" and found["silicon"] == "forsterite"
    assert "kamacite" not in found and "hibonite" not in found and "osmiridium" not in found
    assert list(found) == sorted(found)
    assert av.hosted({}) == {}


def test_a_belt_stores_the_names_its_taxon_mixes_can_carry():
    makeup = sbc.draw_belt(2.7, 3.828e26, 5.0e-5)
    hosted = makeup.detail["hosted"]
    assert hosted and set(hosted) <= set(av.HOSTED)
    assert all(row in makeup.components for row in hosted.values())


def test_an_unknown_name_is_refused():
    with pytest.raises(materials.UnknownMaterial):
        av.canonical("a made-up mineral")
