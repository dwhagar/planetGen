"""GEN.1430: the flavor text tree (`physics/flavor_tree.py`, `physics/data/flavor_tree.json`): the file, the walk, the draw."""
import json
from types import SimpleNamespace

import pytest

from planetgen.physics import flavor_tree, oddities
from planetgen.physics.elements import DataError


def _body(body_type="t", surface_liquid=None, life_stage=None, life_class=None, life_class_id=None, life_animal_grade=0, **extra):
    return SimpleNamespace(body_type=body_type, surface_liquid=surface_liquid, life_class=life_class, life_class_id=life_class_id, life_animal_grade=life_animal_grade, life_stage=life_stage, mass=5.0e24, radius=6400.0, distance=1.0, **extra)


def test_the_file_keeps_every_old_line_under_the_key_it_was_stored_under():
    ids = flavor_tree.load()[1]
    legacy = [i for i in ids if i.split(":")[0] in ("system", "planet", "orbital", "habitable")]
    assert len(legacy) == 82 and len(set(ids)) == len(ids)
    assert {i.split(":")[0] for i in ids} == {"system", "planet", "orbital", "habitable", "life", "belt", "star", "comet", "blackhole", "pulsar", "nebula", "remnant", "quasar", "icomet", "field"}
    assert flavor_tree.text_of("planet:3").startswith("a solar-powered beacon")
    assert flavor_tree.text_of(None) is None


def test_no_gate_names_a_class_letter_or_a_quantity_that_does_not_exist():
    document = json.load(open(flavor_tree.DATA_FILE, encoding="utf-8"))
    known = set(json.load(open(oddities.DATA_FILE, encoding="utf-8"))["quantities"]) | set(flavor_tree.QUANTITIES) | set(flavor_tree.BELT_QUANTITIES)
    known |= {name for tree in document["scopes"].values() for name in tree.get("quantities", {})}     # the kinds that name their own (star_teff_k, ...)

    def visit(node):
        for condition in node.get("gate") or []:
            name = condition["quantity"]
            assert "class" not in name.split("_") and name in known, name
        for child in node.get("children") or []:
            visit(child)
    for tree in document["scopes"].values():
        visit(tree)


def test_the_walk_offers_only_the_lines_whose_gates_pass_and_prunes_a_tagged_family():
    ground = flavor_tree.walk("planet", _body())
    assert ground and all(option.id.startswith("planet:") for option in ground)
    assert flavor_tree.walk("planet", _body(), off=("ancient_works",)) == []
    life = flavor_tree.walk("planet", _body(life_stage="multicellularity", life_class=1, life_class_id="1a", life_animal_grade=2), off=("ancient_works",))
    assert life and all(option.id.startswith(("habitable:", "life:")) for option in life)


def test_the_draw_depends_on_the_body_alone_and_skips_the_ids_already_used():
    body = _body()
    first = flavor_tree.draw_line("planet", body)
    assert first == flavor_tree.draw_line("planet", _body())                    # the same body, however the system made its other bodies
    second = flavor_tree.draw_line("planet", body, used_ids=[first.id])
    assert second.id != first.id
    every = [option.id for option in flavor_tree.walk("planet", body)]
    assert flavor_tree.draw_line("planet", body, used_ids=every) is None         # an empty list is no line, with no retry


def test_a_malformed_file_is_refused(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps({"scopes": {"system": {"id": "system", "lines": [{"id": "a", "text": "x"}, {"id": "a", "text": "y"}]}, "planet": {"id": "planet"}}}))
    with pytest.raises(DataError):
        flavor_tree.load(str(path))
    path.write_text(json.dumps({"scopes": {"system": {"id": "system", "gate": [{"quantity": "mass_earth"}]}, "planet": {"id": "planet"}}}))
    with pytest.raises(DataError):
        flavor_tree.load(str(path))
