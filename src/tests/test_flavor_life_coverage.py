# tests/test_flavor_life_coverage.py

"""GEN.1431, GEN.1438: the life subtree of the flavor tree offers at least three lines for every medium and energy source at every rung from 2 up that the life reaches, offers a vent world
none of the lines of a sunlit ocean, and a world with microbes only now gets a life line."""

from types import SimpleNamespace

import pytest

from planetgen.physics import flavor_tree, life_classes

# The highest rung each solvent medium and energy source reaches, from the `reaches` column of `life_classes.json`: water and light, and water and chemicals,
# reach animals; the acid, carbon dioxide and hydrocarbon classes reach mats at most; ammonia and water reaches microbes (a world past that is still offered the mat lines).
REACH = {("water", "light"): 5, ("water", "chemical"): 5, ("sulfuric", "chemical"): 4, ("co2", "chemical"): 4, ("hydrocarbon", "chemical"): 4, ("ammonia_water", "chemical"): 4}
STAGE_OF_RUNG = {2: "photosynthesis", 3: "complex_cells", 4: "multicellularity", 5: "multicellularity", 6: "technological_civilization"}
OXIDIZED_ID = {"water": {"light": "1a", "chemical": "3a"}}


def _body(class_id, rung):
    row = life_classes.BY_ID[class_id]
    return SimpleNamespace(body_type="t", surface_liquid=None, life_class=row.number, life_class_id=class_id, life_stage=STAGE_OF_RUNG[rung], life_animal_grade=2 if rung >= 5 else 0,
                           mass=5.0e24, radius=6400.0, distance=1.0, atmospheric_pressure=1.0e5)


def _classes_of(medium, energy):
    return [row.id for row in life_classes.CLASSES if row.solvent["medium"] == medium and row.energy["source"] == energy]


def _life_options(class_id, rung):
    return [o for o in flavor_tree.walk("planet", _body(class_id, rung), off=("ancient_works",)) if o.id.startswith(("life:", "habitable:"))]


@pytest.mark.parametrize("group", sorted(REACH))
def test_every_medium_and_energy_has_three_lines_at_every_rung_it_reaches(group):
    medium, energy = group
    classes = _classes_of(medium, energy)
    assert classes, group
    for class_id in classes:
        for rung in range(2, REACH[group] + 1):
            options = _life_options(class_id, rung)
            assert len(options) >= 3, (class_id, rung, len(options))


@pytest.mark.parametrize("class_id", [row.id for row in life_classes.CLASSES])
def test_a_world_with_any_life_of_rung_2_or_more_is_offered_a_life_line(class_id):
    for rung in (2, 3):
        assert _life_options(class_id, rung), (class_id, rung)


def test_a_vent_world_is_not_offered_the_lines_of_a_sunlit_ocean():
    sunlit = {o.id for o in _life_options("1a", 5)}
    vent = {o.id for o in _life_options("4a", 5)}
    assert sunlit and vent and not (sunlit & vent)


def test_a_world_with_no_animals_is_not_offered_the_animal_lines():
    animals = {o.id for o in _life_options("1a", 5)} - {o.id for o in _life_options("1a", 4)}
    assert animals and not ({o.id for o in _life_options("1a", 4)} & animals)


def test_a_sentient_world_is_offered_signs_of_a_civilization():
    assert any(o.id.startswith("life:") and "ships" in o.text for o in _life_options("1a", 6))
