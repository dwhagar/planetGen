# tests/test_life_coexist.py

"""GEN.1209: several life classes on one body (`generation.life`): the dominant class first, any other supported class that takes hold on a schedule of its own, kind scores that are the
highest of each over the classes, and the `planet_life` / `moon_life` rows that keep them."""

import copy

import pytest

from planetgen.db import store
from planetgen.generation import life
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import atmosphere, life_classes
from planetgen.util import draw

SEEDS = range(5000, 5016)


def _config():
    config = SystemConfig()
    config.PLANETS = True
    return config


def _bodies(system):
    for planet in system.all_planets:
        if planet.body_type == "a":
            continue
        yield from (planet, *planet.moons)


@pytest.fixture(scope="module")
def systems():
    made = []
    for seed in SEEDS:
        with draw.bound(seed):
            made.append(StarSystem(system_config=_config()))
    return made


@pytest.fixture(scope="module")
def mixed(systems):
    found = [body for system in systems for body in _bodies(system) if len(body.life_living or ()) > 1]
    if not found:
        pytest.skip("no body with several classes in these seeds")
    return found


def test_the_dominant_class_is_first_and_the_shares_add_up_to_one(systems):
    seen = 0
    for system in systems:
        for body in _bodies(system):
            if not body.life_class:
                assert body.life_living is None
                continue
            living = body.life_living
            assert living[0]["life_class"] == life_classes.class_id(body.life_class, body.life_class_variant)
            assert sum(entry["share"] for entry in living) == pytest.approx(1.0)
            assert len({life_classes.BY_ID[entry["life_class"]].number for entry in living}) == len(living), "one class number, one entry"
            seen += 1
    assert seen


def test_a_class_that_lives_beside_the_dominant_one_began_inside_the_star_age_and_is_not_poisoned(mixed):
    for body in mixed:
        scores = life_classes.score_classes(body)
        for entry in body.life_living[1:]:
            assert entry["ages"][0] <= body.star.age
            assert scores[entry["life_class"]].score > 0, "the air the classes made does not poison a class that lives"


def test_the_kind_scores_of_the_body_are_the_highest_of_each_over_its_classes(mixed):
    for body in mixed:
        by_kind = [max(entry["kinds"][kind] for entry in body.life_living) for kind in range(4)]
        assert by_kind == [body.life_microbial, body.life_flora, body.life_fauna, body.life_sentient]


def test_a_minor_class_has_a_stage_no_later_than_its_schedule_allows(mixed):
    order = {key: index for index, key in enumerate(life_classes.tuning.LIFE_MILESTONES)}
    for body in mixed:
        for entry in body.life_living[1:]:
            ages = entry["ages"]
            if entry["stage"] is not None:
                assert ages[order[entry["stage"]]] <= body.star.age
            assert entry["animal_grade"] in (0, 1, 2, 3)


def test_a_body_with_life_keeps_its_classes_through_an_edit(mixed):
    body = copy.deepcopy(mixed[0])
    before = [(entry["life_class"], entry["ages"]) for entry in body.life_living]
    life.apply_life_data(body, keep_stored=True)
    after = [(entry["life_class"], entry["ages"]) for entry in body.life_living]
    assert [c for c, _ in after] == [c for c, _ in before]
    assert after == before


def test_a_trial_air_leaves_the_body_as_it_was(mixed):
    body = mixed[0]
    clone = life._trial_copy(body)
    cands = [life._Cand(life_classes.BY_ID[entry["life_class"]], life_classes.score_classes(body)[entry["life_class"]], entry["share"], ages=tuple(entry["ages"]))
             for entry in body.life_living]
    before = copy.deepcopy(atmosphere.gas_partials_pa(body))
    life._build_air(clone, cands)
    assert atmosphere.gas_partials_pa(body) == before


def test_the_classes_are_saved_and_loaded_with_the_body(mixed, mysql_config):
    config = _config()   # saving needs the whole system: make one again that has such a body
    for seed in SEEDS:
        with draw.bound(seed):
            made = StarSystem(system_config=config)
        if any(len(b.life_living or ()) > 1 for b in _bodies(made)):
            break
    system_id = store.save_system(made, config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        loaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    stored = {(b.mass, b.distance): b.life_living for b in _bodies(made)}
    for body in _bodies(loaded):
        want = stored[(body.mass, body.distance)]
        got = body.life_living
        assert [e["life_class"] for e in (got or [])] == [e["life_class"] for e in (want or [])]
        for g, w in zip(got or [], want or []):
            assert g["share"] == pytest.approx(w["share"]) and g["stage"] == w["stage"] and g["kinds"] == w["kinds"]
            assert g["ages"] == pytest.approx(w["ages"], rel=1e-9) if all(a is not None for a in w["ages"]) else True
