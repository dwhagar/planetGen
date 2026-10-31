# tests/test_life_standing.py

"""
GEN.873: the standing life tests of the life chemistry review (reports/life-chemistry-review-2026-10-28.md, section 9.3), read in terms of
the life classes now that they have replaced the chemistries (GEN.889, docs/design/habitability/life-classes.md section 11.2).

Tests 2, 7 and 8 of the review already stand in class terms and are not repeated: the pace reaches complex cells inside the age of the
universe (`test_life_timeline.py`), every life constant has a registry row (`test_life_standin_registry.py`,
`test_life_class_standin_rows.py`) and an oxygenic world past photosynthesis has more O2 than before it (`test_biogenic_oxygen.py`). Test 3
(the prevalence band per class) belongs to the prevalence suite (TEST.430). The four here:

1. The star never picks the pigment from a table: the class rules hold no star letter, a light class takes the star's light by its band
   (`life_classes.band_ratio`), every light class gets light from some star of the main sequence and every star letter feeds some light class.
2. A body that has a life class has a life stage, and a body with a stage has a class.
3. A body with Flora is not under an ice-covered ocean (no sunlight reaches it).
4. A class lives only where its solvent is: a water class on a body whose liquid is water, a hydrocarbon class on a hydrocarbon sea, and
   the CO2, sulfuric acid and ammonia classes in their own media.

Every rule is a function of one finished body, run over the shared body sweep (`tests.audit.sweep`) with the life step of the generator
(`generation.life.apply_life_data`) put on a copy of each terrestrial body, and over the bodies of the sweep's generated systems. Each rule is
shown to fail on a body doctored to break it, because a rule that cannot fail is not a test.
"""

import copy
import functools
import json

import pytest

from planetgen.generation import life
from planetgen.physics import habitability_world, life_classes
from planetgen.util import draw
from tests.audit import sweep

STAR_LETTERS = "OBAFGKM"
TEMPERATURE_K = {"O": 35000.0, "B": 20000.0, "A": 8500.0, "F": 6500.0, "G": 5800.0, "K": 4500.0, "M": 3000.0}
"""dict: A representative effective temperature of each star letter of the main sequence."""

LIFE_STREAM_SEED = 873


def _rules():
    with open(life_classes.DATA_FILE, encoding="utf-8") as handle:
        return json.load(handle)


def _keys(node):
    """Every key of a JSON tree."""
    if isinstance(node, dict):
        for key, value in node.items():
            yield key
            yield from _keys(value)
    elif isinstance(node, list):
        for value in node:
            yield from _keys(value)


# --- 1. The star picks no pigment from a table ----------------------------------------------------------------------------------------

def test_the_class_rules_hold_no_star_letter_table():
    keys = set(_keys(_rules()["classes"]))
    assert not keys & set(STAR_LETTERS), f"a class rule is keyed by a star letter: {sorted(keys & set(STAR_LETTERS))}"
    assert not [key for key in keys if "star" in key.lower() and "weight" in key.lower()]


def _light_classes():
    return [row for row in life_classes.CLASSES if row.energy["source"] == "light"]


def test_the_light_a_class_gets_from_a_star_is_the_photons_of_its_band():
    for row in _light_classes():
        ratios = [life_classes.band_ratio(row, {"_teff_k": TEMPERATURE_K[letter]}) for letter in STAR_LETTERS]
        assert all(0.0 <= ratio <= 1.0 for ratio in ratios), (row.id, ratios)
        assert max(ratios) > 0.0, f"class {row.id} gets light from no star of the main sequence: a pigment that can never occur"


def test_every_star_letter_feeds_some_light_class():
    for letter in STAR_LETTERS:
        best = max(life_classes.band_ratio(row, {"_teff_k": TEMPERATURE_K[letter]}) for row in _light_classes())
        assert best > 0.0, f"no light class gets light from a {letter} star"


def test_a_chemical_class_takes_no_light_from_the_star():
    for row in life_classes.CLASSES:
        if row.energy["source"] == "chemical":
            assert life_classes.band_ratio(row, {"_teff_k": 3000.0}) == 1.0, row.id


def test_the_pigment_follows_the_light_of_the_star_not_its_letter():
    """The band ratio of a class moves with the star's temperature alone: two stars of one temperature give one answer whatever their letter."""
    for row in _light_classes():
        assert life_classes.band_ratio(row, {"_teff_k": 5800.0}) == life_classes.band_ratio(row, {"_teff_k": 5800.0, "star_class": "M"}), row.id
    cool = {row.id: life_classes.band_ratio(row, {"_teff_k": 3000.0}) for row in _light_classes()}
    hot = {row.id: life_classes.band_ratio(row, {"_teff_k": 20000.0}) for row in _light_classes()}
    assert cool != hot, "the pigment ratios are the same for a 3,000 K and a 20,000 K star"


# --- The bodies the rules 2 to 4 read ---------------------------------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _lived():
    """`[(label, body)]`: the terrestrial planets and moons of the sweep with the life step run on a copy, and the terrestrial bodies of the
    sweep's generated systems as the generator left them."""
    bodies = []
    for entry in sweep.entries():
        if getattr(entry.body, "body_type", None) != "t":
            continue
        clone = copy.deepcopy(entry.body)
        draw.set_run_seed(LIFE_STREAM_SEED)
        life.apply_life_data(clone)
        bodies.append((entry.label, clone))
    for entry in sweep.system_entries():
        if getattr(entry.body, "body_type", None) == "t":
            bodies.append((entry.label, entry.body))
    return tuple(bodies)


def _class_id(body):
    return life_classes.class_id(getattr(body, "life_class", None), getattr(body, "life_class_variant", None))


def _with_class(body):
    return _class_id(body) is not None


# --- 2. A class and a stage go together -------------------------------------------------------------------------------------------------

def stage_problems(body):
    """The problems of a body's class and stage: a class with no stage, or a stage with no class."""
    return class_without_stage(body) + stage_without_class(body)


def class_without_stage(body):
    if _with_class(body) and getattr(body, "life_stage", None) is None:
        return [f"life class {_class_id(body)} with no life stage"]
    return []


def stage_without_class(body):
    if getattr(body, "life_stage", None) is not None and not _with_class(body):
        return [f"life stage {body.life_stage!r} with no life class"]
    return []


def test_a_body_with_a_life_class_has_a_life_stage():
    problems = [f"{label}: {problem}" for label, body in _lived() for problem in class_without_stage(body)]
    assert not problems, "\n".join(problems[:20])


def test_a_body_with_a_life_stage_has_a_life_class():
    problems = [f"{label}: {problem}" for label, body in _lived() for problem in stage_without_class(body)]
    assert not problems, "\n".join(problems[:20])


def test_the_stage_rule_can_fail():
    body = next(body for _label, body in _lived() if _with_class(body) and body.life_stage is not None)
    assert stage_problems(body) == []
    assert class_without_stage(sweep.broken_copy(body, life_stage=None))
    unclassed = next(body for _label, body in _lived() if not _with_class(body))
    assert stage_without_class(sweep.broken_copy(unclassed, life_stage="abiogenesis"))


def test_the_sweep_gives_life_to_some_bodies_and_not_to_others():
    classed = sum(1 for _label, body in _lived() if _with_class(body))
    assert 0 < classed < len(_lived()), f"{classed} of {len(_lived())} sweep bodies have a life class: the rules above would hold for nothing"


# --- 3. No Flora under ice --------------------------------------------------------------------------------------------------------------

def flora_problems(body):
    """Flora on a body whose ocean is under ice: no light reaches the producers, so only the chemical classes' microbial mats can be there."""
    if getattr(body, "life_flora", 0) and getattr(body, "hydrosphere", None) == "ice-covered ocean":
        return [f"Flora {body.life_flora} under an ice-covered ocean (class {_class_id(body)})"]
    return []


def test_no_body_has_flora_under_an_ice_covered_ocean():
    problems = [f"{label}: {problem}" for label, body in _lived() for problem in flora_problems(body)]
    assert not problems, "\n".join(problems[:20])


def test_no_light_class_is_supported_under_an_ice_covered_ocean():
    """The sweep has few bodies with Flora, so the other half of the rule is read off the model: put every terrestrial body of the sweep under
    an ice-covered ocean and no class that eats light may be supported there (a chemical class's mats can)."""
    light = {row.id for row in life_classes.CLASSES if row.light}
    problems = []
    for entry in sweep.entries():
        if getattr(entry.body, "body_type", None) != "t":
            continue
        scores = life_classes.score_classes(sweep.broken_copy(entry.body, hydrosphere="ice-covered ocean")) or {}
        problems += [f"{entry.label}: light class {class_id} scores {result.score} under ice" for class_id, result in scores.items()
                     if class_id in light and result.score > 0]
    assert not problems, "\n".join(problems[:20])


def test_the_flora_rule_can_fail():
    body = next(body for _label, body in _lived() if getattr(body, "hydrosphere", None) == "ice-covered ocean" or _with_class(body))
    assert flora_problems(sweep.broken_copy(body, hydrosphere="ice-covered ocean", life_flora=3))
    assert not flora_problems(sweep.broken_copy(body, hydrosphere="ice-covered ocean", life_flora=0))
    assert not flora_problems(sweep.broken_copy(body, hydrosphere="surface ocean", life_flora=3))


# --- 4. A class lives only in its own solvent -------------------------------------------------------------------------------------------

def solvent_problems(body):
    """The problems of the solvent under the class that lives on a body: the liquid the class needs is not there."""
    class_id = _class_id(body)
    if class_id is None:
        return []
    row = life_classes.describe(class_id)
    medium = row.solvent["medium"]
    values = habitability_world.values_from_body(body)
    found = life_classes.quantities(values)
    found.update(life_classes.derived_quantities(values))
    if medium in ("water", "hydrocarbon"):
        have = habitability_world.solvent(values)
        return [] if have == medium else [f"class {class_id} needs {medium} and the body's liquid is {have!r}"]
    if medium == "co2":
        phase = found.get("co2_phase")
        return [] if phase in ("liquid", "supercritical") else [f"class {class_id} needs liquid CO2 and the body's is {phase!r}"]
    quantity, line = ("acid_medium", 1.0) if medium == "sulfuric" else ("ammonia_water_mass_fraction", 0.1)
    stored = found.get(quantity)
    return [] if stored is not None and stored >= line else [f"class {class_id} needs {quantity} of {line} or more and the body has {stored!r}"]


def test_no_class_lives_on_a_solvent_it_does_not_allow():
    problems = [f"{label}: {problem}" for label, body in _lived() for problem in solvent_problems(body)]
    assert not problems, "\n".join(problems[:20])


def test_the_solvent_rule_can_fail():
    water = next((body for _label, body in _lived() if _class_id(body) and life_classes.describe(_class_id(body)).solvent["medium"] == "water"), None)
    assert water is not None, "no body of the sweep carries a water class"
    assert solvent_problems(water) == []
    dry = sweep.broken_copy(water, hydrosphere="vapour", aquifer_thickness_km=0.0, surface_temperature=900.0)
    assert solvent_problems(dry)


def test_no_class_is_supported_where_its_solvent_is_missing():
    """The class that lives on a body is the dominant one, so the drawn bodies show the water classes only; every class that scores above 0 on
    a body of the sweep, drawn or not, must have its solvent there (the hydrocarbon, CO2, sulfuric and ammonia classes included)."""
    problems = []
    for entry in sweep.entries():
        body = entry.body
        if getattr(body, "body_type", None) != "t":
            continue
        for class_id, result in (life_classes.score_classes(body) or {}).items():
            if result.score > 0:
                problems += [f"{entry.label}: {problem}" for problem in solvent_problems(sweep.broken_copy(body, life_class=int(class_id[0]), life_class_variant=class_id[1:] or None))]
    assert not problems, "\n".join(problems[:20])


def test_the_solvent_rule_fails_a_class_of_another_medium():
    body = next(entry.body for entry in sweep.entries() if getattr(entry.body, "body_type", None) == "t"
                and habitability_world.solvent(habitability_world.values_from_body(entry.body)) != "hydrocarbon")
    for number, variant in ((7, None), (6, None), (5, "a"), (8, None)):
        assert solvent_problems(sweep.broken_copy(body, life_class=number, life_class_variant=variant)), (number, variant)
