"""
GEN.666: the sentences about a world's life state one pace, one chemistry, one set of milestones and ages.

A planet page once said "life based on melanin would evolve at a normal pace" and, a paragraph later, "on a speculative
timeline, life here would evolve at a hyper-accelerated pace ... technological civilization, reached at 80 Myr": the pace
was drawn twice and written by two writers. The generator now draws it once (`evolution.draw_scale`), and
`describe.life_facts.resolve` reads every fact into one `LifeFacts` from which every sentence is written. These tests
fail if two sentences of a generated page disagree, and show that the checker itself catches a disagreement.
"""

import itertools
import re

import pytest

from planetgen import tuning
from planetgen.describe import build as describe_build, life_facts, model
from planetgen.describe.quantity import inline_text
from planetgen.generation import evolution
from planetgen.generation import life as planet_life
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem

from tests.fuzz_support import deterministic_entropy

_PACE_LABELS = {life_facts.pace_label(scale): scale for scale in tuning.EVOLUTIONARY_TIMELINES}
_PACE_WORDS = re.compile(r"\b(" + "|".join(re.escape(label) for label in _PACE_LABELS) + r") pace\b")
_STAGE_NAMES = {evolution.stage_name(stage, chemical).lower()
                for stage in evolution.MILESTONE_KEYS for chemical in tuning.LIFE_CHEMICALS}


def _bodies(system):
    for planet in system.planets + system.secondary_planets:
        if planet.body_type == "a":
            continue
        yield planet
        yield from getattr(planet, "moons", ())


def _life_text(body):
    """The life and timeline paragraphs of a body's section, as the page shows them."""
    section = describe_build._planet_section(body)
    blocks = [block for block in section["blocks"] if block["type"] == "prose" and block["role"] in ("life", "timeline")]
    return {block["role"]: " ".join(inline_text(sentence) for sentence in block["text"]) for block in blocks}, blocks


def conflicts(body, paragraphs):
    """What the paragraphs `_life_text` returned for `body` say that its stored life facts do not."""
    found = []
    text = " ".join(paragraphs.values())
    if not body.life_chemical:
        return [f"life is described on a body with no chemistry: {text!r}"] if text else []
    said = set(_PACE_WORDS.findall(text))
    expected = set() if body.evolutionary_pace is None else {life_facts.pace_label(body.evolutionary_pace)}
    if said != expected:
        found.append(f"pace: the text says {sorted(said)}, the body has {sorted(expected)}")
    if f"based on {body.life_chemical.lower()}" not in text:
        found.append(f"chemistry: {body.life_chemical!r} is not named")
    for other in tuning.LIFE_CHEMICALS:
        if other != body.life_chemical and f"based on {other.lower()}" in text:
            found.append(f"chemistry: {other!r} is named on a {body.life_chemical!r} world")
    timeline = paragraphs.get("timeline")
    if timeline is None:
        if body.has_timeline and body.planet_class in tuning.HABITABLE_PLANET_CLASSES:
            found.append("timeline: the body has one and the page shows none")
        return found
    if not body.has_timeline:
        return found + ["timeline: the page shows one for a body that has none"]
    stage_names = {name for name in _STAGE_NAMES if f"milestone would be {name}," in timeline}
    if body.life_stage is None:
        if stage_names or "No significant milestone" not in timeline:
            found.append(f"milestone: none was reached and the text says {sorted(stage_names)}")
    else:
        expected_name = evolution.stage_name(body.life_stage, body.life_chemical).lower()
        if stage_names != {expected_name}:
            found.append(f"milestone: the text says {sorted(stage_names)}, the body reached {expected_name!r}")
        if body.life_milestone_age_gyr > body.life_system_age_gyr + 1e-9:
            found.append("milestone: reached after the system's own age")
        if evolution.stage_text(body.life_stage, body.life_chemical) not in timeline:
            found.append("milestone: the biosphere text is another stage's")
    return found


def _generated_bodies(count, **overrides):
    for seed in range(count):
        config = SystemConfig()
        for key, value in overrides.items():
            setattr(config, key, value)
        with deterministic_entropy(7000 + seed):
            system = StarSystem(system_config=config)
        yield from (body for body in _bodies(system) if body.life_chemical)


def test_no_generated_page_contradicts_itself_about_its_life():
    seen = with_timeline = 0
    for overrides in ({"HABITABLE_WORLD": True, "BINARY_SYSTEM": False}, {"HABITABLE_WORLD": True, "INTELLIGENT_LIFE": True},
                      {"INTELLIGENT_LIFE": False}, {"MOONS": True}):
        for body in _generated_bodies(30, **overrides):
            paragraphs, _blocks = _life_text(body)
            assert conflicts(body, paragraphs) == [], (body.name, paragraphs)
            seen += 1
            with_timeline += body.has_timeline
    assert seen > 10 and with_timeline > 5


def test_the_pace_is_one_the_chemistry_allows_whenever_the_star_supports_one():
    seen = 0
    for body in _generated_bodies(40, HABITABLE_WORLD=True):
        allowed = planet_life.get_chemical_scales(body.life_chemical)
        supported = evolution.get_star_evolutionary_profile(body.star)["supported_evolutionary_scales"]
        if set(allowed) & set(supported):
            assert body.evolutionary_pace in allowed, (body.name, body.life_chemical, body.evolutionary_pace)
        assert body.evolutionary_pace in supported
        seen += 1
    assert seen > 5


def test_the_timeline_ages_in_the_page_are_the_stored_ones():
    seen = 0
    for body in _generated_bodies(40, HABITABLE_WORLD=True, INTELLIGENT_LIFE=True):
        if not body.has_timeline:
            continue
        _text, blocks = _life_text(body)
        timeline = next((block for block in blocks if block["role"] == "timeline"), None)
        if timeline is None:
            continue
        quantities = [item for sentence in timeline["text"] for item in sentence if isinstance(item, dict) and "qty" in item]
        expected = [describe_build._gy(body.life_system_age_gyr)]
        if body.life_stage is not None:
            expected.append(describe_build._gy(body.life_milestone_age_gyr))
        assert quantities == expected
        seen += 1
    assert seen > 3


# --- every combination of the facts is worded from them alone ----------------------------------------------------

def _facts(scale, chemical, stage, system_age=13.8, milestone_age=0.08, with_timeline=True):
    timeline = None
    if with_timeline:
        timeline = life_facts.Timeline(
            system_age, stage, None if stage is None else evolution.stage_name(stage, chemical).lower(),
            None if stage is None else milestone_age, None if stage is None else evolution.stage_text(stage, chemical))
    spectra = tuning.LIFE_CHEMICALS[chemical]
    return life_facts.LifeFacts("planet", chemical.lower(), None if scale is None else life_facts.pace_label(scale),
                                tuple(spectra["reflection_spectrum_visible"]),
                                tuple(spectra["reflection_spectrum_non_visible"]), timeline)


@pytest.mark.parametrize("scale, chemical", list(itertools.product(tuning.EVOLUTIONARY_TIMELINES, tuning.LIFE_CHEMICALS)))
def test_the_pace_is_worded_once_and_the_timeline_follows_it(scale, chemical):
    for stage in (None, *evolution.MILESTONE_KEYS):
        facts = _facts(scale, chemical, stage)
        life = " ".join(inline_text(sentence) for sentence in describe_build._life_sentences(facts))
        timeline = " ".join(inline_text(sentence) for sentence in describe_build._timeline_sentences(facts))
        assert _PACE_WORDS.findall(life) == [life_facts.pace_label(scale)]
        assert _PACE_WORDS.findall(timeline) == [], "the timeline says no pace of its own"
        assert timeline.startswith("On a speculative timeline at that pace, the system is about 13.8 Gyr old.")


def test_a_world_with_no_pace_says_none_and_the_timeline_does_not_refer_to_one():
    facts = _facts(None, "Melanin", "photosynthesis")
    life = " ".join(inline_text(sentence) for sentence in describe_build._life_sentences(facts))
    timeline = " ".join(inline_text(sentence) for sentence in describe_build._timeline_sentences(facts))
    assert "pace" not in life and "at that pace" not in timeline


def test_no_life_is_no_sentence():
    assert describe_build._life_sentences(None) == []


def test_a_melanin_world_is_not_told_it_releases_oxygen():
    facts = _facts("normal", "Melanin", "photosynthesis")
    timeline = " ".join(inline_text(sentence) for sentence in describe_build._timeline_sentences(facts))
    assert "radiotrophy" in timeline and "release no oxygen" in timeline and "releasing oxygen" not in timeline


# --- the checker itself ----------------------------------------------------------------------------------------

def _a_world_with_a_timeline():
    for body in _generated_bodies(60, HABITABLE_WORLD=True, INTELLIGENT_LIFE=True):
        paragraphs, _blocks = _life_text(body)
        if body.has_timeline and body.life_stage is not None and "timeline" in paragraphs:
            return body, paragraphs
    pytest.fail("no generated world with a timeline")


def test_the_checker_catches_a_changed_pace_milestone_chemistry_or_missing_timeline():
    body, paragraphs = _a_world_with_a_timeline()
    assert conflicts(body, paragraphs) == []
    label = life_facts.pace_label(body.evolutionary_pace)
    other = next(name for name in _PACE_LABELS if name != label)
    assert conflicts(body, {**paragraphs, "life": paragraphs["life"].replace(f"{label} pace", f"{other} pace")})
    assert conflicts(body, {**paragraphs, "timeline": paragraphs["timeline"] + f" It runs at a {other} pace."})
    name = evolution.stage_name(body.life_stage, body.life_chemical).lower()
    wrong = next(n for n in _STAGE_NAMES if n != name)
    assert conflicts(body, {**paragraphs, "timeline": paragraphs["timeline"].replace(f"would be {name},", f"would be {wrong},")})
    swapped = next(c for c in tuning.LIFE_CHEMICALS if c != body.life_chemical)
    assert conflicts(body, {**paragraphs, "life": paragraphs["life"].replace(
        f"based on {body.life_chemical.lower()}", f"based on {swapped.lower()}")})
    assert conflicts(body, {"life": paragraphs["life"]})


def test_a_describe_document_of_a_generated_system_is_still_a_checked_document():
    with deterministic_entropy(7001):
        system = StarSystem(system_config=SystemConfig())
    from planetgen.describe import describe
    model.check(describe(system))
