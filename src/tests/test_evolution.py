"""
evolution.get_evolutionary_timeline regression tests.

The most advanced milestone a timeline reports is capped by `max_stage` (the habitability index's ceiling,
`generation.life.stage_ceiling`) and by INTELLIGENT_LIFE=False, never by a class letter (GEN.92, GEN.647): a forced
INTELLIGENT_LIFE=True at a very old age must still stop at the cap it is given.

Run with: pytest src/tests/test_evolution.py
"""
import re

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.evolution import (
    MILESTONE_KEYS, get_evolutionary_state, get_evolutionary_timeline, timeline_paragraphs,
)
from planetgen.util import draw
from planetgen.generation.star import Star
from planetgen import tuning as prog_c


def _star(intelligent_life=None, age=None):
    """A real main-sequence G2V Star -- get_evolutionary_timeline only
    reads type/yerkes_class/age/lifespan/system_config off it, and a G2V's
    `yerkes_class == 'V'` short-circuits get_star_evolutionary_profile
    straight to STAR_EVOLUTION['G'] without needing a finite lifespan
    special case (see that function's own docstring)."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.INTELLIGENT_LIFE = intelligent_life
    star = Star(cfg)
    if age is not None:
        star.age = age
    return star


def _milestone_text(evolutionary_data):
    assert len(evolutionary_data) == 1
    return evolutionary_data[0]


@pytest.mark.parametrize("max_stage,forbidden_milestones", [
    ("abiogenesis", ["Photosynthesis", "Complex Cells", "Multicellularity", "Technological Civilization"]),
    ("photosynthesis", ["Complex Cells", "Multicellularity", "Technological Civilization"]),
    ("multicellularity", ["Technological Civilization"]),
])
def test_a_capped_timeline_never_exceeds_its_ceiling_even_when_forced(max_stage, forbidden_milestones):
    """INTELLIGENT_LIFE=True at a very old age (every milestone is reachable) must still never report a milestone
    above the cap."""
    star = _star(intelligent_life=True, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star, max_stage))
    for forbidden in forbidden_milestones:
        assert forbidden not in text, f"{max_stage} reported {forbidden!r} despite its cap: {text}"


def test_an_uncapped_timeline_can_still_reach_technological_civilization():
    star = _star(intelligent_life=True, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star))
    assert "Technological Civilization" in text


def test_intelligent_life_false_downgrades_an_uncapped_timeline():
    star = _star(intelligent_life=False, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star))
    assert "Technological Civilization" not in text


def test_intelligent_life_false_does_not_raise_a_capped_timeline():
    # INTELLIGENT_LIFE=False alone caps at Multicellularity; a tighter cap still wins.
    star = _star(intelligent_life=False, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star, "photosynthesis"))
    assert "Multicellularity" not in text
    assert "Technological Civilization" not in text


def test_forced_milestone_for_a_capped_timeline_is_its_ceiling_not_no_milestone():
    star = _star(intelligent_life=True, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star, "photosynthesis"))
    assert "Photosynthesis" in text


# --- DB.51: the structure, and the paragraph written from it ---------------------------------------------------

_SYSTEM_AGE = re.compile(r"current estimated age of the system is ([0-9.]+) (Billion|Million) Years")
_MILESTONE = re.compile(r"would have been ([A-Za-z ]+) at ([0-9.]+) (Billion|Million) Years")


def _years(amount, unit):
    return float(amount) * {"Billion": 1.0, "Million": 1e-3}[unit]


@pytest.mark.parametrize("age", [0.2, 1.0, 2.0, 3.0, 4.2, 20.0])
@pytest.mark.parametrize("intelligent_life", [None, True, False])
def test_the_state_says_what_the_paragraph_says(age, intelligent_life):
    """The stage, pace and ages that are stored are the ones the old paragraph spelled out (to its rounding)."""
    star = _star(intelligent_life=intelligent_life, age=age)
    with draw.bound(51):
        state = get_evolutionary_state(star)
    with draw.bound(51):
        paragraph = get_evolutionary_timeline(star)
    assert paragraph == timeline_paragraphs(state)
    text = paragraph[0]
    system_age = _SYSTEM_AGE.search(text)
    assert _years(*system_age.groups()) == pytest.approx(state["system_age_gyr"], rel=0.01)
    assert prog_c.EVOLUTIONARY_TIMELINES[state["scale"]]["evolutionary_pace"].lower() in text
    milestone = _MILESTONE.search(text)
    if state["stage"] is None:
        assert milestone is None and state["milestone_age_gyr"] is None
        assert "No significant evolutionary milestones" in text
    else:
        assert milestone.group(1).lower().replace(" ", "_") == state["stage"]
        assert _years(milestone.group(2), milestone.group(3)) == pytest.approx(state["milestone_age_gyr"], rel=0.01)
        assert state["stage"] in MILESTONE_KEYS
