"""
evolution.get_evolutionary_timeline regression tests.

Covers the fix for a bug where a planet's fixed per-class flavor text (e.g.
Class G: "a rocky, barren world with simple life") and its evolutionary/
civilization-tech milestone were rolled by two entirely independent code
paths -- `planet_class` was never passed into `get_evolutionary_timeline`
at all -- so an old/fast-evolving star, or a forced `INTELLIGENT_LIFE=True`,
could still produce "Technological Civilization" for a planet whose own
class text says life there tops out at "simple"/"bacterial". Fixed via
`program_constants.PLANET_CLASS_MAX_LIFE_STAGE`.

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


@pytest.mark.parametrize("capped_class,forbidden_milestones", [
    ("E", ["Photosynthesis", "Complex Cells", "Multicellularity", "Technological Civilization"]),
    ("F", ["Complex Cells", "Multicellularity", "Technological Civilization"]),
    ("G", ["Complex Cells", "Multicellularity", "Technological Civilization"]),
    ("L", ["Technological Civilization"]),
])
def test_capped_class_never_exceeds_its_own_ceiling_even_when_forced(capped_class, forbidden_milestones):
    """
    The direct regression case: INTELLIGENT_LIFE=True at a very old age
    (every milestone, including Technological Civilization, is reachable)
    must still never report a milestone above the class's own cap.
    """
    star = _star(intelligent_life=True, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star, capped_class))
    for forbidden in forbidden_milestones:
        assert forbidden not in text, f"Class {capped_class} reported {forbidden!r} despite its own cap: {text}"


def test_uncapped_class_can_still_reach_technological_civilization():
    # Regression guard for the fix itself over-restricting: a class with no
    # entry in PLANET_CLASS_MAX_LIFE_STAGE (M: "a terrestrial Earth-like
    # world") must behave exactly as before -- INTELLIGENT_LIFE=True at an
    # old age still reports a full civilization.
    star = _star(intelligent_life=True, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star, "M"))
    assert "Technological Civilization" in text


def test_none_planet_class_stays_uncapped():
    # Same as above but exercising the explicit planet_class=None default,
    # for a caller that doesn't know/care about class (matches this
    # function's pre-fix signature and behavior exactly).
    star = _star(intelligent_life=True, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star))
    assert "Technological Civilization" in text


def test_intelligent_life_false_still_downgrades_uncapped_class():
    # Pre-existing behavior, unchanged by this fix: INTELLIGENT_LIFE=False
    # never reports a civilization even for an otherwise-uncapped class.
    star = _star(intelligent_life=False, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star, "M"))
    assert "Technological Civilization" not in text


def test_intelligent_life_false_does_not_raise_a_capped_class_above_its_own_ceiling():
    # A class capped below Multicellularity (F: "photosynthesis") combined
    # with INTELLIGENT_LIFE=False (which alone would only cap at
    # Multicellularity) must still respect the tighter of the two caps.
    star = _star(intelligent_life=False, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star, "F"))
    assert "Multicellularity" not in text
    assert "Technological Civilization" not in text


def test_forced_milestone_for_capped_class_is_its_own_ceiling_not_no_milestone():
    # INTELLIGENT_LIFE=True on a capped class shouldn't silently report
    # "no significant milestone" either -- it should force the highest
    # stage that class's own cap actually allows.
    star = _star(intelligent_life=True, age=20.0)
    text = _milestone_text(get_evolutionary_timeline(star, "G"))
    assert "Photosynthesis" in text


@pytest.mark.parametrize("planet_class", list(prog_c.PLANET_CLASS_MAX_LIFE_STAGE.keys()))
def test_every_capped_class_is_a_real_milestone_key(planet_class):
    # Guards PLANET_CLASS_MAX_LIFE_STAGE itself against a typo'd milestone
    # key that MILESTONE_KEYS.index() would raise ValueError on at
    # generation time.
    assert prog_c.PLANET_CLASS_MAX_LIFE_STAGE[planet_class] in MILESTONE_KEYS


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
        state = get_evolutionary_state(star, "M")
    with draw.bound(51):
        paragraph = get_evolutionary_timeline(star, "M")
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
