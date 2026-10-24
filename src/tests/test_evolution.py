"""
evolution.get_evolutionary_state regression tests.

Covers the fix for a bug where a planet's fixed per-class flavor text (e.g.
Class G: "a rocky, barren world with simple life") and its evolutionary/
civilization-tech milestone were rolled by two entirely independent code
paths -- `planet_class` was never passed into the timeline
at all -- so an old/fast-evolving star, or a forced `INTELLIGENT_LIFE=True`,
could still produce "Technological Civilization" for a planet whose own
class text says life there tops out at "simple"/"bacterial". Fixed via
`program_constants.PLANET_CLASS_MAX_LIFE_STAGE`.

Run with: pytest src/tests/test_evolution.py
"""
import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation import evolution
from planetgen.generation.evolution import MILESTONE_KEYS, draw_scale, get_evolutionary_state
from planetgen.util import draw
from planetgen.generation.star import Star
from planetgen import tuning as prog_c


def _star(intelligent_life=None, age=None):
    """A real main-sequence G2V Star -- get_evolutionary_state only
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


def _stage(star, planet_class=None):
    """The milestone `get_evolutionary_state` reports for a star at the pace drawn for it."""
    return get_evolutionary_state(star, planet_class, draw_scale(star))["stage"]


@pytest.mark.parametrize("capped_class,forbidden_milestones", [
    ("E", ["photosynthesis", "complex_cells", "multicellularity", "technological_civilization"]),
    ("F", ["complex_cells", "multicellularity", "technological_civilization"]),
    ("G", ["complex_cells", "multicellularity", "technological_civilization"]),
    ("L", ["technological_civilization"]),
])
def test_capped_class_never_exceeds_its_own_ceiling_even_when_forced(capped_class, forbidden_milestones):
    """
    The direct regression case: INTELLIGENT_LIFE=True at a very old age
    (every milestone, including Technological Civilization, is reachable)
    must still never report a milestone above the class's own cap.
    """
    star = _star(intelligent_life=True, age=20.0)
    stage = _stage(star, capped_class)
    assert stage not in forbidden_milestones, f"Class {capped_class} reported {stage!r} despite its own cap"


def test_uncapped_class_can_still_reach_technological_civilization():
    # Regression guard for the fix itself over-restricting: a class with no
    # entry in PLANET_CLASS_MAX_LIFE_STAGE (M: "a terrestrial Earth-like
    # world") must behave exactly as before -- INTELLIGENT_LIFE=True at an
    # old age still reports a full civilization.
    star = _star(intelligent_life=True, age=20.0)
    assert _stage(star, "M") == "technological_civilization"


def test_none_planet_class_stays_uncapped():
    # Same as above but exercising the explicit planet_class=None default,
    # for a caller that doesn't know/care about class.
    star = _star(intelligent_life=True, age=20.0)
    assert _stage(star) == "technological_civilization"


def test_intelligent_life_false_still_downgrades_uncapped_class():
    # Pre-existing behavior, unchanged by this fix: INTELLIGENT_LIFE=False
    # never reports a civilization even for an otherwise-uncapped class.
    star = _star(intelligent_life=False, age=20.0)
    assert _stage(star, "M") != "technological_civilization"


def test_intelligent_life_false_does_not_raise_a_capped_class_above_its_own_ceiling():
    # A class capped below Multicellularity (F: "photosynthesis") combined
    # with INTELLIGENT_LIFE=False (which alone would only cap at
    # Multicellularity) must still respect the tighter of the two caps.
    star = _star(intelligent_life=False, age=20.0)
    assert _stage(star, "F") in ("abiogenesis", "photosynthesis")


def test_forced_milestone_for_capped_class_is_its_own_ceiling_not_no_milestone():
    # INTELLIGENT_LIFE=True on a capped class shouldn't silently report
    # "no significant milestone" either -- it should force the highest
    # stage that class's own cap actually allows.
    star = _star(intelligent_life=True, age=20.0)
    assert _stage(star, "G") == "photosynthesis"


@pytest.mark.parametrize("planet_class", list(prog_c.PLANET_CLASS_MAX_LIFE_STAGE.keys()))
def test_every_capped_class_is_a_real_milestone_key(planet_class):
    # Guards PLANET_CLASS_MAX_LIFE_STAGE itself against a typo'd milestone
    # key that MILESTONE_KEYS.index() would raise ValueError on at
    # generation time.
    assert prog_c.PLANET_CLASS_MAX_LIFE_STAGE[planet_class] in MILESTONE_KEYS


# --- GEN.666: one pace -----------------------------------------------------------------------------------------

@pytest.mark.parametrize("age", [0.2, 1.0, 2.0, 3.0, 4.2, 20.0])
@pytest.mark.parametrize("intelligent_life", [None, True, False])
def test_the_state_carries_the_pace_it_was_given_and_draws_nothing(age, intelligent_life):
    """The state is built at the pace `draw_scale` drew; it never draws one of its own."""
    star = _star(intelligent_life=intelligent_life, age=age)
    for scale in prog_c.EVOLUTIONARY_TIMELINES:
        with draw.bound(51):
            state = get_evolutionary_state(star, "M", scale)
            next_draw = draw.random()
        with draw.bound(51):
            assert draw.random() == next_draw  # the state took nothing from the stream
        assert state["scale"] == scale
        assert state["system_age_gyr"] == age
        if state["stage"] is None:
            assert state["milestone_age_gyr"] is None
        else:
            assert state["milestone_age_gyr"] <= age or intelligent_life is True


def test_the_pace_follows_the_chemistry_when_the_star_allows_it():
    """Melanin allows only the standard pace, and a G star supports it, so a melanin world is never "fast"."""
    star = _star(age=5.0)
    for seed in range(40):
        with draw.bound(seed):
            assert draw_scale(star, ["normal"]) == "normal"


def test_a_star_that_supports_none_of_the_chemistry_scales_decides_alone():
    star = _star(age=5.0)
    supported = evolution.get_star_evolutionary_profile(star)["supported_evolutionary_scales"]
    outside = [scale for scale in prog_c.EVOLUTIONARY_TIMELINES if scale not in supported]
    if not outside:
        pytest.skip("this star supports every pace")
    for seed in range(20):
        with draw.bound(seed):
            assert draw_scale(star, outside) in supported


@pytest.mark.parametrize("chemical", sorted(prog_c.LIFE_CHEMICALS))
def test_every_chemistry_has_a_metabolism_with_words_for_its_first_light_stage(chemical):
    metabolism = evolution.metabolism(chemical)
    assert metabolism in prog_c.LIFE_METABOLISMS
    assert evolution.stage_name("photosynthesis", chemical) and evolution.stage_text("photosynthesis", chemical)


def test_only_a_water_splitting_chemistry_releases_oxygen():
    for chemical in prog_c.LIFE_CHEMICALS:
        text = evolution.stage_text("photosynthesis", chemical)
        if evolution.metabolism(chemical) == "oxygenic":
            assert "releasing oxygen" in text
        else:
            assert "release no oxygen" in text and "releasing oxygen" not in text
