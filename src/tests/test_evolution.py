"""
evolution.get_evolutionary_state regression tests.

The most advanced milestone a timeline reports is capped by `max_stage` (the habitability index's ceiling,
`generation.life.stage_ceiling`) and by INTELLIGENT_LIFE=False, never by a class letter (GEN.92, GEN.647): a forced
INTELLIGENT_LIFE=True at a very old age must still stop at the cap it is given. GEN.666 adds the pace: one draw
(`draw_scale`) the state is built at.

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


def _stage(star, max_stage=None):
    """The milestone `get_evolutionary_state` reports for a star at the pace drawn for it."""
    return get_evolutionary_state(star, draw_scale(star), max_stage)["stage"]


@pytest.mark.parametrize("max_stage,forbidden_milestones", [
    ("abiogenesis", ["photosynthesis", "complex_cells", "multicellularity", "technological_civilization"]),
    ("photosynthesis", ["complex_cells", "multicellularity", "technological_civilization"]),
    ("multicellularity", ["technological_civilization"]),
])
def test_a_capped_timeline_never_exceeds_its_ceiling_even_when_forced(max_stage, forbidden_milestones):
    """INTELLIGENT_LIFE=True at a very old age (every milestone is reachable) must still never report a milestone
    above the cap."""
    star = _star(intelligent_life=True, age=20.0)
    stage = _stage(star, max_stage)
    assert stage not in forbidden_milestones, f"{max_stage} reported {stage!r} despite its cap"


def test_an_uncapped_timeline_can_still_reach_technological_civilization():
    star = _star(intelligent_life=True, age=20.0)
    assert _stage(star) == "technological_civilization"


def test_intelligent_life_false_downgrades_an_uncapped_timeline():
    star = _star(intelligent_life=False, age=20.0)
    assert _stage(star) != "technological_civilization"


def test_intelligent_life_false_does_not_raise_a_capped_timeline():
    # INTELLIGENT_LIFE=False alone caps at Multicellularity; a tighter cap still wins.
    star = _star(intelligent_life=False, age=20.0)
    assert _stage(star, "photosynthesis") in ("abiogenesis", "photosynthesis")


def test_forced_milestone_for_a_capped_timeline_is_its_ceiling_not_no_milestone():
    star = _star(intelligent_life=True, age=20.0)
    assert _stage(star, "photosynthesis") == "photosynthesis"


# --- GEN.666: one pace -----------------------------------------------------------------------------------------

@pytest.mark.parametrize("age", [0.2, 1.0, 2.0, 3.0, 4.2, 20.0])
@pytest.mark.parametrize("intelligent_life", [None, True, False])
def test_the_state_carries_the_pace_it_was_given_and_draws_nothing(age, intelligent_life):
    """The state is built at the pace `draw_scale` drew; it never draws one of its own."""
    star = _star(intelligent_life=intelligent_life, age=age)
    for scale in prog_c.EVOLUTIONARY_TIMELINES:
        with draw.bound(51):
            state = get_evolutionary_state(star, scale)
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
