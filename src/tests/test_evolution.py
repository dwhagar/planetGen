"""
evolution.get_evolutionary_state regression tests.

The most advanced milestone a timeline reports is capped by `max_stage` (the habitability index's ceiling,
`life_complexity.top_stage`) and by INTELLIGENT_LIFE=False, never by a class letter (GEN.92, GEN.647). INTELLIGENT_LIFE=True forces no
milestone (GEN.1130): at a very old age the state is the latest stage on the schedule, under the cap it is given. GEN.750: the state reads the schedule the world drew
(`physics.life_timeline`) and draws nothing.

Run with: pytest src/tests/test_evolution.py
"""
import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation import evolution
from planetgen.generation.evolution import MILESTONE_KEYS, get_evolutionary_state
from planetgen.util import draw
from planetgen.generation.star import Star
from planetgen import tuning as prog_c
from planetgen.physics import life_classes


AGES = (0.1, 0.4, 1.0, 2.0, 3.0)
"""tuple: A schedule in Gyr, one age for each stage of `MILESTONE_KEYS`."""


def _star(intelligent_life=None, age=None):
    """A real main-sequence G2V Star -- get_evolutionary_state only
    reads type/yerkes_class/age/lifespan/system_config off it, and a G2V's
    `lifespan` is finite, so no white dwarf special case applies
    (see `supported_evolutionary_scales`)."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.INTELLIGENT_LIFE = intelligent_life
    star = Star(cfg)
    if age is not None:
        star.age = age
    return star


def _stage(star, max_stage=None):
    """The milestone `get_evolutionary_state` reports for a star on the schedule `AGES`."""
    return get_evolutionary_state(star, AGES, max_stage)["stage"]


@pytest.mark.parametrize("max_stage,forbidden_milestones", [
    ("abiogenesis", ["photosynthesis", "complex_cells", "multicellularity", "technological_civilization"]),
    ("photosynthesis", ["complex_cells", "multicellularity", "technological_civilization"]),
    ("multicellularity", ["technological_civilization"]),
])
def test_a_capped_timeline_never_exceeds_its_ceiling_when_a_civilization_is_requested(max_stage, forbidden_milestones):
    """INTELLIGENT_LIFE=True at a very old age (every milestone is reachable) never reports a milestone above the cap."""
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


def test_a_requested_civilization_reports_the_stage_the_schedule_reached():
    star = _star(intelligent_life=True, age=20.0)
    assert _stage(star, "photosynthesis") == "photosynthesis"
    assert get_evolutionary_state(star, (0.1, 0.5, None, None, None))["stage"] == "photosynthesis"      # nothing is forced past it


# --- GEN.750: the state reads a schedule and draws nothing ---------------------------------------------------------

@pytest.mark.parametrize("age", [0.2, 1.0, 2.0, 3.0, 4.2, 20.0])
@pytest.mark.parametrize("intelligent_life", [None, True, False])
def test_the_state_reads_the_schedule_it_was_given_and_draws_nothing(age, intelligent_life):
    star = _star(intelligent_life=intelligent_life, age=age)
    with draw.bound(51):
        state = get_evolutionary_state(star, AGES)
        next_draw = draw.random()
    with draw.bound(51):
        assert draw.random() == next_draw  # the state took nothing from the stream
    assert state["system_age_gyr"] == age
    assert "scale" not in state
    if state["stage"] is None:
        assert state["milestone_age_gyr"] is None
    else:
        assert state["milestone_age_gyr"] <= age


def test_the_stage_is_the_latest_one_whose_age_is_not_past_the_star():
    star = _star(age=1.5)
    state = get_evolutionary_state(star, AGES)
    assert state["stage"] == "complex_cells" and state["milestone_age_gyr"] == AGES[2]


def test_a_schedule_that_never_reaches_a_stage_stops_before_it():
    """A world that cannot make the oxygen has no age for complex cells and the stages after (`None`)."""
    star = _star(age=20.0)
    assert get_evolutionary_state(star, (0.1, 0.5, None, None, None))["stage"] == "photosynthesis"


@pytest.mark.parametrize("class_id", sorted(life_classes.BY_ID))
def test_every_life_class_has_words_for_each_milestone(class_id):
    for stage in MILESTONE_KEYS:
        assert evolution.stage_name(stage, class_id) and evolution.stage_text(stage, class_id)


def test_only_a_water_splitting_class_puts_oxygen_into_its_water():
    for class_id, row in life_classes.BY_ID.items():
        text = evolution.stage_text("photosynthesis", class_id)
        assert evolution.releases_oxygen(class_id) == row.releases_oxygen
        if row.releases_oxygen:
            assert "putting oxygen" in text
        else:
            assert "putting oxygen" not in text and "releasing oxygen" not in text


def test_a_grade_names_the_animals_of_multicellular_life():
    assert evolution.stage_name("multicellularity", "1a", 3) != evolution.stage_name("multicellularity", "1a", 0)


def _fake_star(letter="G", age=5.0, lifespan=10.0, yerkes="V"):
    import types
    return types.SimpleNamespace(type=f"{letter}2{yerkes}", yerkes_class=yerkes, age=age, lifespan=lifespan)


def test_the_paces_a_star_supports_come_from_its_own_lifespan_not_its_letter():
    """GEN.942: the paces whose technological civilization comes within the star's lifespan; a G star that dies at 10 Gy cannot be asked to wait for
    the 35 Gy of the slow pace, and one that lives 40 Gy can."""
    from planetgen.generation.stellar_object import supported_evolutionary_scales as paces
    assert paces(_fake_star("G", lifespan=10.0)) == ["fast", "normal"]
    assert paces(_fake_star("G", lifespan=40.0)) == ["fast", "normal", "slow"]
    assert paces(_fake_star("M", lifespan=2000.0)) == ["fast", "normal", "slow"]
    assert paces(_fake_star("A", lifespan=2.2)) == ["fast"], "the model's A star lives 0.65 to 2.2 Gy, short of the 4.5 Gy of the normal pace"
    assert paces(_fake_star("B", lifespan=0.05)) == ["fast"], "a star that gives less time than the fastest pace needs keeps the fastest"


def test_a_white_dwarf_has_the_time_it_has_had_and_a_remnant_letter_has_no_paces():
    from planetgen.generation.stellar_object import supported_evolutionary_scales as paces
    assert paces(_fake_star("F", age=5.0, lifespan=float("inf"), yerkes="VII")) == ["fast", "normal"]
    assert paces(_fake_star("D", age=5.0, lifespan=float("inf"), yerkes="VII")) == []


def test_the_per_letter_lifespan_table_is_gone():
    from planetgen import tuning
    assert not hasattr(tuning, "STAR_EVOLUTION")
