# planetgen/describe/life_facts.py

"""
The facts of a world's life, resolved once, from which every sentence about that life is written.

The generator draws the chemistry, the pace and the timeline of a world once and stores them (`life.apply_life_data`).
`resolve` reads them into one `LifeFacts`, and `build.py` writes the life paragraph and the timeline paragraph from it
alone, so two sentences cannot name a different pace, milestone or age. A fact a sentence needs is a field here; a
sentence never reads the planet again.
"""

from dataclasses import dataclass

from planetgen import tuning
from planetgen.generation import evolution


@dataclass(frozen=True)
class Timeline:
    """The speculative timeline of a world: the system's age and the latest milestone reached by it."""
    system_age_gyr: float
    stage: str | None
    stage_name: str | None
    milestone_age_gyr: float | None
    stage_text: str | None


@dataclass(frozen=True)
class LifeFacts:
    """Everything the text says about a world's life."""
    noun: str
    chemical: str
    pace: str | None
    visible: tuple
    non_visible: tuple
    timeline: Timeline | None


def pace_label(scale):
    """How a pace is worded, from a `tuning.EVOLUTIONARY_TIMELINES` key: `fast` is "hyper-accelerated", `normal`
    "standard", `slow` "decelerated"."""
    return tuning.EVOLUTIONARY_TIMELINES[scale]["evolutionary_pace"].lower()


def resolve(planet):
    """The `LifeFacts` of a planet or moon, or `None` for a body with no life chemistry."""
    if not planet.life_chemical:
        return None
    timeline = None
    if planet.has_timeline:
        stage = planet.life_stage
        timeline = Timeline(
            planet.life_system_age_gyr, stage,
            None if stage is None else evolution.stage_name(stage, planet.life_chemical).lower(),
            planet.life_milestone_age_gyr, None if stage is None else evolution.stage_text(stage, planet.life_chemical))
    return LifeFacts(
        "moon" if planet.is_moon else "planet", planet.life_chemical.lower(),
        None if planet.evolutionary_pace is None else pace_label(planet.evolutionary_pace),
        tuple(planet.reflection_spectrum_visible or ()), tuple(planet.reflection_spectrum_non_visible or ()), timeline)
