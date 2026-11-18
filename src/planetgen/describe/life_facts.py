# planetgen/describe/life_facts.py

"""
The facts of a world's life, resolved once, from which every sentence about that life is written.

The generator draws the life class, the pace and the timeline of a world once and stores them (`life.apply_life_data`).
`resolve` reads them into one `LifeFacts`, and `build.py` writes the life paragraph and the timeline paragraph from it
alone, so two sentences cannot name a different pace, milestone or age. A fact a sentence needs is a field here; a
sentence never reads the planet again.
"""

from dataclasses import dataclass

from planetgen import tuning
from planetgen.generation import evolution
from planetgen.physics import life_classes


@dataclass(frozen=True)
class Timeline:
    """The timeline of a world: the system's age and the latest milestone reached by it."""
    system_age_gyr: float
    stage: str | None
    stage_name: str | None
    milestone_age_gyr: float | None
    stage_text: str | None


@dataclass(frozen=True)
class LifeFacts:
    """Everything the text says about a world's life."""
    noun: str
    class_id: str
    class_name: str      # the name of the class before its colon (the pigment or habitat of a variant is in its id)
    pace: str | None
    visible: tuple
    non_visible: tuple
    timeline: Timeline | None
    habitat: str | None = None      # where a class that does not live at the surface lives, as a phrase ("aloft in the sulfuric-acid cloud layer, at 366 K and 290 kPa", GEN.1307)


def pace_label(scale):
    """How a pace is worded, from a `tuning.EVOLUTIONARY_TIMELINES` key: `fast` is "hyper-accelerated", `normal`
    "standard", `slow` "decelerated" (`tuning.EVOLUTIONARY_PACE_LABELS`)."""
    return tuning.EVOLUTIONARY_PACE_LABELS[scale]


def resolve(planet):
    """The `LifeFacts` of a planet or moon, or `None` for a body with no life class."""
    class_id = life_classes.class_id(planet.life_class, planet.life_class_variant)
    if class_id is None:
        return None
    grade = planet.life_animal_grade or 0
    timeline = None
    if planet.has_timeline:
        stage = planet.life_stage
        timeline = Timeline(
            planet.life_system_age_gyr, stage,
            None if stage is None else evolution.stage_name(stage, class_id, grade).lower(),
            planet.life_milestone_age_gyr, None if stage is None else evolution.stage_text(stage, class_id, grade))
    row = life_classes.describe(class_id)
    return LifeFacts(
        "moon" if planet.is_moon else "planet", class_id, row.name.split(":")[0],
        None if planet.evolutionary_pace is None else pace_label(planet.evolutionary_pace),
        tuple(planet.reflection_spectrum_visible or ()), tuple(planet.reflection_spectrum_non_visible or ()), timeline,
        life_classes.habitat_words(row, life_classes.values_from_body(planet)) if row.habitat else None)
