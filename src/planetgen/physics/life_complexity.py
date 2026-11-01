# planetgen/physics/life_complexity.py

"""
How complex the life of a class can get: the complexity index and the animal grades (GEN.890, GEN.877, GEN.882 to GEN.884)
===========================================================================================================================

A life class (`physics/life_classes.py`) reaches a rung of the ladder when its complexity index C passes a threshold
(docs/design/habitability/life-classes.md 7.3). C is the class's net primary production over Earth's (`ClassScore.npp`), times the yield of
what it breathes, times the temperature rate factor of its clock (`life_timeline.rate_factor`):

* the yield is 1 for aerobic respiration and about a ninth for none (`tuning.LIFE_COMPLEXITY["anaerobic_yield"]`); an oxygen breather's follows
  the oxygen of the air, from the anaerobic yield with none to 1 at `oxygen_full_yield_kpa`; any other oxidant has the class's own yield;
* `tuning.LIFE_COMPLEXITY["thresholds"]` give the C for complex cells, multicellular bodies and the animal grades A1 (sessile, sponge-like),
  A2 (small mobile) and A3 (large, active). A civilization needs A3.

An oxygen-breathing class also needs oxygen for the grade, whatever its C: `oxygen_floor_kpa`, about 0.1 to 0.8 kPa for sponge-grade animals
rising with body size (GEN.877; [R3084] and [R0689, R0691], recalled). A few anoxic worlds get animals that have lost aerobic
respiration (route R4, GEN.884): the draw is `generation/life.py`'s and the result a grade of 1.

Chemosymbiotic animals (Class 4, route R2) and sponge-grade animals on thin air (route R3) are not special cases: they are what C and the
oxygen floor give those classes. How much Fauna a grade shows is `tuning.LIFE_COMPLEXITY["fauna_cap"]`; `life_scores` applies it.
"""

import dataclasses
import math

from planetgen import tuning
from planetgen.physics import life_timeline

GRADES = ("A1", "A2", "A3")
"""tuple: The animal grades, 1 to 3 in order (grade 0 is no animals)."""


@dataclasses.dataclass(frozen=True)
class Complexity:
    """What a class on a world can reach: the index and its factors, the animal grade and the highest milestone."""
    index: float
    oxidant_yield: float
    rate_factor: float
    grade: int
    top_stage: str


def oxidant_yield(row, o2_kpa):
    """The yield of what the class `row` breathes, 0 to 1: an oxygen breather's follows the oxygen partial pressure `o2_kpa` of the air."""
    if row.oxidant == "o2":
        low = tuning.LIFE_COMPLEXITY["anaerobic_yield"]
        return low + (1.0 - low) * min(1.0, max(0.0, o2_kpa) / tuning.LIFE_COMPLEXITY["oxygen_full_yield_kpa"])
    return row.life["oxidant_yield"]


def index_of(row, npp, o2_kpa, temperature_k):
    """The complexity index C: `npp` times the oxidant yield times the temperature rate factor."""
    return npp * oxidant_yield(row, o2_kpa) * life_timeline.rate_factor(temperature_k)


def grade_by_index(row, index, o2_kpa):
    """The animal grade (0 to 3) the index allows, less where an oxygen breather's air is under the floor of the grade."""
    grade = 0
    for number, key in enumerate(GRADES, start=1):
        if index < tuning.LIFE_COMPLEXITY["thresholds"][key]:
            break
        if row.oxidant == "o2" and o2_kpa < tuning.LIFE_COMPLEXITY["oxygen_floor_kpa"][key]:
            break
        grade = number
    return grade


def top_stage(index, grade):
    """The highest `tuning.LIFE_MILESTONES` milestone the index and the grade allow: a civilization needs grade 3, multicellular bodies and
    complex cells their thresholds, and producers need nothing."""
    thresholds = tuning.LIFE_COMPLEXITY["thresholds"]
    if grade >= 3:
        return "technological_civilization"
    if index >= thresholds["multicellularity"]:
        return "multicellularity"
    if index >= thresholds["complex_cells"]:
        return "complex_cells"
    return "photosynthesis"


def assess(row, npp, o2_kpa, temperature_k, anoxic_animals=False):
    """The `Complexity` of class `row` on a world: its net primary production `npp`, the oxygen partial pressure `o2_kpa` of its air and its
    surface temperature. `anoxic_animals` (route R4) gives a world whose index allows multicellular bodies but no animal a grade of 1."""
    index = index_of(row, npp, o2_kpa, temperature_k)
    grade = grade_by_index(row, index, o2_kpa)
    if anoxic_animals and grade == 0 and index >= tuning.LIFE_COMPLEXITY["thresholds"]["multicellularity"]:
        grade = 1
    return Complexity(index, oxidant_yield(row, o2_kpa), life_timeline.rate_factor(temperature_k), grade, top_stage(index, grade))


def complex_factor(index):
    """The factor, 0 to 1, that the index gives the complex-life score: 0 at the multicellular threshold, 1 from `full_complexity`, on a log scale."""
    low = tuning.LIFE_COMPLEXITY["thresholds"]["multicellularity"]
    if index <= low:
        return 0.0
    return min(1.0, math.log10(index / low) / math.log10(tuning.LIFE_COMPLEXITY["full_complexity"] / low))
