# planetgen/physics/extinction_effects.py

"""
What a mass extinction does to a world's animals and civilization (GEN.1746, slice 4; design: reports/air-history-design-gen1739-2026-11-15.md, section 5.6)
==============================================================================================================================================================

`life_ages` cuts the history of a body at its mass extinctions. This module is the first place where an extinction changes the world: the Boss's rulings of 11-20 are that an extinction **never
kills microbial life** ("life finds a way") and that **a large kill may reset the animals and the civilization**. So the schedule of a class (`life_timeline.schedule`, one age for each of
`tuning.LIFE_MILESTONES`) keeps its first three steps (abiogenesis, photosynthesis, complex cells) whatever happens, and the last two, the animals (`multicellularity`) and the civilization, are drawn again
when a large enough extinction falls after the animals had come:

* *Which extinctions.* The kill curve (`air_events.KILL_CURVE`) gives 0.65 of the species for the events the Earth has once in 100 Myr and 0.90 for the end-Permian's once in 540 Myr, and every extinction of the Earth's animal record is of that size or less:
  the Big Five of the last 540 Myr (the largest impact among them about 10 km) left the animals standing. So an extinction that resets them is one the Earth's animals have not met: an impact the Earth has less
  often than once in `RESET_RECURRENCE_MYR`, 1 Gyr (an impactor of about 37 km on the size law of `air_events`). The flood basalts, which recur every 25 Myr, never do. [Stand-in 6.326 EF1, tier T6, tweak L1.]
* *Animals.* The animals of the world are gone at the extinction and return `air_events.RECOVERY_GYR` (the origination lag [R5321]) plus `REGROWTH_GYR` later. An extinction that falls before the animals came, or
  while they are still returning, changes nothing. [Stand-in 6.326 EF2, tier T7, tweak L1.]
* *Civilization.* One that existed is lost, one that had not come yet starts again: it arises the same hard-step wait (`life_timeline.schedule`, the draw it already had) after the animals are back. [Stand-in 6.326 EF3.]
* *Stored ages.* The ages the world stores are the ones of the lineage that is there now: the age of the animals is when they returned, the age of the civilization when it arose again.

Nothing is drawn here: the extinctions are the events of the body's own stream (`air_events.events_of`, the very call `life_ages.ages_of` makes), so a world made again has the same ones.
"""

from planetgen import tuning
from planetgen.physics import air_events, life_ages

RESET_RECURRENCE_MYR = 1000.0
"""float: How rarely the Earth has an extinction that resets the animals, Myr: one in 1 Gyr, over the 540 Myr of its animal record. [Stand-in 6.326 EF1, tier T6, tweak L1.]"""

REGROWTH_GYR = 0.020
"""float: The time, Gyr, the animals take to return after the origination lag (`air_events.RECOVERY_GYR`) that follows a reset: 20 Myr, so 30 Myr from the extinction. [Stand-in 6.326 EF2, tier T7 educated guess based on inferred physics, tweak L1.]"""

ANIMALS = tuning.LIFE_MILESTONES.index("multicellularity")
CIVILIZATION = tuning.LIFE_MILESTONES.index("technological_civilization")


def reset_impact_km(recurrence_myr=RESET_RECURRENCE_MYR):
    """The smallest impactor, km, the Earth has less often than once in `recurrence_myr` (the size law of `air_events` solved for it; the largest the law draws if none does)."""
    per_gyr = 1000.0 / recurrence_myr
    diameter = 5.0 * (air_events.IMPACT_FIVE_KM_PER_GYR / per_gyr) ** (1.0 / air_events.IMPACT_SIZE_EXPONENT)
    return min(diameter, air_events.IMPACT_LARGEST_KM)


def extinctions_of(body, begin_gyr, end_gyr=None):
    """`[(age_gyr, diameter_km)]` of the impacts of `body`'s history from `begin_gyr` (the age its first producers began) that reset the animals: those of `reset_impact_km` or more, in order of age. They are the events
    `life_ages.ages_of` cuts the life ages at (the same call on the body's own stream), so the extinction of an age in its table is one of these when it was this large."""
    star_age = getattr(getattr(body, "star", None), "age", None)
    end_gyr = end_gyr if end_gyr is not None else star_age
    if begin_gyr is None or not end_gyr or begin_gyr >= end_gyr:
        return []
    events = air_events.events_of(body, begin_gyr, end_gyr, min_impact_km=life_ages.smallest_impact_km())
    return [(event.age_gyr, event.size) for event in events if event.kind == air_events.IMPACT and event.size >= reset_impact_km()]


def regrown(ages, extinctions):
    """
    The schedule `ages` (one age for each `tuning.LIFE_MILESTONES`, `None` for a step never reached) after the `extinctions` (`extinctions_of`): the animals return `air_events.RECOVERY_GYR` plus `REGROWTH_GYR` after each
    extinction that falls once they are there, and the civilization arises again the same wait after them as before. The first three ages are never changed.
    """
    out = list(ages)
    for age, _diameter in sorted(extinctions):
        animals = out[ANIMALS]
        if animals is None or animals > age:
            continue
        wait = None if out[CIVILIZATION] is None else out[CIVILIZATION] - animals
        back = age + air_events.RECOVERY_GYR + REGROWTH_GYR
        out[ANIMALS] = back
        if wait is not None:
            out[CIVILIZATION] = back + wait
    return tuple(out)
