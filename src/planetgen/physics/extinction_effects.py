# planetgen/physics/extinction_effects.py

"""
What a mass extinction does to a world's animals and civilization (GEN.1746, slice 4; design: reports/air-history-design-gen1739-2026-11-15.md, section 5.6)
==============================================================================================================================================================

`life_ages` cuts the history of a body at its mass extinctions. This module is the first place where an extinction changes the world: the Boss's rulings of 11-20 are that an extinction **never
kills microbial life** ("life finds a way") and that **a large kill may reset the animals and the civilization**. So the schedule of a class (`life_timeline.schedule`, one age for each of
`tuning.LIFE_MILESTONES`) keeps its first three steps (abiogenesis, photosynthesis, complex cells) whatever happens, and the last two, the animals (`multicellularity`) and the civilization, are drawn again
when a large enough extinction falls after the animals had come:

* *Which extinctions.* Those that kill `RESET_KILL_SHARE` of the species or more: the top of the kill curve (`air_events.KILL_CURVE`), the events the Earth has once in 100 Myr or less often (the largest
  impacts; the flood basalts, which recur every 25 Myr, kill less). [Stand-in 6.326 EF1, tier T6, tweak L1.]
* *Animals.* The animals of the world are gone at the extinction and return `air_events.RECOVERY_GYR` (the origination lag [R5321]) plus `REGROWTH_GYR` later. An extinction that falls before the animals came, or
  while they are still returning, changes nothing. [Stand-in 6.326 EF2, tier T7, tweak L1.]
* *Civilization.* One that existed is lost, one that had not come yet starts again: it arises the same hard-step wait (`life_timeline.schedule`, the draw it already had) after the animals are back.
* *Stored ages.* The ages the world stores are the ones of the lineage that is there now: the age of the animals is when they returned, the age of the civilization when it arose again. The extinction
  itself stays in the life ages table.

Nothing is drawn here: the extinctions are those of the body's own stream (`life_ages.ages_of`), the same ones its table keeps, so a world made again has the same ones.
"""

from planetgen import tuning
from planetgen.physics import air_events

RESET_KILL_SHARE = 0.65
"""float: The share of species an extinction must kill to reset the animals and the civilization: the top of the kill curve (`air_events.KILL_CURVE`, the events that recur every 100 Myr or less often). [Stand-in 6.326 EF1, tier T6, tweak L1.]"""

REGROWTH_GYR = 0.020
"""float: The time, Gyr, the animals take to return after the origination lag (`air_events.RECOVERY_GYR`) that follows a reset: 20 Myr, so 30 Myr from the extinction. [Stand-in 6.326 EF2, tier T7 educated guess based on inferred physics, tweak L1.]"""

ANIMALS = tuning.LIFE_MILESTONES.index("multicellularity")
CIVILIZATION = tuning.LIFE_MILESTONES.index("technological_civilization")


def reset_extinctions(life_age_entries):
    """`[(age_gyr, kill_share)]` of the extinctions of a body's life ages (`life_ages.entries_of`) that kill at least `RESET_KILL_SHARE`, in order of age."""
    return sorted((entry["end_gyr"], entry["kill_share"]) for entry in life_age_entries or ()
                  if entry["end_gyr"] is not None and entry["kill_share"] is not None and entry["kill_share"] >= RESET_KILL_SHARE)


def regrown(ages, extinctions):
    """
    The schedule `ages` (one age for each `tuning.LIFE_MILESTONES`, `None` for a step never reached) after the `extinctions` (`reset_extinctions`): the animals return `air_events.RECOVERY_GYR` plus `REGROWTH_GYR` after each
    extinction that falls once they are there, and the civilization arises again the same wait after them as before. The first three ages are never changed.
    """
    out = list(ages)
    for age, _share in extinctions:
        animals = out[ANIMALS]
        if animals is None or animals > age:
            continue
        wait = None if out[CIVILIZATION] is None else out[CIVILIZATION] - animals
        back = age + air_events.RECOVERY_GYR + REGROWTH_GYR
        out[ANIMALS] = back
        if wait is not None:
            out[CIVILIZATION] = back + wait
    return tuple(out)
