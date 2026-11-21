# planetgen/physics/life_ages.py

"""
The life ages of a world: the stretches between the mass extinctions of its history (GEN.1746, slice 2; design: reports/air-history-design-gen1739-2026-11-15.md, section 5.6)
==============================================================================================================================================================================

`air_events` draws the shocks of a body's history and the share of species each kills. A **life age** is what lies between two of the kills that are large enough to count as a mass
extinction: it starts when life begins or when the origination that follows an extinction does, and it ends with the next extinction. This module turns the events of a body into that
table, as the geological record divides the history of the Earth into its eras by its mass extinctions. Nothing in the generator reads it yet (the later slices of GEN.1746 store the table,
give the API and the page their rows, and let the extinctions touch the air), so the output of a run is the one it was.

* *What counts.* An event kills at least `AGE_KILL_THRESHOLD` of the species: the 30 percent that the kill curve of the marine record gives once per 10 Myr ([R5320], the size of an extinction of the Big Five
  kind in the proposal of the design); a threshold of 65 percent (once per 100 Myr) gives fewer, longer ages. [Stand-in 6.314 LA1, tier T6, tweak L1: the Boss's pick of 0.30 inside the measured range, 11-19.]
* *No background.* The design proposed drawing a background of kills for whatever the events do not explain, so that a world with no events still has the Earth's turnover. The kill of an event here is read from
  the kill curve at how often the Earth has events as large (`air_events.impact_kill_share`), so the impacts that kill the threshold's share recur exactly as often as the curve says kills of it do: on the
  Earth's record the impacts alone fill the curve, nothing is left over for a background, and the flood basalts are added on top of it. A body with fewer events than the Earth's has fewer extinctions.
* *Recovery.* The origination that follows an extinction lags it by `air_events.RECOVERY_GYR` [R5321], the same whatever its size, so an event that falls within the lag of the one before is part of the
  same extinction (it does not start an age of its own): the extinction keeps the larger kill and its cause. The next age starts one lag after the extinction. [Stand-in 6.314 LA2, tier T7 educated guess
  based on inferred physics, tweak L0.]

Every number that is not read from a source is a design constant with its tier and tweak level in `docs/design/research/gap-filling-plan.md` (GEN.725 rule, section 6.314).
"""

import math
from dataclasses import dataclass

from planetgen.physics import air_events

AGE_KILL_THRESHOLD = 0.30
"""float: The share of species an event must kill to end a life age. [Stand-in 6.314 LA1, tier T6, tweak L1: the Boss's pick of 0.30 inside the measured range, 11-19.]"""


@dataclass(frozen=True)
class LifeAge:
    """
    One stretch of a world's life between two mass extinctions.

    Attributes:
        number (int): Its place in the table, from 1.
        start_gyr (float): When it started, Gyr: the age life began, or the extinction that ended the one before plus the recovery lag.
        end_gyr (float or None): The age of the extinction that ended it, Gyr; `None` for the age the world is in.
        cause (str or None): What ended it: `air_events.FLOOD_BASALT` or `air_events.IMPACT`; `None` for the age the world is in.
        kill_share (float or None): The share of species the extinction killed; `None` for the age the world is in.
    """

    number: int
    start_gyr: float
    end_gyr: float | None
    cause: str | None
    kill_share: float | None


def smallest_impact_km(threshold=AGE_KILL_THRESHOLD):
    """The smallest impactor, km, whose kill reaches `threshold`: the size the Earth has an impact of at the interval the kill curve gives that kill (the size law solved for it); the largest the size law draws if none does."""
    interval_myr = _recurrence_myr(threshold)
    per_gyr = 1000.0 / interval_myr
    diameter = 5.0 * (air_events.IMPACT_FIVE_KM_PER_GYR / per_gyr) ** (1.0 / air_events.IMPACT_SIZE_EXPONENT)
    return min(diameter, air_events.IMPACT_LARGEST_KM)


def _recurrence_myr(threshold):
    """The interval, Myr, at which the kill curve gives a kill of `threshold` (the curve read the other way, straight in the log between its points; the ends of the curve hold their intervals)."""
    curve = air_events.KILL_CURVE
    if threshold <= curve[0][1]:
        return curve[0][0]
    if threshold >= curve[-1][1]:
        return curve[-1][0]
    for (low_myr, low_share), (high_myr, high_share) in zip(curve, curve[1:]):
        if threshold <= high_share:
            weight = (threshold - low_share) / (high_share - low_share)
            return math.exp(math.log(low_myr) + weight * (math.log(high_myr) - math.log(low_myr)))
    return curve[-1][0]


def extinctions_of(events, begin_gyr, end_gyr, threshold=AGE_KILL_THRESHOLD):
    """
    The mass extinctions of a history from `begin_gyr` to `end_gyr`, `(age_gyr, cause, kill_share)` in order: the `events` that kill at least `threshold` and fall inside, the events of one recovery
    lag (`air_events.RECOVERY_GYR`) from the one before merged into it, keeping the larger kill and its cause.
    """
    found = sorted(((event.age_gyr, event.kind, event.kill_share) for event in events if begin_gyr <= event.age_gyr <= end_gyr and event.kill_share >= threshold), key=lambda item: item[0])
    merged = []
    for age, cause, share in found:
        if merged and age - merged[-1][3] < air_events.RECOVERY_GYR:
            first_age, kept_cause, kept_share, _ = merged[-1]
            if share > kept_share:
                kept_cause, kept_share = cause, share
            merged[-1] = (first_age, kept_cause, kept_share, age)
        else:
            merged.append((age, cause, share, age))
    return [(first_age, cause, share) for first_age, cause, share, _ in merged]


def ages_of(body, begin_gyr, end_gyr=None, threshold=AGE_KILL_THRESHOLD):
    """
    The life ages of `body` from the age its life began (`begin_gyr`) to `end_gyr` (the star's age by default): the events of `air_events.events_of` from `begin_gyr` (drawn on the body's own stream), cut into ages at the mass extinctions. Empty for a body whose life never began or whose star has no age. The last age is open (`end_gyr` `None`) unless life began only at the end.
    """
    star_age = getattr(getattr(body, "star", None), "age", None)
    end_gyr = end_gyr if end_gyr is not None else star_age
    if begin_gyr is None or not end_gyr or begin_gyr >= end_gyr:
        return []
    events = air_events.events_of(body, begin_gyr, end_gyr, min_impact_km=smallest_impact_km(threshold))
    ages, start = [], begin_gyr
    for number, (age, cause, share) in enumerate(extinctions_of(events, begin_gyr, end_gyr, threshold), start=1):
        if age < start:
            continue
        ages.append(LifeAge(number, start, age, cause, share))
        start = age + air_events.RECOVERY_GYR
    ages.append(LifeAge(len(ages) + 1, min(start, end_gyr), None, None, None))
    return ages
