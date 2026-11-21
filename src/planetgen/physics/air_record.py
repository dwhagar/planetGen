# planetgen/physics/air_record.py

"""
The record of a world's air history: what is kept of it (GEN.1745 slice 4a; design: reports/air-history-design-gen1739-2026-11-15.md, sections 5.9 and 6.3)
===========================================================================================================================================================

`air_history` follows the clocks of a rocky body (the star's light, the interior's heat, the water lost, the oxygen law) and the temperature of its air at each age, `air_events` draws the shocks
of its history, `life_ages` cuts it into ages at the mass extinctions and `air_history.gas_trajectory` counts the gases the life made. This module puts them together as the record a body keeps:
the recorded steps (`MAX_STEPS` at most) with each one's age, surface temperature, total pressure, the partial pressure of every gas of `PARTIAL_FLOOR_PA` or more and the event that starts it,
and the life ages. `pack` and `unpack` turn it to and from the text the body's column `air_history` holds (JSON, MEDIUMTEXT; slice 4b: `generation/life._history_of` makes it at the end of the life pass from the finished body,
the air as the life found it and the air its oxygen's rise left, which `life._build_air` keeps on the build record). Nothing reads it back yet (the API route and the page are the next slices), so the numbers of a world are the ones they were.

* *Which steps.* The steps of `air_history.schedule` (logarithmic to 1 Gyr, then one every `air_history.MAX_STEP_GYR`), the ages of the oxygen's rise and of the producers of each class, and
  the age of each mass extinction (`life_ages.extinctions_of`); when that is more than `MAX_STEPS`, the steps that are not special (the first, the last, a cut or an extinction) are thinned evenly.
* *The gases.* A world with no life holds its stored air at every age (the air history only moves a gas where life or a shock does; the air effects of the events are a later slice). A world with life holds
  the air `air_history.gas_trajectory` gives at each age (the life as the classes of `living` made it), from the air as the life found it (`found`) and, from the oxygen's rise, the air the rise leaves (`risen`).
* *The temperature.* From `air_history.temperatures` over the inputs of the body's final air solve (`planet.build.air_input`): a full solve at the anchors and the line between them elsewhere. It is the
  temperature of the air as it was first solved, so a world whose life warmed or cooled it has the history's last temperature under its stored one until the life gases reach the path (slice 3 of GEN.1745).
* *Who gets one.* A rocky body (not a giant) with a star that has an age and an air of `MIN_AIR_PA` or more (the Boss's rule of 11-16: under 500 Pa the air is kept as it is built today).
"""

import json
from dataclasses import dataclass

from planetgen.physics import air_history, atmosphere, biogenic_oxygen, life_ages

MAX_STEPS = 60
"""int: The most steps a record keeps: a design size limit of GEN.1745 (the design of GEN.1739 section 5.9 keeps 20 to 60 recorded steps; the size of a packed column goes with it), not a physical constant."""
PARTIAL_FLOOR_PA = 0.1
"""float: The partial pressure, Pa, under which a gas is left out of a step (the design of GEN.1739 section 5.9)."""
MIN_AIR_PA = 500.0
"""float: The surface pressure, Pa, from which a body gets a record: a decision of the Boss (GEN.1745, 5:14 PM EST 2026-11-16: under 500 Pa the air is built as today, with no full history), not a physical constant."""


@dataclass(frozen=True)
class RecordStep:
    """
    One recorded age.

    Attributes:
        age_gyr (float): The age, Gyr.
        temperature_k (float or None): The surface temperature of the air at that age, K; `None` when the history has no clocks for the body.
        pressure_pa (float): The total pressure of the air at that age, Pa.
        partials_pa (dict): Gas formula to partial pressure, Pa, for the gases of `PARTIAL_FLOOR_PA` or more.
        oxygen_law (str or None): `air_history.Step.oxygen_law`.
        event (tuple or None): `(cause, kill_share)` of the mass extinction that starts the step (`cause` an `air_events` kind: `air_events.FLOOD_BASALT` or `air_events.IMPACT`), else `None`.
    """

    age_gyr: float
    temperature_k: float | None
    pressure_pa: float
    partials_pa: dict
    oxygen_law: str | None
    event: tuple | None


@dataclass(frozen=True)
class Record:
    """What a body keeps of its air history: the `steps` (`RecordStep`, ascending age) and the life `ages` (`life_ages.LifeAge`, empty for a world with no life)."""

    steps: tuple
    ages: tuple


def eligible(body):
    """Whether `body` gets a record: a rocky body with a star that has an age and an air of `MIN_AIR_PA` or more."""
    star_age = getattr(getattr(body, "star", None), "age", None)
    return (getattr(body, "body_type", "t") != "g" and bool(star_age) and (getattr(body, "atmospheric_pressure", None) or 0.0) >= MIN_AIR_PA
            and bool(atmosphere.gas_partials_pa(body)))


def record_of(body, living=(), found=None, risen=None, rise_gyr=None, life_begin_gyr=None, ages=None):
    """
    The record of `body` (see the module), or `None` when it gets none (`eligible`). `living` are the classes of life on it with the air `found` as the life found it and the air the oxygen's `risen` at
    `rise_gyr` (as `air_history.gas_trajectory` reads them); `life_begin_gyr` is the age its first producers began, from which the life ages count (none for a world with no life); `ages` are those life ages when the caller has them already (`life_ages.ages_of`, which the record would otherwise make again).
    """
    if not eligible(body):
        return None
    cuts = list(air_history.oxygen_cuts(body, living)) if living else []
    ages_of_life = tuple(ages if ages is not None else life_ages.ages_of(body, life_begin_gyr)) if life_begin_gyr is not None else ()
    kills = {age.end_gyr: (age.cause, age.kill_share) for age in ages_of_life if age.end_gyr is not None}
    cuts.extend(kills)
    steps = air_history.run(body, cuts, getattr(body.build, "air_input", None))
    steps = _thinned(steps, set(cuts))
    ages = [step.age_gyr for step in steps]
    air = atmosphere.gas_partials_pa(body)
    if living and found:
        moles = air_history.gas_trajectory(body, living, ages, found, risen, rise_gyr)
        partials = [biogenic_oxygen.partials_pa(body, state) for _, state in moles]
    else:
        partials = [air] * len(steps)
    recorded = []
    for step, gases in zip(steps, partials):
        kept = {formula: pascal for formula, pascal in gases.items() if pascal >= PARTIAL_FLOOR_PA}
        event = next(((cause, share) for end, (cause, share) in kills.items() if abs(end - step.age_gyr) <= air_history.SAME_AGE_GYR), None)
        recorded.append(RecordStep(step.age_gyr, step.temperature_k, sum(kept.values()), kept, step.oxygen_law, event))
    return Record(tuple(recorded), ages_of_life)


def _thinned(steps, special_ages):
    """`steps` cut to `MAX_STEPS`: the first, the last and the steps at `special_ages` stay, the others are taken out evenly."""
    if len(steps) <= MAX_STEPS:
        return steps
    special = [index for index, step in enumerate(steps) if index in (0, len(steps) - 1) or any(abs(step.age_gyr - age) <= air_history.SAME_AGE_GYR for age in special_ages)]
    others = [index for index in range(len(steps)) if index not in special]
    room = max(MAX_STEPS - len(special), 0)
    kept = {others[int(position * len(others) / room)] for position in range(room)} if room else set()
    return [step for index, step in enumerate(steps) if index in special or index in kept]


def pack(record):
    """The record as the short JSON text a body's column holds: `{"s": [[age, temperature, pressure, event, oxygen law, {gas: pressure}], ...], "a": [[start, end, cause, kill], ...]}`, numbers to five significant digits."""
    def number(value):
        return None if value is None else float(f"{value:.5g}")

    steps = [[number(step.age_gyr), number(step.temperature_k), number(step.pressure_pa), None if step.event is None else [step.event[0], number(step.event[1])],
              step.oxygen_law, {formula: number(pascal) for formula, pascal in step.partials_pa.items()}] for step in record.steps]
    ages = [[number(age.start_gyr), number(age.end_gyr), age.cause, number(age.kill_share)] for age in record.ages]
    return json.dumps({"s": steps, "a": ages}, separators=(",", ":"))


def unpack(text):
    """The `Record` that `pack` made `text` from (to the digits it kept)."""
    data = json.loads(text)
    steps = tuple(RecordStep(age, temperature, pressure, dict(gases), law, None if event is None else tuple(event)) for age, temperature, pressure, event, law, gases in data["s"])
    ages = tuple(life_ages.LifeAge(number, start, end, cause, share) for number, (start, end, cause, share) in enumerate(data["a"], start=1))
    return Record(steps, ages)
