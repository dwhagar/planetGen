# planetgen/physics/life_timeline.py

"""
The time to each stage of life on a world: five hard steps (GEN.750, docs/design/habitability/life-timeline-and-oxygen.md)
==========================================================================================================================

Nothing predicts when life takes each step, so the generator does what the statistics of rare events allow [R4642]: each
stage, in the order of `tuning.LIFE_MILESTONES`, is a step with an exponential waiting time whose mean (`tuning.LIFE_STEP_MEANS_GYR` times the
overall `tuning.LIFE_STEP_SCALE`) is set by hand and is often longer than the world is old. The waiting for a step starts when the one before
it is done, and for the first step at `tuning.LIFE_START_AGE_GYR`, after the surface has settled. So a world has its own schedule, one age for
each stage, and it has reached every stage whose age is not past the age of its star. No Earth history is copied: Earth is a test point the
schedule can land near, not a preset.

Two things change a world's clock:

* **Temperature.** A step goes faster on a warm world and slower on a cold one by the Boltzmann factor of metabolic theory
  (`rate_factor`, activation energy `tuning.LIFE_ACTIVATION_ENERGY_EV`), kept inside `tuning.LIFE_RATE_FACTOR_RANGE`. Its pace
  ("fast", "normal", "slow") is a word for that factor.
* **Oxygen.** Complex cells need free oxygen, so the wait for them starts only when the air has some: at once if it already does, else when
  the world's oxygenic photosynthesis has built it (`biogenic_oxygen`), and never if it cannot (`oxygen_age_gyr` is `None`).

A request for a technological civilization (`INTELLIGENT_LIFE`) raises the odds of a timeline that finishes inside the world's age and forces none
(`raised_waits`; the draw is made on a stream of its own, so the run's stream moves by one draw however many tries it takes).
"""

import math

from planetgen import tuning
from planetgen.physics import constants
from planetgen.util import draw, request_odds

STEPS = len(tuning.LIFE_MILESTONES)
OXYGEN_STEP = tuning.LIFE_MILESTONES.index("complex_cells")
"""int: The step whose wait starts only when the air holds free oxygen."""


def rate_factor(temperature_k):
    """
    How much faster than at the reference temperature a step goes at `temperature_k` (the Boltzmann factor
    `exp(E / k (1 / T_ref - 1 / T))`), inside `tuning.LIFE_RATE_FACTOR_RANGE`; 1 for a world with no temperature.
    """
    if not temperature_k or temperature_k <= 0.0:
        return 1.0
    energy_j = tuning.LIFE_ACTIVATION_ENERGY_EV * constants.ELECTRON_VOLT_J
    exponent = energy_j / constants.BOLTZMANN * (1.0 / tuning.LIFE_REFERENCE_TEMPERATURE_K - 1.0 / temperature_k)
    low, high = tuning.LIFE_RATE_FACTOR_RANGE
    return min(max(math.exp(min(max(exponent, -50.0), 50.0)), low), high)


def pace_of(factor):
    """The word for a rate factor: `"fast"`, `"normal"` or `"slow"` (`tuning.LIFE_PACE_FAST_ABOVE`, `LIFE_PACE_SLOW_BELOW`)."""
    if factor > tuning.LIFE_PACE_FAST_ABOVE:
        return "fast"
    if factor < tuning.LIFE_PACE_SLOW_BELOW:
        return "slow"
    return "normal"


def unit_waits(stream=None):
    """Five waiting times of mean 1 (one exponential draw for each step), from `stream` (a `draw.Stream`) or the run's stream."""
    source = stream if stream is not None else draw
    return tuple(-math.log(1.0 - source.random()) for _ in range(STEPS))


def schedule(waits, factor, oxygen_delay_gyr=0.0, oxygen_at_start=False):
    """
    The age, Gyr, at which a world reaches each stage of `tuning.LIFE_MILESTONES`, from its `unit_waits` and its temperature `factor`
    (`rate_factor`). The wait for complex cells starts when the air holds free oxygen: at once if `oxygen_at_start` (an abiotic oxygen
    air), else `oxygen_delay_gyr` after photosynthesis (`biogenic_oxygen.Profile.delay_gyr`; `None` if the world never makes enough), in
    which case the stages from complex cells on are `None` (never reached). Each age is later than the one before it.
    """
    ages = []
    now = tuning.LIFE_START_AGE_GYR
    for step, (wait, mean) in enumerate(zip(waits, tuning.LIFE_STEP_MEANS_GYR)):
        if step == OXYGEN_STEP and not oxygen_at_start:
            if oxygen_delay_gyr is None:
                ages.extend([None] * (STEPS - step))
                break
            now = max(now, ages[-1] + oxygen_delay_gyr)
        now += wait * mean * tuning.LIFE_STEP_SCALE / factor
        ages.append(now)
    return tuple(ages)


def raised_waits(age_gyr, factor, oxygen_delay_gyr=0.0, oxygen_at_start=False, tries=None):
    """
    The waiting times of a world under a request for a technological civilization (`INTELLIGENT_LIFE`, GEN.1130): timelines are drawn
    from one stream, and one that reaches the last stage inside `age_gyr` is kept, any other one time in `tuning.REQUEST_ODDS`
    (`util.request_odds.keeps`), else another is drawn, up to `tries` (`tuning.REQUEST_MAX_TIMELINE_TRIES`); the last stands, so the
    world may never get there (it never makes the oxygen, or only after the age is gone). The stream is seeded by a single draw of the
    run's, so the run moves by one draw whatever the outcome.
    """
    seed = draw.random()
    stream = draw.Stream(f"life-timeline-request:{seed!r}")
    for _ in range(tries if tries is not None else tuning.REQUEST_MAX_TIMELINE_TRIES):
        waits = unit_waits(stream)
        last = schedule(waits, factor, oxygen_delay_gyr, oxygen_at_start)[-1]
        if request_odds.keeps(last is not None and last <= age_gyr, stream):
            break
    return waits
