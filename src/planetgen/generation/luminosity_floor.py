"""
The luminosity floor presets (GEN.184): the brightness above which every
star is scattered galaxy-wide. Boss (2026-10-10 03:01Z): the user picks from
presets, default 9,000 solar luminosities (3,000, then 5,000, then raised), nothing below 2,500, rising on an
exponential scale to 4 million, with steps of 100 near 2,500 and about
500,000 near the top.

The ladder grows each step geometrically with where it sits on the log
scale (step = 100 * 5000 ** x, x from 0 at the floor to 1 at the ceiling),
and each value is rounded to a round unit for its step so the list reads
cleanly. The default is always one of the presets.
"""
import math

from planetgen import tuning


def _unit(step):
    """The rounding unit for a gap of `step`: 100 up to a gap of 999, then 1,000 ..."""
    return 10.0 ** max(2, math.floor(math.log10(step)))


def _build():
    low, high = tuning.BRIGHT_STAR_FLOOR_MIN_SOL, tuning.BRIGHT_STAR_FLOOR_MAX_SOL
    first, last = tuning.BRIGHT_STAR_FLOOR_MIN_STEP_SOL, tuning.BRIGHT_STAR_FLOOR_MAX_STEP_SOL
    span = math.log(high / low)
    values = [low]
    value = low
    while value < high:
        x = min(1.0, math.log(value / low) / span)
        step = first * (last / first) ** x
        value = max(value + _unit(step), round((value + step) / _unit(step)) * _unit(step))
        values.append(min(value, high))
    values.append(tuning.BRIGHT_STAR_MIN_LUMINOSITY_SOL)
    return tuple(sorted(set(values)))


PRESETS = _build()
"""tuple: Every luminosity floor the Generate page offers, lowest first."""


def nearest(value):
    """The preset closest to `value` (on the log scale)."""
    return min(PRESETS, key=lambda preset: abs(math.log(preset / value)))


def check(value):
    """
    `value` as a float if it is a floor the system accepts: a number from
    `BRIGHT_STAR_FLOOR_MIN_SOL` to `BRIGHT_STAR_FLOOR_MAX_SOL`.

    Raises:
        ValueError: Not a number, or out of range (the message says the range).
    """
    low, high = tuning.BRIGHT_STAR_FLOOR_MIN_SOL, tuning.BRIGHT_STAR_FLOOR_MAX_SOL
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"the luminosity floor must be a number from {low:,.0f} to {high:,.0f} solar luminosities") from None
    if not math.isfinite(number) or not low <= number <= high:
        raise ValueError(f"the luminosity floor must be from {low:,.0f} to {high:,.0f} solar luminosities, got {value}")
    return number
