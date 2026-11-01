# planetgen/web/lib/day_counter.py

"""
The day-and-night counter of an object page (UX.243): what a planet or moon's own day is doing now, counted at the real rate, apart from the
speed its picture turns at (`spin.js` shows every body turning far faster than it does).

The page reads what is stored (GEN.1097, GEN.1108): the state word (`day_night_state`), the solar day (`solar_day_hours`), the hours of night
in it at the equator (`night_hours`) and, for each of the seven latitude bands, the annual mean temperature, the hours of day and night and the
state. `counter` turns a row into the few facts the page and its script need; the clock itself is `static/daycounter.js`, which counts the
object's days from `EPOCH_MS` so that every object page agrees on the day number it is in.

A body with no stored state (made before GEN.1097, a star, a belt) has no counter: nothing is made up for it.
"""

from planetgen.physics import day_night, latitude_bands, regional_temperature
from planetgen.web.lib.fmt import measured_html

EPOCH_MS = 946728000000
"""int: Where every object's days are counted from, milliseconds since 1970: 2000-01-01 12:00 UTC (J2000). A display choice (the day number is
a count of the object's own solar days since then, the same moment for every object, so two pages agree), not a physical one."""

LOCKED_TEXT = {
    "planet": "Tidally locked: one face stays toward its star, so that side has permanent day and the far side permanent night.",
    "moon": "Tidally locked: one face stays toward its planet, so there is no day and night cycle of its own; what the star does to that face "
            "is on the bands below.",
}
STATE_TEXT = {
    day_night.CYCLE: "Day and night alternate.",
    day_night.POLAR: "Day and night alternate, but a polar zone covers a large part of the surface, so much of it has a polar day or night.",
    day_night.RESONANT: "Locked in a spin-orbit resonance: the star crosses the sky slowly, and a solar day is longer or shorter than the spin.",
    day_night.NON_ROTATING: "Not rotating relative to its star: the star goes round the sky once a year.",
    day_night.NO_SURFACE: "No solid surface: the day below counts the turn of its clouds.",
    day_night.NO_NIGHT: "No night.",
}
"""dict: What the page says for each stored state word (`day_night.STATES`), apart from the locked one, which depends on what the body orbits."""

DISPLAY_NOTE = "The picture of the body is shown turning much faster than the real rotation; the clock below counts at the real rate."
"""str: The sentence that keeps the spin of the picture and the clock apart (spin.js `DISPLAY_NOTE` says the same on the picture)."""

BAND_COLUMNS = (("mean_k", "Annual mean temperature", "K"), ("day_hours", "Day", "h"), ("night_hours", "Night", "h"))
"""tuple: The figures shown beside each band's name: the stored field, its heading and its unit."""


def _positive(value):
    return value is not None and value > 0


def _band_rows(row):
    """The latitude bands (or, for a world that turns once an orbit, its rings) as `{"label", "cells", "state"}`, or `[]` for a body that stores none."""
    from planetgen.describe.sheet import BAND_LABELS
    if row.get(regional_temperature.column(1, "mean_k")) is None:
        return []
    locked = row.get(regional_temperature.column(1, "state")) in latitude_bands.RING_STATES
    names = latitude_bands.RING_NAMES if locked else latitude_bands.NAMES
    rows = []
    for number, name in enumerate(names, start=1):
        cells = []
        for field, _heading, unit in BAND_COLUMNS:
            value = row.get(regional_temperature.column(number, field))
            cells.append(measured_html(unit, value) if value is not None else None)
        state = row.get(regional_temperature.column(number, "state"))
        rows.append({"label": BAND_LABELS[name], "cells": cells, "state": str(state).capitalize() if state else None})
    return rows


def counter(row, kind):
    """
    The counter of the planet or moon whose stored row is `row`, or `None` when it has no stored state.

    Returns:
        dict: `state` (the stored word), `text` (what the page says of it), `ticks` (whether a day is being counted), `day_hours`, `night_hours`
            (the solar day and the night in it, hours, or `None`), `days_per_year`, `epoch_ms`, `note` (`DISPLAY_NOTE`), `bands` (`_band_rows`) and
            `band_headings` (the headings of `BAND_COLUMNS` and the state).
    """
    state = row.get("day_night_state") if row else None
    if state is None:
        return None
    day_hours = row.get("solar_day_hours")
    locked = state == day_night.LOCKED
    ticks = not locked and _positive(day_hours)
    return {
        "state": state,
        "text": LOCKED_TEXT.get(kind, LOCKED_TEXT["planet"]) if locked else STATE_TEXT.get(state, state.capitalize() + "."),
        "ticks": ticks,
        "day_hours": day_hours if ticks else None,
        "night_hours": row.get("night_hours") if ticks else None,
        "days_per_year": row.get("days_per_year"),
        "epoch_ms": EPOCH_MS,
        "note": DISPLAY_NOTE,
        "bands": _band_rows(row),
        "band_headings": [heading for _field, heading, _unit in BAND_COLUMNS] + ["State"],
    }
