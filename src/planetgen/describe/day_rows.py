# planetgen/describe/day_rows.py

"""
The derived day-and-night row of the object sheet for a star and a brown dwarf with no star (GEN.1153)
=======================================================================================================

GEN.1097 stores the state of the day and night of every planet and moon (`physics/day_night`). A star and a free brown dwarf store nothing, because the answer is the same word every time: a
star is `no night` (it shines, it does not stand under one) and a brown dwarf with no star is a giant, `no solid surface`, with no solar day (`daylight.compute` finds no sun for it, GEN.1496).
They are worked out when the sheet is made and are not columns, the derive-on-demand route (no migration): a pure function of what the object is. A rogue planet the planet chain built keeps the
state it stores; only a row with none (a brown dwarf) gets the derived one. Each row is `(group, key, label, value, raw)`.
"""

from planetgen.physics import day_night

STAR_KEY = "star_day_night_state"
"""str: The key of a star's derived row."""

FREE_KEY = "free_day_night_state"
"""str: The key of a free brown dwarf's derived row."""

LABEL = "State of its day and night"

KEYS = (STAR_KEY, FREE_KEY)


def is_derived_key(key):
    """Whether `key` is one of the rows this module gives."""
    return key in KEYS


def rows(conn, kind, object_id):
    """The derived day-and-night row of a star, or of a rogue planet that stores none (a brown dwarf): `[(group, key, label, value, raw)]`, empty for any other object."""
    if kind == "star":
        return [("rotation", STAR_KEY, LABEL, [day_night.state(None, None, star=True)], None)]
    if kind == "rogue_planet":
        row = conn.execute("SELECT day_night_state, planet_type FROM rogue_planets WHERE id = ?", (object_id,)).fetchone()
        if row is not None and row["day_night_state"] is None:
            return [("rotation", FREE_KEY, LABEL, [day_night.state(None, None, solid=False, no_sun=True)], None)]
    return []
