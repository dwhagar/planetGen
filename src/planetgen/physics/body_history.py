# planetgen/physics/body_history.py

"""
The numbers a body's history leaves that cannot be worked out again from what is stored (GEN.1497, the "kind 4" of GEN.1465)
=========================================================================================================================

The gravity effects of the effect ledger are recomputed on demand from stored masses and orbits (Boss 8:37 PM EST 11-08). A few numbers cannot be: they were
drawn, or they are the state a run started from, and the stored body keeps only where the run ended. This module names those numbers and reads them off a
generated body when the system is saved. It draws nothing and changes no result: the draws are kept on the body by the generator (`planets.draw_birth_spin`
keeps `_birth_spin`, `planets.generate_orbit_shape` keeps `_orbit_birth`) for the passes that run the history again, and this table keeps them for the stored galaxy.

| Column | Unit | From |
|---|---|---|
| `birth_eccentricity` | fraction | the eccentricity a planet or captured moon was born on, before the tides and the secular chain settled it (`_orbit_birth`) |
| `birth_spin_rate` | rad/s | the spin rate it was born with (`_birth_spin[0]`) |
| `birth_obliquity_deg` | deg | the tilt it was born with, before the tides, the Cassini state or a Venus-type spin changed it (`_birth_spin[1]`) |
| `capture_draw_1`, `capture_draw_2` | fraction | one uniform for each resonance the spin can pass without being caught (`spin_torques.CAPTURE_RATIOS`, in that order; `_birth_spin[2]`): they decide the tide history of the spin |
| `formula_version` | count | `FORMULA_VERSION`: which set of formulas wrote the row, so a later change to them can tell the rows it must not refill from the rows it can |

The Kozai cycle count and start age (GEN.1447) and the age a tidal lock was reached (GEN.1454) are the numbers still to come: each lane that builds one adds its
column here (a migration, a unit in `db/column_units.py`, a reader below) and a test checks the columns. Nothing writes them yet.

**One row a body** (`body_history`: system, body kind, body, then the columns), only for a planet or moon that has any; a column is NULL when the body does not have that
number (a captured moon has a birth orbit, a moon born with its planet does not). The census of `reports/gen1465-ledger-store-vs-compute-2026-11-09.md` section 4
budgeted about 63 bytes a body; a row is about 70 bytes.
"""

from planetgen.physics import effect_ledger, spin_torques

FORMULA_VERSION = 1
"""int: The version of the formulas the numbers were drawn under; raise it when a formula these numbers feed changes in a way that a stored row would no longer match."""

COLUMNS = (
    ("birth_eccentricity", None, "the eccentricity the body was born on, before its orbit was settled"),
    ("birth_spin_rate", "rad_s", "the spin rate the body was born with"),
    ("birth_obliquity_deg", "deg", "the tilt the body was born with, before its history moved it"),
    ("capture_draw_1", None, "the uniform that decides whether the spin is caught in the first resonance (`spin_torques.CAPTURE_RATIOS[0]`)"),
    ("capture_draw_2", None, "the uniform that decides whether the spin is caught in the second resonance (`spin_torques.CAPTURE_RATIOS[1]`)"),
    ("formula_version", None, "the version of the formulas that wrote the row"),
)
"""tuple: `(column, unit key of the unit registry or None for a plain number, what it is)` for each number after the three key columns of `body_history`."""

COLUMN_NAMES = tuple(name for name, *_rest in COLUMNS)
"""tuple: The numbers' column names, in table order."""

ROW_FIELDS = ("star_system_id", "body_type", "body_id", *COLUMN_NAMES)
"""tuple: The columns of a `body_history` row, in table order."""

assert len(spin_torques.CAPTURE_RATIOS) == 2, "capture_draw_1 and capture_draw_2 name the two resonances of spin_torques.CAPTURE_RATIOS"


def numbers_of(body):
    """
    The numbers of `body`'s history that it keeps, as `{column: value}`; empty for a body with none.

    Args:
        body: A generated planet or moon (any object: one with no `_orbit_birth` or `_birth_spin` has none of those).
    """
    out = {}
    birth_orbit = getattr(body, "_orbit_birth", None)
    if birth_orbit is not None:
        out["birth_eccentricity"] = float(birth_orbit)
    birth_spin = getattr(body, "_birth_spin", None)
    if birth_spin is not None:
        rate, tilt, captures = birth_spin
        out["birth_spin_rate"] = float(rate)
        out["birth_obliquity_deg"] = float(tilt)
        for index, unit in enumerate(captures[:2], start=1):
            out[f"capture_draw_{index}"] = float(unit)
    if out:
        out["formula_version"] = FORMULA_VERSION
    return out


def rows(star_system_id, bodies, ids):
    """
    The `body_history` rows of a saved system.

    Args:
        star_system_id (int): The system's `star_systems.id`.
        bodies (iterable): The system's generated planets and moons.
        ids (dict): `id(object)` to `(table, row id)` of the saved bodies (`object_rows.logged_rows`); a body with no entry was not saved and is skipped.

    Returns:
        list: One dict of `ROW_FIELDS` for each body that has any number; a number the body lacks is `None`.
    """
    out = []
    for body in bodies:
        where = ids.get(id(body))
        numbers = numbers_of(body) if where is not None else None
        if not numbers:
            continue
        out.append({"star_system_id": star_system_id, "body_type": effect_ledger.BODY_CODES[where[0]], "body_id": where[1],
                    **{name: numbers.get(name) for name in COLUMN_NAMES}})
    return out
