# planetgen/physics/body_history.py

"""
The numbers a body's history leaves that cannot be worked out again from what is stored (GEN.1497, the "kind 4" of GEN.1465)
=========================================================================================================================

The gravity effects of the effect ledger are recomputed on demand from stored masses and orbits (Boss 8:37 PM EST 11-08). A few numbers cannot be: they were
drawn, or they are the state a run started from, and the stored body keeps only where the run ended. This module names those numbers and reads them off a
generated body when the system is saved. It draws nothing and changes no result: the draws are kept on the body by the generator (`planets.draw_birth_spin`
keeps `build.birth_spin`, `planets.generate_orbit_shape` keeps `build.orbit_birth`) for the passes that run the history again, and this table keeps them for the stored galaxy.

| Column | Unit | From |
|---|---|---|
| `birth_eccentricity` | fraction | the eccentricity a planet or captured moon was born on, before the tides and the secular chain settled it (`build.orbit_birth`) |
| `birth_spin_rate` | rad/s | the spin rate it was born with (`build.birth_spin[0]`) |
| `birth_obliquity_deg` | deg | the tilt it was born with, before the tides, the Cassini state or a Venus-type spin changed it (`build.birth_spin[1]`) |
| `capture_draw_1`, `capture_draw_2` | fraction | one uniform for each resonance the spin can pass without being caught (`spin_torques.CAPTURE_RATIOS`, in that order; `build.birth_spin[2]`): they decide the tide history of the spin |
| `formation_distance_au` | AU | where the body formed (`Planet.formation_distance_au`; a moon's is its planet's): the key of the end-stage draw streams (`draw.body_stream_of`) and the line the water's origin is read against, so a body that moves keeps drawing what it drew |
| `birth_distance_au` | AU | a moon's orbit when it was born, before the planet's tide moved it (`_birth_distance_au`); the moon history of the planet starts from it |
| `born_index` | count | the place of a moon in its planet's list of born moons (`_born_moons`), which the lost moons share; the place of a planet in its star's plan (`Planet.plan_slot`), which no trim of the planets moves. Both name the body's slot, so its ID follows it (DB.75b) |
| `held_eccentricity` | fraction | the eccentricity a resonance holds the body on (`build.resonance_e`: a planet of a resonant chain, oddity W5; a moon of a first-order pair or caught in the evection resonance), the floor `settle_orbit` keeps under the tide (DB.132) |
| `forced_eccentricity` | fraction | the eccentricity the planets of a resonant chain force on a planet (`oddity_quantities["forced_eccentricity"]`, `resonant_chains.give_forced_eccentricities`), the quantity the oddity of the chain reads (DB.132) |
| `formula_version` | count | `FORMULA_VERSION`: which set of formulas wrote the row, so a later change to them can tell the rows it must not refill from the rows it can |

**Read back (GEN.1533, GEN.1753, DB.132):** `attach`, called by `store.load_star_system`, gives a restored planet or moon its `build.orbit_birth`, `build.birth_spin`, `formation_distance_au`, (a moon) `_birth_distance_au`, the eccentricity a resonance holds it on (`build.resonance_e`) and a planet's forced one (`oddity_quantities`) from its row, so
the stored galaxy carries what a run would (the air history of GEN.1739 starts from them, and a moved body's history starts from its own birth). The passes that run on a moved or edited body read them: `planets.generate_spin`
takes the `build.birth_spin` it is given (`validation.reconcile_moved_planet`), the draw streams are keyed by the formation distance, and the moon history (`moon_history.apply`) starts from the born moons. A planet's born moons are
rebuilt by `store.attach_body_history`: the moons it still has, each at its birth distance, and the ones the history lost, kept in the `lost_moons` table (`lost_rows`: mass, birth distance, density, direction) and given back as
`moon_history.LostMoon` stand-ins that take part in the history and are never a body of the system. The capture draws are a `FLOAT` column, so a restored draw is the stored one to six digits.

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
    ("formation_distance_au", "au", "where the body formed, which keys its end-stage draw streams"),
    ("birth_distance_au", "au", "the orbit a moon was born on, before the planet's tide moved it"),
    ("born_index", None, "the place of a moon in its planet's list of born moons, or of a planet in its star's plan"),
    ("held_eccentricity", None, "the eccentricity a resonance holds the body on, the floor its tide cannot wear below"),
    ("forced_eccentricity", None, "the eccentricity the planets of a resonant chain force on a planet"),
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
        body: A generated planet or moon (any object: one with no `build.orbit_birth` or `build.birth_spin` has none of those).
    """
    out = {}
    record = getattr(body, "build", None)
    if record is None:      # a belt builds nothing and carries no record
        return out
    birth_orbit = record.orbit_birth
    if birth_orbit is not None:
        out["birth_eccentricity"] = float(birth_orbit)
    birth_spin = record.birth_spin
    if birth_spin is not None:
        rate, tilt, captures = birth_spin
        out["birth_spin_rate"] = float(rate)
        out["birth_obliquity_deg"] = float(tilt)
        for index, unit in enumerate(captures[:2], start=1):
            out[f"capture_draw_{index}"] = float(unit)
    formed = getattr(body, "formation_distance_au", None)
    if formed is not None:
        out["formation_distance_au"] = float(formed)
    born_at = getattr(body, "_birth_distance_au", None)
    if born_at is not None:
        out["birth_distance_au"] = float(born_at)
    held = record.resonance_e
    if held:
        out["held_eccentricity"] = float(held)
    forced = (getattr(body, "oddity_quantities", None) or {}).get("forced_eccentricity")
    if forced:
        out["forced_eccentricity"] = float(forced)
    if out:
        out["formula_version"] = FORMULA_VERSION
    return out


def attach(body, row):
    """
    Gives a planet or moon restored from the database the numbers of its history that `body_history` holds (GEN.1533, the reader of the table): the birth eccentricity as `build.orbit_birth` and the
    birth spin, tilt and capture draws as `build.birth_spin`, the two fields a generated body keeps them in (`BuildRecord`), so a pass that runs a body's history again reads the body it was made as and not a fresh draw.
    The spin is given whole or not at all (a body that has any of it has all of it). A row written under another `FORMULA_VERSION` is not read: its numbers belong to formulas this code no longer runs.
    Nothing changes for a body with no row (a galaxy made before v243, or a body with no history).

    Args:
        body: A restored planet or moon.
        row (dict): The body's `body_history` row (`ROW_FIELDS`).

    Returns:
        bool: Whether any number was attached.
    """
    if row["formula_version"] != FORMULA_VERSION:
        return False
    attached = False
    if row["formation_distance_au"] is not None:      # GEN.1753: where it formed keys the draw streams, so a body that moves or is edited keeps drawing what it drew
        body.formation_distance_au = float(row["formation_distance_au"])
        attached = True
    if row["birth_distance_au"] is not None:
        body._birth_distance_au = float(row["birth_distance_au"])
        attached = True
    if row["birth_eccentricity"] is not None:
        body.build.orbit_birth = float(row["birth_eccentricity"])
        attached = True
    if row["held_eccentricity"] is not None:      # DB.132: the floor a resonance keeps under the tide, so a body settled again is held on what it was held on
        body.build.resonance_e = float(row["held_eccentricity"])
        attached = True
    if row["forced_eccentricity"] is not None:
        body.oddity_quantities = {**(getattr(body, "oddity_quantities", None) or {}), "forced_eccentricity": float(row["forced_eccentricity"])}
        attached = True
    draws = (row["capture_draw_1"], row["capture_draw_2"])
    if row["birth_spin_rate"] is not None and row["birth_obliquity_deg"] is not None and all(draw is not None for draw in draws):
        body.build.birth_spin = (float(row["birth_spin_rate"]), float(row["birth_obliquity_deg"]), tuple(float(draw) for draw in draws))
        attached = True
    return attached


def rows(star_system_id, bodies, ids, born=None):
    """
    The `body_history` rows of a saved system.

    Args:
        star_system_id (int): The system's `star_systems.id`.
        bodies (iterable): The system's generated planets and moons.
        ids (dict): `id(object)` to `(table, row id)` of the saved bodies (`object_rows.logged_rows`); a body with no entry was not saved and is skipped.
        born (dict, optional): `id(body)` to its `born_index`: a moon's place in its planet's born moons (`born_places`), a planet's slot in its star's plan (`Planet.plan_slot`); a body with none has no `born_index`.

    Returns:
        list: One dict of `ROW_FIELDS` for each body that has any number; a number the body lacks is `None`.
    """
    out = []
    born = born or {}
    for body in bodies:
        where = ids.get(id(body))
        numbers = numbers_of(body) if where is not None else None
        if not numbers:
            continue
        if id(body) in born:
            numbers["born_index"] = born[id(body)]
        out.append({"star_system_id": star_system_id, "body_type": effect_ledger.BODY_CODES[where[0]], "body_id": where[1],
                    **{name: numbers.get(name) for name in COLUMN_NAMES}})
    return out


LOST_FIELDS = ("star_system_id", "planet_id", "born_index", "mass_kg", "birth_distance_au", "density_kg_m3", "retrograde")
"""tuple: The columns of a `lost_moons` row, in table order."""


def born_places(planet):
    """`{id(moon): place}` of the moons `planet` was born with (`moon_history.born_moons`) that it still has, the place being their index in that list; empty for a planet with no born list."""
    born = planet.build.born_moons or ()
    return {id(moon): place for place, moon in enumerate(born) if any(moon is kept for kept in getattr(planet, "moons", ()))}


def lost_rows(star_system_id, planet, planet_id):
    """
    The `lost_moons` rows of one saved planet (GEN.1753): the moons it was born with that its history took away (fell below the Roche limit, were pushed out of the stable orbits or dragged in), each as the
    four numbers `moon_orbits.Moon` reads (its mass, birth orbit, density and direction) and its place in the born list, so that the history runs again as it ran.

    Args:
        star_system_id (int): The system's `star_systems.id`.
        planet: A generated planet.
        planet_id (int): Its `planets.id`.

    Returns:
        list: One dict of `LOST_FIELDS` for each born moon the planet no longer has; empty when it has no born list.
    """
    from planetgen.physics import moon_history
    kept = born_places(planet)
    out = []
    for place, moon in enumerate(planet.build.born_moons or ()):
        if id(moon) in kept:
            continue
        out.append({"star_system_id": star_system_id, "planet_id": planet_id, "born_index": place, "mass_kg": float(moon.mass),
                    "birth_distance_au": float(getattr(moon, "_birth_distance_au", moon.distance)), "density_kg_m3": moon_history.density_kg_m3(moon),
                    "retrograde": int(moon_history.retrograde(moon))})
    return out
