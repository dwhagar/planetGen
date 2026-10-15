# planetgen/db/species_fit.py

"""
Which species could live where (POP.12).

`limits_for_species` is the one function that says what a species' environment limits are. Until DB.30 adds the
columns that store them, it generates them from each species' homeworld (`species_limits.generated`); when the
columns exist, a stored value replaces the generated one there and nowhere else. Everything below reads the
limits through it.

`species_for_body` lists the species that could live on a planet or moon, and `bodies_for_species` the worlds of a
system or sector that one species could live on. Both score the body against the species' limits with
`species_limits.assess`, the same six domains as the human index, and keep two lists: able to live there unaided
(no Red domain) and able to live there with the equipment its tech level can make.
"""

from planetgen.db import store
from planetgen.physics import habitability_world, species_limits
from planetgen.physics.protection_limits import LIMIT_ROWS

CHUNK = 500
"""int: Homeworlds read per query."""

DOMAINS = habitability_world.DOMAINS
"""tuple: The six domains as the API names them (lower case)."""

SCOPE_LIMIT = 2000
"""int: The most bodies `bodies_for_species` scores in one call (a sector's, or a system's)."""


def stored_overrides(conn, species_ids):
    """
    The limits stored for `species_ids`, `{species id: overrides}` for those that have them. Nothing is stored yet:
    DB.30 adds the columns, and this is the one place that then reads them.
    """
    return {}


def limits_for_species(conn, species_ids=None):
    """
    The environment limits of species, `{species id: overrides}` (`species_limits`: the rows that differ from the
    human table; `{}` for a species with the human limits). `species_ids` None is every species. A stored limit wins;
    otherwise they are generated from the species' homeworld, and a species whose homeworld has no habitability score
    has the human limits.
    """
    if species_ids is None:
        wanted = [row["id"] for row in conn.execute("SELECT id FROM species ORDER BY id").fetchall()]
    else:
        wanted = list(species_ids)
    out = dict(stored_overrides(conn, wanted))
    missing = [species_id for species_id in wanted if species_id not in out]
    for start in range(0, len(missing), CHUNK):
        chunk = missing[start:start + CHUNK]
        marks = ", ".join("?" for _ in chunk)
        for row in conn.execute(
            f"SELECT s.id AS species_id, p.* FROM species s JOIN planets p ON p.id = s.homeworld_planet_id "
            f"WHERE s.id IN ({marks})", chunk
        ).fetchall():
            world = store.stored_world(conn, "planets", row, with_life=False)
            out[row["species_id"]] = species_limits.generated(world) if world is not None else {}
    for species_id in missing:
        out.setdefault(species_id, {})
    return out


def limits_table(overrides):
    """
    A species' limits as rows for the API and the species page: `name`, `unit`, `meaning`, `human` (the human
    value), `value` (the species' own) and `changed`.
    """
    return [{"name": row.name, "unit": row.unit, "meaning": row.meaning, "human": row.value,
             "value": overrides.get(row.name, row.value), "changed": row.name in overrides} for row in LIMIT_ROWS]


class PlacementError(ValueError):
    """A species put where it could not live (POP.14). `str(error)` says why in plain words; `failing` lists the
    domains that are Red, lower case."""

    def __init__(self, message, failing=()):
        super().__init__(message)
        self.failing = list(failing)


def check_placement(conn, species_id, kind, body_id, habitat=False, overrides=None, native=False):
    """
    POP.14: the one check every way of putting a species on a planet or moon calls before it writes. Passes when the
    species could live on the body unaided (no Red habitability domain against its own limits); with `habitat` true
    (the placement is a habitat or a vessel, which brings its own kit) it also passes when the species could live there
    with the equipment its tech level can make. Otherwise raises `PlacementError` saying which conditions fail.

    A body with no habitability score (saved before v76) cannot be judged and passes; a gas giant has no surface to
    live on and fails. `overrides` are the species' limits if the caller already has them (a species being created).
    `native` is true for the species' own homeworld, where its Suitability is not judged (`species_limits.assess`).

    Returns:
        dict: The verdict fields of `_verdict`, or `{}` for an unscored body.
    """
    species = conn.execute("SELECT id, name, tech_level FROM species WHERE id = ?", (species_id,)).fetchone()
    if species is None:
        raise PlacementError(f"There is no species {species_id}.")
    table, row = _body_row(conn, kind, body_id)
    if row is None:
        raise PlacementError(f"There is no {kind} {body_id}.")
    if row["body_type"] == "g":
        raise PlacementError(f"The {species['name']} could not live on {row['name']}: it is a gas giant with no surface.")
    world = store.stored_world(conn, table, row)
    if world is None:
        return {}
    if overrides is None:
        overrides = limits_for_species(conn, [species_id])[species_id]
    tech = None if species["tech_level"] is None else float(species["tech_level"])
    assessment = species_limits.assess(world, overrides, tech, native)
    if assessment["unaided"] or (habitat and assessment["equipped"]):
        return _verdict(assessment)
    raise _refusal(f"The {species['name']}", row["name"], world, overrides, assessment)


def _refusal(who, body_name, world, overrides, assessment):
    reasons = species_limits.explain(world, overrides, assessment["failing"]) or ["it has no drinkable water"]
    return PlacementError(f"{who} could not live on {body_name}: " + "; ".join(reasons) + ".",
                          [name.lower() for name in assessment["failing"]])


def check_homeworld(conn, planet_id):
    """
    POP.14: the check a species gets before it is made on planet `planet_id` as its homeworld. Its limits are the ones
    its homeworld would give it (`species_limits.generated`), so this fails only for what a homeworld cannot fix in a
    species' limits (undrinkable water). A planet with no score passes. Raises `PlacementError`.
    """
    row = conn.execute("SELECT * FROM planets WHERE id = ?", (planet_id,)).fetchone()
    world = None if row is None else store.stored_world(conn, "planets", row)
    if world is None:
        return
    overrides = species_limits.generated(world)
    assessment = species_limits.assess(world, overrides, None, native=True)
    if not assessment["unaided"]:
        raise _refusal("A species", row["name"], world, overrides, assessment)


def inconsistencies(conn):
    """
    POP.14: every species placed outside its limits, `[{species_id, name, kind, id, body_name, star_system_id,
    reasons, failing}]`. The only placement a species has today is its homeworld; a later settlement (POP.11) adds
    its own and is listed the same way. A homeworld with no score is not listed.
    """
    out = []
    species = conn.execute(
        "SELECT s.id, s.name, s.homeworld_planet_id, s.star_system_id FROM species s ORDER BY s.id").fetchall()
    limits = limits_for_species(conn, [s["id"] for s in species])
    for s in species:
        try:
            check_placement(conn, s["id"], "planet", s["homeworld_planet_id"], overrides=limits[s["id"]], native=True)
        except PlacementError as error:
            body = conn.execute("SELECT name FROM planets WHERE id = ?", (s["homeworld_planet_id"],)).fetchone()
            out.append({"species_id": s["id"], "name": s["name"], "kind": "planet", "id": s["homeworld_planet_id"],
                        "body_name": body["name"] if body else None, "star_system_id": s["star_system_id"],
                        "reasons": str(error), "failing": error.failing})
    return out


def _verdict(assessment):
    """One body's `species_limits.assess` as the API's fields."""
    return {
        "unaided": assessment["unaided"], "equipped": assessment["equipped"], "equipment": assessment["equipment"],
        "failing": [name.lower() for name in assessment["failing"]],
        "phi4": assessment["domains"]["PHI-4"],
        "domains": {name: {"score": assessment["domains"][name.capitalize()][0],
                           "tier": assessment["domains"][name.capitalize()][1]} for name in DOMAINS},
    }


def _body_row(conn, kind, body_id):
    table = "moons" if kind == "moon" else "planets"
    return table, conn.execute(f"SELECT * FROM {table} WHERE id = ?", (body_id,)).fetchone()


def _all_species(conn):
    """Every species with its limits: `[(row, overrides)]`."""
    species = conn.execute("SELECT id, name, homeworld_planet_id, tech_level FROM species ORDER BY id").fetchall()
    limits = limits_for_species(conn, [s["id"] for s in species])
    return [(s, limits[s["id"]]) for s in species]


def _fits(world, species, planet_id=None):
    """`(unaided, equipped)` lists of the `species` (`_all_species`) that could live on `world`, best PHI-4 first.
    `planet_id` is the row id of the world if it is a planet: its own species are judged as natives there."""
    unaided, equipped = [], []
    for s, overrides in species:
        tech = None if s["tech_level"] is None else float(s["tech_level"])
        verdict = _verdict(species_limits.assess(world, overrides, tech, s["homeworld_planet_id"] == planet_id))
        entry = {"id": s["id"], "name": s["name"], "homeworld_planet_id": s["homeworld_planet_id"],
                 "tech_level": tech, **verdict}
        if verdict["unaided"]:
            unaided.append(entry)
        elif verdict["equipped"]:
            equipped.append(entry)
    for found in (unaided, equipped):
        found.sort(key=lambda e: (-e["phi4"], e["name"]))
    return unaided, equipped


def species_for_body(conn, kind, body_id):
    """
    The species that could live on a planet or moon (`kind` "planet" or "moon", `body_id` its row id), or `None`
    for a body that does not exist.

    Returns:
        dict: `scored` (false for a gas giant or a body saved without a score: both lists are then empty), `unaided`
            (species with no Red domain there, best PHI-4 first) and `equipped` (those that could with the kit their
            tech level can make and are not already in `unaided`). Each entry has `id`, `name`, `tech_level`, the
            verdict fields of `_verdict` and `homeworld_planet_id`.
    """
    table, row = _body_row(conn, kind, body_id)
    if row is None:
        return None
    world = store.stored_world(conn, table, row)
    out = {"kind": kind, "id": body_id, "scored": world is not None, "unaided": [], "equipped": []}
    if world is not None:
        out["unaided"], out["equipped"] = _fits(world, _all_species(conn), body_id if kind == "planet" else None)
    return out


def species_for_system(conn, system_id):
    """
    The species that could live on each rocky planet and moon of one system, or `None` for a system that does not
    exist: `{system_id, bodies}`, `bodies` the scored ones that have a species that could live there, each as
    `{kind, id, name, unaided, equipped}` (the lists of `species_for_body`). The limits are worked out once.
    """
    if conn.execute("SELECT 1 FROM star_systems WHERE id = ?", (system_id,)).fetchone() is None:
        return None
    species = _all_species(conn)
    bodies = []
    for table, kind in (("planets", "planet"), ("moons", "moon")):
        for row in conn.execute(
            f"SELECT * FROM {table} WHERE star_system_id = ? AND body_type = 't' AND energy_flux_w_m2 IS NOT NULL "
            f"AND surface_temperature_k IS NOT NULL ORDER BY id", (system_id,)
        ).fetchall():
            world = store.stored_world(conn, table, row)
            if world is None:
                continue
            unaided, equipped = _fits(world, species, row["id"] if kind == "planet" else None)
            if unaided or equipped:
                bodies.append({"kind": kind, "id": row["id"], "name": row["name"], "unaided": unaided,
                               "equipped": equipped})
    return {"system_id": system_id, "bodies": bodies}


def bodies_for_species(conn, species_id, system_id=None, sector_id=None):
    """
    The rocky planets and moons of one system (`system_id`) or of every system of one sector (`sector_id`, its row
    id) that species `species_id` could live on, or `None` for a species that does not exist.

    Returns:
        dict: `species_id`, `scored` (bodies scored), `truncated` (more than `SCOPE_LIMIT` bodies were in scope; the
            first were scored), `unaided` and `equipped` lists of `{kind, id, name, star_system_id}` and the verdict
            fields of `_verdict`, best PHI-4 first.
    """
    species = conn.execute("SELECT id, tech_level, homeworld_planet_id FROM species WHERE id = ?",
                           (species_id,)).fetchone()
    if species is None:
        return None
    home = species["homeworld_planet_id"]
    overrides = limits_for_species(conn, [species_id])[species_id]
    tech = None if species["tech_level"] is None else float(species["tech_level"])
    where, params = ("star_system_id = ?", (system_id,)) if system_id is not None else (
        "star_system_id IN (SELECT id FROM star_systems WHERE sector_id = ?)", (sector_id,))
    out = {"species_id": species_id, "scored": 0, "truncated": False, "unaided": [], "equipped": []}
    budget = SCOPE_LIMIT
    for table, kind in (("planets", "planet"), ("moons", "moon")):
        rows = conn.execute(
            f"SELECT * FROM {table} WHERE {where} AND body_type = 't' AND energy_flux_w_m2 IS NOT NULL "
            f"AND surface_temperature_k IS NOT NULL ORDER BY id LIMIT ?", (*params, budget + 1)).fetchall()
        if len(rows) > budget:
            out["truncated"] = True
            rows = rows[:budget]
        budget -= len(rows)
        for row in rows:
            world = store.stored_world(conn, table, row)
            if world is None:
                continue
            out["scored"] += 1
            verdict = _verdict(species_limits.assess(world, overrides, tech, kind == "planet" and row["id"] == home))
            entry = {"kind": kind, "id": row["id"], "name": row["name"], "star_system_id": row["star_system_id"],
                     **verdict}
            if verdict["unaided"]:
                out["unaided"].append(entry)
            elif verdict["equipped"]:
                out["equipped"].append(entry)
    for key in ("unaided", "equipped"):
        out[key].sort(key=lambda e: (-e["phi4"], e["name"] or "", e["id"]))
    return out
