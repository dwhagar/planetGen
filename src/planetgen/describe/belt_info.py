# planetgen/describe/belt_info.py

"""
The information page of an asteroid belt (UX.241): the bodies it keeps one by one, the rest of its mass by size, and how crowded it is
================================================================================================================================

`belt_info(conn, belt_id)` reads a stored belt and its star, rebuilds the ledger the belt's mass and its stored members make (the members are the ledger's stored bodies, the rest is the belt's
mass less theirs, under the size law of a belt of that place), and answers the structure `physics/belt_page.build` makes: the members by designation (as they were saved, not drawn again), the
size table of everything not stored, the density factor and the totals. `GET /api/objects/<ref>/belt-page` serves it and the belt's information page shows it. Nothing here is stored.
"""

from planetgen.db import store
from planetgen.galaxy import object_uid
from planetgen.physics import belt_ledger, radiogenic
from planetgen.physics import unit_registry as units
from planetgen.physics import belt_page as page_data


class BeltPageNotFound(ValueError):
    """No belt has that id."""


def host_star(system, star_id):
    """The star a belt of the system goes round: the one with that stored id, else the system's primary."""
    stars = [system.star, getattr(system, "secondary_star", None), *getattr(system, "member_stars", [])]
    for star in stars:
        if star is not None and getattr(star, "db_id", None) == star_id:
            return star
    return system.star


def belt_of(system, belt_id):
    bodies = [*system.planets, *getattr(system, "secondary_planets", []), *(body for group in (getattr(system, "member_planets", None) or []) for body in group)]
    for body in bodies:
        if getattr(body, "body_type", None) == "a" and getattr(body, "db_id", None) == belt_id:
            return body
    raise BeltPageNotFound(f"no belt has the id {belt_id}")


def _al26_ratio(system_config, star):
    """The aluminium-26 the star was born with; a star with no [Fe/H] or age takes its galaxy's mean (as the generator does for the same star)."""
    shape = getattr(system_config, "GALAXY_SHAPE", None)
    ratio = radiogenic.star_al26_ratio(star, shape)
    return radiogenic.al26_mean(shape, getattr(star, "metallicity_fe_h", None)) if ratio is None else ratio


def belt_page_record(belt, star, object_id):
    """The `belt_page.BeltRecord` of a stored `belt` going round `star`: its members' masses are the ledger's stored bodies and its remainder is the belt's mass less theirs."""
    masses = tuple(float(member["mass_kg"]) for member in belt.members)
    belt_mass = float(belt.mass_kg)
    law = belt_ledger.ICY_LAW if belt.distance >= star.snow_line_au else belt_ledger.ROCKY_LAW
    ledger = belt_ledger.Ledger(masses, belt_mass, max(belt_mass - sum(masses), 0.0), len(masses), False, law)
    return page_data.BeltRecord(object_id, belt.lower_limit, belt.upper_limit, units.convert(star.mass, "kg", "Msun"), star.luminosity,
                                _al26_ratio(belt.system_config, star), ledger, ())


def belt_info(conn, belt_id):
    """
    The page data of the stored belt `belt_id` (`physics/belt_page.build`) with the name of the star it goes round under `belt.star`.

    Raises:
        BeltPageNotFound: No belt has that id.
    """
    row = conn.execute("SELECT uid, star_system_id, star_id FROM asteroid_belts WHERE id = ?", (belt_id,)).fetchone()
    if row is None:
        raise BeltPageNotFound(f"no belt has the id {belt_id}")
    system = store.load_star_system(conn, row["star_system_id"])
    belt = belt_of(system, belt_id)
    star = host_star(system, row["star_id"])
    object_id = object_uid.from_bytes(row["uid"])
    page = page_data.build(belt_page_record(belt, star, object_id), page_data.stored_bodies(object_id, belt.members))
    page["belt"]["star"] = star.name
    return page
