# planetgen/describe/belt_member_info.py

"""
The page of one body a belt keeps (GEN.1070): its designation, status, mass, size, orbit, makeup and the history of the belt it belongs to
================================================================================================================================

`member_info(conn, belt_id, rank)` reads a stored belt and answers the structure of the member of that rank: the entry the belt's information page lists for it (`physics/belt_page.stored_bodies`,
so the two pages never differ), its orbit (the elements as they were saved, with the closest and farthest distance and the period the star's mass gives), its share of the belt's mass, the
belt it is a member of (its designation, the star it goes round, where it came from and what shaped it) and the designations of the heavier and the lighter neighbour. A rank above the number of
bodies the belt keeps individually names a body that is drawn on demand and never stored: it has no page here yet. `GET /api/objects/<ref>/belt-members/<rank>` serves it and the member's page shows it.
Nothing here is stored.
"""

import math

from planetgen.db import store
from planetgen.describe import belt_info
from planetgen.galaxy import object_uid
from planetgen.generation import belt_ondemand, irregular_moons
from planetgen.physics import planet_class_text
from planetgen.physics import belt_page as page_data
from planetgen.physics import unit_registry as units


class BeltMemberNotFound(ValueError):
    """No belt has that id, or the belt keeps no body of that rank."""


def orbit(member, star_mass_msun):
    """The orbit of a stored member: its six elements and the closest distance, the farthest distance and the period they and the star's mass give (Kepler's third law)."""
    a, e = member["distance_au"], member["eccentricity"]
    q = page_data._q
    return {
        "semi_major_axis": q(a, "au"),
        "eccentricity": e,
        "inclination": q(member["inclination_deg"], "deg"),
        "ascending_node": q(member["node_deg"], "deg"),
        "argument_of_periapsis": q(member["periapsis_deg"], "deg"),
        "mean_anomaly": q(member["anomaly_deg"], "deg"),
        "closest_distance": q(a * (1.0 - e), "au"),
        "farthest_distance": q(a * (1.0 + e), "au"),
        "period": q(math.sqrt(a ** 3 / star_mass_msun), "yr"),
    }


def _belt_block(belt, star, belt_id, class_d_count):
    """What the member page says of its belt: its designation, the star it goes round and the sentences of where it came from and what shaped it."""
    return {
        "designation": object_uid.format_id(belt_id),
        "star": star.name,
        "origin": getattr(belt, "origin_sentence", None),
        "history": list(getattr(belt, "origin_history", ()) or ()),
        "history_text": getattr(belt, "history", None) or None,
        "individual_count": len(belt.members),
        "class_d_count": class_d_count,
    }


def planet_block(system, star, belt, rank):
    """
    The planet the body of `rank` is when asked for (GEN.1674): `belt_ondemand.planet_of_rank` makes it with the one planet generator, and this is its class and a few of what the pass gave it. `None` for a
    belt whose body the generator cannot make a planet of (a body under `FLOOR_KG` is not asked for).
    """
    order = next((place for place, held in irregular_moons.belts_of(system.planets).items() if held is belt), 1)
    planet = belt_ondemand.planet_of_rank(system, star, belt, order, rank)
    q = page_data._q
    return {
        "class": planet.planet_class,
        "class_text": planet_class_text.description(planet.planet_class, False),
        "mass": q(planet.mass, "kg"),
        "radius": q(planet.radius, "km"),
        "density": q(planet.density * 1000.0, "kg_m3"),
        "gravity": q(planet.gravity, "g0"),
        "surface_temperature": q(planet.surface_temperature, "K"),
        "surface_pressure": q(planet.atmospheric_pressure or 0.0, "Pa"),
        "zone": planet.zone,
    }


def member_info(conn, belt_id, rank):
    """
    The page data of the member of `rank` (1 is the heaviest) of the stored belt `belt_id`: `member` (the belt page's entry), `orbit`, `share` of the belt's mass, `belt`, the `previous` and `next`
    designations (`None` at the ends) and `planet` (the body made into a planet on demand, `planet_block`). A rank past the stored ones is drawn from the belt's size law (`belt_ondemand`).

    Raises:
        BeltMemberNotFound: No belt has that id, or it holds no body of that rank at the Class D size or larger.
    """
    row = conn.execute("SELECT uid, star_system_id, star_id FROM asteroid_belts WHERE id = ?", (belt_id,)).fetchone()
    if row is None:
        raise BeltMemberNotFound(f"no belt has the id {belt_id}")
    system = store.load_star_system(conn, row["star_system_id"])
    belt = belt_info.belt_of(system, belt_id)
    star = belt_info.host_star(system, row["star_id"])
    members = list(belt.members)
    found = belt_ondemand.member_of_rank(belt, star, system.system_config, rank)
    if found is None:
        raise BeltMemberNotFound(f"the belt holds no body of rank {rank} at or over the Class D size (it keeps {len(members)} one by one and {belt_ondemand.individual_count(belt_ondemand.settled_of(belt, star)):,} in all)")
    object_id = object_uid.from_bytes(row["uid"])
    count = belt_ondemand.individual_count(belt_ondemand.settled_of(belt, star))
    ranks = [member["rank"] for member in members]
    if rank in ranks:
        place = ranks.index(rank)
        neighbour = lambda i: object_uid.designation(object_id, ranks[i]) if 0 <= i < len(ranks) else None
        previous, following = neighbour(place - 1), neighbour(place + 1)
    else:
        previous = object_uid.designation(object_id, rank - 1) if rank > 1 else None
        following = object_uid.designation(object_id, rank + 1) if rank < count else None
    shown = page_data.stored_bodies(object_id, [found])[0]
    return {
        "designation": shown["designation"],
        "rank": rank,
        "member": shown,
        "orbit": orbit(found, units.convert(star.mass, "kg", "Msun")),
        "share": found["mass_kg"] / float(belt.mass_kg),
        "belt": _belt_block(belt, star, object_id, count),
        "previous": previous,
        "next": following,
        "planet": planet_block(system, star, belt, rank),
    }
