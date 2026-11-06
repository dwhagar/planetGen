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


def _belt_block(belt, star, belt_id):
    """What the member page says of its belt: its designation, the star it goes round and the sentences of where it came from and what shaped it."""
    return {
        "designation": object_uid.format_id(belt_id),
        "star": star.name,
        "origin": getattr(belt, "origin_sentence", None),
        "history": list(getattr(belt, "origin_history", ()) or ()),
        "history_text": getattr(belt, "history", None) or None,
        "individual_count": len(belt.members),
    }


def member_info(conn, belt_id, rank):
    """
    The page data of the member of `rank` (1 is the heaviest) of the stored belt `belt_id`: `member` (the belt page's entry), `orbit`, `share` of the belt's mass, `belt` and the `previous` and `next`
    designations (`None` at the ends).

    Raises:
        BeltMemberNotFound: No belt has that id, or it keeps no stored body of that rank.
    """
    row = conn.execute("SELECT uid, star_system_id, star_id FROM asteroid_belts WHERE id = ?", (belt_id,)).fetchone()
    if row is None:
        raise BeltMemberNotFound(f"no belt has the id {belt_id}")
    system = store.load_star_system(conn, row["star_system_id"])
    belt = belt_info.belt_of(system, belt_id)
    star = belt_info.host_star(system, row["star_id"])
    members = list(belt.members)
    found = next((member for member in members if member["rank"] == rank), None)
    if found is None:
        raise BeltMemberNotFound(f"the belt keeps {len(members)} bodies one by one; the body of rank {rank} is drawn when asked for and has no page")
    object_id = object_uid.from_bytes(row["uid"])
    ranks = [member["rank"] for member in members]
    place = ranks.index(rank)
    neighbour = lambda i: object_uid.designation(object_id, ranks[i]) if 0 <= i < len(ranks) else None
    shown = page_data.stored_bodies(object_id, [found])[0]
    return {
        "designation": shown["designation"],
        "rank": rank,
        "member": shown,
        "orbit": orbit(found, units.convert(star.mass, "kg", "Msun")),
        "share": found["mass_kg"] / float(belt.mass_kg),
        "belt": _belt_block(belt, star, object_id),
        "previous": neighbour(place - 1),
        "next": neighbour(place + 1),
    }
