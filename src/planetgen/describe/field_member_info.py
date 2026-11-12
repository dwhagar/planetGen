# planetgen/describe/field_member_info.py

"""
The bodies an interstellar asteroid field keeps one by one, and the page of one of them (GEN.1105 part 2)
=========================================================================================================

`members_of(conn, field_id)` reads the `field_members` rows of a stored field and answers the list its page shows (heaviest first, each under `<field ID>.<rank hex>`, with its mass, diameter, density,
makeup and status word); `member_info(conn, field_id, rank)` answers one of them with its share of the field's mass, the field it belongs to (its name, where it came from and its history) and the
designations of the heavier and the lighter neighbour. A field has no star, so a body has no orbit of its own here: where it is at any age follows from the field's seed (`field_dispersal.member_position_m`).
`GET /api/objects/<ref>/field-members` and `GET /api/objects/<ref>/field-members/<rank>` serve them and the pages show them. Nothing here is stored.
"""

import math

from planetgen.galaxy import object_uid
from planetgen.generation.phenomena import asteroid_field
from planetgen.physics import belt_page as page_data, constants, field_dispersal


class FieldMemberNotFound(ValueError):
    """No field has that id, or the field keeps no body of that rank."""


def _field_row(conn, field_id):
    row = conn.execute("SELECT id, uid, name, mass_kg, origin, history, shape, rms_speed_kms, birth_age_yr FROM asteroid_fields WHERE id = ?", (field_id,)).fetchone()
    if row is None or row["uid"] is None:
        raise FieldMemberNotFound(f"no field has the id {field_id}")
    return row


def _entry(object_id, row):
    """One body as the page lists it: the same quantities as a belt body, less the orbit a free body has none of."""
    q = page_data._q
    return {
        "rank": row["member_rank"],
        "designation": object_uid.field_designation(object_id, row["member_rank"]),
        "mass": q(row["mass_kg"], "kg"),
        "diameter": q(2.0 * row["radius_km"], "km"),
        "density": q(row["density_g_cm3"] * 1000.0, "kg_m3"),
        "composition": {"rock": q(row["rock_fraction"], "frac"), "metal": q(row["metal_fraction"], "frac"),
                        "water": q(row["water_fraction"], "frac"), "ice": q(row["ice_fraction"], "frac")},
        "status": {"word": row["status"]},
    }


def _rows(conn, field_id):
    return conn.execute("SELECT member_rank, mass_kg, radius_km, density_g_cm3, rock_fraction, metal_fraction, water_fraction, ice_fraction, status"
                        " FROM field_members WHERE field_id = ? ORDER BY member_rank", (field_id,)).fetchall()


def _position(field, object_id, rank):
    """Where the body of `rank` is now, light-years from the center of its field, and how fast it moves, from the field's designation as the seed (`field_dispersal.member_position_m`, GEN.1559): none for a field stored without its speed or age."""
    if not field["rms_speed_kms"] or field["birth_age_yr"] is None:
        return None
    annulus = field["shape"] == "annulus"
    seed = object_uid.format_id(object_id)
    speed_ms = field["rms_speed_kms"] * 1000.0
    shape = "annulus" if annulus else "ball"
    ring_m = asteroid_field.RING_RADIUS_AU * constants.AU_M if annulus else 0.0
    x, y, z = field_dispersal.member_position_m(seed, rank, field["birth_age_yr"], speed_ms, shape, ring_m)
    vx, vy, vz = field_dispersal.member_velocity_ms(seed, rank, speed_ms, shape)
    return {"from_center": page_data._q(math.sqrt(x * x + y * y + z * z) / constants.LIGHTYEAR_M, "ly"),
            "speed": page_data._q(math.sqrt(vx * vx + vy * vy + vz * vz) / 1000.0, "km_s")}


def members_of(conn, field_id):
    """
    The bodies the stored field `field_id` keeps: `field` (its designation, name and mass), `count` and `bodies`, heaviest first (empty for a field with no ledger).

    Raises:
        FieldMemberNotFound: No field has that id.
    """
    field = _field_row(conn, field_id)
    object_id = object_uid.from_bytes(field["uid"])
    rows = _rows(conn, field_id)
    return {
        "field": {"designation": object_uid.format_id(object_id), "name": field["name"], "mass": page_data._q(field["mass_kg"], "kg") if field["mass_kg"] else None},
        "count": len(rows),
        "bodies": [_entry(object_id, row) for row in rows],
    }


def member_info(conn, field_id, rank):
    """
    The page data of the body of `rank` (1 is the heaviest) of the stored field `field_id`: `designation`, `rank`, `member` (the entry the field's list holds), `share` of the field's mass, `field`, its `position` (light-years from the field's center and its speed, none for a field stored without them) and the
    `previous` and `next` designations (`None` at the ends).

    Raises:
        FieldMemberNotFound: No field has that id, or it keeps no stored body of that rank.
    """
    field = _field_row(conn, field_id)
    object_id = object_uid.from_bytes(field["uid"])
    rows = _rows(conn, field_id)
    ranks = [row["member_rank"] for row in rows]
    if rank not in ranks:
        raise FieldMemberNotFound(f"the field keeps {len(rows)} bodies one by one and none of rank {rank}")
    place = ranks.index(rank)
    neighbour = lambda i: object_uid.field_designation(object_id, ranks[i]) if 0 <= i < len(ranks) else None
    shown = _entry(object_id, rows[place])
    mass = field["mass_kg"]
    return {
        "designation": shown["designation"],
        "rank": rank,
        "member": shown,
        "share": rows[place]["mass_kg"] / float(mass) if mass else None,
        "field": {"designation": object_uid.format_id(object_id), "name": field["name"], "origin": field["origin"], "history": field["history"] or None, "individual_count": len(rows)},
        "position": _position(field, object_id, rank),
        "previous": neighbour(place - 1),
        "next": neighbour(place + 1),
    }
