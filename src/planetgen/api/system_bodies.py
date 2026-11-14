# planetgen/api/system_bodies.py

"""
The body ledger of one system and the barycenter that includes its belts (GEN.1241, the wiring of GEN.1066; `docs/design/systems/system-bodies-ledger.md`).

`GET /api/systems/<id>/bodies` builds the rows of `physics/system_bodies` as a view over the typed tables, with no new storage: the stars, planets, moons and comets of the system
scene (`web/maps/systemscene.build_scene`, which places them from the stored rows), the belts with the mass they store, the bodies each belt keeps (`belt_members`), the remainder of the
belt below the stored bodies (one extended mass at the belt's center) and the rings of the planets (`planets.ring_mass_kg`). It answers the total mass, the barycenter, where the central
star sits against it, and how far the belts move that star (`belt_wobble_au`), with the mass each belt's rows add up to against the mass the belt stores.

The frame is the primary star's center at the origin (the star the scene calls the origin of a wide binary, else the first star of the system). A body that goes round the center of
mass of a close pair is placed from that center. Positions are those the typed tables hold: a planet's, moon's and comet's at the epoch the orbits were last advanced to, a belt body's
at the mean anomaly it was drawn with (its orbit is not advanced), so the belts' wobble is the size of the effect, not a place on a date.

    payload = build(conn, system_id)       # None for a system that does not exist
"""

import math

from planetgen.physics import system_bodies as ledger
from planetgen.physics import unit_registry as units
from planetgen.util import detmath
from planetgen.web.maps.systemscene import build_scene

_km_to_au = units.converter("km", "au")
_au_to_km = units.converter("au", "km")

PAIR_ROLES = ("primary", "secondary")
"""tuple: The roles of the stars of a close pair, the center of mass of which its own bodies go round."""


def _position_au(km):
    return tuple(_km_to_au(value) for value in km)


def _shifted(position, offset):
    """`position` moved by `offset`, both AU triples."""
    return tuple(position[axis] + offset[axis] for axis in range(3))


def _au(km):
    """A distance in km as AU, or `None` for a body that has none stored."""
    return None if km is None else _km_to_au(km)


def _density(mass_kg, radius_km):
    if not mass_kg or not radius_km or radius_km <= 0:
        return None
    return mass_kg / (4.0 / 3.0 * math.pi * (radius_km * 1000.0) ** 3)


class _Anchors:
    """Where a body that goes round `around` (a ref of the scene) is anchored: the id of the parent row and the place of that center from the parent's, in AU."""

    def __init__(self, stars):
        self.root = stars[0]
        self.absolute = {s["ref"]: _position_au(s["position_km"]) for s in stars}
        pair = [s for s in stars if s.get("role") in PAIR_ROLES]
        mass = math.fsum(s["mass_kg"] or 0.0 for s in pair)
        self.pair_center = tuple(math.fsum((s["mass_kg"] or 0.0) * self.absolute[s["ref"]][axis] for s in pair) / mass for axis in range(3)) if mass > 0 else None
        self.origin = self.absolute[self.root["ref"]]

    def of(self, around):
        """`(parent row id, offset of the center from the parent, AU)`."""
        if around in self.absolute:
            return around, (0.0, 0.0, 0.0)
        center = self.pair_center if self.pair_center is not None else (0.0, 0.0, 0.0)
        return self.root["ref"], tuple(center[axis] - self.origin[axis] for axis in range(3))


def _rows(conn, system_id, scene):
    """The ledger rows and the stored mass of each belt, from the scene and the typed tables."""
    stars = scene["stars"]
    if not stars:
        return [], {}, {}
    root = stars[0]
    if scene["tree"]["origin"]:
        root = next(s for s in stars if s["ref"] == scene["tree"]["origin"])
        stars = [root] + [s for s in stars if s is not root]
    anchors = _Anchors(stars)
    rows, names = [], {}

    def add(identifier, parent, kind, name, mass_kg, position, radius_km=None, a_au=None, e=None, i_deg=None, status=None, tier="typed", density=None):
        rows.append(ledger.Body(identifier, parent, kind, mass_kg, radius_km, density if density is not None else _density(mass_kg, radius_km), a_au, e, i_deg, status, tier,
                                position))
        names[identifier] = name

    for star in stars:
        if star is root:
            add(star["ref"], None, "star", star["name"], star["mass_kg"] or 0.0, None, star["radius_km"], status="star")
        else:
            offset = tuple(anchors.absolute[star["ref"]][axis] - anchors.origin[axis] for axis in range(3))
            orbit = star.get("orbit") or {}
            add(star["ref"], root["ref"], "star", star["name"], star["mass_kg"] or 0.0, offset, star["radius_km"], a_au=_au(orbit.get("distance_km")),
                i_deg=orbit.get("inclination_deg"), status="star")

    rings = {r["id"]: r["ring_mass_kg"] for r in conn.execute("SELECT id, ring_mass_kg FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()}
    for planet in scene["planets"]:
        parent, offset = anchors.of(planet["orbit"]["around"])
        position = _shifted(_position_au(planet["position_km"]), offset)
        orbit = planet["orbit"]
        add(planet["ref"], parent, "planet", planet["name"], planet["mass_kg"] or 0.0, position, planet["radius_km"], _au(orbit["distance_km"]), orbit["eccentricity"],
            orbit["inclination_deg"], status="planet")
        ring_mass = rings.get(planet["id"])
        if ring_mass and ring_mass > 0:
            add(f"ring:{planet['id']}", planet["ref"], "ring", f"Rings of {planet['name']}", ring_mass, None, status="ring")
        for moon in planet["moons"]:
            orbit = moon["orbit"]
            add(moon["ref"], planet["ref"], "moon", moon["name"], moon["mass_kg"] or 0.0, _position_au(moon["position_km"]), moon["radius_km"], _au(orbit["distance_km"]),
                orbit["eccentricity"], orbit["inclination_deg"], status="moon")

    belt_mass = {}
    belts = {r["id"]: dict(r) for r in conn.execute("SELECT id, mass_kg FROM asteroid_belts WHERE star_system_id = ?", (system_id,)).fetchall()}
    for belt in scene["belts"]:
        parent, offset = anchors.of(belt["around"])
        add(belt["ref"], parent, "belt", belt["name"], 0.0, offset, status="belt")
        stored = (belts.get(belt["id"]) or {}).get("mass_kg")
        belt_mass[belt["ref"]] = stored if stored and stored > 0 else None
        kept = 0.0
        members = conn.execute(
            "SELECT member_rank, mass_kg, radius_km, distance_km, eccentricity, inclination_deg, ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, density_g_cm3, status, taxon "
            "FROM belt_members WHERE belt_id = ? ORDER BY member_rank", (belt["id"],)).fetchall()
        for member in members:
            a_au = _km_to_au(member["distance_km"])
            position = ledger.relative_position_au(a_au, min(max(member["eccentricity"], 0.0), 0.999999), member["inclination_deg"], member["ascending_node_deg"],
                                                   member["arg_periapsis_deg"], member["mean_anomaly_deg"])
            status = member["status"] if member["status"] in ledger.STATUSES else "minor body"
            add(f"belt_member:{belt['id']}:{member['member_rank']}", belt["ref"], "belt_member", f"{belt['name']} {member['member_rank']:04X}", member["mass_kg"], position,
                member["radius_km"], a_au, member["eccentricity"], member["inclination_deg"], status, "L1", member["density_g_cm3"] * 1000.0)
            kept += member["mass_kg"]
        total = belt_mass[belt["ref"]]
        if total is not None and total > kept:
            add(f"extended:{belt['id']}", belt["ref"], "extended", f"{belt['name']} (the rest)", total - kept, None, status="extended", tier="L3")

    comets = {r["id"]: r["mass_kg"] for r in conn.execute("SELECT id, mass_kg FROM comets WHERE star_system_id = ?", (system_id,)).fetchall()}
    for comet in scene["comets"]:
        parent, offset = anchors.of(comet["orbit"]["around"])
        add(comet["ref"], parent, "comet", comet["name"], comets.get(comet["id"]) or 0.0, _shifted(_position_au(comet["position_km"]), offset), comet["radius_km"],
            status="comet")
    return rows, belt_mass, names


def _entry(row, names):
    """One row as the answer prints it: a typed body carries its `kind:id` ref (the API prints the object ID), a body of a belt its rank or its belt."""
    kind = row.kind
    entry = {"kind": kind, "name": names[row.id], "mass_kg": row.mass_kg, "radius_km": row.radius_km, "density_kg_m3": row.density_kg_m3, "a_au": row.a_au, "e": row.e,
             "i_deg": row.i_deg, "status": row.status, "tier": row.tier, "parent": row.parent,
             "position_au": list(row.position_au) if row.position_au is not None else None}
    if kind in ("belt_member", "extended"):
        parts = row.id.split(":")
        entry["belt"] = f"belt:{parts[1]}"
        entry["rank"] = f"{int(parts[2]):04X}" if kind == "belt_member" else None
        entry["ref"] = None
        entry["parent"] = entry["belt"]
    elif kind == "ring":
        entry["ref"] = None
    else:
        entry["ref"] = row.id
        entry["id"] = int(row.id.split(":")[1])
    return entry


def build(conn, system_id):
    """
    The ledger payload of one system, or `None` for a system that does not exist: `system_id`, `bodies` (each `{ref, kind, name, mass_kg, radius_km, density_kg_m3, a_au, e, i_deg, status, tier,
    parent, position_au, distance_from_barycenter_au}`; a body of a belt has `belt` and `rank` in place of an `id`), `total_mass_kg`, `barycenter_au`, `central_offset_au`, `central_offset_km`,
    `belt_wobble_au`, `belt_wobble_km`, `central_radius_km` and `belts` (`{belt, name, mass_kg, members, stored_members_kg, extended_kg, residual}`; `residual` is `None` for a belt that stores no mass).
    """
    try:
        scene = build_scene(conn, system_id)
    except ValueError:
        return None
    rows, belt_mass, names = _rows(conn, system_id, scene)
    payload = {"system_id": system_id, "bodies": [], "belts": [], "total_mass_kg": 0.0, "barycenter_au": None, "central_offset_au": None, "central_offset_km": None,
               "belt_wobble_au": 0.0, "belt_wobble_km": 0.0, "central_radius_km": None, "frame": "the center of the primary star at the origin"}
    if not rows:
        return payload
    order = ledger.check(rows)
    by_id = {row.id: row for row in rows}
    places = ledger.absolute_positions_au(rows)
    barycenter = ledger.barycenter_au(rows)
    offset = ledger.central_offset_au(rows)
    for identifier in order:
        entry = _entry(by_id[identifier], names)
        entry["distance_from_barycenter_au"] = detmath.dist(places[identifier], barycenter)
        payload["bodies"].append(entry)
    wobble = ledger.belt_wobble_au(rows)
    for belt in scene["belts"]:
        members = [r for r in rows if r.parent == belt["ref"] and r.kind == "belt_member"]
        extended = [r for r in rows if r.parent == belt["ref"] and r.kind == "extended"]
        stored = belt_mass.get(belt["ref"])
        payload["belts"].append({
            "belt": belt["ref"], "name": belt["name"], "mass_kg": stored, "members": len(members), "stored_members_kg": math.fsum(r.mass_kg for r in members),
            "extended_kg": math.fsum(r.mass_kg for r in extended),
            "residual": ledger.belt_mass_residual(rows, belt["ref"], stored) if stored else None})
    payload.update(total_mass_kg=ledger.total_mass_kg(rows), barycenter_au=list(barycenter), central_offset_au=list(offset), central_offset_km=_au_to_km(detmath.hypot(*offset)),
                   belt_wobble_au=wobble, belt_wobble_km=_au_to_km(wobble), central_radius_km=by_id[order[0]].radius_km)
    return payload
