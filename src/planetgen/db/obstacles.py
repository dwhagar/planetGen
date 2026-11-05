# planetgen/db/obstacles.py

"""
The obstacles along a path (NAV.25)
===================================

Course planning (NAV.26) bends a straight course around the keep-out spheres it crosses, so it first needs every object whose
sphere reaches a straight segment of the galaxy. `obstacles_along_segment` answers that for a segment in the galaxy frame (between
systems and inside a sector; a course inside a system works from the system's bodies at a time, NAV.27):

    every placed object whose keep-out sphere (NAV.54, `db.query.keep_out_radius`) comes within the sphere's own radius of the segment,
    that is, whose sphere the segment touches, in order along the segment.

The radii differ by orders of magnitude (a rogue planet's is under 0.1 pc, a quasar's can be 200 pc), so the search is made one kind at
a time with that kind's own widest radius (`widest_reach_pc`, from aggregates of the kind's table). Stars are not listed apart: a system's
sphere is its stars' widest, the one the star-level rule gives, so a system stands for them. Nebulae, remnants, belts and comets have
no sphere and are not obstacles (their notes belong to the course's `passed_through` list, NAV.28).

Only placed, generated objects count; nothing is generated here. Parsecs, galaxy frame.
"""

import math

from planetgen.db import corridor
from planetgen.galaxy import avoid, keepout
from planetgen.physics import constants
from planetgen.physics.units import ly_to_pc, pc_to_ly

KM_PER_PC = constants.PARSEC_M / 1000.0
PC_CM = constants.PARSEC_M * 100.0

CONTAIN_MARGIN = avoid.CONTAIN_MARGIN
"""float: A sphere contains a point when the center is no farther than `radius * (1 + CONTAIN_MARGIN)` (the closed test of the design,
course-avoidance.md section 2.1)."""

OBSTACLE_KINDS = ("system", "black_hole", "neutron_star", "rogue_planet", "quasar", "asteroid_field")
"""tuple: The kinds that have a keep-out sphere between systems."""

_CHUNK = 500


def _scalar(conn, sql):
    row = conn.execute(sql).fetchone()
    return float(row["v"]) if row is not None and row["v"] is not None else 0.0


def widest_reach_pc(conn, kind):
    """
    The widest keep-out radius, parsecs, any placed object of `kind` can have, from aggregates of its table (a bound, not an exact
    maximum: the heaviest object is taken at the farthest distance from the center, and the galactic Hill radius only grows outward).
    `0.0` when the kind has no placed object.
    """
    from planetgen.db import query  # not at the top: query imports this module's neighbor, corridor

    def hill_pc(mass_kg, radius_pc, size_km):
        if mass_kg <= 0.0 or radius_pc <= 0.0:
            return size_km / KM_PER_PC
        hill = keepout.galactic_hill_radius_km(mass_kg, radius_pc, query._enclosed_mass_kg(conn, radius_pc))
        return max(hill, size_km) / KM_PER_PC

    sun = constants.SOLAR_MASS_TO_KG
    if kind == "system":
        stars = _scalar(conn, "SELECT MAX(system_perimeter_km) AS v FROM stars")
        binary = _scalar(conn, "SELECT MAX(binary_system_perimeter_km) AS v FROM star_systems")
        return max(stars, binary) / KM_PER_PC
    if kind == "asteroid_field":
        return ly_to_pc(_scalar(conn, "SELECT MAX(radius_ly) AS v FROM asteroid_fields WHERE center_x_pc IS NOT NULL"))
    if kind == "quasar":
        if not conn.execute("SELECT 1 FROM quasars WHERE center_x_pc IS NOT NULL LIMIT 1").fetchone():
            return 0.0
        found = keepout.quasar_keep_out(_scalar(conn, "SELECT MAX(black_hole_mass_solar) AS v FROM quasars WHERE center_x_pc IS NOT NULL"),
                                        _scalar(conn, "SELECT MAX(luminosity_w) AS v FROM quasars WHERE center_x_pc IS NOT NULL"))
        return found.radius_km / KM_PER_PC
    if kind == "black_hole":
        reach = 0.0
        row = conn.execute("SELECT MAX(mass_solar) AS m FROM black_holes WHERE center_x_pc IS NOT NULL AND mass_class = 'supermassive'"
                           ).fetchone()
        if row["m"]:
            reach = keepout.sphere_of_influence_km(float(row["m"])) / KM_PER_PC
        row = conn.execute("SELECT MAX(mass_solar) AS m, MAX(galactic_radius_pc) AS r, MAX(event_horizon_radius_km) AS s FROM black_holes "
                           "WHERE star_id IS NULL AND center_x_pc IS NOT NULL AND mass_class <> 'supermassive'").fetchone()
        if row["m"]:
            reach = max(reach, hill_pc(float(row["m"]) * sun, float(row["r"] or 0.0), float(row["s"] or 0.0)))
        return reach
    if kind == "neutron_star":
        row = conn.execute("SELECT MAX(mass_solar) AS m, MAX(galactic_radius_pc) AS r, MAX(radius_km) AS s FROM neutron_stars "
                           "WHERE star_id IS NULL AND center_x_pc IS NOT NULL").fetchone()
        return hill_pc(float(row["m"]) * sun, float(row["r"] or 0.0), float(row["s"] or 0.0)) if row["m"] else 0.0
    if kind == "rogue_planet":
        row = conn.execute("SELECT MAX(mass_kg) AS m, MAX(galactic_radius_pc) AS r, MAX(radius_km) AS s FROM rogue_planets "
                           "WHERE center_x_pc IS NOT NULL").fetchone()
        return hill_pc(float(row["m"]), float(row["r"] or 0.0), float(row["s"] or 0.0)) if row["m"] else 0.0
    raise ValueError(f"{kind!r} has no keep-out sphere between systems")


def _system_radii_km(conn, system_ids):
    """`{system id: keep-out radius km}`: the widest stored perimeter of the system's stars and of the binary."""
    found = {}
    ids = list(system_ids)
    for offset in range(0, len(ids), _CHUNK):
        chunk = ids[offset:offset + _CHUNK]
        marks = ", ".join("?" * len(chunk))
        for row in conn.execute(
                "SELECT s.id, GREATEST(COALESCE(s.binary_system_perimeter_km, 0), COALESCE(MAX(st.system_perimeter_km), 0)) AS km "
                f"FROM star_systems s LEFT JOIN stars st ON st.star_system_id = s.id WHERE s.id IN ({marks}) GROUP BY s.id",
                tuple(chunk)).fetchall():
            found[row["id"]] = float(row["km"] or 0.0)
    return found


def obstacles_along_segment(conn, start_pc, end_pc, kinds=None):
    """
    The keep-out spheres a straight segment touches, in order along it.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        start_pc, end_pc (tuple): The ends of the segment, galaxy frame, parsecs.
        kinds (iterable, optional): Any of `OBSTACLE_KINDS`; all by default.

    Returns:
        list[dict]: One per sphere: `ref`, `kind`, `id`, `name`, `center_pc`, `radius_pc`, `radius_km`, `basis` and `note` (the keep-out
            of NAV.54), `distance_pc` (center to the segment), `along_pc` (where along it the nearest point lies) and
            `contains_start` / `contains_end` (an end lies inside the sphere, so a course that starts or ends there is exempt from it).
            A sphere is listed when `distance_pc <= radius_pc`.

    Raises:
        ValueError: For a kind that has no sphere.
    """
    from planetgen.db import query

    wanted = list(OBSTACLE_KINDS if kinds is None else kinds)
    unknown = [kind for kind in wanted if kind not in OBSTACLE_KINDS]
    if unknown:
        raise ValueError(f"no keep-out sphere for: {', '.join(sorted(unknown))}")
    start, end = tuple(float(v) for v in start_pc), tuple(float(v) for v in end_pc)
    found = []
    for kind in wanted:
        reach = widest_reach_pc(conn, kind)
        if reach <= 0.0:
            continue
        rows = corridor.objects_near_segment(conn, start, end, reach, kinds=[kind])
        radii_km = _system_radii_km(conn, [row["id"] for row in rows]) if kind == "system" else {}
        for row in rows:
            if kind == "system":
                sphere = keepout.KeepOut(radii_km.get(row["id"]) or None, "perimeter", None)
            else:
                sphere = query.keep_out_radius(conn, kind, row["id"])
            if not sphere.radius_km or sphere.radius_km <= 0.0:
                continue
            radius_pc = sphere.radius_km / KM_PER_PC
            if row["distance_pc"] > radius_pc:
                continue
            center = tuple(row["point_pc"])
            limit = radius_pc * (1.0 + CONTAIN_MARGIN)
            found.append({
                "ref": row["ref"], "kind": kind, "id": row["id"], "name": row["name"], "center_pc": center,
                "radius_pc": radius_pc, "radius_km": sphere.radius_km, "basis": sphere.basis, "note": sphere.note,
                "distance_pc": row["distance_pc"], "along_pc": row["along_pc"],
                "contains_start": math.dist(center, start) <= limit, "contains_end": math.dist(center, end) <= limit,
            })
    found.sort(key=lambda row: (row["along_pc"], row["distance_pc"], row["kind"], str(row["id"])))
    return found


PASS_THROUGH_KINDS = ("nebula", "supernova_remnant")
"""tuple: The kinds a course passes through with a note and no keep-out (NAV.54, course-avoidance.md section 3.4)."""

_PASS_THROUGH_TABLES = {"nebula": "nebulae", "supernova_remnant": "supernova_remnants"}


def passed_through_along_segment(conn, start_pc, end_pc):
    """
    NAV.28: the clouds and remnants a straight segment crosses, in order along it. A nebula or a supernova remnant is no obstacle (the gas
    is thin and the drawn shape fills a fraction of its bounding sphere), so the course goes through it and says so.

    The length inside is the chord of the bounding sphere clipped to the segment, an upper bound: the drawn shape fills about 13% of the
    sphere. A nebula's column density is its stored density times that length, also an upper bound.

    Returns:
        list[dict]: `ref`, `kind`, `id`, `name`, `note` (`keepout.PASS_THROUGH_NOTE`), `inside_length` and `inside_length_ly` (the length
            inside, parsecs and light-years; never more than the segment) and, for a nebula with a density, `column_cm2` (hydrogen
            atoms per cm^2 along that length). Only the ones the segment actually enters.
    """
    start, end = tuple(float(v) for v in start_pc), tuple(float(v) for v in end_pc)
    length = math.dist(start, end)
    if length <= 0.0:
        return []
    direction = tuple((e - s) / length for s, e in zip(start, end))
    found = []
    for kind in PASS_THROUGH_KINDS:
        table = _PASS_THROUGH_TABLES[kind]
        widest = ly_to_pc(_scalar(conn, f"SELECT MAX(radius_ly) AS v FROM {table} WHERE center_x_pc IS NOT NULL"))
        if widest <= 0.0:
            continue
        density = "density_cm3" if kind == "nebula" else "NULL"
        for row in corridor.objects_near_segment(conn, start, end, widest, kinds=[kind]):
            extra = conn.execute(f"SELECT radius_ly, {density} AS density FROM {table} WHERE id = ?", (row["id"],)).fetchone()
            radius = ly_to_pc(float(extra["radius_ly"] or 0.0))
            center = tuple(row["point_pc"])
            along = sum((c - s) * d for c, s, d in zip(center, start, direction))
            off = math.sqrt(sum((c - s - along * d) ** 2 for c, s, d in zip(center, start, direction)))
            if radius <= 0.0 or off >= radius:
                continue
            half = math.sqrt(radius * radius - off * off)
            inside = max(0.0, min(length, along + half) - max(0.0, along - half))
            if inside <= 0.0:
                continue
            entry = {"ref": row["ref"], "kind": kind, "id": row["id"], "name": row["name"], "note": keepout.PASS_THROUGH_NOTE,
                     "inside_length": inside, "inside_length_ly": pc_to_ly(inside), "_along": max(0.0, along - half)}
            if extra["density"]:
                entry["column_cm2"] = float(extra["density"]) * inside * PC_CM
            found.append(entry)
    found.sort(key=lambda entry: (entry["_along"], entry["ref"]))
    for entry in found:
        del entry["_along"]
    return found


def adjusted_course(conn, start_pc, end_pc, kinds=None):
    """
    The course from `start_pc` to `end_pc` bent around the keep-out spheres between them (NAV.6): the obstacles of
    `obstacles_along_segment` handed to the planner (`galaxy.avoid.plan_course`). The `adjusted` object of
    course-avoidance.md section 7, in the galactic frame, parsecs, with `names` (`{ref: name}` of the objects it goes round) and `passed_through`
    (`passed_through_along_segment`). Nothing is generated: space with no placed object has none to avoid.
    """
    found = obstacles_along_segment(conn, start_pc, end_pc, kinds=kinds)
    course = avoid.plan_course(start_pc, end_pc, avoid.spheres_from_obstacles(found), frame="galactic", unit="pc",
                               passed_through=passed_through_along_segment(conn, start_pc, end_pc))
    names = {row["ref"]: row["name"] for row in found}
    course["names"] = {ref: names[ref] for ref in course["avoided"] + course["ignored_below_resolution"] if ref in names}
    return course


def system_keep_outs(conn, scene):
    """
    The keep-out radius, in km, of every body of a system's scene (NAV.27): a star its radiation radius, a planet or moon its Hill radius
    (`db.query.keep_out_radius` at system scale). `{ref: radius_km}`; a comet and a belt keep nothing out (a comet is too small and a belt
    is thin: `keepout.PASS_THROUGH_NOTE`). The scene is `web.maps.systemscene.build_scene`'s.
    """
    from planetgen.db import query

    found = {}
    for star in scene["stars"]:
        found[star["ref"]] = query.keep_out_radius(conn, "star", star["id"], scale="system").radius_km
    for planet in scene["planets"]:
        found[planet["ref"]] = query.keep_out_radius(conn, "planet", planet["id"], scale="system").radius_km
        for moon in planet["moons"]:
            found[moon["ref"]] = query.keep_out_radius(conn, "moon", moon["id"], scale="system").radius_km
    return {ref: radius for ref, radius in found.items() if radius}
