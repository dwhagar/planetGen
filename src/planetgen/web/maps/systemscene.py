# planetgen/web/maps/systemscene.py

"""
The 3D system view's data (MAP.69): one star system as JSON that a client
can draw with orbits in three dimensions.

`GET /api/systems/<id>/scene` returns `build_scene`'s dict. Where the flat
System Map (`systemmap.py`) is drawn in Python as SVG with log-scaled
distances and no z, this keeps every number real: kilometres, the orbit
elements the database stores, and the position each body had at `epoch`
(the time `planetgen.cli.orbits` last advanced the phases). A client moves
bodies to any other time from the elements (`period_years`, `phase_deg`);
that is MAP.70.

Every object carries a `ref` of the form `kind:id` (`star:3`, `planet:12`,
`moon:40`, `belt:5`, `comet:8`) and an `orbit` naming what it goes round:

    {"around": ref or "barycenter", "distance_km", "period_years",
     "inclination_deg", "ascending_node_deg", "phase_deg"}

for circular orbits (planets, moons, a binary pair), or for a comet the
Kepler elements (`kepler`: perihelion, eccentricity, argument of
periapsis, mean anomaly). `position_km` is the object's position relative
to what it goes round, as stored; the scene's frame has the system's
barycenter (a single star: the star) at the origin, except a wide binary,
whose second star orbits the first and whose own bodies orbit each star.
"""

from datetime import datetime, timezone

from planetgen.db import query
from planetgen.web.maps.starmap import star_color
from planetgen.web.maps.systemmap import class_color


def _xyz(row, prefix="position_", suffix="_km"):
    values = [row.get(f"{prefix}{axis}{suffix}") for axis in "xyz"]
    return [0.0 if v is None else float(v) for v in values]


def _orbit(row, around):
    """A circular orbit's elements from a `planets`/`moons` row."""
    return {
        "around": around,
        "distance_km": row["distance_km"], "period_years": row["period_years"],
        "inclination_deg": row["orbital_inclination_deg"],
        "ascending_node_deg": row["orbital_ascending_node_deg"],
        "phase_deg": row["orbital_phase_deg"],
    }


def _star_around(planet_or_belt, stars, system):
    """What a planet or belt goes round: its own star in a wide binary or a
    single-star system, else the barycenter."""
    if system["binary_configuration"] == "close":
        return "barycenter"
    star_id = planet_or_belt.get("star_id")
    if star_id is not None and any(s["id"] == star_id for s in stars):
        return f"star:{star_id}"
    return f"star:{stars[0]['id']}" if stars else "barycenter"


def _stars(stars, row):
    """The scene's stars, positioned in the system frame."""
    primary_xyz = _xyz(row, "binary_primary_position_")
    secondary_xyz = _xyz(row, "binary_secondary_position_")
    mutual_xyz = _xyz(row, "binary_mutual_position_")
    out = []
    for star in stars:
        ref = f"star:{star['id']}"
        fill, _stroke = star_color(star["star_type"], star["temperature_k"], star["luminosity_w"])
        item = {
            "ref": ref, "id": star["id"], "kind": "star", "name": star["name"], "role": star["role"],
            "star_type": star["star_type"], "radius_km": star["radius_km"], "mass_kg": star["mass_kg"],
            "temperature_k": star["temperature_k"], "luminosity_w": star["luminosity_w"], "color": fill,
            "position_km": [0.0, 0.0, 0.0], "orbit": None,
        }
        if row["is_binary"] and star["role"] == "secondary":
            if row["binary_configuration"] == "close":
                item["position_km"] = secondary_xyz
            else:
                item["position_km"] = mutual_xyz
            item["orbit"] = {
                "around": "barycenter" if row["binary_configuration"] == "close" else f"star:{stars[0]['id']}",
                "distance_km": row["binary_separation_km"],
                "period_years": row["binary_mutual_orbital_period_years"],
                "inclination_deg": row["binary_mutual_orbital_inclination_deg"],
                "ascending_node_deg": row["binary_mutual_orbital_ascending_node_deg"],
                "phase_deg": row["binary_mutual_orbital_phase_deg"],
            }
        elif row["is_binary"] and row["binary_configuration"] == "close":
            item["position_km"] = primary_xyz
        out.append(item)
    return out


def _planet(planet, around, kind="planet", parent=None):
    item = {
        "ref": f"{kind}:{planet['id']}", "id": planet["id"], "kind": kind, "name": planet["name"],
        "planet_class": planet.get("planet_class"), "radius_km": planet["radius_km"], "mass_kg": planet.get("mass_kg"),
        "color": class_color(planet.get("planet_class")), "orbit": _orbit(planet, around),
        "position_km": _xyz(planet),
        "habitable": bool(planet.get("habitable")), "inhabited": bool(planet.get("inhabited")),
    }
    if parent is not None:
        item["parent"] = parent
    return item


def _belt(belt, around):
    return {
        "ref": f"belt:{belt['id']}", "id": belt["id"], "kind": "belt", "name": belt.get("name"),
        "body_type": belt.get("body_type"), "around": around,
        "inner_km": belt["lower_limit_km"], "outer_km": belt["upper_limit_km"], "distance_km": belt["distance_km"],
    }


def _comet(comet, around):
    return {
        "ref": f"comet:{comet['id']}", "id": comet["id"], "kind": "comet", "name": comet["name"],
        "radius_km": comet["nucleus_diameter_km"] / 2, "color": "#b9d6e8", "is_active": bool(comet["is_active"]),
        "orbit": {
            "around": around, "type": comet["orbit_type"],
            "kepler": {
                "perihelion_distance_km": comet["perihelion_distance_km"], "eccentricity": comet["eccentricity"],
                "inclination_deg": comet["inclination_deg"], "arg_periapsis_deg": comet["arg_periapsis_deg"],
                "ascending_node_deg": comet["ascending_node_deg"],
                "period_years": comet["orbital_period_years"], "mean_anomaly_deg": comet["mean_anomaly_deg"],
                "parabolic_mean_anomaly": comet["parabolic_mean_anomaly"],
            },
        },
        "position_km": _xyz(comet),
    }


def _epoch(conn):
    """When the stored phases were last advanced, as (ISO UTC text, Unix
    seconds), or (None, None) when `planetgen.cli.orbits` has never run."""
    row = conn.execute(
        "SELECT UNIX_TIMESTAMP(last_updated_at) AS unix FROM orbit_simulation_state WHERE id = 1"
    ).fetchone()
    if row is None or row["unix"] is None:
        return None, None
    unix = int(row["unix"])
    return datetime.fromtimestamp(unix, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), unix


def build_scene(conn, system_id):
    """
    One system's 3D scene.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        system_id (int): The `star_systems.id`.

    Returns:
        dict: `system` (id, name, ref, binary facts, heliopause), `epoch`
            and `epoch_unix`, then `stars`, `planets` (each with its
            `moons`), `belts` and `comets`, as described in the module
            docstring.

    Raises:
        ValueError: If no such system exists.
    """
    detail = query.system_detail(conn, system_id)
    row = dict(conn.execute("SELECT * FROM star_systems WHERE id = ?", (system_id,)).fetchone())
    stars = detail["stars"]
    planets = []
    for planet in detail["planets"]:
        item = _planet(planet, _star_around(planet, stars, row))
        item["moons"] = [_planet(moon, item["ref"], kind="moon", parent=item["ref"]) for moon in planet["moons"]]
        planets.append(item)
    epoch, epoch_unix = _epoch(conn)
    return {
        "system": {
            "id": detail["id"], "ref": f"system:{detail['id']}", "name": detail["name"],
            "is_binary": bool(detail["is_binary"]), "binary_configuration": detail["binary_configuration"],
            "heliopause_au": detail["heliopause_au"],
        },
        "epoch": epoch, "epoch_unix": epoch_unix,
        "stars": _stars(stars, row),
        "planets": planets,
        "belts": [_belt(b, _star_around(b, stars, row)) for b in detail["belts"]],
        "comets": [_comet(c, _star_around(c, stars, row)) for c in detail["comets"]],
    }
