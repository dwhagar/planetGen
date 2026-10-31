# planetgen/web/maps/systemscene.py

"""
The 3D system view's data (MAP.69): one star system as JSON that a client
can draw with orbits in three dimensions.

`GET /api/systems/<id>/scene` returns `build_scene`'s dict. Where the flat
System Map (`systemmap.py`) is drawn in Python as SVG with log-scaled
distances and no z, this keeps every number real: kilometres, the orbit
elements the database stores, and the position each body had at `epoch`
(the time `planetgen.cli.orbits` last advanced the orbits). A client moves
bodies to any other time: a planet or moon from its state vector
(`position_km`, `velocity_kms`, `epoch_unix`), a pair's second star from its
elements (`period_years`, `phase_deg`); that is MAP.70.

Every object carries a `ref` of the form `kind:id` (`star:3`, `planet:12`,
`moon:40`, `belt:5`, `comet:8`) and an `orbit` naming what it goes round:

    {"around": ref or "barycenter", "distance_km", "period_years",
     "inclination_deg", "ascending_node_deg", "phase_deg"}

for a binary pair's circular orbit, the same without `phase_deg` and with
`position_km`, `velocity_kms` and `epoch_unix` for a planet or moon (its
state), or for a comet the
Kepler elements (`kepler`: perihelion, eccentricity, argument of
periapsis, mean anomaly). `position_km` is the object's position relative
to what it goes round, as stored; the scene's frame has the system's
barycenter (a single star: the star) at the origin, except a wide binary,
whose second star orbits the first and whose own bodies orbit each star. A close pair that hangs below another orbit
of the tree (a third star far out) has its bodies go round `node:N`, the pair's centre of mass, not the barycenter.

`tree` (GEN.311) is the stored tree of orbits the stars hang in: `nodes` (`physics.system_tree.OrbitNode` as dicts,
each the relative orbit of its slot-1 child about its slot-0 child, with its own epoch), `members` (one per star: `ref`,
`member_no`, `parent_node_no`, `slot`, `mass_solar`) and `origin`, the ref of the star at the origin (a wide binary's
first star) or `None` for the barycenter. `physics.body_positions.positions_at` places the stars from it.

A star, planet or moon carries `look` (MAP.287), the same body look record the System Map's `data-look` holds
(`bodylook.star_look`, `bodylook.planet_look`), which the client colours and paints the sphere from; a comet carries
`color`, its nucleus colour.
"""

import dataclasses
import math
from datetime import datetime, timezone

from planetgen.db import orbit_trees, query
from planetgen.physics import body_state, habitability, state_vectors
from planetgen.web.maps.bodylook import comet_nucleus_color, planet_look, star_look
from planetgen.web.maps.systemmap import finite_row, put_star_temperatures
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")

_au_to_km = units.converter("au", "km")


def _xyz(row, prefix="position_", suffix="_km"):
    values = [row.get(f"{prefix}{axis}{suffix}") for axis in "xyz"]
    return [0.0 if v is None else float(v) for v in values]


def _orbit(row, around):
    """A planet's or moon's orbit from its `planets`/`moons` row (GEN.289): the elements the row keeps (the orbit's
    semi-major axis, period and plane) and the state vector that is the body, in km and km/s relative to what it goes
    round, with the epoch (Unix seconds, or None before the first orbit update) the state is true at. A client moves the
    state with `kepler.universal_step` (`physics/body_state.py`, `static/orbitpositions.js`).

    The path the browser draws is the ellipse of the state (GEN.727): `eccentricity` and `arg_periapsis_deg` are read off the state
    (`state_vectors.elements_from_state`), so a pull that changed the orbit changes the drawn path with it; a circle has 0 and 0."""
    position, velocity = _xyz(row), _xyz(row, "velocity_", "_kms")
    eccentricity = arg_periapsis_deg = 0.0
    try:
        elements = state_vectors.elements_from_state(position, velocity, body_state.mu_km3_s2(row["distance_km"], row["period_years"]))
        if elements["kind"] == "elliptical":
            eccentricity, arg_periapsis_deg = elements["eccentricity"], math.degrees(elements["arg_periapsis"])
    except ValueError:
        pass   # a body with no motion yet has no ellipse to draw: its orbit is the circle of its distance
    return {
        "around": around,
        "distance_km": row["distance_km"], "period_years": row["period_years"],
        "inclination_deg": row["orbital_inclination_deg"],
        "ascending_node_deg": row["orbital_ascending_node_deg"],
        "eccentricity": eccentricity, "arg_periapsis_deg": arg_periapsis_deg,
        "position_km": position, "velocity_kms": velocity,
        "epoch_unix": row.get("epoch_unix"),
    }


def _pair_around(tree, stars):
    """What a close pair's own bodies go round: the pair's centre of mass. That is the barycenter when the pair is the root
    of the tree of orbits, else the node it hangs in (`node:N`), which a further star holds the pair's centre away from the
    barycenter of."""
    roles = {f"star:{star['id']}": star["role"] for star in stars}
    nodes = {member["parent_node_no"] for member in tree["members"] if roles.get(member["ref"]) in ("primary", "secondary")}
    if len(nodes) != 1:
        return "barycenter"
    node_no = nodes.pop()
    parent = next((node["parent_node_no"] for node in tree["nodes"] if node["node_no"] == node_no), None)
    return "barycenter" if parent is None else f"node:{node_no}"


def _star_around(planet_or_belt, stars, system, pair_around="barycenter"):
    """What a planet or belt goes round: its own star in a wide binary, a
    single-star system or a further star of the tree (GEN.332), else the close pair's centre of mass."""
    star_id = planet_or_belt.get("star_id")
    if system["binary_configuration"] == "close" and not any(
            s["id"] == star_id and s["role"] == "member" for s in stars):    # a member star's own planets (GEN.332)
        return pair_around
    if star_id is not None and any(s["id"] == star_id for s in stars):
        return f"star:{star_id}"
    return f"star:{stars[0]['id']}" if stars else "barycenter"


def _offset_km(row):
    """A star row's stored place in the whole tree (`offset_*_au`), in km."""
    return [_au_to_km(row[f"offset_{axis}_au"] or 0.0) for axis in "xyz"]


def _stars(conn, system_id, stars, row):
    """The scene's stars, positioned in the system frame."""
    offsets = {r["id"]: r for r in conn.execute(
        "SELECT id, offset_x_au, offset_y_au, offset_z_au FROM stars WHERE star_system_id = ?", (system_id,)).fetchall()}
    in_tree = any(star["role"] == "member" for star in stars)    # a pair inside a larger tree sits where the tree puts it
    primary_xyz = _xyz(row, "binary_primary_position_")
    secondary_xyz = _xyz(row, "binary_secondary_position_")
    mutual_xyz = _xyz(row, "binary_mutual_position_")
    out = []
    for star in stars:
        ref = f"star:{star['id']}"
        item = {
            "ref": ref, "id": star["id"], "kind": "star", "name": star["name"], "role": star["role"],
            "star_type": star["star_type"], "radius_km": star["radius_km"], "mass_kg": star["mass_kg"],
            "temperature_k": star["temperature_k"], "luminosity_w": star["luminosity_w"], "look": star_look(finite_row(star)),
            "position_km": [0.0, 0.0, 0.0], "orbit": None,
        }
        if row["is_binary"] and star["role"] == "secondary":
            if row["binary_configuration"] == "close":
                item["position_km"] = _offset_km(offsets[star["id"]]) if in_tree else secondary_xyz
            else:
                item["position_km"] = mutual_xyz
            item["orbit"] = {
                "around": "barycenter" if row["binary_configuration"] == "close" else f"star:{stars[0]['id']}",
                "distance_km": row["binary_separation_km"],
                "period_years": row["binary_mutual_orbital_period_years"],
                "inclination_deg": row["binary_mutual_orbital_inclination_deg"],
                "ascending_node_deg": row["binary_mutual_orbital_ascending_node_deg"],
                "phase_deg": row["binary_mutual_orbital_phase_deg"],
                "secondary_mass_fraction": row["binary_secondary_mass_fraction"],
            }
        elif star["role"] == "member":    # GEN.630: a star beyond the pair, at the offset its row stores
            item["position_km"] = _offset_km(offsets[star["id"]])
        elif row["is_binary"] and row["binary_configuration"] == "close":
            item["position_km"] = _offset_km(offsets[star["id"]]) if in_tree else primary_xyz
        out.append(item)
    return out


def _tree(conn, system_id, stars, row):
    """The scene's `tree`: the system's `orbit_nodes` and where each star hangs in them, the node and the slot its row
    stores (a pair's primary is slot 0, its secondary slot 1), in the order of its id."""
    nodes = orbit_trees.read_nodes(conn, system_id)
    places = {r["id"]: r for r in conn.execute(
        "SELECT id, parent_node_id, slot FROM stars WHERE star_system_id = ?", (system_id,)).fetchall()}
    members = []
    for number, star in enumerate(sorted(stars, key=lambda s: s["id"])):
        parent, slot = places[star["id"]]["parent_node_id"], places[star["id"]]["slot"] or 0
        members.append({"ref": f"star:{star['id']}", "member_no": number, "parent_node_no": parent, "slot": slot,
                        "mass_solar": _kg_to_msun(star["mass_kg"])})
    origin = members[0]["ref"] if row["binary_configuration"] == "wide" else None
    return {"nodes": [dataclasses.asdict(node) for node in nodes], "members": members, "origin": origin}


def _planet(planet, around, kind="planet", parent=None):
    item = {
        "ref": f"{kind}:{planet['id']}", "id": planet["id"], "kind": kind, "name": planet["name"],
        "planet_class": planet.get("planet_class"), "radius_km": planet["radius_km"], "mass_kg": planet.get("mass_kg"),
        "look": planet_look(finite_row(planet), kind), "orbit": _orbit(planet, around),
        "position_km": _xyz(planet),
        "habitable": bool(planet.get("habitable")), "inhabited": bool(planet.get("inhabited")),
        # The symbol a body carries on the Galaxy Map (MAP.321): terrestrial or gas giant, its overall habitability and the colour of it.
        "body_type": planet.get("body_type"), "phi4": planet.get("phi4"),
        "tier": habitability.tier(planet["phi4"]) if planet.get("phi4") is not None else None,
        # What the map's reasons to look read (MAP.343, as `db.interests.REASONS` does for the Galaxy Map): the stored life scores (0 none to 9
        # abundant, None for no score) and the stage of the timeline, kept for a class that can carry life only.
        "life_microbial": planet.get("life_microbial"), "life_flora": planet.get("life_flora"),
        "life_fauna": planet.get("life_fauna"), "life_sentient": planet.get("life_sentient"),
        "life_stage": planet.get("life_stage"),
        # The life class that lives on the body, as one id such as "1a" or None (MAP.351, `interests.LIFE_CLASS_PREFIX`).
        "life_class_id": planet.get("life_class_id"),
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
        "radius_km": comet["nucleus_diameter_km"] / 2, "color": comet_nucleus_color(comet), "is_active": bool(comet["is_active"]),
        "orbit": {
            "around": around, "type": comet["orbit_type"],
            "kepler": {
                "perihelion_distance_km": comet["perihelion_distance_km"], "eccentricity": comet["eccentricity"],
                "inclination_deg": comet["inclination_deg"], "arg_periapsis_deg": comet["arg_periapsis_deg"],
                "ascending_node_deg": comet["ascending_node_deg"],
                "period_years": comet["orbital_period_years"], "mean_anomaly_deg": comet["mean_anomaly_deg"],
                "parabolic_mean_anomaly": comet["parabolic_mean_anomaly"],
                "primary_mass_solar": comet["primary_mass_solar"],
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
    put_star_temperatures(stars, detail["planets"])
    tree = _tree(conn, system_id, stars, row)
    pair_around = _pair_around(tree, stars)
    planets = []
    for planet in detail["planets"]:
        item = _planet(planet, _star_around(planet, stars, row, pair_around))
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
        "stars": _stars(conn, system_id, stars, row),
        "tree": tree,
        "planets": planets,
        "belts": [_belt(b, _star_around(b, stars, row, pair_around)) for b in detail["belts"]],
        "comets": [_comet(c, _star_around(c, stars, row, pair_around)) for c in detail["comets"]],
    }
