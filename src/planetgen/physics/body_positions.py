# planetgen/physics/body_positions.py

"""
Where every body of a system is at any time (MAP.70), from the scene
`GET /api/systems/<id>/scene` returns (`web/maps/systemscene.py`).

`positions_at(scene, years)` moves each body `years` past the scene's
`epoch` along its stored orbit and returns `{ref: [x, y, z]}` in
kilometres in the scene's frame (the system barycenter, or the first star
for a wide binary, at the origin):

- A planet or moon goes round on the orbit its stored state vector says
  (GEN.289): the state is moved from its own epoch to the requested time
  by `kepler.universal_step` (`body_state.advance`).
- The stars sit where the scene's tree of orbits puts them (GEN.311): each node is the relative Kepler orbit of its
  slot-1 child about its slot-0 child, advanced from the node's own epoch (`system_tree.member_offsets_au`, the
  function that places the members of a scatter row).
- A comet goes round on its Kepler orbit: an elliptical one's mean anomaly
  advances by 360 degrees a period (`kepler.comet_orbital_state`), a
  parabolic one's parabolic mean anomaly by `sqrt(mu / 2q^3)` a year.
- A close pair's stars sit either side of the barycenter in proportion to
  the other's mass, a wide pair's second star orbits the first (the tree's `origin` says which star is at the origin,
  or `None` for the barycenter).

`static/orbitpositions.js` is the browser twin; `tests/test_js_unit.py`
checks the two against each other. NAV.6 plans courses on this one.
"""

import math

from planetgen.physics import body_state, unit_registry as units
from planetgen.physics.kepler import comet_orbital_state, gravitational_parameter_au3_yr2, wrapped_turn_deg
from planetgen.physics.system_tree import OrbitNode, TreeMember, member_offsets_au

_au_to_km = units.converter("au", "km")
_km_to_au = units.converter("km", "au")
_yr_to_s = units.converter("yr", "s")


def _star_positions_km(scene, years):
    """Every star's position in the scene's frame, km, from the scene's tree of orbits: `{ref: [x, y, z]}`. The tree's
    nodes hold their phase at their own epoch; `years` count from the scene's epoch (the first node's epoch when the
    scene has none: the database has never been updated)."""
    tree = scene["tree"]
    nodes = [OrbitNode(**node) for node in tree["nodes"]]
    members = [TreeMember(**{key: member[key] for key in ("member_no", "parent_node_no", "slot", "mass_solar")})
               for member in tree["members"]]
    start = scene.get("epoch_unix")
    if start is None:
        start = nodes[0].epoch_unix if nodes else 0.0
    offsets = member_offsets_au(members, nodes, start + _yr_to_s(years))
    refs = {member["member_no"]: member["ref"] for member in tree["members"]}
    origin = next((number for number, ref in refs.items() if ref == tree.get("origin")), None)
    shift = offsets[origin] if origin is not None else (0.0, 0.0, 0.0)
    return {refs[number]: [_au_to_km(offset[axis] - shift[axis]) for axis in range(3)] for number, offset in offsets.items()}


def _node_centre_km(tree, positions, node_no):
    """The centre of mass of tree node `node_no` (a pair's centre, which bodies can go round: `node:N`), km: the mean of the
    stars below it, each by its mass."""
    below = {node_no}
    grew = True
    while grew:
        grew = False
        for node in tree["nodes"]:
            if node["parent_node_no"] in below and node["node_no"] not in below:
                below.add(node["node_no"])
                grew = True
    held = [member for member in tree["members"] if member["parent_node_no"] in below]
    total = sum(member["mass_solar"] for member in held)
    return [sum(member["mass_solar"] * positions[member["ref"]][axis] for member in held) / total for axis in range(3)]


def _relative_state_km(orbit, years, scene_epoch_unix):
    """A planet's or moon's position relative to what it goes round, in km: its stored state moved from the epoch it
    is true at to `years` past the scene's epoch."""
    epoch = orbit.get("epoch_unix")
    offset = 0.0 if epoch is None or scene_epoch_unix is None else scene_epoch_unix - epoch
    position, _velocity = body_state.advance(
        orbit["position_km"], orbit["velocity_kms"], orbit["distance_km"], orbit["period_years"],
        _yr_to_s(years) + offset)
    return list(position)


def _comet_relative_km(orbit, years):
    """A comet's position relative to its star, in km."""
    k = orbit["kepler"]
    q_au = _km_to_au(k["perihelion_distance_km"])
    if orbit["type"] == "elliptical":
        anomaly = math.radians(k["mean_anomaly_deg"] + wrapped_turn_deg(years, k["period_years"]))
        state = comet_orbital_state(
            "elliptical", q_au, k["eccentricity"], k["inclination_deg"], k["arg_periapsis_deg"],
            k["ascending_node_deg"], k["primary_mass_solar"], mean_anomaly_rad=anomaly,
            orbital_period_years=k["period_years"])
    else:
        mu = gravitational_parameter_au3_yr2(k["primary_mass_solar"])
        anomaly = k["parabolic_mean_anomaly"] + math.sqrt(mu / (2 * q_au ** 3)) * years
        state = comet_orbital_state(
            "parabolic", q_au, 1.0, k["inclination_deg"], k["arg_periapsis_deg"],
            k["ascending_node_deg"], k["primary_mass_solar"], parabolic_mean_anomaly_value=anomaly)
    return [_au_to_km(state["position_x_au"]), _au_to_km(state["position_y_au"]), _au_to_km(state["position_z_au"])]


def _add(a, b):
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]


def positions_at(scene, years):
    """
    Every star, planet, moon and comet's position `years` after the
    scene's epoch (negative for before), in km, as `{ref: [x, y, z]}`.
    Belts have no single point and are left out.
    """
    out = _star_positions_km(scene, years)

    def origin(around):
        if around == "barycenter":
            return [0.0, 0.0, 0.0]
        return _node_centre_km(scene["tree"], out, int(around[5:])) if around.startswith("node:") else out[around]

    epoch_unix = scene.get("epoch_unix")
    for planet in scene["planets"]:
        out[planet["ref"]] = _add(origin(planet["orbit"]["around"]), _relative_state_km(planet["orbit"], years, epoch_unix))
        for moon in planet["moons"]:
            out[moon["ref"]] = _add(out[planet["ref"]], _relative_state_km(moon["orbit"], years, epoch_unix))
    for comet in scene["comets"]:
        out[comet["ref"]] = _add(origin(comet["orbit"]["around"]), _comet_relative_km(comet["orbit"], years))
    return out
