# planetgen/web/system_orbits_view.py

"""
The orbits of a built system's stars as its page lists them (MAP.213): one row per `orbit_nodes` row of
`db.query.system_detail`, naming the two things the node joins (a star, or the pair another node makes), its mass,
separation, period, eccentricity and how close the pair is (`interaction`). The page's tree of bodies lists a star's
planets and moons; this is the tree of the stars themselves, the one the 3D system view places them by.
"""

from planetgen.physics import unit_registry as units
from planetgen.web.lib.fmt import format_distance_au, format_mass_solar, format_number, format_period_years
from planetgen.web.scatter_members_view import ORBIT_LABELS

_km_to_au = units.converter("km", "au")


def _joined(detail, node_no):
    """What fills the two places of a node, slot 0 (the centre) first: a star's name, or `orbit N` for a pair hung below
    the node. `None` for a place nothing hangs in."""
    places = {}
    for star in detail.get("stars") or []:
        if star.get("parent_node_no") == node_no:
            places[star.get("slot") or 0] = star["name"]
    for node in detail.get("orbits") or []:
        if node["parent_node_no"] == node_no:
            places[node["slot"] or 0] = f"orbit {node['node_no'] + 1}"
    return [places.get(0), places.get(1)]


def orbit_rows(detail):
    """The rows of the "Orbits of the stars" panel, each `[number, joins, kind, mass, separation, period, eccentricity,
    interaction]` as strings; empty for a single star or a system built before it stored a tree."""
    rows = []
    for node in detail.get("orbits") or []:
        inner, outer = _joined(detail, node["node_no"])
        interaction = node.get("interaction")
        rows.append([
            str(node["node_no"] + 1), " and ".join(name for name in (inner, outer) if name) or "–",
            ORBIT_LABELS.get(node["kind"], node["kind"].capitalize()), format_mass_solar(node["mass_solar"]),
            format_distance_au(_km_to_au(node["separation_km"])), format_period_years(node["period_years"]),
            format_number(node["eccentricity"], ".2f"),
            interaction.replace("-", " ").capitalize() if interaction else "–"])
    return rows
