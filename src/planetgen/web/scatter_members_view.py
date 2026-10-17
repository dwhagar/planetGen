# planetgen/web/scatter_members_view.py

"""
The members of an unbuilt scatter row as a page shows them (UX.114): the stars and compact objects the scatter drew
around one point and the orbits that bind them (`db.query.scatter_member_tree`). Their positions are not shown:
the fill places them when it builds the row's system, and the page says so.
"""

from planetgen.web.lib.fmt import format_distance_au, format_mass_solar, format_number, format_period_years

_KIND_LABELS = {
    "main-sequence": "Main-sequence star", "giant": "Giant", "white-dwarf": "White dwarf",
    "neutron-star": "Neutron star", "black-hole": "Black hole", "brown-dwarf": "Brown dwarf",
    "pre-main-sequence": "Pre-main-sequence star", "protostar": "Protostar", "stripped-helium": "Stripped helium star",
    "thorne-zytkow": "Thorne-Zytkow object", "wolf-rayet": "Wolf-Rayet star",
}
"""dict: `scatter_members.kind` -> its label."""

_ORBIT_LABELS = {"close": "Close pair", "wide": "Wide pair", "compact": "Compact pair", "group": "Young group"}
"""dict: `scatter_nodes.kind` -> its label."""


def members_view(tree):
    """
    The page rows for a `scatter_member_tree` (or a detail that carries one): `members` (number, kind, mass) and
    `orbits` (number, what it joins, separation, period, eccentricity), as lists of strings. Both are empty for a
    single object.
    """
    members = [[str(row["member_no"] + 1),
                _KIND_LABELS.get(row["kind"], row["kind"].replace("-", " ").capitalize()),
                format_mass_solar(row["mass_solar"])] for row in tree.get("members") or []]
    orbits = [[str(row["node_no"] + 1),
               _ORBIT_LABELS.get(row["kind"], row["kind"].capitalize()),
               format_distance_au(row["separation_au"]), format_period_years(row["period_years"]),
               format_number(row["eccentricity"], ".2f")] for row in tree.get("nodes") or []]
    return {"members": members, "orbits": orbits}
