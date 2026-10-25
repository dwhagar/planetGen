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


def _member_label(row):
    """The kind's label, with a brown dwarf's L, T or Y type after it (GEN.252): "Brown dwarf (T5)"."""
    label = _KIND_LABELS.get(row["kind"], row["kind"].replace("-", " ").capitalize())
    if row["kind"] == "brown-dwarf" and row.get("state"):
        label += f" ({row['state']})"
    return label


def _orbit_label(row, members, age_gy):
    """The node kind's label; a young group's also says how long it has left (GEN.327): its members' decay age less the
    star's age, "Young group, about 21 thousand years left" (a group is only stored while that is positive)."""
    label = _ORBIT_LABELS.get(row["kind"], row["kind"].capitalize())
    if row["kind"] != "group" or age_gy is None:
        return label
    decay = next((m["decay_age_yr"] for m in members if m["parent_node_no"] == row["node_no"]
                  and m.get("decay_age_yr") is not None), None)
    if decay is None:
        return label
    left = decay - age_gy * 1.0e9
    return f"{label}, {format_period_years(left)} left" if left > 0 else f"{label}, past its decay age"


def members_view(tree):
    """
    The page rows for a `scatter_member_tree` (or a detail that carries one): `members` (number, kind, mass) and
    `orbits` (number, what it joins, separation, period, eccentricity), as lists of strings. Both are empty for a
    single object. A young group (GEN.327) is listed with the life it has left, from the detail's `age_gy`; its period is
    the crossing time and it has no eccentricity.
    """
    members = [[str(row["member_no"] + 1), _member_label(row), format_mass_solar(row["mass_solar"])]
               for row in tree.get("members") or []]
    stored = tree.get("members") or []
    orbits = [[str(row["node_no"] + 1), _orbit_label(row, stored, tree.get("age_gy")),
               format_distance_au(row["separation_au"]), format_period_years(row["period_years"]),
               "\u2013" if row["kind"] == "group" else format_number(row["eccentricity"], ".2f")]
              for row in tree.get("nodes") or []]
    return {"members": members, "orbits": orbits}
