# planetgen/describe/heavy_rows.py

"""
The derived rows of the object sheet for the heavy elements of a planet or moon (GEN.1371)
==========================================================================================

`rows(conn, kind, object_id)` reads the stored planet or moon, its rock and metal shares and the star it circles ([Fe/H], age, and the sector's distance from the galactic center, which key the star's
abundance draws as `radiogenic` keys the thorium) and gives the rows `physics/heavy_elements` works out: the rare earths in the rock, the platinum-group metals and the gold in the rock and in the
metal core. They are worked out when the sheet is made and are not columns (the derive-on-demand route Boss's TODO approved for GEN.1371). Each row is `(group, key, label, value, raw)`; none for a
body with no metal or rock share, a star with no [Fe/H] or age, or a kind that is not a planet or moon.
"""

from planetgen.physics import constants, heavy_elements, stellar_composition

ROWS = (
    ("heavy_ree_rock_ppm", "Rare earths in the rock", "rare earths", 0),
    ("heavy_pge_rock_ppm", "Platinum-group metals in the rock", "platinum group", 0),
    ("heavy_pge_core_ppm", "Platinum-group metals in the metal core", "platinum group", 1),
    ("heavy_gold_rock_ppm", "Gold in the rock", "gold", 0),
    ("heavy_gold_core_ppm", "Gold in the metal core", "gold", 1),
)
"""tuple: `(key, label, group, side)` of each row: `side` 0 is the silicate, 1 the metal."""

KEYS = tuple(row[0] for row in ROWS)


def is_derived_key(key):
    """Whether `key` is one of the rows this module gives."""
    return key in KEYS


def _star_and_body(conn, kind, object_id):
    if kind == "planet":
        body = conn.execute("SELECT star_system_id, star_id, metal_mass_fraction, rock_mass_fraction FROM planets WHERE id = ?", (object_id,)).fetchone()
        host = body
    elif kind == "moon":
        body = conn.execute("SELECT parent_type, parent_id, metal_mass_fraction, rock_mass_fraction FROM moons WHERE id = ?", (object_id,)).fetchone()
        if body is None or body["parent_type"] != "planet":
            return None, None
        host = conn.execute("SELECT star_system_id, star_id FROM planets WHERE id = ?", (body["parent_id"],)).fetchone()
    else:
        return None, None
    if body is None or host is None:
        return None, None
    columns = "mass_kg, metallicity_fe_h, age_gy, star_system_id"
    star = conn.execute(f"SELECT {columns} FROM stars WHERE id = ?", (host["star_id"],)).fetchone() if host["star_id"] else None
    if star is None and host["star_system_id"]:
        star = conn.execute(f"SELECT {columns} FROM stars WHERE star_system_id = ? ORDER BY id LIMIT 1", (host["star_system_id"],)).fetchone()
    return star, body


def rows(conn, kind, object_id):
    """The derived heavy-element rows of a planet or moon, `(group, key, label, value, raw)` each (see the module text)."""
    from planetgen.describe import model, sheet            # the sheet imports this module's neighbours, so it is read late
    star, body = _star_and_body(conn, kind, object_id)
    if star is None or body is None or star["metallicity_fe_h"] is None or star["age_gy"] is None or not star["mass_kg"]:
        return []
    metal, rock = body["metal_mass_fraction"], body["rock_mass_fraction"]
    if not metal or not rock or metal <= 0.0 or rock <= 0.0:
        return []
    key = stellar_composition.abundance_key(mass=star["mass_kg"], distance_ly=sheet._star_distance_ly(conn, star["star_system_id"]))
    ratios = heavy_elements.star_to_mg(float(star["metallicity_fe_h"]), float(star["age_gy"]), key)
    groups = heavy_elements.group_ppm(heavy_elements.bulk_ppm(ratios, float(rock), float(metal)))
    found = []
    for row_key, label, group, side in ROWS:
        ppm = groups[group][side]
        value = [sheet._plain_number(ppm, "ppm")]
        found.append(("composition", row_key, label, value, ppm))
    return found
