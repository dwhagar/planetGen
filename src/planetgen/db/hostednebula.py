# planetgen.db.hostednebula

"""
The nebula a pre-placed star lights, read from its `bright_stars` row (GEN.530).

GEN.320 rolls a star's hosted nebula from the star's own seed (`nebula_contents.hosted_nebula_of`), so the row of a
pre-placed O, B, A or Wolf-Rayet star says which nebula its sector will build before the sector exists. This module
rolls it for rows read by the pages and maps, with the sector profile's `"nebula"` factor of the row's cell, the same
factor its sector rolls with, so the nebula shown before the sector is built is the one it holds after.
"""

import re

from planetgen import tuning
from planetgen.generation import nebula_contents
from planetgen.galaxy.profile import sector_profile

HOST_STAR_PATTERN = "^(" + "|".join(sorted({rule[0] for rule in tuning.NEBULA_HOST_RULES if rule[0] != "W"})) \
    + ")[0-9]V|^W[NCO][0-9]"
"""str: The `star_type`s that can host a nebula (`tuning.NEBULA_HOST_RULES`), as the regular expression the database
filters a cell's rows with."""

ROW_COLUMNS = "id, shell_index, band_index, wedge_index, star_type, age_gy, mass_kg, temperature_k"
"""str: The `bright_stars` columns `hosted_by_rows` reads."""

_HOST_STAR = re.compile(HOST_STAR_PATTERN)


def summary(nebula):
    """The facts of a hosted nebula a page or a map shows, as a dict: `name`, `nebula_class`, `nebula_type` (its family: emission, reflection, ...),
    `class_name`, `radius_ly` and `age_myr`."""
    return {"name": str(nebula.name), "nebula_class": nebula.nebula_class, "nebula_type": nebula.nebula_type,
            "class_name": nebula.class_name, "radius_ly": nebula.radius_ly, "age_myr": nebula.age_myr}


def hosted_by_rows(conn, rows):
    """
    The hosted nebulae of `rows` (`bright_stars` rows with `ROW_COLUMNS`), as `{row id: summary}`; a row that lights
    none is left out. Each is rolled with the `"nebula"` factor of the profile of the row's cell, as its sector rolls
    it.

    Args:
        conn (planetgen.db.store.Connection): An open, read-only connection.
        rows (iterable): The rows.
    """
    from planetgen.db import store
    rows = [row for row in rows if _HOST_STAR.match(row["star_type"] or "")]
    if not rows:
        return {}
    skeleton = store.get_galaxy_shape(conn)
    factors, found = {}, {}
    for row in rows:
        address = (row["shell_index"], row["band_index"], row["wedge_index"])
        if address not in factors:
            factors[address] = (sector_profile(address, skeleton).kind_factors.get("nebula", 1.0)
                                if skeleton is not None else 1.0)
        nebula = nebula_contents.row_nebula(row, factors[address])
        if nebula is not None:
            found[row["id"]] = summary(nebula)
    return found


def hosted_in_cell(conn, address):
    """
    The nebulae the unbuilt pre-placed stars of one cell light (GEN.530): one `{"star_id", "x", "y", "z",
    **summary}` for each, the point being the star's (galaxy-frame parsecs), in the order of the star ids.
    """
    rows = conn.execute(
        f"SELECT {ROW_COLUMNS}, position_x_mpc, position_y_mpc, position_z_mpc FROM bright_stars"
        " WHERE shell_index = ? AND band_index = ? AND wedge_index = ? AND star_system_id IS NULL"
        " AND star_type REGEXP ? ORDER BY id", (*address, HOST_STAR_PATTERN)).fetchall()
    found = hosted_by_rows(conn, rows)
    return [{"star_id": row["id"], **{axis: row[f"position_{axis}_mpc"] / 1000.0 for axis in "xyz"},
             **found[row["id"]]} for row in rows if row["id"] in found]
