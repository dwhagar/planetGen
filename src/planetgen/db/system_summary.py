# planetgen/db/system_summary.py

"""
What each star system adds up to, stored beside it (GEN.1478, audit F1 to F4; `system_summary` in schema.sql)
=============================================================================================================

Four places used to add up the stars, planets and moons of a system each time they needed a system's mass (the gravity sources, the sector
paths' point masses, the object index, the sector totals), and the sky read counted a system's stars with a correlated query for every
system. The sums are kept in `system_summary`, a narrow row per system (`star_systems` is wide, with a FULLTEXT key), worked out in the
statement `REFRESH_SQL` from the rows themselves:

* `stars_mass_kg`: the summed mass of the stars (the point mass the sector paths use);
* `system_mass_kg`: the stars, planets and moons together (the point mass of the gravity map: each table's `SUM`, added in that order);
* `brightest_star_luminosity_w`, `total_luminosity_w`: the brightest star's and all the stars' luminosity, W;
* `star_count`, `planet_count`, `moon_count`, `belt_count`, `comet_count`.

Whatever writes a system's bodies writes the row after: the save of a new system (`store.insert_star_system`) and the editor
(`edits.save_system_edits`); no other code changes a body's mass, luminosity or whether it exists. A system's row goes with it. A galaxy saved
before the columns existed is filled by migration 0237, in the same statement.
"""

COLUMNS = ("stars_mass_kg", "system_mass_kg", "brightest_star_luminosity_w", "total_luminosity_w",
           "star_count", "planet_count", "moon_count", "belt_count", "comet_count")
"""tuple: The summary columns of `system_summary`, in the order `REFRESH_SQL` writes them."""

SELECT_SQL = (
    "SELECT ss.id AS star_system_id,"
    " (SELECT SUM(b.mass_kg) FROM stars b WHERE b.star_system_id = ss.id) AS stars_mass_kg,"
    " COALESCE((SELECT SUM(b.mass_kg) FROM stars b WHERE b.star_system_id = ss.id), 0)"
    " + COALESCE((SELECT SUM(b.mass_kg) FROM planets b WHERE b.star_system_id = ss.id), 0)"
    " + COALESCE((SELECT SUM(b.mass_kg) FROM moons b WHERE b.star_system_id = ss.id), 0) AS system_mass_kg,"
    " (SELECT MAX(b.luminosity_w) FROM stars b WHERE b.star_system_id = ss.id) AS brightest_star_luminosity_w,"
    " (SELECT SUM(b.luminosity_w) FROM stars b WHERE b.star_system_id = ss.id) AS total_luminosity_w,"
    " (SELECT COUNT(*) FROM stars b WHERE b.star_system_id = ss.id) AS star_count,"
    " (SELECT COUNT(*) FROM planets b WHERE b.star_system_id = ss.id) AS planet_count,"
    " (SELECT COUNT(*) FROM moons b WHERE b.star_system_id = ss.id) AS moon_count,"
    " (SELECT COUNT(*) FROM asteroid_belts b WHERE b.star_system_id = ss.id) AS belt_count,"
    " (SELECT COUNT(*) FROM comets b WHERE b.star_system_id = ss.id) AS comet_count"
    " FROM star_systems ss")
"""str: The summary of every system, `FROM star_systems ss`, ready for a `WHERE`."""

REFRESH_SQL = ("INSERT INTO system_summary (star_system_id, " + ", ".join(COLUMNS) + ") {select} WHERE {where}"
               " ON DUPLICATE KEY UPDATE " + ", ".join(f"{name} = VALUES({name})" for name in COLUMNS))
"""str: Writes the summary of the systems matching `{where}` (SQL on `star_systems` aliased `ss`)."""


def refresh(conn, system_ids):
    """Works the summary of the systems `system_ids` out again from their rows (see the module doc). Part of the caller's transaction."""
    ids = sorted(set(system_ids))
    for start in range(0, len(ids), 500):
        chunk = ids[start:start + 500]
        conn.execute(REFRESH_SQL.format(select=SELECT_SQL, where=f"ss.id IN ({', '.join('?' for _ in chunk)})"), tuple(chunk))
