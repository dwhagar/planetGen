# planetgen/db/sector_mass.py

"""
The total mass of each generated sector (DB.108, `sector_stats.total_mass_kg`).

The Galaxy Map names and draws the sectors of a region densest first (design section 0.4 of the region names), which needs
what each holds: the sum of the masses of its stars, planets, moons, asteroid belts and comets (the bodies of its star
systems) and of its black holes, neutron stars, quasars, rogue planets, interstellar comets and asteroid fields. Nebulae and
supernova remnants are diffuse gas, counted by volume under their own label slice, and are not in the sum.

Every body's mass is in kilograms (a black hole's, neutron star's and quasar's black hole's are stored in solar masses and are
multiplied by `SOLAR_MASS_TO_KG`). A sector's sum is made exactly (the masses are added as decimals and rounded once to a double), so
the same bodies give the same bits whatever order they were generated or stored in.

A sector holding a body with no stored mass cannot be summed and is not guessed: its `total_mass_kg` stays NULL, it counts as
nothing in the block and region sums, and `fill_all` and `refresh_sectors` return it so the caller lists it. A generated sector
with no body at all holds 0 kg.

`refresh_sectors` puts the rows of some sectors right (the hook `block_stats.refresh_chain` and `refresh_sectors` run first, so the
blocks above read fresh sums), `fill_all` makes them all, a shell at a time (`block_stats.rebuild`, `planetgen fill-block-totals`, the
migration).
"""

import decimal

from planetgen.physics import constants

SOLAR = repr(constants.SOLAR_MASS_TO_KG)
"""str: Kilograms in a solar mass, as SQL text."""

SOURCES = (
    ("stars", "system", "b.mass_kg"),
    ("planets", "system", "b.mass_kg"),
    ("moons", "system", "b.mass_kg"),
    ("asteroid_belts", "system", "b.mass_kg"),
    ("comets", "system", "b.mass_kg"),
    ("black_holes", "sector", f"b.mass_solar * {SOLAR}"),
    ("neutron_stars", "sector", f"b.mass_solar * {SOLAR}"),
    ("quasars", "sector", f"b.black_hole_mass_solar * {SOLAR}"),
    ("rogue_planets", "sector", "b.mass_kg"),
    ("interstellar_comets", "sector", "b.mass_kg"),
    ("asteroid_fields", "sector", "b.mass_kg"),
)
"""tuple: `(table, how it reaches its sector, SQL of its mass in kilograms on the table aliased `b`)`: a body of a star
system is reached through `star_systems.sector_id`, a phenomenon holds its `sector_id`."""

CHUNK = 200
"""int: Sector addresses per statement."""

_TOTAL = "CAST(COALESCE({mass}, 0) AS DECIMAL(65, 0))"


def _sums_sql(table, via, mass, where):
    join = ("JOIN star_systems ss ON ss.sector_id = s.id JOIN {table} b ON b.star_system_id = ss.id" if via == "system"
            else "JOIN {table} b ON b.sector_id = s.id").format(table=table)
    return (f"SELECT s.shell_index AS shell, s.band_index AS band, s.wedge_index AS wedge,"
            f" SUM({_TOTAL.format(mass=mass)}) AS total, SUM(({mass}) IS NULL) AS unknown"
            f" FROM sectors s {join} WHERE {where} GROUP BY s.id")


def masses(conn, where, params):
    """
    The summed mass of the sectors matching `where` (SQL on `sectors` aliased `s`, with `params`) that hold a body.

    Returns:
        dict: `(shell, band, wedge)` -> `(kilograms as a decimal, bodies with no stored mass)`.
    """
    found = {}
    for table, via, mass in SOURCES:
        for row in conn.execute(_sums_sql(table, via, mass, where), tuple(params)).fetchall():
            key = (int(row["shell"]), int(row["band"]), int(row["wedge"]))
            total, unknown = found.get(key, (decimal.Decimal(0), 0))
            found[key] = (total + decimal.Decimal(row["total"] or 0), unknown + int(row["unknown"] or 0))
    return found


def sector_mass_kg(conn, sector_id):
    """The total mass of the sector `sector_id` in kilograms (0.0 with no body), or `None` when a body of it has no stored mass."""
    found = masses(conn, "s.id = ?", (sector_id,))
    if not found:
        return 0.0
    total, unknown = next(iter(found.values()))
    return None if unknown else float(total)


def _write(conn, addresses, found):
    """Sets `total_mass_kg` of each generated sector in `addresses`; returns the addresses it could not sum."""
    skipped = []
    for address in addresses:
        total, unknown = found.get(address, (decimal.Decimal(0), 0))
        if unknown:
            skipped.append(address)
        conn.execute(
            "UPDATE sector_stats SET total_mass_kg = ? WHERE shell_index = ? AND band_index = ? AND wedge_index = ?"
            " AND bright_level_sol = 0", (None if unknown else float(total), *address))
    return skipped


def refresh_sectors(conn, addresses):
    """
    Puts `sector_stats.total_mass_kg` right for the generated sectors at `addresses`, each `(shell, band, wedge)`.

    Returns:
        list: The addresses that hold a body with no stored mass (their value is NULL).
    """
    addresses = sorted({tuple(int(value) for value in address) for address in addresses})
    skipped = []
    for start in range(0, len(addresses), CHUNK):
        chunk = addresses[start:start + CHUNK]
        where = "(s.shell_index, s.band_index, s.wedge_index) IN (" + ", ".join("(?, ?, ?)" for _ in chunk) + ")"
        skipped += _write(conn, chunk, masses(conn, where, [value for address in chunk for value in address]))
    return skipped


def fill_all(conn, on_progress=None):
    """
    Makes `sector_stats.total_mass_kg` of every generated sector again, a shell at a time.

    Returns:
        list: The addresses that hold a body with no stored mass, sorted.
    """
    shells = [row["shell_index"] for row in conn.execute(
        "SELECT DISTINCT shell_index FROM sector_stats WHERE bright_level_sol = 0 ORDER BY shell_index").fetchall()]
    skipped = []
    for done, shell in enumerate(shells, start=1):
        addresses = [(int(r["shell_index"]), int(r["band_index"]), int(r["wedge_index"])) for r in conn.execute(
            "SELECT shell_index, band_index, wedge_index FROM sector_stats WHERE shell_index = ? AND bright_level_sol = 0",
            (shell,)).fetchall()]
        skipped += _write(conn, addresses, masses(conn, "s.shell_index = ?", (shell,)))
        if on_progress is not None:
            on_progress(done, len(shells))
    return sorted(skipped)
