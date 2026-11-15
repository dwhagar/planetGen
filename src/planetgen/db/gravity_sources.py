# planetgen/db/gravity_sources.py

"""
The point masses of the gravity map: one per system (MAP.386; docs/design/orbits/gravity-map.md, `reports/perf72-gravity-grid-cost-2026-11-05.md`)
=====================================================================================================================================================

The sector grid of the gravity map sums point masses over zones 0.25 pc wide. A star system is smaller than that by a factor of a thousand or more (a
planet at 50 AU is 2.4e-4 pc from its star), so for the grid a system is **one point mass at the system's place with the mass of everything in it**: its
stars, its planets and its moons. Folding the planets in matters (their mass is about 0.1 to 0.3 percent of the system's, more than the 0.1 percent the
folded field may be off by at the zone size) and the offsets inside the system do not (the error of putting the mass at the center instead of in its
places is the offset over the distance, 1e-3 at 0.25 pc, times the planets' share, so under one part in a million). The System Map's orbital-plane version
reads the bodies themselves (MAP.169) and does not come through here.

`system_point_masses(conn, sector_ids)` reads, for the sectors given:

- every `star_systems` row with a place: its stars' mass plus its planets' plus its moons', at the system's position;
- every `black_holes` and `neutron_stars` row with a place that is not filed on a star (`star_id` null). A remnant filed on a star is that star's row in
  `stars`, whose mass the system already holds.

Rogue planets, interstellar comets and asteroid fields are left out: the heaviest is a thousandth of a solar mass, below the folded field's tolerance
beside any star within the sector grid's reach (a stand-in; they enter when the grid has a reason to show them).

Positions are the galaxy frame in parsecs, masses in solar masses, `keys` the `(table, id)` of the row; rows come in key order, so the same catalogue
gives the same arrays. Nothing is stored or cached here.
"""

import dataclasses

import numpy as np

from planetgen.db import kind_tables
from planetgen.galaxy import spherical as geometry
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")
STAND_ALONE_TABLES = kind_tables.tables("star_filed")
"""tuple: The remnant tables whose rows can stand alone; one with a `star_id` is already a row of `stars`."""


@dataclasses.dataclass(frozen=True)
class PointMasses:
    """The point masses of some sectors: `positions_pc` `(n, 3)`, `masses_solar` `(n,)`, and `keys`, `n` pairs of table name and row id."""

    positions_pc: np.ndarray
    masses_solar: np.ndarray
    keys: tuple


def system_point_masses(conn, sector_ids):
    """The `PointMasses` of the sectors `sector_ids` (see the module doc); empty when none is given or none holds anything."""
    ids = sorted(set(sector_ids))
    if not ids:
        return PointMasses(np.zeros((0, 3)), np.zeros(0), ())
    marks = ", ".join("?" for _ in ids)
    centers = {row["id"]: (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]) for row in conn.execute(
        f"SELECT id, center_x_pc, center_y_pc, center_z_pc FROM sectors WHERE id IN ({marks})", ids).fetchall()}
    rows = []
    for row in conn.execute(
            f"SELECT ss.id, ss.sector_id, ss.position_x_mpc AS x, ss.position_y_mpc AS y, ss.position_z_mpc AS z, sm.system_mass_kg AS mass_kg"
            f" FROM star_systems ss JOIN system_summary sm ON sm.star_system_id = ss.id"
            f" WHERE ss.sector_id IN ({marks}) AND ss.position_x_mpc IS NOT NULL ORDER BY ss.id", ids).fetchall():
        mass_kg = row["mass_kg"]       # GEN.1478: the stars', planets' and moons' mass, kept with the system (`db/system_summary.py`)
        if row["sector_id"] in centers and mass_kg is not None and mass_kg > 0.0:
            rows.append((("star_systems", row["id"]), geometry.local_to_galaxy_pc(centers[row["sector_id"]], (row["x"] / 1000.0, row["y"] / 1000.0, row["z"] / 1000.0)),
                         _kg_to_msun(mass_kg)))
    for table in STAND_ALONE_TABLES:
        for row in conn.execute(
                f"SELECT id, center_x_pc, center_y_pc, center_z_pc, mass_solar FROM {table}"
                f" WHERE sector_id IN ({marks}) AND star_id IS NULL AND center_x_pc IS NOT NULL ORDER BY id", ids).fetchall():
            if row["mass_solar"] and row["mass_solar"] > 0.0:
                rows.append(((table, row["id"]), (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"]), float(row["mass_solar"])))
    rows.sort(key=lambda entry: entry[0])
    return PointMasses(np.array([entry[1] for entry in rows], dtype=np.float64).reshape(-1, 3), np.array([entry[2] for entry in rows], dtype=np.float64),
                       tuple(entry[0] for entry in rows))
