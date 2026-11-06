# planetgen/db/sky_tier_ledger.py

"""
The spheres the sky tiers have drawn (VIEW.25)
===============================================

Charting draws the stars a system's sky needs (the tier table of `generation/sky_tiers.py`) over spheres, a luminosity band or a mass ring a sphere, with the star
backfill (`generation/sky_charting.py`). A sphere of a 7.5 sky holds billions of cells, so the draw cannot give each cell a `sector_stats` row; it records the sphere here
(`sky_tier_draws`) instead, and a cell inside it is complete down to the tier's floor.

`Ledger` holds the rows. `Ledger.luminosity_reached(addresses, edge_pc, starts)` and `Ledger.mass_reached(...)` say how deep a cell is complete: from the level the cell
started at (the galaxy scatter's level, or the mass the scatter placed down to) the covering spheres are followed downward while their ranges touch, so a band that no
sphere drew leaves everything under it where it was. A cell is covered by a sphere when its center is (`spherical.cell_in_sphere`, the test the draw used).

The backfill and the fill read this through `settle_levels` before they draw, which writes the depth into the `sector_stats` row of each sector they are about to touch, so the
readers of a sector's level (`store.bright_star_fill_level`, `store.bright_star_fill_mass_limit`) see it without knowing about the spheres.
"""

import dataclasses
import math

import numpy as np

from planetgen.galaxy import spherical

LUMINOSITY, MASS = "luminosity", "mass"
KINDS = (LUMINOSITY, MASS)

EDGE_TOLERANCE = 1e-6
"""float: How far apart two band edges may be and still touch (the same number written twice by `canonical_bands` is exactly equal; this is for the ones read back from the database)."""


@dataclasses.dataclass(frozen=True)
class Draw:
    """One sphere a tier was drawn over: `kind`, the band's range (`low`, `high`, `None` for no upper edge), the sphere's `center` and `radius_pc`, the limiting magnitude and what it cost."""

    kind: str
    low: float
    high: object
    center: tuple
    radius_pc: float
    limit_mag: float = 0.0
    stars: int = 0
    seconds: float = 0.0

    def covers(self, address, edge_pc):
        """Whether the cell `address` lies in the sphere."""
        return spherical.cell_in_sphere(address, self.center, self.radius_pc, edge_pc)

    def contains_sphere(self, center, radius_pc):
        """Whether the whole sphere `(center, radius_pc)` lies inside this one."""
        gap = math.dist(self.center, center)
        return gap + radius_pc <= self.radius_pc + 1e-9


class Ledger:
    """The recorded draws of one kind or both, read with `read`. Build one directly from `Draw`s for a test."""

    def __init__(self, draws=()):
        self.draws = tuple(draws)

    def __bool__(self):
        return bool(self.draws)

    def of_kind(self, kind):
        """The draws of one kind, the highest band first."""
        found = [draw for draw in self.draws if draw.kind == kind]
        found.sort(key=lambda draw: (-draw.low, 0 if draw.high is None else 1))
        return found

    def drawn(self, kind, low, high, center, radius_pc):
        """Whether a recorded sphere of the same band holds the whole sphere `(center, radius_pc)`, so drawing it again would write nothing new."""
        for draw in self.of_kind(kind):
            if abs(draw.low - low) <= EDGE_TOLERANCE and (draw.high is None) == (high is None) \
                    and (high is None or abs(draw.high - high) <= EDGE_TOLERANCE) and draw.contains_sphere(center, radius_pc):
                return True
        return False

    def _reached(self, kind, addresses, edge_pc, starts):
        """
        How low each cell is complete: from `starts[address]` downward through the covering draws of `kind` while their ranges touch. The distances are worked out on arrays; a cell
        within a hair of a sphere's edge is judged by `spherical.cell_in_sphere`, the test the draw used.
        """
        draws = self.of_kind(kind)
        addresses = [a for a in addresses if a in starts]
        if not draws or not addresses:
            return {}
        positions = np.array([spherical.sector_position_pc(*a, edge_pc) for a in addresses], dtype=float)
        start = np.array([starts[a] for a in addresses], dtype=float)
        current = start.copy()
        for draw in draws:
            gap = np.linalg.norm(positions - np.array(draw.center, dtype=float), axis=1) - draw.radius_pc
            covered = gap < -1e-6
            for k in np.nonzero(np.abs(gap) <= 1e-6)[0]:
                covered[k] = draw.covers(addresses[k], edge_pc)
            usable = covered & (draw.low < current - EDGE_TOLERANCE)
            if draw.high is not None:
                usable &= draw.high >= current - EDGE_TOLERANCE
            current = np.where(usable, draw.low, current)
        return {addresses[k]: float(current[k]) for k in np.nonzero(current < start - EDGE_TOLERANCE)[0]}

    def sphere_reached(self, kind, center, radius_pc, start):
        """
        How low a whole sphere `(center, radius_pc)` is complete when its cells started at `start` (a level or mass): the draws of `kind` that hold the whole sphere, followed downward
        while their ranges touch. The answer for a sphere too large to list cell by cell.
        """
        current = float(start)
        for draw in self.of_kind(kind):
            if draw.low >= current - EDGE_TOLERANCE or (draw.high is not None and draw.high < current - EDGE_TOLERANCE):
                continue
            if draw.contains_sphere(center, radius_pc):
                current = draw.low
        return current

    def luminosity_reached(self, addresses, edge_pc, starts):
        """
        `{address: level}` for the cells whose luminosity band draws reach below the level they started at (`starts`, `{address: Lsun}`; a cell without an entry has no scatter level and
        is left out). A cell no sphere reaches below its start is not in the result.
        """
        return self._reached(LUMINOSITY, addresses, edge_pc, starts)

    def mass_reached(self, addresses, edge_pc, starts):
        """`{address: mass}` as `luminosity_reached`, for the mass rings, from the mass the cell was placed down to (`starts`, `{address: Msun}`)."""
        return self._reached(MASS, addresses, edge_pc, starts)


@dataclasses.dataclass(frozen=True)
class Settled:
    """What `settle` found a sector's rows lacking: `levels` (`{address: Lsun}`, for `store.set_sector_bright_levels`) and `masses` (`{address: (Msun, level or None)}`, for `store.set_sector_bright_masses`)."""

    levels: dict
    masses: dict


NOTHING = Settled({}, {})


def settle(ledger, addresses, edge_pc, limits, levels, masses):
    """
    What the sky spheres say about some sectors that their `sector_stats` rows do not yet: the depth the spheres reached is written into the rows before the backfill or the fill reads
    them, so nothing is drawn twice and the readers of a sector's level need not know about the spheres.

    Args:
        ledger (Ledger): `read(conn)`.
        addresses (iterable): The sectors about to be touched.
        edge_pc (float): The sector edge.
        limits (ScatterLimits): The scatter's limits by place.
        levels (dict): `{address: bright_level_sol}` as stored (-1 untouched, 0 filled), for the sectors that have a row.
        masses (dict): `{address: bright_mass_sol}` as stored, for the sectors that have one.

    Returns:
        Settled: Only what the spheres reach below what is stored.
    """
    if not ledger:
        return NOTHING
    addresses = list(addresses)
    lum_start, mass_start, own = {}, {}, {}
    for address in addresses:
        stored = levels.get(address)
        if stored is not None and stored == 0.0:
            continue                                          # filled: its systems were drawn below its level
        galaxy = limits.level_at(address[0], address[1])
        level = stored if stored is not None and stored > 0 else galaxy
        if level is not None:
            lum_start[address] = float(level)
        limit = limits.mass_limit_at(address[0], address[1])
        held = masses.get(address)
        held = limit if held is None else held if limit is None else min(held, limit)
        if held is not None:
            mass_start[address] = float(held)
            own[address] = galaxy if galaxy is not None and (stored is None or stored < 0) else None
    new_levels = ledger.luminosity_reached(addresses, edge_pc, lum_start)
    new_masses = {address: (mass, own[address]) for address, mass in ledger.mass_reached(addresses, edge_pc, mass_start).items()}
    return Settled(new_levels, new_masses)


def read(conn):
    """Every complete draw, as a `Ledger` (a sphere still being drawn does not count)."""
    rows = conn.execute(
        "SELECT kind, band_low, band_high, center_x_pc, center_y_pc, center_z_pc, radius_pc, limit_mag, stars, seconds FROM sky_tier_draws"
        " WHERE complete = 1 ORDER BY kind, band_low, id").fetchall()
    return Ledger(Draw(row["kind"], float(row["band_low"]), None if row["band_high"] is None else float(row["band_high"]),
                       (float(row["center_x_pc"]), float(row["center_y_pc"]), float(row["center_z_pc"])), float(row["radius_pc"]),
                       float(row["limit_mag"]), int(row["stars"]), float(row["seconds"])) for row in rows)


def begin(conn, draw, generator_epoch=None):
    """
    Starts (or finds the start of) a sphere's draw: a row that does not count yet (`complete` 0). A draw stopped part way is found again by its kind, band, center and radius, and goes
    on from the batch it had reached.

    Returns:
        tuple: `(id, batches_done)`.
    """
    key = (draw.kind, draw.low, *draw.center, draw.radius_pc)
    conn.execute(
        "INSERT INTO sky_tier_draws (kind, band_low, band_high, center_x_pc, center_y_pc, center_z_pc, radius_pc, limit_mag, generator_epoch, complete)"
        " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, 0) ON DUPLICATE KEY UPDATE id = LAST_INSERT_ID(id)",
        (draw.kind, draw.low, draw.high, *draw.center, draw.radius_pc, draw.limit_mag, generator_epoch))
    row = conn.execute(
        "SELECT id, batches_done, complete FROM sky_tier_draws WHERE kind = ? AND band_low = ? AND center_x_pc = ? AND center_y_pc = ? AND center_z_pc = ? AND radius_pc = ?",
        key).fetchone()
    return int(row["id"]), int(row["batches_done"])


def advance(conn, draw_id, batches_done, stars, seconds):
    """Notes that the draw has written `batches_done` batches in all (inside the transaction that wrote the last one's stars); `stars` and `seconds` are the totals so far."""
    conn.execute("UPDATE sky_tier_draws SET batches_done = ?, stars = ?, seconds = ? WHERE id = ?", (int(batches_done), int(stars), float(seconds), int(draw_id)))


def finish(conn, draw_id, stars, seconds):
    """Makes the sphere count: every batch is in. The caller commits."""
    conn.execute("UPDATE sky_tier_draws SET complete = 1, stars = ?, seconds = ?, drawn_at = CURRENT_TIMESTAMP(3) WHERE id = ?", (int(stars), float(seconds), int(draw_id)))


def forget_all(conn):
    """Removes every record (the stars they drew are gone, `store.clear_bright_stars`)."""
    conn.execute("DELETE FROM sky_tier_draws")
