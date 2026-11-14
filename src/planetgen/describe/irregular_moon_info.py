# planetgen/describe/irregular_moon_info.py

"""
The irregular moons of a planet, as a list and as pages of one body each (GEN.1086, docs/design/systems/captured-asteroid-moons.md sections 6 and 7)
=====================================================================================================================================================

`list_info(conn, planet_id)` reads a stored planet's rows (no system is built) and answers the irregular moons it keeps (`generation/irregular_moons.IrregularSet`, saved by GEN.1084): each body by its designation with its
mass, size, makeup and the orbit about the planet (the six elements as saved, the closest and farthest distance, the period the planet's mass gives, and the share of the planet's Hill radius it
reaches), and what was not stored (the bodies below the storage floor, the ones lost, the flag for temporary captures). `moon_info(conn, planet_id, rank)` answers one of them with the heavier and
lighter neighbour. `GET /api/objects/<ref>/irregular-moons` and `.../irregular-moons/<rank>` serve them and the planet's irregular-moon pages show them. Nothing here is stored.

**Designation.** A body caught from the outer disk or from a planet's own zone has no earlier name, so it is `<planet ID>.<rank hex>` (`object_uid.designation`, the same shape as a belt member's), its
rank the stored `moon_rank` plus one (a designation's rank starts at 1; the stored rank starts at 0 for the heaviest). A body that was a stored member of an inner belt keeps the designation it had
there, `<belt ID>.<rank hex>`, and the belt's page for that member (`belt_url`), because the belt still keeps the row (status `moon`).
"""

import math

from planetgen.db import store
from planetgen.generation import irregular_moons
from planetgen.galaxy import object_uid
from planetgen.physics import belt_page as page_data
from planetgen.physics import constants, unit_registry as units
from planetgen.physics import small_body_composition as sbc


class IrregularMoonNotFound(ValueError):
    """No planet has that id, or it keeps no irregular moon of that rank."""


_KM_PER_AU = units.convert(1.0, "au", "km")


def designation_rank(moon):
    """The rank a stored irregular moon's own designation carries: its stored rank plus one (a designation's rank starts at 1)."""
    return moon.rank + 1


def _period_yr(distance_km, mass_kg):
    mu = constants.G * mass_kg
    seconds = 2.0 * math.pi * math.sqrt((distance_km * 1000.0) ** 3 / mu)
    return units.convert(seconds, "s", "yr")


def _entry(planet_id, moon, belts, hill_km):
    """One body's entry, in the units the page shows."""
    q = page_data._q
    a, e = moon.orbit.semi_major_axis_au * _KM_PER_AU, moon.orbit.eccentricity
    belt_id = belts.get(moon.belt_order) if moon.belt_rank is not None else None
    text, belt_designation = object_uid.irregular_designation(planet_id, moon.rank, belt_id, moon.belt_rank), belt_id is not None
    return {
        "rank": designation_rank(moon),
        "designation": text,
        "belt_designation": belt_designation,
        "window": moon.window,
        "retrograde": moon.retrograde,
        "mass": q(moon.mass_kg, "kg"),
        "diameter": q(2.0 * moon.radius_km, "km"),
        "density": q(moon.density_kg_m3, "kg_m3"),
        "taxon": moon.taxon,
        "analogue": sbc.TAXON_DATA[moon.taxon].analogue,
        "composition": {"rock": q(moon.rock, "frac"), "metal": q(moon.metal, "frac"), "water": q(moon.water, "frac"), "ice": q(moon.ice, "frac")},
        "distance": q(a, "km"),
        "eccentricity": e,
        "inclination": q(moon.orbit.inclination_deg, "deg"),
        "hill_fraction": None if not hill_km else a / hill_km,
        "status": {"word": None, "soter_mu": None},
        "population": None,
    }


def _planet_block(planet, planet_id):
    return {"designation": object_uid.format_id(planet_id), "name": planet["name"], "star": planet["star_name"], "mass": page_data._q(planet["mass_kg"], "kg"),
            "hill_radius": None if not planet["hill_radius_km"] else page_data._q(planet["hill_radius_km"], "km")}


def _load(conn, planet_id):
    """`(planet object ID, the planet's row, {belt order: belt object ID}, IrregularSet)` read from the rows alone, with no system built."""
    row = conn.execute("SELECT p.uid, p.name, p.mass_kg, p.hill_radius_km, p.star_id, s.name AS star_name FROM planets p LEFT JOIN stars s ON s.id = p.star_id WHERE p.id = ?", (planet_id,)).fetchone()
    if row is None:
        raise IrregularMoonNotFound(f"no planet has the id {planet_id}")
    summary = conn.execute("SELECT remainder_kg, above_floor, lost, `temporary` FROM irregular_sets WHERE planet_id = ?", (planet_id,)).fetchone()
    rows = conn.execute(f"SELECT {', '.join(store.IRREGULAR_MOON_COLUMNS)} FROM irregular_moons WHERE planet_id = ? ORDER BY moon_rank", (planet_id,)).fetchall()
    held = irregular_moons.IrregularSet.from_dict(store.irregular_set_dict(summary, rows))
    belts = {order: object_uid.from_bytes(found["uid"]) for order, found in enumerate(conn.execute(
        "SELECT uid FROM asteroid_belts WHERE star_id = ? ORDER BY orbital_index, id", (row["star_id"],)).fetchall())}
    return object_uid.from_bytes(row["uid"]), row, belts, held


def list_info(conn, planet_id):
    """
    The irregular moons of the stored planet `planet_id`: `planet`, `moons` (the entries, heaviest first), and `not_stored` (`remainder`, the mass of the bodies below the storage floor; `above_floor`, the
    bodies above it; `lost`, the stored ones given up; `temporary`, whether the planet's zone is wide enough to hold temporary captures).

    Raises:
        IrregularMoonNotFound: No planet has that id.
    """
    uid, planet, belts, held = _load(conn, planet_id)
    return {
        "planet": _planet_block(planet, uid),
        "moons": [_entry(uid, moon, belts, planet["hill_radius_km"]) for moon in held.moons],
        "not_stored": {"remainder": page_data._q(held.remainder_kg, "kg"), "above_floor": held.above_floor, "lost": held.lost, "temporary": held.temporary},
    }


def moon_info(conn, planet_id, rank):
    """
    The page data of the irregular moon of `rank` (the rank its designation writes, from 1) of the stored planet `planet_id`: `moon` (the entry), `orbit` (the six elements, the closest and farthest
    distance and the period), `planet`, and the `previous` and `next` neighbours (`{rank, designation}`, `None` at the ends; the heavier is previous).

    Raises:
        IrregularMoonNotFound: No planet has that id, or it keeps no irregular moon of that rank.
    """
    uid, planet, belts, held = _load(conn, planet_id)
    moons = list(held.moons)
    place = next((i for i, moon in enumerate(moons) if designation_rank(moon) == rank), None)
    if place is None:
        raise IrregularMoonNotFound(f"the planet keeps {len(moons)} irregular moons; none has the rank {rank}")
    q = page_data._q
    moon = moons[place]
    shown = _entry(uid, moon, belts, planet["hill_radius_km"])
    a, e = moon.orbit.semi_major_axis_au * _KM_PER_AU, moon.orbit.eccentricity
    neighbour = lambda i: None if not 0 <= i < len(moons) else {"rank": designation_rank(moons[i]), "designation": _entry(uid, moons[i], belts, planet["hill_radius_km"])["designation"]}
    return {
        "designation": shown["designation"],
        "rank": rank,
        "moon": shown,
        "orbit": {
            "semi_major_axis": q(a, "km"), "eccentricity": e, "inclination": q(moon.orbit.inclination_deg, "deg"),
            "ascending_node": q(moon.orbit.ascending_node_deg, "deg"), "argument_of_periapsis": q(moon.orbit.arg_periapsis_deg, "deg"),
            "mean_anomaly": q(moon.orbit.mean_anomaly_deg, "deg"), "closest_distance": q(a * (1.0 - e), "km"), "farthest_distance": q(a * (1.0 + e), "km"),
            "period": q(_period_yr(a, planet["mass_kg"] + moon.mass_kg), "yr"), "retrograde": moon.retrograde,
            "hill_fraction": shown["hill_fraction"],
        },
        "planet": _planet_block(planet, uid),
        "previous": neighbour(place - 1),
        "next": neighbour(place + 1),
    }
