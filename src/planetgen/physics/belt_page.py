# planetgen/physics/belt_page.py

"""
The data of an asteroid belt's information page (UX.240; docs/design/systems/belt-evolution.md section 12)
===========================================================================================================

Boss asked (4:37 PM EDT 2026-10-30) for a page per belt that shows "the bodies in the asteroid belt that are named and the contents and sizes of the rest of the mass of the asteroid belt, as well as a density factor". Belt members have a
designation, not a name (`galaxy/object_uid.designation`), so the first block is the individually tracked bodies by designation. This module is the pure structure that the page and the API both read; the route and the template are UX.241.
It is a plain dict of strings, numbers and lists (so already JSON), and **every number is a quantity**: `{"value", "unit", "text"}` with `unit` a key of the unit registry (or `count`) and `text` the number written with its unit.

The three blocks, over a `BeltRecord` (the belt's limits and star, its `belt_ledger.Ledger`, and where each stored body is):

1. **`bodies`**: one entry for each stored body of the ledger, heaviest first: its designation, mass, diameter, orbit distance, taxon and composition shares (`belt_member.member` at its own distance, the taxon the caller's uniform number picks), density, and
   status word (`body_status.status` over the zone mass the other stored bodies, the planets and the belt's leftover mass give). The word and Soter's mu are both kept (decision 2 of the design).
2. **`bins`**: the rest of the belt mass, everything the ledger did not store, as bins of diameter that end at the smallest stored body: the decades 1, 10, 100, 1,000 km below it, the first running down to `BIN_LOW_KM`. Each bin has the expected count, the
   mass and the share of the belt's mass, the tier (bodies above `L2_FLOOR_KG` are drawn on demand from the belt's seed, the rest are a total only), and the typical composition: the taxon shares the kernel gives over the belt, weighted by where the belt's
   solids are, and the density of a typical body of that bin. The masses follow the belt's size law and are **scaled together so that they add up to the ledger's `remainder_kg`**; the counts are scaled by the same factor, so the table agrees with the
   stored belt mass whatever the draw put in the stored bodies.
3. **`density_factor`**: the objects per square astronomical unit above each diameter of `DENSITY_DIAMETERS_KM` (1, 10 and 100 km) over the belt's area `pi (upper^2 - lower^2)`, and the mass per square astronomical unit. The count above a diameter is the stored bodies
   above it plus the size law's count between it and the smallest stored body (scaled as in 2).

`totals` gives the belt's mass, the listed bodies' mass, the bins' mass and the difference, which is zero to rounding.

It sits in `physics`, not in `describe` as the item's text had it: the wiring guard counts the whole `describe` package as the running app, so every module this one imports (the ledger, the kernel, the status module and the rows waiting on the
generator) would count as read before anything serves the page. UX.241 imports it from the route.

Stand-ins BP1 to BP4 of `docs/design/research/gap-filling-plan.md` section 6.120. Nothing serves it yet (UX.241).
"""

import dataclasses
import math
from dataclasses import dataclass

from planetgen.galaxy import object_uid
from planetgen.physics import belt_history, belt_ledger, belt_member, body_status, compaction, unit_registry
from planetgen.physics import small_body_composition as sbc
from planetgen.util.checks import finite_domain
from planetgen.util.format import format_number

BIN_EDGES_KM = (1.0, 10.0, 100.0, 1000.0)
"""tuple: The decades where the size table splits, km (presentation). [Stand-in BP1, tier T9, tweak L0.]"""
BIN_LOW_KM = 1.0e-3
"""float: The size the lowest bin runs down to, km (a metre; the mass below it is a few parts in 10^4 of a belt whose size law has a cumulative index under 2). [Stand-in BP1, tier T9, tweak L0.]"""
DENSITY_DIAMETERS_KM = (1.0, 10.0, 100.0)
"""tuple: The diameters the density factor counts objects above, km (Boss's 1, 10 and 100 km). [Stand-in BP2, tier T9, tweak L0.]"""
L2_FLOOR_KG = 1.0e15
"""float: Bodies from this mass up to the stored ones are drawn on demand from the belt's seed (the design's L2 tier, about 10 km). [Stand-in BP3, tier T9, tweak L1.]"""
SHARE_SAMPLES = 24
"""int: The distances across the belt at which the kernel's taxon shares are averaged for the typical composition (numerical). [Stand-in BP4, tier T9, tweak L0.]"""
INTEGRAL_STEPS = 200
"""int: The log-spaced steps of the sum that gives a bin's mass (numerical)."""
MIN_SHARE = 0.005
"""float: A taxon with a share under this is left out of a composition list (presentation). [Stand-in BP4, tier T9, tweak L0.]"""


class BeltPageError(ValueError):
    """An input the page cannot be built from."""


@dataclass(frozen=True)
class Placement:
    """Where a stored body is and the two uniform numbers that pick its taxon and, if light, its porosity (`belt_member.member`)."""

    distance_au: float
    u_taxon: float
    u_porosity: float = 0.5


@dataclass(frozen=True)
class BeltRecord:
    """
    What the page needs of a belt.

    Attributes:
        belt_id (int): The belt's object ID (its body in the system).
        inner_au, outer_au (float): The belt's limits.
        star_mass_solar, luminosity_w, aluminium_26_ratio (float): The star's.
        ledger (belt_ledger.Ledger): The stored bodies and the rest (the bodies are `ledger.masses_kg`, heaviest first).
        placements (tuple): A `Placement` for each stored body, in the ledger's order.
        planets (tuple): `(a_au, mass_kg)` of the system's planets, for the zone of each body (Soter's mu).
    """

    belt_id: int
    inner_au: float
    outer_au: float
    star_mass_solar: float
    luminosity_w: float
    aluminium_26_ratio: float
    ledger: object
    placements: tuple
    planets: tuple = ()


def _q(value, unit):
    """A quantity: the number, the unit's registry key and the text of the number with its unit symbol."""
    if unit == "frac":
        text = f"{format_number(value * 100.0, '.3g')} percent"
    elif unit == "count":
        text = f"{format_number(value, ',.0f') if value >= 100 else format_number(value, '.3g')} bodies"
    else:
        text = f"{format_number(value, '.3g')} {unit_registry.unit(unit).symbol}"
    return {"value": value, "unit": unit, "text": text}


def _diameter_km(mass_kg, density):
    return belt_ledger.diameter_km(mass_kg, density)


def _mass_of_diameter(diameter_km, density):
    return math.pi / 6.0 * (diameter_km * 1000.0) ** 3 * density


def _shares(record):
    """The taxon shares over the belt: the kernel at `SHARE_SAMPLES` distances taken at equal steps of the belt's mass (its solids go as `a^-SURFACE_DENSITY_SLOPE`), averaged."""
    k = 2.0 - belt_history.SURFACE_DENSITY_SLOPE
    total = {}
    for i in range(SHARE_SAMPLES):
        u = (i + 0.5) / SHARE_SAMPLES
        distance = (record.inner_au ** k + u * (record.outer_au ** k - record.inner_au ** k)) ** (1.0 / k)
        weights = belt_member.shares(distance, record.luminosity_w, record.aluminium_26_ratio)
        norm = math.fsum(weights.values())
        for taxon, weight in weights.items():
            total[taxon] = total.get(taxon, 0.0) + weight / norm / SHARE_SAMPLES
    return total


def _composition_list(shares):
    return [{"taxon": t, "analogue": sbc.TAXON_DATA[t].analogue, "share": _q(s, "frac")}
            for t, s in sorted(shares.items(), key=lambda item: -item[1]) if s >= MIN_SHARE]


def _check(record):
    ledger = record.ledger
    if not 0.0 < record.inner_au < record.outer_au:
        raise BeltPageError("the belt needs 0 < inner < outer")
    if len(record.placements) != len(ledger.masses_kg):
        raise BeltPageError("one placement for each stored body")
    if not (record.star_mass_solar > 0.0 and record.luminosity_w > 0.0):
        raise BeltPageError("the star's mass and luminosity must be positive")
    for placement in record.placements:
        if not record.inner_au <= placement.distance_au <= record.outer_au:
            raise BeltPageError("a stored body lies outside its belt")


def _member(record, mass_kg, placement):
    """`belt_member.member`, except that a body over the compaction table (`compaction.TAPER_END_KG`, for the structure solver) takes the table's top density and its own radius, and is marked."""
    clamped = min(mass_kg, compaction.TAPER_END_KG)
    member = belt_member.member(placement.u_taxon, placement.u_porosity, clamped, placement.distance_au, record.luminosity_w, record.aluminium_26_ratio)
    if clamped == mass_kg:
        return member, False
    radius_km = (3.0 * mass_kg / (4.0 * math.pi * member.density_kg_m3)) ** (1.0 / 3.0) / 1.0e3
    return dataclasses.replace(member, mass_kg=mass_kg, radius_km=radius_km), True


def _bodies(record):
    ledger = record.ledger
    pairs = [_member(record, m, p) for m, p in zip(ledger.masses_kg, record.placements)]
    neighbours = tuple(record.planets) + tuple((p.distance_au, m) for m, p in zip(ledger.masses_kg, record.placements))
    belts = ((record.inner_au, record.outer_au, ledger.remainder_kg),)
    out = []
    for rank, ((member, held), placement) in enumerate(zip(pairs, record.placements), start=1):
        others = body_status.zone_mass_kg(placement.distance_au, member.mass_kg, record.star_mass_solar, neighbours, belts)
        ice = member.makeup.water + member.makeup.ice
        state = body_status.status(member.mass_kg, ice, placement.distance_au, record.star_mass_solar, others)
        out.append({
            "rank": rank,
            "designation": object_uid.designation(record.belt_id, rank),
            "mass": _q(member.mass_kg, "kg"),
            "diameter": _q(2.0 * member.radius_km, "km"),
            "distance": _q(placement.distance_au, "au"),
            "density": _q(member.density_kg_m3, "kg_m3"),
            "density_held": held,
            "taxon": member.taxon,
            "analogue": sbc.TAXON_DATA[member.taxon].analogue,
            "differentiated": member.differentiated,
            "composition": {"rock": _q(member.makeup.rock, "frac"), "metal": _q(member.makeup.metal, "frac"),
                            "water": _q(member.makeup.water, "frac"), "ice": _q(member.makeup.ice, "frac")},
            "status": {"word": state.word, "soter_mu": None if state.mu is None or math.isinf(state.mu) else state.mu, "roundness": state.roundness,
                       "margot_pi": state.pi, "agrees": state.pi_agrees},
        })
    return out


def _raw_mass_and_count(law, belt_mass_kg, low_kg, high_kg):
    """The size law's expected count and mass between two masses (a log-spaced sum of `m dN`)."""
    step = (high_kg / low_kg) ** (1.0 / INTEGRAL_STEPS)
    mass, r = 0.0, low_kg
    n_low = belt_ledger.count_above(law, belt_mass_kg, low_kg)
    for _ in range(INTEGRAL_STEPS):
        upper = r * step
        n_high = belt_ledger.count_above(law, belt_mass_kg, upper)
        mass += max(n_low - n_high, 0.0) * math.sqrt(r * upper)
        r, n_low = upper, n_high
    return belt_ledger.count_above(law, belt_mass_kg, low_kg) - belt_ledger.count_above(law, belt_mass_kg, high_kg), mass


def _smallest_stored(ledger, bodies):
    """`(diameter_km, mass_kg)` of the smallest stored body, the diameter as the page lists it (from the body's own density); with none stored, the size law's typical largest body."""
    law = ledger.law
    if ledger.masses_kg:
        return bodies[-1]["diameter"]["value"], ledger.masses_kg[-1]
    mass = belt_ledger.mass_at_count(law, ledger.belt_mass_kg, 1.0)
    return _diameter_km(mass, law.density_kg_m3), mass


def _edges(smallest_km, smallest_kg, density):
    """The bins' diameters: the bottom, the decades whose mass (at the law's density) is under the smallest stored body's, and that body's own diameter."""
    return [BIN_LOW_KM] + [e for e in BIN_EDGES_KM if BIN_LOW_KM < e and _mass_of_diameter(e, density) < smallest_kg] + [smallest_km]


def _raw_bins(ledger, smallest_km, smallest_kg):
    """`(low_km, high_km, count, mass_kg)` of each bin by the size law alone; the diameters between edges are turned into masses at the law's density, the top edge is the smallest stored body's mass."""
    law, density = ledger.law, ledger.law.density_kg_m3
    edges = _edges(smallest_km, smallest_kg, density)
    masses = [_mass_of_diameter(e, density) for e in edges[:-1]] + [smallest_kg]
    return [(low, high) + _raw_mass_and_count(law, ledger.belt_mass_kg, m_low, m_high)
            for low, high, m_low, m_high in zip(edges, edges[1:], masses, masses[1:])]


def _scale(ledger, raw_bins):
    """The factor that makes the bins' masses add up to `remainder_kg`."""
    raw = math.fsum(mass for _, _, _, mass in raw_bins)
    return ledger.remainder_kg / raw if raw > 0.0 else 0.0


def _bins(record, raw_bins, scale, shares):
    ledger = record.ledger
    density = ledger.law.density_kg_m3
    l2_km = _diameter_km(L2_FLOOR_KG, density)
    out = []
    for low, high, count, mass in raw_bins:
        count, mass = count * scale, mass * scale
        middle_kg = _mass_of_diameter(math.sqrt(low * high), density)
        typical = None
        try:
            typical = math.fsum(share * belt_member.typical_density_kg_m3(t, middle_kg) for t, share in shares.items()) / math.fsum(shares.values())
        except compaction.OutOfRange:
            pass
        out.append({
            "from": None if low == BIN_LOW_KM else _q(low, "km"),
            "to": _q(high, "km"),
            "tier": "on demand" if low >= l2_km else "total only",
            "count": _q(count, "count"),
            "mass": _q(mass, "kg"),
            "mass_share": _q(mass / ledger.belt_mass_kg, "frac"),
            "composition": _composition_list(shares),
            "typical_density": None if typical is None else _q(typical, "kg_m3"),
        })
    return out


def _count_above(record, bodies, diameter_km, scale, smallest_km, smallest_kg):
    ledger = record.ledger
    law, density = ledger.law, ledger.law.density_kg_m3
    stored = sum(1 for b in bodies if b["diameter"]["value"] >= diameter_km)
    if diameter_km >= smallest_km:
        return float(stored)
    low, high = _mass_of_diameter(diameter_km, density), smallest_kg
    return stored + scale * max(belt_ledger.count_above(law, ledger.belt_mass_kg, low) - belt_ledger.count_above(law, ledger.belt_mass_kg, high), 0.0)


@finite_domain()
def build(record):
    """
    The page data of a belt (see the module text): `{"belt", "bodies", "bins", "density_factor", "totals"}`.

    Raises:
        BeltPageError: A record that does not hang together.
        compaction.OutOfRange: A stored body over `compaction.TAPER_END_KG` (for the structure solver).
    """
    _check(record)
    ledger = record.ledger
    bodies = _bodies(record)
    smallest, smallest_kg = _smallest_stored(ledger, bodies)
    raw_bins = _raw_bins(ledger, smallest, smallest_kg)
    scale = _scale(ledger, raw_bins)
    shares = _shares(record)
    bins = _bins(record, raw_bins, scale, shares)
    area = math.pi * (record.outer_au ** 2 - record.inner_au ** 2)
    factor = [{"above": _q(d, "km"), "objects": _q(_count_above(record, bodies, d, scale, smallest, smallest_kg) / area, "au-2")} for d in DENSITY_DIAMETERS_KM]
    listed = math.fsum(m for m in ledger.masses_kg)
    in_bins = math.fsum(b["mass"]["value"] for b in bins)
    return {
        "belt": {"designation": object_uid.format_id(record.belt_id), "inner": _q(record.inner_au, "au"), "outer": _q(record.outer_au, "au"),
                 "mass": _q(ledger.belt_mass_kg, "kg"), "composition": _composition_list(shares)},
        "bodies": bodies,
        "bins": bins,
        "density_factor": {"by_diameter": factor, "mass": _q(ledger.belt_mass_kg / area, "kg_au2")},
        "totals": {"belt_mass": _q(ledger.belt_mass_kg, "kg"), "listed": _q(listed, "kg"), "bins": _q(in_bins, "kg"), "difference": _q(ledger.belt_mass_kg - listed - in_bins, "kg"),
                   "scale": scale},
    }
