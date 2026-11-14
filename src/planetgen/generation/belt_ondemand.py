# planetgen/generation/belt_ondemand.py

"""
Every body of a belt at the Class D size or larger, on demand (GEN.1674; Boss, 10:07 PM EST on 2026-11-12: every body in the Class D range or larger gets a full planet class, and, 12:30 PM EST on
2026-11-13, "on demand, from the one generator")
==================================================================================================================================================================================

A belt stores at most `belt_ledger.LEDGER_MAX_BODIES` (300) bodies, but its size law holds many more over `belt_ledger.L1_MASS_KG` (1e19 kg, about 200 km across, the floor the ledger already calls "its own
object"): a median outer disk holds about 230,000 (reports/gen1674-class-d-bodies-slice0-2026-11-13.md). The planet pass takes 0.05 to 0.15 s for any body, so classing all of them at generation would take hours a
belt and could not be stored. A body is therefore never stored for its class: its class is made when someone asks for the body, by the same planet generator that makes every planet (`belt_planets.planet_of`,
one path, no second rule).

* **The bodies.** Rank 1 to `stored_count` are the ledger's (`AsteroidBelt.members`; the heaviest ten and any over the oligarch mass are planets already, `belt_planets`). Rank `stored_count + 1` onward are the
  tail of the same size law: the gap from the smallest stored body (`belt_ledger.count_above`) plus the exponential gaps of a stream of the belt's own (the ledger's own rule, `belt_ledger.ledger_from`), so a
  rank's mass falls as its rank grows and a body is the same whenever it is asked for. The tail ends where the mass falls under `L1_MASS_KG`; `individual_count` is the number the law expects over it.
* **Its orbit and makeup** come from a stream keyed on the belt and the rank (`belt_roster.entry_at`, so a body never depends on another), made into a belt member the way a stored one is
  (`belt_page.drawn_body`, `belt.member_row`), and judged for its status against the stored bodies and the rest of the belt.
* **Its class** is the planet the generator makes of that member (`planet_of_rank`): the one stream keyed on the belt and the rank that GEN.1510's promoted bodies use.

A belt read back from a store holds neither the giant's distance nor the stretches the giants cleared (they were not stored), so a tail body of a stored outer disk is drawn cold (`belt_roster`: giant `None`).
The tail does not take the companions' forcing of the stored members (`AsteroidBelt._force_by_companions`, GEN.1492), which works on the stored rows alone. Nothing here is stored, so no output changes.
"""

from planetgen.generation import belt as belt_module, belt_planets
from planetgen.physics import belt_ledger, belt_page, belt_populations, belt_roster, bulk_draw
from planetgen.physics import unit_registry as units
from planetgen.util import draw

FLOOR_KG = belt_ledger.L1_MASS_KG
"""float: The mass from which a body is of the Class D size or larger: the ledger's own object floor."""


class OnDemandError(ValueError):
    """A body that the belt does not hold."""


def settled_of(belt, star):
    """The `belt_ledger.BeltLedger` the belt's stored members were drawn from (`belt.present_settled` over the mass the giants left), or `None` for a belt whose solids the giants took."""
    left = None if belt.origin is None else belt.primordial_mass_earth - belt.cleared_mass_earth - belt.scattered_mass_earth
    return belt_module.present_settled(star, belt.distance, belt.lower_limit, belt.upper_limit, left)


def stored_count(settled):
    """How many bodies the ledger stores (the ranks 1 to this are the stored ones, whether they are members or planets now)."""
    return 0 if settled is None else len(settled.ledger.masses_kg)


def individual_count(settled):
    """How many bodies of the Class D size or larger the belt holds: the size law's count over `FLOOR_KG`, never fewer than the stored ones over it."""
    if settled is None:
        return 0
    ledger = settled.ledger
    expected = belt_ledger.count_above(ledger.law, ledger.belt_mass_kg, FLOOR_KG)
    return max(int(expected), sum(1 for mass in ledger.masses_kg if mass >= FLOOR_KG))


def tail_mass_kg(belt, star, settled, rank):
    """
    The mass of the body of `rank` past the stored ones, from the gap of the smallest stored body and one exponential for each rank since (`belt_ledger.ledger_from`'s rule), or `None` when it is
    under `FLOOR_KG`. A pure function of the belt: the stream is keyed on it, and the ranks are in order of falling mass.
    """
    ledger = settled.ledger
    first = len(ledger.masses_kg)
    if rank <= first:
        raise OnDemandError(f"rank {rank} is a stored body, not a tail body")
    law, belt_mass = ledger.law, ledger.belt_mass_kg
    gap = belt_ledger.count_above(law, belt_mass, ledger.masses_kg[-1])
    stream = draw.Stream(f"belt-tail:{belt.distance!r}:{belt.lower_limit!r}:{star.mass!r}")
    for _ in range(rank - first):
        gap += belt_ledger.next_exponential(stream)
    mass = belt_ledger.mass_at_count(law, belt_mass, gap)
    return mass if mass >= FLOOR_KG else None


def _record(belt, star, settled, system_config):
    return belt_page.BeltRecord(None, belt.lower_limit, belt.upper_limit, units.convert(star.mass, "kg", "Msun"), belt_module.lit_luminosity_w(star, belt.lower_limit, belt.upper_limit),
                                belt_module.star_al26_ratio(system_config, star), settled.ledger, ())


def tail_member(belt, star, system_config, rank):
    """The member row (`AsteroidBelt.members`'s shape) of the tail body of `rank`, or `None` for a rank past where the law's bodies fall under `FLOOR_KG`."""
    settled = settled_of(belt, star)
    if settled is None:
        return None
    mass = tail_mass_kg(belt, star, settled, rank)
    if mass is None:
        return None
    populations = None if belt.scatter_fates is None else belt_populations.member_weights(belt.scatter_fates)
    stream = draw.Stream(f"belt-body:{belt.distance!r}:{belt.lower_limit!r}:{star.mass!r}:{rank}")
    entry = belt_roster.entry_at(stream, belt.lower_limit, belt.upper_limit, star.mass, star.snow_line_au, bulk_draw.metal_ratio_to_solar(star), populations, belt.scatter_giant_au,
                                 belt.cleared_spans)
    ledger = settled.ledger
    neighbours = tuple((member["distance_au"], member["mass_kg"]) for member in belt.members)
    belts = ((belt.lower_limit, belt.upper_limit, ledger.remainder_kg),)
    body = belt_page.drawn_body(_record(belt, star, settled, system_config), rank, mass, entry.placement, neighbours, belts)
    return belt_module.member_row(body, entry)


def member_of_rank(belt, star, system_config, rank):
    """The member row of the body of `rank`, stored or drawn, or `None` for a body the belt does not hold (a planet already, or past the tail)."""
    stored = next((member for member in belt.members if member["rank"] == rank), None)
    if stored is not None:
        return stored
    settled = settled_of(belt, star)
    if settled is None or rank <= stored_count(settled):
        return None
    return tail_member(belt, star, system_config, rank)


def planet_of_rank(system, star, belt, order, rank):
    """
    The full planet (`Planet`, with its class, air, interior, surface and the rest) of the belt body of `rank`: the member made into a planet by the one planet generator, from the stream keyed on
    the belt and the rank (`belt_planets.planet_of`). `order` is the belt's place among its star's belts (`irregular_moons.belts_of`). Raises `OnDemandError` for a body the belt does not hold.
    """
    member = member_of_rank(belt, star, system.system_config, rank)
    if member is None:
        raise OnDemandError(f"the belt holds no body of rank {rank} at or over {FLOOR_KG:.0e} kg")
    return belt_planets.planet_of(system, star, belt, order, member)
