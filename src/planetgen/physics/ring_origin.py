# planetgen/physics/ring_origin.py

"""
Where a planet's ring comes from (GEN.1054; docs/design/systems/weird-worlds.md 9.2, R1, R2, R9)
==================================================================================================

`rings.draw_ring` gives every giant a ring and a terrestrial planet one in a hundred, whatever happened to the planet. A ring is a result of history, and this is the
history: a ring exists when something came apart inside the planet's Roche limit and has not yet spread away, and a giant with nothing in its past has none.

* **R1, a moon that fell** ([R0637], recalled; [R1927], recalled; [R0414] for Phobos, recalled). A moon the tide drew below the Roche limit
  (`moon_orbits.Fate`, `time_s`) was torn apart: the share `RING_RETAINED_FRACTION` of its mass stays in orbit as a ring and the rest rains onto the planet. That ring spreads
  and feeds new moons beyond the limit, which is the ring-moon cycle: its mass falls by `1/e` every `RING_LIFETIME_YR`.
* **R2, a comet that passed inside the limit** ([R1140], recalled). A body that passes inside the Roche limit of a giant is torn apart in the same way and
  `COMET_CAPTURE_FRACTION` of it is left in orbit. The passages come at a rate `COMET_PASSAGES_PER_GYR` for Saturn's situation, scaled to the planet by the cross-section of its Roche
  sphere with gravitational focusing, `pi q^2 (1 + v_esc(q)^2 / v_inf^2)`, with `v_inf` a share (`COMET_ENCOUNTER_FRACTION`) of the circular speed at the planet's distance. Each is a
  body of `COMET_MASS_KG_RANGE` drawn log-uniform, at a time drawn uniformly over the planet's age.
* **R9, close to the star.** A ring inside `RING_MIN_DISTANCE_AU` is gone at once ([R3957], recalled): no ring.

The ring that results is the sum of the rings of every event, each spread for the time since it came to be, capped at `rings.MAX_MASS_OVER_PLANET` of the planet; it exists above
`MIN_RING_MASS_KG`, its origin is that of the event that holds most of it, and its age that event's. The Roche limit is `rings.roche_limit_km` for the ring's own material.
The module draws only through `draw_comet_passages`, from the unit's stream; the rest is arithmetic on its arguments. `rings.apply` reads it (GEN.1057): the moons that fell come from `moon_history`
(the planet's `fallen_moons`), the comets from `draw_comet_passages`, and the chances of the old `rings.draw_ring` are gone. Stand-ins are RO1 to RO8 of `docs/design/research/gap-filling-plan.md` section 6.70; the draws of `rings` that this
replaces (BC20 to BC23) are T4, L3 until the ring age and disruption rates of the Cassini end of mission replace them.
"""

import math
from dataclasses import dataclass

from planetgen.physics import constants, rings
from planetgen.util import draw
from planetgen.util.random import log_uniform

RING_RETAINED_FRACTION = 0.3
"""float: The share of a moon torn apart inside the Roche limit that stays in orbit as a ring (the rest falls on the planet); [R0414] keep most of Phobos in orbit and [R0637] only the
ice mantle of a differentiated moon, so the share is argued between the two. Source tier T4, tweak level L3."""
RING_LIFETIME_YR = 5.0e8
"""float: The time in which a ring's mass falls by 1/e as it spreads and feeds moons beyond the Roche limit (the ring-moon cycle of [R1927] runs over 1e8 to 1e9 years; the age of Saturn's
rings is argued between 1e8 and the age of the planet). Source tier T4, tweak level L3."""
COMET_PASSAGES_PER_GYR = 0.02
"""float: How many ring-making comets (`COMET_MASS_KG_RANGE`) pass inside the Roche limit of a body in Saturn's situation per Gyr; no survey measures it (about one in 50 Gyr would leave Saturn with a comet's ring
only once in the age of the Solar System). Source tier T4, tweak level L3."""
COMET_MASS_KG_RANGE = (1.0e14, 1.0e20)
"""tuple: The masses, kg, of the comets that make rings, drawn log-uniform (a 5 km comet to a 300 km Centaur). Source tier T6, tweak level L1."""
COMET_CAPTURE_FRACTION = 0.1
"""float: The share of a torn comet left in orbit ([R1140] find a few percent to tens of percent captured). Source tier T4, tweak level L3."""
COMET_ENCOUNTER_FRACTION = 0.2
"""float: A comet's speed at infinity as a share of the circular speed at the planet's distance (the Jupiter-family comets meet Jupiter at about 5 km/s of its 13). Source tier T6, tweak level L1."""
RING_MIN_DISTANCE_AU = 0.1
"""float: A planet closer to its star than this holds no ring ([R3957], recalled). Source tier T2, tweak level L0."""
MIN_RING_MASS_KG = 1.0e7
"""float: The least ring that counts (the dust of Jupiter's ring is 1e7 to 1e9 kg). Source tier T6, tweak level L0."""

SATURN = dict(mass_kg=5.683e26, radius_km=58232.0, distance_au=9.58, star_mass_kg=1.989e30)
"""dict: The reference situation of `COMET_PASSAGES_PER_GYR`."""


@dataclass(frozen=True)
class Planet:
    """A planet as the rule reads it: mass, radius, distance from its star and the star's mass."""
    mass_kg: float
    radius_km: float
    distance_au: float
    star_mass_kg: float


@dataclass(frozen=True)
class FallenMoon:
    """A moon that fell inside the Roche limit: its mass and when, years after the planet formed."""
    mass_kg: float
    time_yr: float


@dataclass(frozen=True)
class CometPassage:
    """A comet that passed inside the Roche limit: its mass and when, years after the planet formed."""
    mass_kg: float
    time_yr: float


@dataclass(frozen=True)
class RingOrigin:
    """The result: whether a ring exists now, its origin (`rings.ORIGINS`) and mass, the years since its main event, and the mass each source left (kg, now)."""
    exists: bool
    origin: str
    mass_kg: float
    age_yr: float
    from_moons_kg: float
    from_comets_kg: float


def _density_kg_m3(planet):
    return planet.mass_kg / (4.0 / 3.0 * math.pi * (planet.radius_km * 1000.0) ** 3)


def roche_limit_km(planet, cold=True):
    """The Roche limit, km, of the planet for the material of a ring beyond the snow line (`cold`, with the ice of a ring) or inside it."""
    ice_fraction = rings.RANGES["ice_fraction"]["broad"][0] if cold else 0.0
    return rings.roche_limit_km(planet.radius_km, _density_kg_m3(planet), rings.mixture_density_kg_m3(ice_fraction))


def cross_section_m2(planet, cold=True):
    """The cross-section, m2, for a comet to pass inside the Roche limit, with gravitational focusing: `pi q^2 (1 + v_esc(q)^2 / v_inf^2)`."""
    q = roche_limit_km(planet, cold) * 1000.0
    v_inf = COMET_ENCOUNTER_FRACTION * math.sqrt(constants.G * planet.star_mass_kg / (planet.distance_au * constants.AU_TO_M))
    v_escape_squared = 2.0 * constants.G * planet.mass_kg / q
    return math.pi * q * q * (1.0 + v_escape_squared / v_inf ** 2)


def comet_passage_rate_per_gyr(planet, cold=True):
    """The rate, per Gyr, of ring-making comets passing inside the Roche limit: `COMET_PASSAGES_PER_GYR` scaled from Saturn by the cross-section; 0 inside `RING_MIN_DISTANCE_AU`."""
    if planet.distance_au < RING_MIN_DISTANCE_AU:
        return 0.0
    return COMET_PASSAGES_PER_GYR * cross_section_m2(planet, cold) / cross_section_m2(Planet(**SATURN), True)


def draw_comet_passages(planet, age_yr, cold=True):
    """The comets that passed inside the Roche limit over `age_yr`: their number is Poisson at the rate above, each with a time uniform over the age and a log-uniform mass. Draws from the unit's stream."""
    expected = comet_passage_rate_per_gyr(planet, cold) * age_yr / 1.0e9
    count, product, limit = 0, draw.random(), math.exp(-expected)
    while product > limit:
        count += 1
        product *= draw.random()
    return tuple(CometPassage(log_uniform(*COMET_MASS_KG_RANGE), draw.uniform(0.0, age_yr)) for _ in range(count))


def spread(mass_kg, since_yr):
    """A ring's mass `since_yr` after it formed: it falls by 1/e every `RING_LIFETIME_YR`."""
    return mass_kg * math.exp(-max(since_yr, 0.0) / RING_LIFETIME_YR)


def ring_origin(planet, age_yr, fallen_moons=(), comet_passages=()):
    """
    The ring of `planet` at `age_yr` from the moons that fell (`FallenMoon`s, from `moon_orbits`) and the comets that passed inside the Roche limit (`CometPassage`s, from `draw_comet_passages`):
    the sum of each event's ring spread for the time since, capped at `rings.MAX_MASS_OVER_PLANET` of the planet, none inside `RING_MIN_DISTANCE_AU` or below `MIN_RING_MASS_KG`.
    """
    if planet.distance_au < RING_MIN_DISTANCE_AU:
        return RingOrigin(False, None, 0.0, 0.0, 0.0, 0.0)
    events = [("disrupted_moon", spread(RING_RETAINED_FRACTION * m.mass_kg, age_yr - m.time_yr), age_yr - m.time_yr) for m in fallen_moons]
    events += [("captured_body", spread(COMET_CAPTURE_FRACTION * c.mass_kg, age_yr - c.time_yr), age_yr - c.time_yr) for c in comet_passages]
    from_moons = sum(mass for origin, mass, _ in events if origin == "disrupted_moon")
    from_comets = sum(mass for origin, mass, _ in events if origin == "captured_body")
    total = min(from_moons + from_comets, rings.MAX_MASS_OVER_PLANET * planet.mass_kg)
    if total < MIN_RING_MASS_KG:
        return RingOrigin(False, None, 0.0, 0.0, from_moons, from_comets)
    origin, _, since = max(events, key=lambda event: event[1])
    return RingOrigin(True, origin, total, since, from_moons, from_comets)
