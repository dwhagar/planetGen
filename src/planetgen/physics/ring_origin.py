# planetgen/physics/ring_origin.py

"""
Where a planet's ring comes from (GEN.1054; docs/design/systems/weird-worlds.md 9.2, R1, R2, R9)
==================================================================================================

`rings.draw_ring` gives every giant a ring and a terrestrial planet one in a hundred, whatever happened to the planet. A ring is a result of history, and this is the
history: a ring exists when something came apart inside the planet's Roche limit and has not yet spread away, and a giant with nothing in its past has none.

* **R1, a moon that fell** ([R0637], recalled; [R1927], recalled; [R0414] for Phobos, recalled). A moon the tide drew below the Roche limit
  (`moon_orbits.Fate`, `time_s`) was torn apart: the share `RING_RETAINED_FRACTION` of its mass stays in orbit as a ring and the rest rains onto the planet. That ring spreads
  and feeds new moons beyond the limit, which is the ring-moon cycle, and the meteoroids that fall on it wear it away (GEN.1171, below).
* **R2, a comet that passed inside the limit** ([R1140], recalled). A body that passes inside the Roche limit of a giant is torn apart in the same way and
  `COMET_CAPTURE_FRACTION` of it is left in orbit. The passages come at a rate `COMET_PASSAGES_PER_GYR` for Saturn's situation, scaled to the planet by the cross-section of its Roche
  sphere with gravitational focusing, `pi q^2 (1 + v_esc(q)^2 / v_inf^2)`, with `v_inf` a share (`COMET_ENCOUNTER_FRACTION`) of the circular speed at the planet's distance. Each is a
  body of `COMET_MASS_KG_RANGE` drawn log-uniform, at a time drawn uniformly over the planet's age.
* **R9, close to the star.** A ring inside `RING_MIN_DISTANCE_AU` is gone at once ([R3957], recalled): no ring.

**How long a ring lasts (GEN.1171).** Two losses act on a ring together, and the heavier the ring the faster the first and the lighter the ring the faster the second.
*Viscous spreading* ([R0920], recalled): a dense ring's viscosity grows with its surface density `Sigma`, so the time `r^2 / nu` falls as the ring gets denser, `tau = RING_SPREAD_TIME_YR (r / r_S)^2 (Sigma_S / Sigma)`
against Saturn's rings (`r` the Roche edge, `Sigma_S` the surface density of `SATURN_RING_MASS_KG` over Saturn's annulus), and the mass that spreads past the Roche limit is gone (it feeds moons or falls on the planet).
*Meteoroid erosion* (Cuzzi and Estrada 1998, recalled): the mass flux on the ring, `METEOROID_FLUX_KG_M2_S` at Saturn's rings, scales with the distance from the star as `(9.58 AU / a)^METEOROID_FLUX_EXPONENT` and with the gravitational focusing at the ring
(the same encounter speed as the comets), and takes away `EROSION_YIELD` times the impactor's mass over the annulus area. The mass obeys `dM/dt = -k M^2 - B`, which has the closed form used by `spread`, so a ring of 1e13 kg lives some 1e5 years and Saturn's of 1.5e19 kg most of the age of the
Solar System. The ring that results is the sum of the rings of every event, each spread for the time since it came to be, capped at `rings.MAX_MASS_OVER_PLANET` of the planet; it exists above
`MIN_RING_MASS_KG`, its origin is that of the event that holds most of it, and its age that event's. The Roche limit is `rings.roche_limit_km` for the ring's own material.
The module draws only through `draw_comet_passages`, from the unit's stream; the rest is arithmetic on its arguments. `rings.apply` reads it (GEN.1057): the moons that fell come from `moon_history`
(the planet's `fallen_moons`), the comets from `draw_comet_passages`, and the chances of the old `rings.draw_ring` are gone. Stand-ins are RO1 to RO8 of `docs/design/research/gap-filling-plan.md` section 6.70; the draws of `rings` that this
replaces (BC20 to BC23) are T4, L3 until the ring age and disruption rates of the Cassini end of mission replace them.
"""

import functools
import math
from dataclasses import dataclass

from planetgen.physics import constants, rings
from planetgen.util import draw
from planetgen.util.random import log_uniform

RING_RETAINED_FRACTION = 0.3
"""float: The share of a moon torn apart inside the Roche limit that stays in orbit as a ring (the rest falls on the planet); [R0414] keep most of Phobos in orbit and [R0637] only the
ice mantle of a differentiated moon, so the share is argued between the two. Source tier T4, tweak level L3."""
RING_SPREAD_TIME_YR = 2.0e9
"""float: The time in which Saturn's rings (`SATURN_RING_MASS_KG` over their annulus) would lose a share 1/e of their mass to viscous spreading; the spreading time of the main rings is argued from 5e8 to 1e10 years
([R0920] for the spreading and re-accretion; the age of Saturn's rings between 1e8 years and the planet's). Source tier T4, tweak level L3."""
SATURN_RING_MASS_KG = 1.54e19
"""float: The mass of Saturn's rings, kg ([R2077], from the Cassini Grand Finale gravity, recalled). Source tier T2, tweak level L0."""
RING_INNER_RADII = 1.24
"""float: The inner edge of a ring in planet radii (the C ring of Saturn begins at 1.24; the D ring at 1.11). Source tier T3, tweak level L0."""
METEOROID_FLUX_KG_M2_S = 1.0e-16
"""float: The meteoroid mass that falls on Saturn's rings, kg per square metre per second, with the focusing of Saturn's gravity (Cuzzi and Estrada 1998 give 1e-17 g per square cm per second, recalled; the Cassini dust analyser
measures the same order). Source tier T3, tweak level L0."""
METEOROID_FLUX_EXPONENT = 1.5
"""float: How the interplanetary meteoroid flux falls with the distance from the star, `(a_S / a)^exponent` (the dust flux of the outer Solar System falls faster than the first power of the distance, Poppe 2016, recalled). Source tier T6, tweak level L1."""
EROSION_YIELD = 1.0
"""float: The ring mass an impactor takes away, in impactor masses (the ejecta of a hypervelocity impact on ice are tens to thousands of impactor masses, most of which falls back; 1 is the least). Source tier T4, tweak level L3."""
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


def annulus_m(planet, cold=True):
    """The inner and outer edges, m, of the planet's ring: `RING_INNER_RADII` planet radii out to the Roche limit (at least 10 percent wider than the inner edge)."""
    inner = RING_INNER_RADII * planet.radius_km * 1000.0
    return inner, max(roche_limit_km(planet, cold) * 1000.0, 1.1 * inner)


def annulus_area_m2(planet, cold=True):
    """The area, m2, of the ring annulus."""
    inner, outer = annulus_m(planet, cold)
    return math.pi * (outer * outer - inner * inner)


def meteoroid_flux_kg_m2_s(planet, cold=True):
    """
    The meteoroid mass that falls on the planet's ring, kg per m2 per second: `METEOROID_FLUX_KG_M2_S` scaled from Saturn's distance by `METEOROID_FLUX_EXPONENT` and
    by the gravitational focusing `1 + v_esc^2 / v_inf^2` at the middle of the annulus against Saturn's, the encounter speed a share `COMET_ENCOUNTER_FRACTION` of the circular speed.
    """
    def focusing(body):
        inner, outer = annulus_m(body, cold)
        v_inf = COMET_ENCOUNTER_FRACTION * math.sqrt(constants.G * body.star_mass_kg / (body.distance_au * constants.AU_TO_M))
        return 1.0 + 2.0 * constants.G * body.mass_kg / (0.5 * (inner + outer)) / v_inf ** 2
    saturn = Planet(**SATURN)
    return (METEOROID_FLUX_KG_M2_S * (saturn.distance_au / planet.distance_au) ** METEOROID_FLUX_EXPONENT
            * focusing(planet) / focusing(saturn))


@functools.lru_cache(maxsize=4096)
def _loss_terms(planet, cold):
    """`(k, b)` of `dM/dt = -k M^2 - b`: the spreading coefficient (per kg per year) and the erosion rate (kg per year) of the planet's ring."""
    saturn = Planet(**SATURN)
    sigma_saturn = SATURN_RING_MASS_KG / annulus_area_m2(saturn, True)
    r_ratio_squared = (annulus_m(planet, cold)[1] / annulus_m(saturn, True)[1]) ** 2
    k = 1.0 / (RING_SPREAD_TIME_YR * r_ratio_squared * sigma_saturn * annulus_area_m2(planet, cold))
    b = EROSION_YIELD * meteoroid_flux_kg_m2_s(planet, cold) * annulus_area_m2(planet, cold) * constants.SECONDS_PER_YEAR
    return k, b


def spread(mass_kg, since_yr, planet, cold=True):
    """
    A ring's mass `since_yr` after it formed with `mass_kg` around `planet`: `dM/dt = -k M^2 - b` (viscous spreading, which quickens with the ring's density, and meteoroid erosion, which
    does not), solved as `M = sqrt(b / k) tan(atan(M0 sqrt(k / b)) - sqrt(k b) t)`; 0 once the angle reaches 0.
    """
    since = max(since_yr, 0.0)
    if mass_kg <= 0.0 or since == 0.0:
        return max(mass_kg, 0.0)
    k, b = _loss_terms(planet, cold)
    scale = math.sqrt(b / k)
    angle = math.atan(mass_kg / scale) - math.sqrt(k * b) * since
    return scale * math.tan(angle) if angle > 0.0 else 0.0


def ring_origin(planet, age_yr, fallen_moons=(), comet_passages=()):
    """
    The ring of `planet` at `age_yr` from the moons that fell (`FallenMoon`s, from `moon_orbits`) and the comets that passed inside the Roche limit (`CometPassage`s, from `draw_comet_passages`):
    the sum of each event's ring spread for the time since, capped at `rings.MAX_MASS_OVER_PLANET` of the planet, none inside `RING_MIN_DISTANCE_AU` or below `MIN_RING_MASS_KG`.
    """
    if planet.distance_au < RING_MIN_DISTANCE_AU:
        return RingOrigin(False, None, 0.0, 0.0, 0.0, 0.0)
    events = [("disrupted_moon", spread(RING_RETAINED_FRACTION * m.mass_kg, age_yr - m.time_yr, planet), age_yr - m.time_yr) for m in fallen_moons]
    events += [("captured_body", spread(COMET_CAPTURE_FRACTION * c.mass_kg, age_yr - c.time_yr, planet), age_yr - c.time_yr) for c in comet_passages]
    from_moons = sum(mass for origin, mass, _ in events if origin == "disrupted_moon")
    from_comets = sum(mass for origin, mass, _ in events if origin == "captured_body")
    total = min(from_moons + from_comets, rings.MAX_MASS_OVER_PLANET * planet.mass_kg)
    if total < MIN_RING_MASS_KG:
        return RingOrigin(False, None, 0.0, 0.0, from_moons, from_comets)
    origin, _, since = max(events, key=lambda event: event[1])
    return RingOrigin(True, origin, total, since, from_moons, from_comets)
