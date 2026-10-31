# planetgen/physics/irregular_population.py

"""
How many captured small bodies a planet holds, and which are its own objects (GEN.1077; docs/design/systems/captured-asteroid-moons.md sections 4.1, 4.2 and 6)
================================================================================================================================================================

A planet does not pick up the bodies of a belt that happens to lie beside it: a body that comes into the Hill sphere from outside leaves again unless something takes its energy away. The irregular moons of the giants came
in **capture windows**, short stretches of the planet's life in which something could (Nesvorny, Vokrouhlicky and Deienno 2014, arXiv:1401.0253, READ; `docs/design/research/captured-moon-sources.md`). Boss decided (4:34 PM EDT 2026-10-30)
that the generator keeps the windows and has no continuous pick-up. This module gives the count that a planet's windows put in orbit about it, how the collisions of the next Gyr grind it down, and the draw of the bodies that
are stored as objects. It does not place them (`irregular_orbits`, `irregular_zone`) or give them a composition (GEN.1080).

1. **The source.** Every window draws from a reservoir of small bodies: a planetesimal disc with a piecewise size law, the cumulative number above a radius `R` falling as `R^-1.8` below `R = 70 km` and `R^-5` above it, `3.5e7` bodies above `R = 40 km`
   in a disc of 30 Earth masses (`source_count_above`; NDD14 section 3.3, NVM07 section 3.5 and Bottke et al. 2010, rows CM3). A heavier or lighter disc scales every count (`Pool.scale`). Until the belt ledger (GEN.1061) hands the generator its own
   reservoirs, this disc is the only source, and so is a disc the generator does not have: a system without a planetesimal disc at the window passes `reservoir_mass_earth` 0 and gets no captures.
2. **The encounter window** is the planetary instability: a giant passes close to another giant while the disc is there, and three-body captures put a share of the disc's bodies about it. NDD14 count, for Jupiter in three histories, 1.3e-10 to 2.6e-10
   captures **per disc particle per encounter** (`ENCOUNTER_YIELD`, row CM1), times 49 to 280 encounters, which is 1.3e-8 to 3.6e-8 per particle (CM2). The caller gives the planet's own encounter count (`Pool.encounters`; the architecture records it,
   GEN.1082). The yield is for a target of the size of Jupiter's Hill sphere and falls with the Hill volume for a smaller one (`geometry_factor`, `min(1, (r_H / 0.35 AU)^3)`): an Earth at 1 AU gets 2e-5 of it (the design's figure),
   and the ice giants, whose Hill spheres are larger, get the same as Jupiter (NDD14's per-particle probabilities for them, 1e-8 to 3e-8, are the same order; the Hill volume would give Neptune ten times Jupiter's).
3. **The gas window** is the last stage of a giant's formation, while its envelope was extended and the disc still there: gas drag slows a body that meets about its own mass of gas, and the envelope then collapses (Pollack et al. 1979, abstract only;
   second-hand in Jewitt and Haghighipour 2007 section 6.2, CM18). It needs a gas envelope (`gas_giant`) and a disc that lasted `GAS_MIN_DISC_LIFETIME_MYR`, and takes only bodies of `GAS_SIZE_WINDOW_KM` (larger ones barely feel the drag, smaller ones
   spiral in). No paper gives its strength against the encounter channel: `GAS_WINDOW_STRENGTH` 0.2 of what the reference instability gives is a guess (CM14, L3), so a giant that had no instability still holds a fraction of the Solar System's number
   and an ice giant none, which is the design's rule.
4. **The exchange window** is a planetesimal binary flying past: one component keeps the energy, the other stays (Vokrouhlicky, Nesvorny and Levison 2008; CM11). It is the only route found for a massive moon such as Triton, far too weak for small
   irregulars, and none at a gas giant (`exchange_chance`, `draw_exchange`). The chance is `EXCHANGE_CHANCE` for a disc in which every planetesimal is a binary, times the binary fraction of the reservoir.
5. **The temporary window** is a flag, not an object: a body that passes through the Hill sphere slowly stays for a year or so, and the Earth always holds one a metre across (Granvik, Vaubaillon and Jedicke 2012, CM13).

**What the windows leave is counted before the stable zone is applied.** NDD14 count the captures that stay bound to the end of their integrations, and the stable limit of the zone (`irregular_zone`) is what makes them stay. Applying the
zone's survival to the counts a second time would take Jupiter from 24 bodies above 5 km to under one (the draw of `irregular_orbits` puts 97 percent of Jupiter's orbits outside the limits): the orbits of the bodies this module returns are
drawn from the stable part of the orbit distribution, never dropped.

**Grinding.** The population then loses bodies to collisions with itself and with the disc: `N(t) = N0 / (1 + t / t_c)` for bodies of one size (Wyatt et al. 2007, the law `belt_mass` uses), with `t_c` the time in which a body of radius `R` is destroyed by
the projectiles above `DISRUPTION_SIZE_RATIO` times its radius, `t_c = 1 / (P_i (R (1 + ratio))^2 N0(> ratio R))`, and `P_i` the intrinsic collision probability of the zone, per square kilometre per year: the relative speed (`sqrt(2)` times the circular
speed at the zone's radius) over the zone's volume (`intrinsic_probability`). It comes out 1.4 to 6 times below the 3.8e-15 to 1.6e-14 that Bottke et al. 2010 compute for the irregulars (CM9), who follow the orbits' eccentricity and
inclination (which concentrate the bodies); no factor is added to match them, so the grinding is a lower bound. The projectile count is the initial one, which overstates the loss of the larger bodies, and the two errors act the other way round.
The saturation the design expected (the four giants ending with the same population whatever they captured) is the late-time limit `N -> 1 / (k t)` of this law; at 4.5 Gyr the loss is a factor 1.2 to 2.4 at 5 km and the giants are only partly
saturated, so the counts above 5 km come out 24, 31, 16 and 33 against the observed 8, 11, 9 and 8 (CM6; the survey of the ice giants is complete only above 7 and 16 km), and above 10.6 km 5, 7, 4 and 8 against 8, 5, 6 and 8. A planet of
another mass or distance is not tuned to them.

**Which are objects.** Boss's rule (4:34 PM EDT 2026-10-30): every irregular above `STORE_FLOOR_KG` (1e16 kg, 10.6 km across 2000 kg/m^3) is stored; if those hold less than `belt_ledger.LEDGER_COVERAGE` of the population's mass, the next heaviest are
added, up to `STORE_MAX_BODIES` (CM15). The rest is the L2 tier (reproducible from the planet's seed by the same draw) and the L3 tier (`remainder_kg`, a total): the same rule as the belt ledger, so one function serves both when the belt storage
merges. The heaviest bodies are the order statistics of a Poisson process, as in `belt_ledger`: rank `k` has `N(> R) = G_k`, the sum of `k` unit exponentials, and the same stream gives the same bodies in the same order whatever the planet.
The Poisson count at the top means a Jupiter-like planet has none of the larger irregulars in some systems and three in others (expected 0.4 above 40 km after grinding).

Pure functions of the planet's numbers; the random numbers come from the caller's stream (`STORE_MAX_BODIES` draws for the sizes and as many for the windows; the exchange takes two more from the caller). Nothing calls it yet; GEN.1083 wires it into the moons. Stand-ins are IP1 to IP14
of `docs/design/research/gap-filling-plan.md` section 6.118.
"""

import math
from dataclasses import dataclass

from planetgen.physics import belt_ledger, constants, irregular_orbits, irregular_zone
from planetgen.util.checks import finite_domain

REFERENCE_DISC_MASS_EARTH = 30.0
"""float: The mass of the disc whose size law the counts are, Earth masses (NDD14 and NVM07 use 20 to 35; Bottke et al. 2010 20.4 and 40.8). [Stand-in IP1, tier T1, tweak L0; row CM3.]"""
DISC_REFERENCE_RADIUS_KM = 40.0
"""float: The radius at which the disc's count is given. [Stand-in IP1, tier T1, tweak L0.]"""
DISC_COUNT_AT_REFERENCE = 3.5e7
"""float: The bodies above `DISC_REFERENCE_RADIUS_KM` in the reference disc ((3 to 4)e7, calibrated to the Jupiter Trojans by NDD14). [Stand-in IP1, tier T1, tweak L0.]"""
DISC_BREAK_RADIUS_KM = 70.0
"""float: The radius where the disc's size law steepens (NDD14; NVM07 and Bottke et al. 2010 break at 100 km across). [Stand-in IP1, tier T1, tweak L0.]"""
DISC_INDEX_BELOW = 1.8
"""float: The cumulative index below the break (1.75 to 2.75 over the three papers). [Stand-in IP1, tier T1, tweak L0.]"""
DISC_INDEX_ABOVE = 5.0
"""float: The cumulative index above the break (4.25 to 5.5 over the three papers). [Stand-in IP1, tier T1, tweak L0.]"""
CAPTURED_DENSITY_KG_M3 = 2000.0
"""float: The bulk density that turns a captured body's radius into a mass (the density of row CM15's 10.6 km at 1e16 kg; the composition of GEN.1080 replaces it body by body). [Stand-in IP2, tier T9, tweak L0.]"""

ENCOUNTER_YIELD = 2.0e-10
"""float: Captures per disc particle per planetary encounter at a planet with Jupiter's Hill sphere or more (1.3e-10 to 2.6e-10 in NDD14's three histories). [Stand-in IP3, tier T2, tweak L1; row CM1.]"""
REFERENCE_HILL_AU = 0.35
"""float: The Hill radius at which the encounter yield applies in full, AU (Jupiter's, Jewitt and Haghighipour 2007 Table 1); a smaller one scales as its cube, a larger one gets the same. [Stand-in IP4, tier T9, tweak L1.]"""
REFERENCE_ENCOUNTERS = 130.0
"""float: The encounters of the reference instability the gas channel is scaled by, the mean of NDD14's 59, 280 and 49. [Stand-in IP5, tier T2, tweak L0.]"""
GAS_WINDOW_STRENGTH = 0.2
"""float: The gas window's captures as a share of the reference instability's (0.1 to 0.3); zero for a planet with no gas envelope. [Stand-in IP5, tier T6, tweak L3; row CM14.]"""
GAS_SIZE_WINDOW_KM = (1.0, 500.0)
"""tuple: The radii, km, that gas drag can capture. [Stand-in IP6, tier T1, tweak L0; row CM18.]"""
GAS_MIN_DISC_LIFETIME_MYR = 1.0
"""float: The shortest gas disc that gives a gas window. [Stand-in IP6, tier T6, tweak L0; row CM14.]"""
EXCHANGE_CHANCE = 0.008
"""float: The chance of a massive capture by binary exchange at an ice giant when every planetesimal is a binary (0.002 to 0.035, geometric mean 0.008). [Stand-in IP7, tier T1, tweak L1; row CM11.]"""
EXCHANGE_MIN_RADIUS_KM = 300.0
"""float: The smallest body the exchange route captures (a few hundred km; the design's reading of the binaries that give Triton). [Stand-in IP8, tier T6, tweak L3.]"""
TEMPORARY_MIN_PLANET_RADII = 3.0
"""float: A planet whose prograde stable limit is under this many of its own radii holds no temporary captures. [Stand-in IP9, tier T9, tweak L0.]"""

DISRUPTION_SIZE_RATIO = 0.1
"""float: A projectile of this fraction of a target's radius (a thousandth of its mass) destroys it at the 1 to 7 km/s of irregulars colliding (specific energy of a few 1e3 J/kg, of the order of a weak body's strength; recalled, not read). [Stand-in IP10, tier T6, tweak L3.]"""
MODEL_MIN_RADIUS_KM = 0.4
"""float: The smallest body counted when a population's mass is added up (Sheppard, Jewitt and Kleyna 2005's survey goes to 0.4 km; the mass below is about a hundredth of the total for a cumulative index under 2). [Stand-in IP11, tier T1, tweak L0.]"""
MODEL_MAX_RADIUS_KM = 3000.0
"""float: The upper end of the search for a body's radius. [Stand-in IP11, tier T9, tweak L0.]"""
INTEGRAL_STEPS = 120
"""int: The steps of the log-spaced sum that gives a population's mass (numerical)."""

STORE_FLOOR_KG = 1.0e16
"""float: Every irregular above this mass is an object (Boss 4:34 PM EDT 2026-10-30; R = 10.6 km at 2000 kg/m^3). [Stand-in IP12, tier T9, tweak L1; row CM15.]"""
STORE_MAX_BODIES = 100
"""int: The most irregulars stored for one planet. [Stand-in IP12, tier T9, tweak L1; row CM15.]"""

_KM = 1.0e3


class IrregularPopulationError(ValueError):
    """An input the model cannot use."""


def _break_count():
    return DISC_COUNT_AT_REFERENCE * (DISC_BREAK_RADIUS_KM / DISC_REFERENCE_RADIUS_KM) ** -DISC_INDEX_BELOW


@finite_domain()
def source_count_above(radius_km, scale=1.0):
    """The reservoir's bodies of radius over `radius_km`: the disc's size law times `scale` (the reservoir's mass over `REFERENCE_DISC_MASS_EARTH`)."""
    if not radius_km > 0.0 or not scale >= 0.0:
        raise IrregularPopulationError("the radius must be positive and the scale not negative")
    if radius_km <= DISC_BREAK_RADIUS_KM:
        return scale * DISC_COUNT_AT_REFERENCE * (radius_km / DISC_REFERENCE_RADIUS_KM) ** -DISC_INDEX_BELOW
    return scale * _break_count() * (radius_km / DISC_BREAK_RADIUS_KM) ** -DISC_INDEX_ABOVE


@finite_domain()
def source_radius_at_count(count, scale=1.0):
    """The inverse of `source_count_above`: the radius with `count` bodies expected above it, km."""
    if not count > 0.0 or not scale > 0.0:
        raise IrregularPopulationError("the count and the scale must be positive")
    at_break = scale * _break_count()
    if count >= at_break:
        return DISC_REFERENCE_RADIUS_KM * (count / (scale * DISC_COUNT_AT_REFERENCE)) ** (-1.0 / DISC_INDEX_BELOW)
    return DISC_BREAK_RADIUS_KM * (count / at_break) ** (-1.0 / DISC_INDEX_ABOVE)


@finite_domain()
def mass_kg(radius_km):
    """The mass of a captured body of this radius, kg, at `CAPTURED_DENSITY_KG_M3`."""
    if not radius_km > 0.0:
        raise IrregularPopulationError("the radius must be positive")
    return 4.0 / 3.0 * math.pi * (radius_km * _KM) ** 3 * CAPTURED_DENSITY_KG_M3


@finite_domain()
def radius_km(mass_kg_):
    """The radius of a captured body of this mass, km."""
    if not mass_kg_ > 0.0:
        raise IrregularPopulationError("the mass must be positive")
    return (3.0 * mass_kg_ / (4.0 * math.pi * CAPTURED_DENSITY_KG_M3)) ** (1.0 / 3.0) / _KM


@finite_domain()
def geometry_factor(hill_radius_au):
    """The share of the reference yield a target of this Hill radius gets: `min(1, (r_H / REFERENCE_HILL_AU)^3)`."""
    if not hill_radius_au > 0.0:
        raise IrregularPopulationError("the Hill radius must be positive")
    return min(1.0, (hill_radius_au / REFERENCE_HILL_AU) ** 3)


@dataclass(frozen=True)
class Pool:
    """
    What a planet's windows captured, before grinding.

    Attributes:
        encounter_fraction (float): Bodies captured by the encounter window per body of the reservoir.
        gas_fraction (float): The same for the gas window (for the bodies of `GAS_SIZE_WINDOW_KM` only).
        scale (float): The reservoir's mass over `REFERENCE_DISC_MASS_EARTH`.
    """

    encounter_fraction: float
    gas_fraction: float
    scale: float

    def count_above(self, radius_km_):
        """The bodies of radius over `radius_km_` that stayed bound, before grinding."""
        n = self.encounter_fraction * source_count_above(radius_km_, self.scale)
        top = GAS_SIZE_WINDOW_KM[1]
        if self.gas_fraction > 0.0 and radius_km_ < top:
            n += self.gas_fraction * (source_count_above(max(radius_km_, GAS_SIZE_WINDOW_KM[0]), self.scale) - source_count_above(top, self.scale))
        return n

    def gas_share(self, radius_km_):
        """The share of the captured bodies of this radius that the gas window brought, 0 to 1."""
        total = self.encounter_fraction + (self.gas_fraction if GAS_SIZE_WINDOW_KM[0] <= radius_km_ < GAS_SIZE_WINDOW_KM[1] else 0.0)
        if total <= 0.0:
            return 0.0
        return (total - self.encounter_fraction) / total


@finite_domain()
def make_pool(hill_radius_au, encounters, reservoir_mass_earth, gas_giant=False, gas_disc_myr=0.0):
    """
    The `Pool` of a planet.

    Args:
        hill_radius_au (float): The planet's Hill radius (`irregular_zone.hill_radius_au`).
        encounters (float): The giant-giant close encounters the planet had while the reservoir was there (0 for a system with no instability).
        reservoir_mass_earth (float): The mass of the small-body reservoir at the window, Earth masses (0 for none).
        gas_giant (bool): Whether the planet had a gas envelope (`envelope_mass_fraction` at or above the giant fraction) in the last stage of its formation.
        gas_disc_myr (float): How long the gas disc lasted, Myr.
    """
    if not hill_radius_au > 0.0 or not encounters >= 0.0 or not reservoir_mass_earth >= 0.0 or not gas_disc_myr >= 0.0:
        raise IrregularPopulationError("the Hill radius must be positive and the encounters, reservoir mass and gas lifetime not negative")
    factor = geometry_factor(hill_radius_au)
    encounter = ENCOUNTER_YIELD * encounters * factor
    gas = GAS_WINDOW_STRENGTH * ENCOUNTER_YIELD * REFERENCE_ENCOUNTERS * factor if gas_giant and gas_disc_myr >= GAS_MIN_DISC_LIFETIME_MYR else 0.0
    return Pool(encounter, gas, reservoir_mass_earth / REFERENCE_DISC_MASS_EARTH)


@finite_domain()
def exchange_chance(binary_fraction, hill_radius_au, gas_giant):
    """The chance of one massive capture by binary exchange: `EXCHANGE_CHANCE` times the binary fraction of the reservoir times the geometry factor; none at a gas giant."""
    if not 0.0 <= binary_fraction <= 1.0:
        raise IrregularPopulationError("the binary fraction must be from 0 to 1")
    return 0.0 if gas_giant else EXCHANGE_CHANCE * binary_fraction * geometry_factor(hill_radius_au)


@finite_domain()
def draw_exchange(u_chance, u_size, chance):
    """The radius of the body an exchange captured, km, or `None` when it did not: `u_size` places it on the reservoir's size law above `EXCHANGE_MIN_RADIUS_KM` (index `DISC_INDEX_ABOVE`)."""
    if not (0.0 <= u_chance <= 1.0 and 0.0 <= u_size < 1.0):
        raise IrregularPopulationError("the draws must be uniform numbers")
    if not u_chance < chance:
        return None
    return EXCHANGE_MIN_RADIUS_KM * (1.0 - u_size) ** (-1.0 / DISC_INDEX_ABOVE)


@finite_domain()
def holds_temporary_captures(prograde_limit_au, planet_radius_km):
    """Whether the planet's irregular zone is wide enough (`TEMPORARY_MIN_PLANET_RADII` planet radii) for bodies to pass through and stay a while: a flag, no object."""
    if not prograde_limit_au > 0.0 or not planet_radius_km > 0.0:
        raise IrregularPopulationError("the limit and the radius must be positive")
    return prograde_limit_au * constants.AU_M > TEMPORARY_MIN_PLANET_RADII * planet_radius_km * _KM


@finite_domain()
def zone_radius_au(hill_radius_au, planet_distance_au):
    """The radius of the sphere the captured bodies fill, AU: the median of their orbits (`irregular_orbits.median_before_survival`) held inside the prograde stable limit."""
    if not hill_radius_au > 0.0 or not planet_distance_au > 0.0:
        raise IrregularPopulationError("the Hill radius and the distance must be positive")
    return min(irregular_orbits.median_before_survival(planet_distance_au), irregular_zone.PROGRADE_STABLE_FRACTION) * hill_radius_au


@finite_domain()
def intrinsic_probability(planet_mass_kg, zone_radius_au_):
    """The intrinsic collision probability of the zone, per square kilometre per year: the relative speed (`sqrt(2)` times the circular speed at the zone's radius) over the zone's volume."""
    if not planet_mass_kg > 0.0 or not zone_radius_au_ > 0.0:
        raise IrregularPopulationError("the mass and the radius must be positive")
    r_m = zone_radius_au_ * constants.AU_M
    speed_km_yr = math.sqrt(2.0 * constants.G * planet_mass_kg / r_m) / _KM * constants.SECONDS_PER_YEAR
    return speed_km_yr / (4.0 / 3.0 * math.pi * (r_m / _KM) ** 3)


@finite_domain()
def collision_time_yr(pool, radius_km_, probability):
    """`t_c` for bodies of this radius: the time in which the projectiles above `DISRUPTION_SIZE_RATIO` times it destroy one, years (infinity where the pool has none)."""
    if not radius_km_ > 0.0 or not probability > 0.0:
        raise IrregularPopulationError("the radius and the probability must be positive")
    projectiles = pool.count_above(DISRUPTION_SIZE_RATIO * radius_km_)
    if projectiles <= 0.0:
        return math.inf
    return 1.0 / (probability * (radius_km_ * (1.0 + DISRUPTION_SIZE_RATIO)) ** 2 * projectiles)


@finite_domain()
def surviving_fraction(age_yr, collision_time_yr_):
    """`1 / (1 + t / t_c)`, the share of the bodies left after `age_yr`."""
    if not age_yr >= 0.0 or not collision_time_yr_ > 0.0:
        raise IrregularPopulationError("the age must not be negative and the collision time must be positive")
    return 1.0 / (1.0 + age_yr / collision_time_yr_)


@finite_domain()
def ground_count_above(pool, radius_km_, probability, age_yr):
    """The bodies of radius over `radius_km_` left at `age_yr`: the pool's count times the surviving fraction of that radius."""
    return pool.count_above(radius_km_) * surviving_fraction(age_yr, collision_time_yr(pool, radius_km_, probability))


@finite_domain()
def radius_at_count(pool, probability, age_yr, count):
    """The radius with `count` ground bodies expected above it, km, or `None` when the whole population has fewer (found by bisection in the log of the radius)."""
    if not count > 0.0:
        raise IrregularPopulationError("the count must be positive")
    if not ground_count_above(pool, MODEL_MIN_RADIUS_KM, probability, age_yr) > count:
        return None
    low, high = MODEL_MIN_RADIUS_KM, MODEL_MAX_RADIUS_KM
    if ground_count_above(pool, high, probability, age_yr) > count:
        return high
    for _ in range(40):
        mid = math.sqrt(low * high)
        if ground_count_above(pool, mid, probability, age_yr) > count:
            low = mid
        else:
            high = mid
    return math.sqrt(low * high)


@finite_domain()
def mass_below_kg(pool, probability, age_yr, radius_km_):
    """The expected mass of the ground bodies smaller than `radius_km_` (down to `MODEL_MIN_RADIUS_KM`), kg, from a log-spaced sum of `m dN`."""
    if not radius_km_ > MODEL_MIN_RADIUS_KM:
        return 0.0
    step = (radius_km_ / MODEL_MIN_RADIUS_KM) ** (1.0 / INTEGRAL_STEPS)
    total, r, n_low = 0.0, MODEL_MIN_RADIUS_KM, ground_count_above(pool, MODEL_MIN_RADIUS_KM, probability, age_yr)
    for _ in range(INTEGRAL_STEPS):
        upper = r * step
        n_high = ground_count_above(pool, upper, probability, age_yr)
        total += max(n_low - n_high, 0.0) * mass_kg(math.sqrt(r * upper))
        r, n_low = upper, n_high
    return total


@dataclass(frozen=True)
class Irregulars:
    """
    The stored irregulars of a planet.

    Attributes:
        radii_km (tuple): The stored bodies' radii, largest first; the index is the rank, stable for a stream whatever the planet.
        masses_kg (tuple): Their masses.
        gas_window (tuple): For each, whether the gas window brought it (otherwise the encounter window).
        remainder_kg (float): The expected mass of everything not stored (the L2 and L3 tiers).
        above_floor (int): How many bodies are above `STORE_FLOOR_KG`.
    """

    radii_km: tuple
    masses_kg: tuple
    gas_window: tuple
    remainder_kg: float
    above_floor: int

    @property
    def covered(self):
        """The share of the population's mass in the stored bodies (1 for a population with no mass)."""
        stored = sum(self.masses_kg)
        total = stored + self.remainder_kg
        return stored / total if total > 0.0 else 1.0


def draw_unit_numbers(stream):
    """The planet's `2 * STORE_MAX_BODIES` draws, the same number from `stream` whatever the planet: the unit exponentials of the size ranks, then one uniform number per rank to name its window."""
    exponentials = [-math.log(1.0 - stream.random()) for _ in range(STORE_MAX_BODIES)]
    return exponentials, [stream.random() for _ in range(STORE_MAX_BODIES)]


@finite_domain()
def irregulars_from(pool, probability, age_yr, exponentials, window_numbers):
    """The `Irregulars` of a pool from `draw_unit_numbers` (see the module text for the rule that decides how many are stored)."""
    if len(exponentials) != STORE_MAX_BODIES or len(window_numbers) != STORE_MAX_BODIES:
        raise IrregularPopulationError("one exponential and one window number per possible rank")
    radii, gap = [], 0.0
    for e in exponentials:
        gap += e
        r = radius_at_count(pool, probability, age_yr, gap)
        if r is None:
            break
        radii.append(r)
    masses = [mass_kg(r) for r in radii]
    above = sum(1 for m in masses if m >= STORE_FLOOR_KG)
    n = min(above, STORE_MAX_BODIES)             # a planet with nothing above the floor stores nothing: its bodies are the L2 and L3 tiers

    def remainder(count):
        return mass_below_kg(pool, probability, age_yr, radii[count - 1] if count else MODEL_MAX_RADIUS_KM)

    left = remainder(n)
    held = sum(masses[:n])
    while 0 < n < len(radii) and held < belt_ledger.LEDGER_COVERAGE * (held + left):
        held += masses[n]
        n += 1
        left = remainder(n)
    windows = tuple(window_numbers[k] < pool.gas_share(radii[k]) for k in range(n))
    return Irregulars(tuple(radii[:n]), tuple(masses[:n]), windows, left, above)


def draw_irregulars(pool, probability, age_yr, stream):
    """The `Irregulars` of a pool drawn from `stream` (a `draw.Stream`; `2 * STORE_MAX_BODIES` draws of `random()`)."""
    exponentials, window_numbers = draw_unit_numbers(stream)
    return irregulars_from(pool, probability, age_yr, exponentials, window_numbers)
