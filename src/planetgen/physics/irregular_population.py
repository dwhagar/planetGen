# planetgen/physics/irregular_population.py

"""
How many irregular moons a planet holds, and which are heavy enough to store (GEN.1077; docs/design/systems/captured-asteroid-moons.md 4.1, 4.2 and 4.5, rows CM1 to CM20 of docs/design/research/gap-filling-plan.md section 6.100)
===========================================================================================================================================================================================================================

An irregular moon is a small body captured from a reservoir of planetesimals and ground down since. The unit answers four questions, none of which draws an orbit (`irregular_orbits`) or a composition:

1. **Which planets capture.** Only a giant (`moon_composition.parent_kind`, an envelope of at least `GIANT_ENVELOPE_FRACTION`): small near-circular moons of a rocky planet are made from an impact disc, not captured (CM16), and a
   dwarf planet has no reservoir to draw from here. A planet whose irregular zone is empty (`irregular_zone.Zone.empty`: a hot planet, its stable limit inside its regular moons) captures nothing that stays.
2. **How much is captured, in windows (CM1, CM2, CM14, CM18; Boss 4:34 PM EDT 2026-10-30, no continuous capture).** The capture probability per particle of the reservoir is a sum of channels, each a window in the planet's history:
   * the **encounter** channel, `encounters` close approaches of a second giant in an instability (the architecture's number, 49 to 280 in the three histories read) times `ENCOUNTER_YIELD` stable captures per source particle per
     encounter (Nesvorny, Vokrouhlicky and Deienno 2014, Table 1: 1.3e-10 to 2.6e-10). It reproduces the per-planet probabilities of CM2, 1.3e-8 to 3.6e-8 for Jupiter. The four giants of that work capture about the same
     per particle although their Hill volumes differ by a factor 100, so the yield is not scaled up for a bigger Hill sphere, only down for one smaller than Jupiter's (`CAPTURE_HILL_AU`, as the Hill volume is the target);
   * the **gas** channel, only for a planet with a massive primordial envelope (`GAS_MIN_ENVELOPE_EARTH`, CM18: Jupiter about 300 and Saturn 100 earth masses of gas, the ice giants about 1), a fraction `GAS_WINDOW_SHARE` of what
     `GAS_REFERENCE_ENCOUNTERS` encounters would give (CM14), so a gas giant that had no instability still holds a fraction of the Solar System's number and an ice giant that had none holds none;
   * the **exchange** channel, the capture of a Triton-class body from a binary planetesimal (CM11), is not a swarm of small moons: `Windows.triton_class` is its expected number for a reservoir whose every body is a binary, for an ice
     giant only (Jupiter and Saturn, 2e-4 or fewer); the wiring multiplies it by the reservoir's binary fraction and the moon draw makes the body;
   * the **temporary** capture (CM13) leaves nothing and has no channel.
3. **What is left now (CM3, CM5, CM9).** The captured bodies follow the reservoir's size law (cumulative 1.8 below a radius of 70 km, 5 above); collisions among the irregulars remove 99 percent of the mass in the first hundreds of
   Myr and leave a steady state with a shallower size law, cumulative 1.0 (CM5, CM17). The design's `N(t) = N0 / (1 + t / t_c)` is used in its two limits: a body of radius `R` with `t_c(R) = GRIND_REFERENCE_YR (R / GRIND_REFERENCE_KM)^GRIND_SIZE_EXPONENT`
   below `t_c = age` keeps the fraction `t_c / age`, and above it all of it. The size where `t_c = age` is `Survivors.full_radius_km` (32 km at 4.5 Gyr); the survivors' size law is the reservoir's with the exponent lowered
   by `GRIND_SIZE_EXPONENT` below it, so it is a `belt_ledger.SizeLaw` again and the belt ledger's machinery serves unchanged. The grinding takes about a tenth of the mass, not the 99 percent of CM9: with a cumulative index of 1.8 the mass is in the biggest bodies, which it does not reach, and the count above 5 km is what is fitted (`GRIND_REFERENCE_YR`). `collision_time_scale` multiplies `t_c`: the check that the irregular zones of the four giants grind alike
   (GEN.1079) sets it per planet.
4. **Which are stored (CM15, Boss's rule).** The bodies of the survivors are an order statistic drawn as a belt's (`belt_ledger.ledger_from` with `realised`, as a handful of bodies has no budget to be scaled to; stable ranks), stored from the top: every body of `STORAGE_FLOOR_KG` (1e16 kg, a radius of
   10.6 km) or more, then the next heaviest until they hold `STORAGE_COVERAGE` of the mass, at most `STORAGE_MAX_BODIES` and none lighter than `STORAGE_MIN_MASS_KG`. Those are the individual objects. The rest of the population is the belt design's
   other two tiers: the bodies after the stored ones are the later ranks of the same exponentials (drawn on demand from the planet's seed, addressable by rank) and the whole remainder is `Ledger.remainder_kg` (the law's mass below the lightest stored body) with the law.

The counts are of bodies that **stay**: the yields of CM1 and CM2 are stable captures, so the wiring draws an orbit for each stored body until one stays (`irregular_orbits.draw_orbit`, `irregular_zone.verdict`) and does not cut a
captured count a second time. The reservoir is the outer planetesimal disc of the instability (`OUTER_DISC`, 3.5e7 bodies above a radius of 40 km); a system whose disc is lighter or heavier scales it (`outer_disc`); a belt of the
system as a reservoir is the belt design's work (GEN.1083 wires it). The giants' own head counts (107, 268, 10 and 8 irregulars today; 8 above 5 km at Jupiter; 8, 5, 6 and 8 above the storage floor; Triton aside) are targets of the tests,
never inputs. Stand-ins are CP1 to CP10 of `docs/design/research/gap-filling-plan.md` section 6.118. Nothing reads the module yet; GEN.1083 and the moon wiring will.

Use:

    population = irregular_population(planet_mass_kg, envelope_fraction, distance_au, star_mass_kg, encounters, age_yr, stream)     # a Population (nothing for a rocky planet)
"""

import bisect
import math
from dataclasses import dataclass

from planetgen.physics import belt_ledger, irregular_zone, moon_composition
from planetgen.util.checks import finite_domain

EARTH_MASS_KG = 5.9722e24
"""float: The mass of the Earth, kg, the unit of the envelope and disc masses here."""

ENCOUNTER_YIELD = 1.8e-10
"""float: Stable captures per reservoir particle per planetary encounter at a gas giant, the geometric middle of 1.3e-10 to 2.6e-10 (Nesvorny, Vokrouhlicky and Deienno 2014 Table 1, READ). [Stand-in CP1 (CM1), tier T2, tweak L1.]"""

GAS_WINDOW_SHARE = 0.2
"""float: The gas-drag channel as a share of the encounter channel of `GAS_REFERENCE_ENCOUNTERS` encounters, the middle of 0.1 to 0.3 (CM14; no paper gives it). [Stand-in CP2, tier T6, tweak L3.]"""

GAS_REFERENCE_ENCOUNTERS = 129.0
"""float: The encounters the gas share is a fraction of, the mean of the three histories read (59, 280 and 49). [Stand-in CP2, tier T2, tweak L1.]"""

GAS_MIN_ENVELOPE_EARTH = 10.0
"""float: The least primordial envelope, earth masses, for the gas channel: between the ice giants' about 1 and Saturn's about 100 (CM18; the gate is the size of the envelope, the threshold is argued). [Stand-in CP3, tier T6, tweak L1.]"""

CAPTURE_HILL_AU = 0.35
"""float: The Hill radius, AU, above which the yield per encounter is Jupiter's; below it falls as the Hill volume, `(r_H / this)^3` (Jupiter's is 0.35 AU, the planet CM1 was run for; CM2 gives the bigger giants no more). [Stand-in CP4, tier T6, tweak L1.]"""

TRITON_EXCHANGE_PROBABILITY = 0.0084
"""float: The chance that one ice giant captures a Triton-class body from a binary when every planetesimal is a binary, the geometric middle of 0.002 to 0.035 (CM11, Vokrouhlicky and Nesvorny). [Stand-in CP5, tier T2, tweak L1.]"""

GRIND_REFERENCE_KM = 5.0
"""float: The radius, km, at which `GRIND_REFERENCE_YR` is the time a body is ground down in. [Stand-in CP6, tier T6, tweak L3.]"""

GRIND_REFERENCE_YR = 1.0e9
"""float: The collision time of a body of radius `GRIND_REFERENCE_KM` in the irregular zone, years: chosen so that a Jupiter with the mean history holds 8 irregulars above a radius of 5 km at 4.5 Gyr (CM4, CM6), within the 9e7 to 1.3e10 yr of CM9. [Stand-in CP6, tier T6, tweak L3.]"""

GRIND_SIZE_EXPONENT = 0.8
"""float: How the collision time grows with radius, `t_c ~ R^this`: the reservoir's cumulative index 1.8 less the survivors' 1.0 (CM3, CM5, CM17), which is what makes a ground population settle on the observed law. [Stand-in CP7, tier T1, tweak L0.]"""

STORAGE_FLOOR_KG = 1.0e16
"""float: Every irregular moon of at least this mass is its own stored object (a radius of 10.6 km at 2000 kg/m3). [Stand-in CP8 (CM15), tier T9, tweak L1.]"""

STORAGE_COVERAGE = 0.9
"""float: The share of the population's mass the stored bodies must hold. [Stand-in CP8 (CM15), tier T9, tweak L1.]"""

STORAGE_MAX_BODIES = 100
"""int: The most irregular moons stored for one planet. [Stand-in CP8 (CM15), tier T9, tweak L1.]"""

STORAGE_MIN_MASS_KG = 1.0e12
"""float: No stored body is lighter than this (a radius of 0.6 km at 1000 kg/m3, near the smallest a survey has found, 0.4 km at Jupiter, CM5), so a population with a tiny mass is not a list of dust. [Stand-in CP9, tier T2, tweak L1.]"""

REFERENCE_DISC_MASS_EARTH = 27.5
"""float: The mass of the outer disc of the instability the reservoir's count belongs to, earth masses, the middle of 20 to 35 (CM3). [Stand-in CP10, tier T1, tweak L0.]"""

OUTER_DISC_ANCHOR_RADIUS_KM = 40.0
OUTER_DISC_ANCHOR_COUNT = 3.5e7
"""float: The reservoir has 3.5e7 bodies above a radius of 40 km, the middle of 3 to 4e7 (CM3, Nesvorny, Vokrouhlicky and Deienno 2014). [Stand-in CP10, tier T1, tweak L0.]"""


def _radius_mass_kg(radius_km, density_kg_m3):
    """The mass of a sphere of this radius and density."""
    return 4.0 / 3.0 * math.pi * (radius_km * 1000.0) ** 3 * density_kg_m3


@dataclass(frozen=True)
class Reservoir:
    """
    A population of small bodies a planet can capture from.

    Attributes:
        law (belt_ledger.SizeLaw): Its size law.
        mass_kg (float): Its mass, the one the law is normalised to.
    """

    law: belt_ledger.SizeLaw
    mass_kg: float


OUTER_DISC_LAW = belt_ledger.SizeLaw(breaks_km=(140.0,), exponents=(1.8, 5.0), density_kg_m3=2000.0)
"""SizeLaw: The size law of the outer disc: cumulative 1.8 below a radius of 70 km (a diameter of 140 km) and 5 above (CM3, Nesvorny, Vokrouhlicky and Deienno 2014, READ); the density is the 2000 kg/m3 of the radius table (CM6, CM15), which gives the disc about 25 earth masses, inside the 20 to 35 of CM3. [Stand-in CP10, tier T1, tweak L0.]"""


def outer_disc(disc_mass_earth=REFERENCE_DISC_MASS_EARTH):
    """
    The outer disc of an instability as a `Reservoir`: `OUTER_DISC_ANCHOR_COUNT` bodies above a radius of `OUTER_DISC_ANCHOR_RADIUS_KM` for a disc of `REFERENCE_DISC_MASS_EARTH`, in proportion to the disc's mass.

    The count above the anchor is the datum; the mass it implies for the law (about 25 earth masses at 2000 kg/m3) is inside the 20 to 35 of the disc.
    """
    if not disc_mass_earth > 0.0:
        raise ValueError("the disc's mass must be positive")
    anchor = _radius_mass_kg(OUTER_DISC_ANCHOR_RADIUS_KM, OUTER_DISC_LAW.density_kg_m3)
    mass = OUTER_DISC_ANCHOR_COUNT / belt_ledger.count_above(OUTER_DISC_LAW, 1.0, anchor)
    return Reservoir(OUTER_DISC_LAW, mass * disc_mass_earth / REFERENCE_DISC_MASS_EARTH)


@dataclass(frozen=True)
class Windows:
    """
    The captures of a planet, per particle of the reservoir.

    Attributes:
        encounter (float): The probability that a reservoir particle is captured by the instability's encounters and stays.
        gas (float): The same by gas drag in the primordial envelope.
        triton_class (float): The expected number of Triton-class bodies captured by binary exchange if every body of the reservoir is a binary (multiply by the binary fraction).
    """

    encounter: float
    gas: float
    triton_class: float

    @property
    def total(self):
        """The capture probability per particle of the small bodies, both channels."""
        return self.encounter + self.gas


NO_WINDOWS = Windows(0.0, 0.0, 0.0)


@finite_domain()
def capture_windows(planet_mass_kg, envelope_fraction, distance_au, star_mass_kg, encounters, outermost_regular_moon_au=0.0):
    """
    The `Windows` of a planet: all zero for a planet that is not a giant or whose irregular zone is empty (CM16).

    Args:
        planet_mass_kg (float): The planet's mass.
        envelope_fraction (float): Its envelope's share of the mass (`moon_composition.parent_kind` decides what is a giant).
        distance_au (float): Its semi-major axis.
        star_mass_kg (float): Its star's mass.
        encounters (float): The close approaches of another giant it has had in the instability (the architecture's number; 0 for none).
        outermost_regular_moon_au (float): The orbit of its outermost regular moon, 0 for none (the inner edge of the zone).
    """
    if encounters < 0.0 or not (0.0 <= envelope_fraction <= 1.0):
        raise ValueError("the encounters must not be negative and the envelope fraction must be between 0 and 1")
    if moon_composition.parent_kind(envelope_fraction, planet_mass_kg) != moon_composition.GIANT:
        return NO_WINDOWS
    zone = irregular_zone.zone(distance_au, planet_mass_kg, star_mass_kg, outermost_regular_moon_au)
    if zone.empty:
        return NO_WINDOWS
    geometry = min(1.0, (zone.hill_radius_au / CAPTURE_HILL_AU) ** 3)
    gas_giant = envelope_fraction * planet_mass_kg / EARTH_MASS_KG >= GAS_MIN_ENVELOPE_EARTH
    encounter = encounters * ENCOUNTER_YIELD * geometry
    gas = GAS_WINDOW_SHARE * GAS_REFERENCE_ENCOUNTERS * ENCOUNTER_YIELD * geometry if gas_giant else 0.0
    triton = TRITON_EXCHANGE_PROBABILITY * geometry if (encounters > 0.0 and not gas_giant) else 0.0
    return Windows(encounter, gas, triton)


@dataclass(frozen=True)
class Survivors:
    """
    The captured bodies that are left after grinding, as a size law.

    Attributes:
        law (belt_ledger.SizeLaw): Their cumulative size law.
        mass_kg (float): Their total mass, the one the law is normalised to.
        full_radius_km (float): The radius from which up, none has been ground (`t_c = age`).
    """

    law: belt_ledger.SizeLaw
    mass_kg: float
    full_radius_km: float

    def count_above(self, radius_km):
        """The expected number of survivors of at least this radius."""
        return belt_ledger.count_above(self.law, self.mass_kg, _radius_mass_kg(radius_km, self.law.density_kg_m3))


@finite_domain()
def full_radius_km(age_yr, collision_time_scale=1.0):
    """The radius at which the collision time is `age_yr`: smaller bodies have been ground to the fraction `t_c / age`, larger ones not at all (0 for no age)."""
    if not collision_time_scale > 0.0:
        raise ValueError("the collision time scale must be positive")
    if age_yr <= 0.0:
        return 0.0
    return GRIND_REFERENCE_KM * (age_yr / (GRIND_REFERENCE_YR * collision_time_scale)) ** (1.0 / GRIND_SIZE_EXPONENT)


@finite_domain()
def survivors(reservoir, capture_probability, age_yr, collision_time_scale=1.0):
    """
    The `Survivors` of a capture probability `capture_probability` per particle of `reservoir` after `age_yr` years of grinding, or `None` for a probability of zero.

    The reservoir's law is cut at `full_radius_km` and below it every stretch has its exponent lowered by `GRIND_SIZE_EXPONENT`, so the count is continuous there and the law mass is finite as before.
    """
    if capture_probability < 0.0:
        raise ValueError("the capture probability must not be negative")
    if capture_probability == 0.0:
        return None
    law = reservoir.law
    r_full = full_radius_km(age_yr, collision_time_scale)
    breaks = list(law.breaks_km)
    if r_full > 0.0:
        d_full = 2.0 * r_full
        if not any(math.isclose(d_full, b, rel_tol=1e-9) for b in breaks):
            bisect.insort(breaks, d_full)
        else:
            d_full = next(b for b in breaks if math.isclose(d_full, b, rel_tol=1e-9))
    else:
        d_full = 0.0
    exponents = []
    for i in range(len(breaks) + 1):
        if i == 0:
            probe = breaks[0] / 2.0
        elif i == len(breaks):
            probe = breaks[-1] * 2.0
        else:
            probe = math.sqrt(breaks[i - 1] * breaks[i])
        q = law.exponents[bisect.bisect_right(law.breaks_km, probe)]
        exponents.append(q - GRIND_SIZE_EXPONENT if i < len(breaks) and breaks[i] <= d_full * (1.0 + 1e-9) else q)
    ground = belt_ledger.SizeLaw(breaks_km=tuple(breaks), exponents=tuple(exponents), density_kg_m3=law.density_kg_m3)
    top_mass = math.pi / 6.0 * (breaks[-1] * 1000.0) ** 3 * law.density_kg_m3
    top_count = capture_probability * belt_ledger.count_above(law, reservoir.mass_kg, top_mass)
    return Survivors(ground, belt_ledger.mass_for_top_count(ground, top_count), r_full)


@dataclass(frozen=True)
class Population:
    """
    The irregular moons of a planet.

    Attributes:
        windows (Windows): The captures that made it.
        survivors (Survivors): What is left, `None` for a planet with no captures.
        ledger (belt_ledger.Ledger): The stored bodies (`masses_kg`, heaviest first, rank stable) and the mass of the rest (`remainder_kg`), `None` for a planet with no captures.
    """

    windows: Windows
    survivors: Survivors
    ledger: belt_ledger.Ledger

    def count_above(self, radius_km):
        """The expected number of irregulars of at least this radius (0 for none)."""
        return self.survivors.count_above(radius_km) if self.survivors else 0.0

    def stored_above(self, mass_kg):
        """How many of the stored bodies are of at least this mass."""
        return sum(1 for m in self.ledger.masses_kg if m >= mass_kg) if self.ledger else 0


def draw_stored(survivors_, stream):
    """The stored `belt_ledger.Ledger` of survivors from `stream`: `STORAGE_MAX_BODIES` draws of `random()`, the same number whatever the planet."""
    exponentials = belt_ledger.draw_exponentials(stream, STORAGE_MAX_BODIES)
    if survivors_ is None:
        return None
    return belt_ledger.ledger_from(survivors_.law, survivors_.mass_kg, exponentials, STORAGE_FLOOR_KG, STORAGE_COVERAGE, STORAGE_MAX_BODIES, STORAGE_MIN_MASS_KG, realised=True)


def irregular_population(planet_mass_kg, envelope_fraction, distance_au, star_mass_kg, encounters, age_yr, stream, reservoir=None,
                         outermost_regular_moon_au=0.0, collision_time_scale=1.0):
    """
    The `Population` of a planet: the windows, the survivors of the grinding and the stored bodies drawn from `stream` (`STORAGE_MAX_BODIES` draws, whatever the planet).

    `reservoir` is the outer disc of the reference instability unless the system's own is given (`outer_disc`).
    """
    windows = capture_windows(planet_mass_kg, envelope_fraction, distance_au, star_mass_kg, encounters, outermost_regular_moon_au)
    left = survivors(reservoir or outer_disc(), windows.total, age_yr, collision_time_scale)
    return Population(windows, left, draw_stored(left, stream))
