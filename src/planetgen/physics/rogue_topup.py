# planetgen/physics/rogue_topup.py

"""
The rogue top-up of a sector: how many rogues make up a mass shortfall, and what each is (GEN.1125)
===================================================================================================

Stage 3 of a sector's generation (`docs/design/systems/region-object-store.md` section 8.5) runs once per region run, after every sector has been built and the rogues that
landed from other sectors have arrived. It compares the mass the sector should hold with what it does, and makes up the shortfall with rogue planets, never more than the
Rogues level allows (None, Few, Normal, Many: `tuning.ROGUE_LEVEL_CEILINGS` rogues per star) and the rest, if the level stops it, left as a remainder for the lighter asteroid
fields (GEN.1126). This module is that planning step as a pure function; nothing calls it yet (the region run and the Rogues dial do).

**How many and how heavy.** The ceiling is `ceiling * stars` rogues for the sector less the rogues already there. Rogues are drawn one at a time, a bin from the existing mass bins
(`tuning.ROGUE_PLANET_MASS_BINS`: 96 percent terrestrial, 4 percent gas giant) by their per-star rates and a mass within the bin from the mass function `dN/dlnM` proportional to
`M^-slope`, until the shortfall is used or the ceiling reached. A rogue that would overshoot the shortfall is cut to what is left, if that is at least the lightest rogue
(`MIN_MASS_EARTH`), and otherwise the drawing stops, so **the placed mass plus the remainder is the shortfall exactly**.

**Where, and how fast.** Each rogue is placed uniformly in the sector inside `tuning.ROGUE_TOPUP_PLACEMENT` of its edge (the one stand-in Boss wants to play with). It is then
given what it would have if it had come from another sector: an age and a time since ejection from the star-formation history (`rogue_age`), an ejection speed from one of the
sector's own host stars thrown by a giant (`rogue_kick`; the scatterer is a Jupiter at a distance drawn log-uniformly from `SCATTERER_AU_RANGE`), an isotropic direction, the
velocity through the Galaxy as the local host's (drawn from the velocity dispersion) plus the drift, and an **apparent origin**: the distance and direction back along the drifted
path to where the donor would be now (`rogue_drift.ring_drift`, the along-orbit slide taken as an arc of the sector's galactocentric orbit, so a Gyr-old rogue is some kpc from its apparent
donor and never farther than the orbit's diameter). The origin is marked `apparent`, there is no donor system, so no link to one can exist: the sector's own systems
are recorded as the donors of the *mass* only.

The random numbers come from the caller's stream (`rng.random()`), so the draw is on the caller's own `draw.py` stream and the same stream gives the same plan. Stand-ins are in
`docs/design/research/gap-filling-plan.md` section 6.86.
"""

import math
from dataclasses import dataclass
from statistics import NormalDist

from planetgen import tuning
from planetgen.galaxy.age_law import AgeLaw
from planetgen.physics import constants, rogue_age, rogue_drift, rogue_kick, unit_registry as units

_kg_to_mearth = units.converter("kg", "Mearth")
_mearth_to_kg = units.converter("Mearth", "kg")

LEVELS = tuple(tuning.ROGUE_LEVEL_CEILINGS)
MIN_MASS_EARTH = min(low for low, _high, _rate in tuning.ROGUE_PLANET_MASS_BINS.values())
"""float: The lightest rogue, the lower edge of the lightest mass bin."""
SCATTERER_AU_RANGE = (1.0, 30.0)
"""tuple: Where the scattering giant of a top-up rogue's apparent origin sits, AU, drawn log-uniformly. Source tier T8, tweak level L3."""
VELOCITY_DISPERSION_KM_S = (35.0, 25.0, 20.0)
"""tuple: The thin disc's velocity dispersion `(sigma_U, sigma_V, sigma_W)` of the local hosts, km/s, where the caller gives none. Source tier T6, tweak level L1."""
DEFAULT_AGE_LAW = AgeLaw.uniform(tuning.ROGUE_PLANET_AGE_RANGE_GY)
GALACTOCENTRIC_RADIUS_PC = 8000.0
"""float: The galactocentric radius used where the caller gives none (about the Sun's). Source tier T6, tweak level L0."""
ORIGIN = "apparent"
"""str: The origin every top-up rogue carries; a donor link never exists for it."""

_JUPITER_KG = constants.JUPITER_MASS_TO_KG
_NORMAL = NormalDist()


class RogueTopUpError(ValueError):
    """An input the planner cannot use."""


@dataclass(frozen=True)
class ApparentOrigin:
    """
    Where a top-up rogue seems to have come from.

    Attributes:
        distance_pc (float): How far away the donor would be now.
        direction (tuple): The unit vector to it `(outward, along the rotation, up)`; `(0, 0, 0)` for a rogue that has not moved.
        age_gy, since_ejection_gy (float): The rogue's age and the time since the ejection.
        speed_km_s (float): The ejection speed.
        origin (str): Always `ORIGIN`.
    """

    distance_pc: float
    direction: tuple
    age_gy: float
    since_ejection_gy: float
    speed_km_s: float
    origin: str = ORIGIN


@dataclass(frozen=True)
class TopUpRogue:
    """
    One rogue the top-up places.

    Attributes:
        mass_earth (float): Its mass, Earth masses.
        mass_bin (str): The mass bin.
        position_pc (tuple): Where it is in the sector, pc from the sector's corner.
        velocity_km_s (tuple): Its velocity relative to the local standard of rest `(U, V, W)`.
        apparent_origin (ApparentOrigin): Where it seems to have come from.
        donor (None): Always `None`: an apparent origin is not a system.
    """

    mass_earth: float
    mass_bin: str
    position_pc: tuple
    velocity_km_s: tuple
    apparent_origin: ApparentOrigin
    donor: None = None


@dataclass(frozen=True)
class TopUp:
    """
    The plan for a sector.

    Attributes:
        rogues (tuple): The `TopUpRogue`s.
        placed_mass_kg (float): Their mass.
        remainder_kg (float): The shortfall they did not make up (the level stopped them), for the asteroid fields.
        ceiling (int): The most rogues the level allowed in all.
        room (int): The room the sector had under it after the rogues already there.
        stopped_by_level (bool): The ceiling, not the shortfall, ended the drawing.
    """

    rogues: tuple
    placed_mass_kg: float
    remainder_kg: float
    ceiling: int
    room: int
    stopped_by_level: bool


def ceiling_count(level, stars):
    """The most rogues a sector of `stars` stars may hold at a level: `floor(ceiling * stars)`."""
    if level not in tuning.ROGUE_LEVEL_CEILINGS:
        raise RogueTopUpError(f"the Rogues level is one of {LEVELS}, not {level!r}")
    if not stars >= 0:
        raise RogueTopUpError("the number of stars must not be negative")
    return int(math.floor(tuning.ROGUE_LEVEL_CEILINGS[level] * stars + 1e-9))


def mass_bin_of(mass_earth):
    """The mass bin a mass falls in (the heaviest bin takes anything above its edge)."""
    found = None
    for name, (low, high, _rate) in tuning.ROGUE_PLANET_MASS_BINS.items():
        if mass_earth >= low:
            found = name
        if mass_earth < high:
            break
    if found is None:
        raise RogueTopUpError("a rogue is at least the lightest bin's mass")
    return found


def draw_mass_earth(rng):
    """A rogue's mass in Earth masses: a bin by the per-star rates, then `dN/dlnM` proportional to `M^-slope` inside it."""
    bins = list(tuning.ROGUE_PLANET_MASS_BINS.items())
    total = math.fsum(rate for _n, (_l, _h, rate) in bins)
    pick = rng.random() * total
    chosen = bins[-1]
    for item in bins:
        pick -= item[1][2]
        if pick <= 0.0:
            chosen = item
            break
    _name, (low, high, _rate) = chosen
    slope = tuning.ROGUE_PLANET_MASS_FUNCTION_SLOPE
    return (low ** -slope + rng.random() * (high ** -slope - low ** -slope)) ** (-1.0 / slope)


def _gauss(rng):
    u = min(max(rng.random(), 1e-12), 1.0 - 1e-12)
    return _NORMAL.inv_cdf(u)


def place_in_sector(rng, sector_edge_pc, placement=None):
    """A uniform position in the sector, inside the central `placement` share of its edge (pc from the corner)."""
    share = tuning.ROGUE_TOPUP_PLACEMENT if placement is None else placement
    if not 0.0 < share <= 1.0:
        raise RogueTopUpError("the placement share must be above 0 and at most 1")
    start = 0.5 * (1.0 - share) * sector_edge_pc
    span = share * sector_edge_pc
    return tuple(start + rng.random() * span for _ in range(3))


def apparent_donor(kick, time_myr, galactocentric_radius_pc, galaxy=rogue_drift.SOLAR_NEIGHBOURHOOD):
    """
    The donor's offset from the rogue now `(outward, along the rotation, up)`, pc, where the along-orbit slide is an arc of the orbit of radius `galactocentric_radius_pc`:
    the straight line of the linear solution would put a rogue that has slid for Gyr tens of kiloparsecs from a star on the same orbit.
    """
    ring = rogue_drift.ring_drift(kick, time_myr, galactocentric_radius_pc, galaxy)
    r = galactocentric_radius_pc
    return (r * (math.cos(ring.azimuth_rad) - 1.0) - ring.radial_pc, -r * math.sin(ring.azimuth_rad), -ring.vertical_pc)


def apparent_history(rng, age_law, dispersion, radius_pc):
    """
    What an object placed in a sector would have if it had come from another one: its velocity relative to the local standard of rest `(U, V, W)`, km/s, and its `ApparentOrigin`.
    The draws, in order: the age and time since ejection, the scatterer's distance, the ejection speed, the direction and the host's velocity (the
    rogue top-up and the lighter asteroid-field top-up, GEN.1126, share it).
    """
    ages = rogue_age.rogue_age(age_law, rng.random(), rng.random())
    low, high = SCATTERER_AU_RANGE
    scatter_au = low * (high / low) ** rng.random()
    speed = rogue_kick.ejection_speed_km_s(rng.random(), _JUPITER_KG, scatter_au)
    direction = rogue_kick.isotropic_direction(rng.random(), rng.random())
    kick = rogue_kick.kick(speed, direction)
    host_velocity = tuple(sigma * _gauss(rng) for sigma in dispersion)
    velocity = rogue_kick.rogue_velocity(host_velocity, kick, ages.since_ejection_gy * 1.0e3)
    offset = apparent_donor(kick, ages.since_ejection_gy * 1.0e3, radius_pc)
    distance = math.sqrt(sum(c * c for c in offset))
    unit = tuple(c / distance for c in offset) if distance > 0.0 else (0.0, 0.0, 0.0)
    return velocity, ApparentOrigin(distance, unit, ages.age_gy, ages.since_ejection_gy, speed)


def _rogue(mass_earth, rng, sector_edge_pc, age_law, dispersion, placement, radius_pc):
    position = place_in_sector(rng, sector_edge_pc, placement)
    velocity, origin = apparent_history(rng, age_law, dispersion, radius_pc)
    return TopUpRogue(mass_earth, mass_bin_of(mass_earth), position, velocity, origin)


def plan(rng, level, stars, shortfall_kg, host_masses_kg, rogues_present=0, sector_edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC, age_law=DEFAULT_AGE_LAW,
         dispersion_km_s=VELOCITY_DISPERSION_KM_S, placement=None, galactocentric_radius_pc=GALACTOCENTRIC_RADIUS_PC):
    """
    The top-up of one sector.

    Args:
        rng: A stream with `random()`, 0 up to 1 (the caller's `draw.py` stream).
        level (str): The Rogues level, one of `LEVELS`.
        stars (int): The stars in the sector (the ceiling is per star).
        shortfall_kg (float): The mass the sector is missing.
        host_masses_kg (sequence): The masses of the sector's own host stars, whose systems donate the mass (at least one when anything is to be placed).
        rogues_present (int): The rogues the sector already holds (arrived or made), counted against the ceiling.
        sector_edge_pc (float): The sector's edge.
        age_law (AgeLaw): The star-formation history's age law for the sector's population.
        dispersion_km_s (tuple): The local hosts' velocity dispersion.
        placement (float, optional): Override of `tuning.ROGUE_TOPUP_PLACEMENT`.
        galactocentric_radius_pc (float): The sector's distance from the galactic centre, for the arc the apparent origin lies along.

    Returns:
        TopUp: The rogues, the mass placed and the remainder.
    """
    if not shortfall_kg >= 0.0 or not math.isfinite(shortfall_kg):
        raise RogueTopUpError("the shortfall must be a finite mass, not negative")
    if rogues_present < 0:
        raise RogueTopUpError("the rogues already present must not be negative")
    ceiling = ceiling_count(level, stars)
    room = max(0, ceiling - rogues_present)
    rogues, placed = [], 0.0
    remaining = shortfall_kg
    lightest_kg = _mearth_to_kg(MIN_MASS_EARTH)
    if room and remaining >= lightest_kg and not host_masses_kg:
        raise RogueTopUpError("a top-up needs at least one host star to donate the mass")
    while len(rogues) < room and remaining >= lightest_kg:
        mass = draw_mass_earth(rng)
        mass_kg = _mearth_to_kg(mass)
        last = mass_kg >= remaining
        if last:
            mass_kg, mass = remaining, _kg_to_mearth(remaining)
        rogues.append(_rogue(mass, rng, sector_edge_pc, age_law, dispersion_km_s, placement, galactocentric_radius_pc))
        placed += mass_kg
        remaining -= mass_kg
        if last:
            remaining = 0.0
    remainder = max(shortfall_kg - placed, 0.0)
    stopped = len(rogues) >= room and remainder >= lightest_kg
    return TopUp(tuple(rogues), placed, remainder, ceiling, room, stopped)
