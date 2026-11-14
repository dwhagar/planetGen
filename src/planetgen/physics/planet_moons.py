# planetgen/physics/planet_moons.py

"""
Moons of planet size around the giants (GEN.1157 slice 0, stand-ins PSM1 to PSM4 of docs/design/research/gap-filling-plan.md section 6.276)
=========================================================================================================================================

How often a giant planet holds a moon as big as a planet (`PLANET_MOON_MIN_MASS_KG`, a tenth of an Earth mass or more), and how big, from the two ways a giant can get one. A pure module:
nothing is stored and nothing is wired, so it changes no generated body. Design: docs/design/systems/planet-size-moons.md. Report: reports/gen1157-planet-size-moons-2026-11-12.md.

* **The circumplanetary disk** makes regular moons whose total is about 1e-4 to 6e-4 of the planet (`moon_composition.total_moon_mass_ratio` at the ends of `ALPHA_OVER_F`, [R0642], [R0643:
  Eq. 4.15]). No disk moon can pass that share of its parent, so `disk_moon_ceiling_kg` is the ceiling the generator should hold the disk's moons to; a planet-size moon of the disk needs a parent of
  a few Jupiter masses (`disk_parent_mass_kg`), the super-Jovians of [R1899]. A disk moon is drawn with the disk, not here.
* **Capture** (a giant that meets a binary, keeps one and flings the other away, [R0037, R5229]; or a planet pulled in during a migration instability, [R5230]; the candidate
  Kepler-1625 b-i of [R4322] is the case studied, and [R1900] read its signal another way) is rare and not tied to the disk: `exchange_chance` of a giant has one, `exchange_ceiling_kg` is the most it can hold (a share of the parent, `EXCHANGE_MAX_RATIO`) and
  `draw_exchange_moon` draws its mass from the generator's own moon-mass law (`tuning.MOON_MASS_FUNCTION_SLOPE`) between `PLANET_MOON_MIN_MASS_KG` and that ceiling.
* **Giant impact** is the rocky parents' route (`tuning.ROCKY_IMPACT_MOON_CHANCE`, GEN.962) and is not repeated here.

The moon itself is then made by the ordinary body generator (no special case, same classes and life rules as a planet); this module gives only its mass and origin.

Use:

    disk_moon_ceiling_kg(10 * JUPITER_MASS_KG)       # 1.2e25 kg, about two Earth masses
    exchange_chance(JUPITER_MASS_KG)                 # 0.014
    draw_exchange_moon(JUPITER_MASS_KG, rng)         # PlanetMoon("exchange", mass) or None
"""

import collections

from planetgen import tuning
from planetgen.physics import moon_composition
from planetgen.physics.constants import EARTH_MASS_TO_KG

JUPITER_MASS_KG = moon_composition.JUPITER_MASS_KG

PLANET_MOON_MIN_MASS_KG = 0.1 * EARTH_MASS_TO_KG
"""float: The lightest moon that counts as a moon of planet size: a tenth of an Earth mass (Mars is 0.107). [Stand-in PSM1, tier "from needed output", tweak L0.]"""

DISK_TOTAL_RATIO_MAX = moon_composition.total_moon_mass_ratio(moon_composition.ALPHA_OVER_F[1])
"""float: The most the moons of one circumplanetary disk can hold of their parent: Canup and Ward's total at the top of the alpha/f range, 6.4e-4 ([R0643: Eq. 4.15]). [Stand-in PSM2, tier T1, tweak L0.]"""

CAPTURE_CHANCE = 1.0 / 70.0
"""float: The share of giants that hold a moon of planet size from exchange capture: one candidate, Kepler-1708 b-i, in the 70 cool giants of the survey [R2365]. A candidate, not a confirmed moon, and
the survey found no more, so this is a loose upper figure. [Stand-in PSM3, tier "marginal unverified", tweak L2.]"""

EXCHANGE_MAX_RATIO = 1.0e-2
"""float: The most an exchange-captured moon holds of its parent: the candidate Kepler-1625 b-i is about 5.4e-3 of its planet (17 Earth masses on 10 Jupiter masses, as restated in [R5229]), the
Earth and Moon 1.2e-2; the binary the planet meets sets the real limit (the bound of Eq. 1 in [R5229], not rebuilt here). [Stand-in PSM4, tier "educated guess from physics", tweak L3.]"""

PlanetMoon = collections.namedtuple("PlanetMoon", ("origin", "mass_kg"))
"""A moon of planet size: its `origin` ("exchange") and its mass in kg."""


def disk_moon_ceiling_kg(parent_mass_kg):
    """The most mass the moons of a circumplanetary disk can hold between them, kg: `DISK_TOTAL_RATIO_MAX` of the parent. No single disk moon is heavier.

    Raises:
        ValueError: The mass is not positive.
    """
    if not parent_mass_kg > 0.0:
        raise ValueError(f"the parent mass {parent_mass_kg!r} is not positive")
    return DISK_TOTAL_RATIO_MAX * parent_mass_kg


def disk_parent_mass_kg(moon_mass_kg=PLANET_MOON_MIN_MASS_KG):
    """The lightest parent whose disk can make a moon of `moon_mass_kg`: the mass at which `disk_moon_ceiling_kg` reaches it (about 0.5 Jupiter masses for a tenth of an Earth mass, and
    in practice several, since the moons of one disk share the total). Pure."""
    if not moon_mass_kg > 0.0:
        raise ValueError(f"the moon mass {moon_mass_kg!r} is not positive")
    return moon_mass_kg / DISK_TOTAL_RATIO_MAX


def exchange_chance(parent_mass_kg):
    """The chance that a giant of this mass holds a moon of planet size from exchange capture: `CAPTURE_CHANCE`, or 0 where the ceiling (`exchange_ceiling_kg`) is under
    `PLANET_MOON_MIN_MASS_KG`, a parent too light to hold one (under about 6 Earth masses)."""
    if not parent_mass_kg > 0.0:
        raise ValueError(f"the parent mass {parent_mass_kg!r} is not positive")
    return CAPTURE_CHANCE if exchange_ceiling_kg(parent_mass_kg) > PLANET_MOON_MIN_MASS_KG else 0.0


def exchange_ceiling_kg(parent_mass_kg):
    """The heaviest exchange-captured moon a parent holds, kg: `EXCHANGE_MAX_RATIO` of the parent."""
    if not parent_mass_kg > 0.0:
        raise ValueError(f"the parent mass {parent_mass_kg!r} is not positive")
    return EXCHANGE_MAX_RATIO * parent_mass_kg


def exchange_mass_kg(parent_mass_kg, u):
    """
    The mass of an exchange-captured moon, kg, from the uniform `u` in 0..1: `dN/dlogM ~ M^-s` (`tuning.MOON_MASS_FUNCTION_SLOPE`, the generator's own moon-mass law) between
    `PLANET_MOON_MIN_MASS_KG` and `exchange_ceiling_kg`, by its inverse distribution.

    Raises:
        ValueError: The parent holds no moon that heavy, or `u` is outside 0..1.
    """
    ceiling = exchange_ceiling_kg(parent_mass_kg)
    if not 0.0 <= u <= 1.0:
        raise ValueError(f"the uniform {u!r} is outside 0..1")
    if ceiling <= PLANET_MOON_MIN_MASS_KG:
        raise ValueError(f"a parent of {parent_mass_kg!r} kg holds no moon of {PLANET_MOON_MIN_MASS_KG!r} kg")
    slope = tuning.MOON_MASS_FUNCTION_SLOPE
    low_term, high_term = PLANET_MOON_MIN_MASS_KG ** -slope, ceiling ** -slope
    return (low_term + u * (high_term - low_term)) ** (-1.0 / slope)


def draw_exchange_moon(parent_mass_kg, rng):
    """
    The moon of planet size a giant of `parent_mass_kg` holds by exchange capture, or `None`. Always two numbers from `rng` (anything with `random()`): the roll against `exchange_chance`,
    then the mass, so the stream never moves with the parent.

    Returns:
        PlanetMoon or None: The moon, or `None` for no capture (including a parent too light to hold one).
    """
    roll, mass_roll = rng.random(), rng.random()
    if roll >= exchange_chance(parent_mass_kg):
        return None
    return PlanetMoon("exchange", exchange_mass_kg(parent_mass_kg, mass_roll))
