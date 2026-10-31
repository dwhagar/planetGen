# planetgen/physics/dwarf_population.py

"""
The dwarf planets of a belt or an outer disc: how many there are, which are round, which have cleared their orbit, and which are Class D (GEN.1090; docs/design/systems/belt-origins.md section 5, belt-evolution.md section 3.2)
============================================================================================================================================================================================================================

The census of GEN.1087 (`reports/dwarf-census-2026-10-30.md`) found no free-orbiting dwarf planet in any generated system: a belt was an annulus with nothing in it, and the planet list stops at 0.1 Earth masses. The
belt ledger (`belt_ledger.py`, GEN.1061) draws the heaviest bodies of a belt from its mass and size law; this module says what they are. Pure functions, no wiring (GEN.1069 stores the ledger and GEN.1071 and GEN.1092 put it
on the outer disc; the system page lists the dwarf planets as members of their belt).

1. **How many.** `count_above_km` is `N(>D)` of the belt's size law for the belt's mass, `typical_largest` the body with `N(>m) = 1` (the order statistic) and `chance_of_body_above` the chance that a belt holds any body above a mass
   (`1 - exp(-N)`, a Poisson process). Nothing is typed in per belt: the count follows from the mass.
2. **Round.** A body is a dwarf planet candidate when its mass is above the mass of a fully relaxed hydrostatic body of its composition: `ROUND_MASS_ICE_KG` for ice and `ROUND_MASS_ROCK_KG` for rock, the geometric centres of the
   dwarf-planet bands (1e20 to 1e21 kg for ice, 1e21 to 1e22 kg for rock, `belt-sources.md` BR9), mixed log-linearly in the ice fraction up to `ICE_BODY_FRACTION` (`round_mass_kg`).
3. **Cleared.** The body is a planet if it has cleared its orbit, a dwarf planet if not. Soter's `mu = M / m` (the body over the mass of everything else in its zone, here the rest of the belt) above `SOTER_MU_PLANET` (100) is a
   planet (Soter 2006, READ); Margot's `Pi = M / M_clear` with `M_clear = 1.9e-4 C^(3/2) (a / au)^(9/8) (M* / Msun)^(5/2)` Earth masses and `C = 2 sqrt 3` (Margot 2015, READ) is the cross-check and is reported with every
   member (`Member.margot_pi`; Jupiter 4.0e4, Mars 54, as the paper prints). `status` is `minor body` below the round mass, `dwarf planet` if round and not cleared, `planet` if round and cleared.
4. **Class D.** The class is the conditions rule of `physics/data/planet_classes.json` run on the body, not a table here: a body of a belt that holds no air (`AIRLESS_MAX_MASS_KG`, the mass under which the generator already
   treats a body as dwarf-sized, 5e22 kg) and at least 1 percent ice is Class D (`class_letter`). Bodies below `belt_ledger.L1_MASS_KG` are not stored and have no letter; the other letters a rocky member could take
   (C, A) need a temperature and a gravity the wiring step has.

The Solar System's main belt (Ceres a dwarf planet, Soter's mu about 0.6) and the Kuiper belt (Pluto, Eris, Haumea, Makemake, Quaoar and Gonggong as dwarf candidates, 15 to 20 in the design note) are test targets, never inputs.
In `tests/test_dwarf_population.py` the model's median counts are within a factor of three of them; Ceres (9.4e20 kg, 39 percent of the belt) is a rare draw, not the median, because a power law's order statistic puts
the top body at a tenth of Ceres at the main belt's mass.

Use:

    count_above_km(belt_ledger.ROCKY_LAW, belt_mass_kg, 100.0)                           # about how many bodies of 100 km and up
    members(ledger, ice_fraction, distance_au, star_mass_solar)                          # every stored body: its rank, size, status and class
    dwarf_planets(members(...))                                                          # only the dwarf planets
"""

import math
from dataclasses import dataclass

from planetgen.physics import belt_ledger, moon_composition, planet_classes, unit_registry as units
from planetgen.util.checks import finite_domain

_kg_to_mearth = units.converter("kg", "Mearth")

ROUND_MASS_ICE_KG = 3.2e20
"""float: The mass above which an ice body is in hydrostatic equilibrium: the geometric centre of the band 1e20 to 1e21 kg. [Stand-in DWP1, tier T4, tweak L1.]"""

ROUND_MASS_ROCK_KG = 3.2e21
"""float: The same for a rock body: the geometric centre of 1e21 to 1e22 kg. [Stand-in DWP1, tier T4, tweak L1.]"""

ICE_BODY_FRACTION = 0.5
"""float: The ice mass fraction at and above which a body is as weak as pure ice (the round mass is the ice value); between zero and this it is mixed log-linearly. [Stand-in DWP2, tier T6, tweak L1.]"""

SOTER_MU_PLANET = 100.0
"""float: Soter's discriminant above which a body has cleared its zone and is a planet (Soter 2006, arXiv:astro-ph/0608359, READ; the gap between planets and the rest is four to five orders of magnitude). [Stand-in BR6, tier T1, tweak L0.]"""

MARGOT_COEFFICIENT_EARTH = 1.9e-4
"""float: Margot's orbit-clearing mass at 1 au round a Sun-like star for `C = 1`, Earth masses (Margot 2015, arXiv:1507.06300, READ). [Stand-in BR7, tier T1, tweak L0.]"""

MARGOT_C = 2.0 * math.sqrt(3.0)
"""float: Margot's constant `C = 2 sqrt 3` (the clearing zone of 2 sqrt 3 Hill radii over the main-sequence lifetime, 10 Gyr; printed in the paper). [Stand-in BR7, tier T1, tweak L0.]"""

AIRLESS_MAX_MASS_KG = moon_composition.DWARF_PLANET_MAX_MASS_KG
"""float: A belt member under this mass holds no air; the generator's existing upper bound of a dwarf body. [Registered with `moon_composition.DWARF_PLANET_MAX_MASS_KG`.]"""

MINOR = "minor body"
DWARF = "dwarf planet"
PLANET = "planet"


@finite_domain()
def count_above_km(law, belt_mass_kg, diameter_km):
    """`N(>D)`: the expected number of bodies of this diameter and up in a belt of this mass whose bodies follow `law`."""
    if not diameter_km > 0.0:
        raise ValueError("the diameter must be positive")
    mass = math.pi / 6.0 * (diameter_km * 1000.0) ** 3 * law.density_kg_m3
    return belt_ledger.count_above(law, belt_mass_kg, mass)


@finite_domain()
def typical_largest(law, belt_mass_kg):
    """The body with `N(>m) = 1` (the order statistic of the belt's size law): its mass in kg and its diameter in km."""
    mass = belt_ledger.mass_at_count(law, belt_mass_kg, 1.0)
    return mass, belt_ledger.diameter_km(mass, law.density_kg_m3)


@finite_domain()
def chance_of_body_above(law, belt_mass_kg, mass_kg):
    """The chance that a belt holds at least one body above `mass_kg`: the bodies are a Poisson process in the size law, so `1 - exp(-N(>m))`."""
    return 1.0 - math.exp(-belt_ledger.count_above(law, belt_mass_kg, mass_kg))


@finite_domain()
def round_mass_kg(ice_fraction):
    """The mass above which a body of this ice mass fraction is in hydrostatic equilibrium (ice and rock mixed log-linearly in the fraction, up to `ICE_BODY_FRACTION`)."""
    if not 0.0 <= ice_fraction <= 1.0:
        raise ValueError("the ice fraction is a share between 0 and 1")
    weight = min(1.0, ice_fraction / ICE_BODY_FRACTION)
    return math.exp((1.0 - weight) * math.log(ROUND_MASS_ROCK_KG) + weight * math.log(ROUND_MASS_ICE_KG))


@finite_domain()
def soter_mu(mass_kg, others_kg):
    """Soter's `mu`: the body's mass over the mass of everything else in its zone."""
    if not (mass_kg > 0.0 and others_kg > 0.0):
        raise ValueError("both masses must be positive")
    return mass_kg / others_kg


@finite_domain()
def margot_clearing_mass_earth(distance_au, star_mass_solar):
    """Margot's `M_clear`, Earth masses: the mass that clears its orbit in a main-sequence lifetime."""
    if not (distance_au > 0.0 and star_mass_solar > 0.0):
        raise ValueError("the distance and the star's mass must be positive")
    return MARGOT_COEFFICIENT_EARTH * MARGOT_C ** 1.5 * distance_au ** (9.0 / 8.0) * star_mass_solar ** 2.5


@finite_domain()
def margot_pi(mass_kg, distance_au, star_mass_solar):
    """Margot's discriminant `Pi = M / M_clear`; a planet has `Pi >= 1`."""
    return _kg_to_mearth(mass_kg) / margot_clearing_mass_earth(distance_au, star_mass_solar)


@finite_domain()
def status(mass_kg, ice_fraction, others_kg):
    """`minor body` below the round mass, else `planet` if it has cleared its zone (Soter's mu over `SOTER_MU_PLANET`) and `dwarf planet` if not."""
    if mass_kg < round_mass_kg(ice_fraction):
        return MINOR
    return PLANET if soter_mu(mass_kg, others_kg) > SOTER_MU_PLANET else DWARF


def class_letter(mass_kg, ice_fraction):
    """
    `D` for a stored belt member that holds no air and is at least 1 percent ice, by the conditions rule of the planet classes; `None` for a body under the ledger's cut (not stored) or over
    `AIRLESS_MAX_MASS_KG` (it has its own pipeline), or for a rocky one (its letter needs the temperature and gravity the wiring step has).
    """
    if not (belt_ledger.L1_MASS_KG <= mass_kg < AIRLESS_MAX_MASS_KG):
        return None
    letter = planet_classes.classify_values({"body_type": planet_classes.TERRESTRIAL, "mass_kg": mass_kg, "atmospheric_pressure_pa": 0.0, "water_mass_fraction": ice_fraction}).letter
    return letter if letter == "D" else None


@dataclass(frozen=True)
class Member:
    """
    A stored body of a belt, said what it is.

    Attributes:
        rank (int): The ledger rank (0 is the heaviest), stable for the belt.
        mass_kg, diameter_km (float): Its mass and diameter.
        status (str): `minor body`, `dwarf planet` or `planet`.
        class_letter (str or None): `D` if it holds no air and is at least 1 percent ice, else `None`.
        soter_mu, margot_pi (float): The two clearing numbers.
    """

    rank: int
    mass_kg: float
    diameter_km: float
    status: str
    class_letter: object
    soter_mu: float
    margot_pi: float


def members(ledger, ice_fraction, distance_au, star_mass_solar):
    """
    Every stored body of a `belt_ledger.Ledger`, with its status and class.

    Args:
        ledger (belt_ledger.Ledger): The belt's stored bodies.
        ice_fraction (float or sequence): The ice mass fraction of the belt, or of each body by rank.
        distance_au (float): The belt's distance from the star.
        star_mass_solar (float): The star's mass.

    Returns:
        tuple of Member
    """
    out = []
    for rank, mass in enumerate(ledger.masses_kg):
        ice = ice_fraction[rank] if hasattr(ice_fraction, "__getitem__") else ice_fraction
        others = ledger.belt_mass_kg - mass
        out.append(Member(rank, mass, belt_ledger.diameter_km(mass, ledger.law.density_kg_m3), status(mass, ice, others), class_letter(mass, ice),
                          soter_mu(mass, others), margot_pi(mass, distance_au, star_mass_solar)))
    return tuple(out)


def dwarf_planets(listed):
    """The dwarf planets among some `members`."""
    return tuple(m for m in listed if m.status == DWARF)
