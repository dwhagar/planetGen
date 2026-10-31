# planetgen/physics/rogue_mass_function.py

"""
The masses of the bodies a system throws out (GEN.1112, RP2)
============================================================

The rogue planets drew their mass from one fitted power law, `dN/dlogM ~ M^-0.65`, from 0.1 Earth masses to 13 Jupiter masses, cut into four bins by hand (`rogue-planet-review.md` rows 3 and 4). A
system does not eject a fitted law; it ejects what it has, and only what the scatterer can lift out. This module gives the ejected masses as the two things that decide them.

**What there is to eject.** The small bodies of a disc are a collisional cascade, and a steady-state cascade has a power-law size spectrum `dN/dD ~ D^-q` with `q = 3.5` for bodies held
by their strength (`n(m) ~ m^-11/6` as quoted by Belyaev and Rafikov 2011, arXiv:1012.3220 eq. 1, and `q = 3.5` by Pan and Schlichting 2012, arXiv:1111.0667 abstract), steepening to `q` up to 4 when the bodies' velocities evolve with the cascade (Pan and Schlichting 2012).
For a body of density that does not depend on size, `M ~ D^3` and the number per logarithm of mass is

    dN/dlogM ~ M^-(q - 1)/3,

so 0.83 at `q = 3.5` and 1.0 at `q = 4`. Sumi et al. 2023 measure `0.96 (+0.47 -0.27)` for the free-floating planets (arXiv:2303.08280, the short events of the MOA and OGLE surveys, 0.33 to 6,660 Earth masses): the supply
of a cascade sits inside the measurement, and no number was tuned to it.

**What the scatterer lifts out.** A body is thrown out by the planet it scatters off only if that planet has the Safronov number to eject (`rogue_kick.ejects`; below it the body is stirred), and
only if it is the lighter of the pair: the energy that unbinds the body comes out of the scatterer's orbit, so a body heavier than the scatterer is the one that stays (the architecture already throws out
the lighter of a crowded pair, `planet_architecture._scatter_crowded_giants`; Barclay et al. 2017 eject only embryos of at most a Mars mass to an Earth mass in the same runs). The ejected spectrum is
the supply's, cut at the scatterer's mass.

The normalisation (how many, how much mass) comes from the ejection record (RP1, GEN.1119), not from here; this module gives the shape, the share of the ejected bodies in a mass range, and a draw from it.
Nothing calls it yet (GEN.1111 creates the rogue rows). Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.110.
"""

import math

from planetgen.physics import rogue_kick
from planetgen.util.checks import finite_domain

SUPPLY_DIFFERENTIAL_INDEX = 3.5
"""float: `q` of the supply, `dN/dD ~ D^-q`: Dohnanyi's steady-state cascade of bodies held by their strength. Source tier T1, tweak level L0 (the velocity-evolved cascade of Pan and Schlichting 2012 gives up to 4)."""
MAX_EJECTED_OVER_SCATTERER = 1.0
"""float: The heaviest body a scatterer ejects, as a multiple of its own mass: the lighter of a pair is the one thrown out. Source tier T6, tweak level L0."""


class RogueMassFunctionError(ValueError):
    """An input the model cannot use."""


def supply_slope(q=SUPPLY_DIFFERENTIAL_INDEX):
    """The index `s` of `dN/dlogM ~ M^-s` for a cascade of differential size index `q` (bodies of one density): `(q - 1) / 3`."""
    if not q > 1.0:
        raise RogueMassFunctionError("a size spectrum index must be over 1")
    return (q - 1.0) / 3.0


def heaviest_ejected_kg(scatterer_mass_kg):
    """The mass above which a scatterer ejects nothing, kg."""
    if not scatterer_mass_kg > 0.0:
        raise RogueMassFunctionError("the scatterer's mass must be positive")
    return MAX_EJECTED_OVER_SCATTERER * scatterer_mass_kg


def can_eject(scatterer_mass_kg, scatterer_radius_m, host_mass_kg, distance_au):
    """Whether the scatterer ejects what it scatters at all (Safronov number over `rogue_kick.EJECTION_SAFRONOV_NUMBER`)."""
    return rogue_kick.ejects(scatterer_mass_kg, scatterer_radius_m, host_mass_kg, distance_au)


@finite_domain()
def mass_function(mass_kg, scatterer_mass_kg, q=SUPPLY_DIFFERENTIAL_INDEX):
    """
    The ejected bodies' `dN/dlogM` in units of the supply's at 1 kg (a shape, not a count): `M^-s` up to the scatterer's mass and 0 above it.

    Args:
        mass_kg (float): The body's mass.
        scatterer_mass_kg (float): The mass of the planet that ejects it.
        q (float): The supply's differential size index.
    """
    if not mass_kg > 0.0:
        raise RogueMassFunctionError("the body's mass must be positive")
    if mass_kg > heaviest_ejected_kg(scatterer_mass_kg):
        return 0.0
    return mass_kg ** -supply_slope(q)


def _range(mass_min_kg, scatterer_mass_kg):
    top = heaviest_ejected_kg(scatterer_mass_kg)
    if not 0.0 < mass_min_kg < top:
        raise RogueMassFunctionError("the lightest ejected mass must be positive and under the heaviest")
    return mass_min_kg, top


@finite_domain()
def number_share(low_kg, high_kg, mass_min_kg, scatterer_mass_kg, q=SUPPLY_DIFFERENTIAL_INDEX):
    """
    The share of the bodies a scatterer ejects that lie between `low_kg` and `high_kg`, 0 to 1.

    Args:
        low_kg, high_kg (float): The mass range asked about (clipped to what is ejected).
        mass_min_kg (float): The lightest body counted: the floor of the debris the system carries (the architecture's `PLANET_MIN_MASS_EARTH` for planets).
        scatterer_mass_kg (float): The scatterer's mass.
        q (float): The supply's differential size index.
    """
    lo, hi = _range(mass_min_kg, scatterer_mass_kg)
    a, b = max(low_kg, lo), min(high_kg, hi)
    if not a < b:
        return 0.0
    s = supply_slope(q)
    return (a ** -s - b ** -s) / (lo ** -s - hi ** -s)


@finite_domain()
def sample_mass_kg(u, mass_min_kg, scatterer_mass_kg, q=SUPPLY_DIFFERENTIAL_INDEX):
    """
    A mass drawn from the ejected spectrum by inverting its cumulative number, kg.

    Args:
        u (float): A uniform random number, 0 to 1.
    """
    if not 0.0 <= u <= 1.0:
        raise RogueMassFunctionError("u must be from 0 to 1")
    lo, hi = _range(mass_min_kg, scatterer_mass_kg)
    s = supply_slope(q)
    return (lo ** -s - u * (lo ** -s - hi ** -s)) ** (-1.0 / s)


@finite_domain()
def log_slope(mass_low_kg, mass_high_kg, scatterer_mass_kg, q=SUPPLY_DIFFERENTIAL_INDEX):
    """The index `alpha` of `dN/dlogM ~ M^-alpha` the ejected spectrum shows between two masses (compare Sumi et al. 2023's 0.96), 0 when either is above the scatterer's mass."""
    if not 0.0 < mass_low_kg < mass_high_kg:
        raise RogueMassFunctionError("the masses must be positive and in order")
    f_low = mass_function(mass_low_kg, scatterer_mass_kg, q)
    f_high = mass_function(mass_high_kg, scatterer_mass_kg, q)
    if f_low <= 0.0 or f_high <= 0.0:
        return 0.0
    return -math.log(f_high / f_low) / math.log(mass_high_kg / mass_low_kg)
