# planetgen/physics/belt_member.py

"""
What one stored member of a belt is made of (GEN.1067; docs/design/systems/belt-evolution.md section 3.1 step 5)
=================================================================================================================

The belt ledger (`belt_ledger`, GEN.1061) says how many bodies there are above each mass; this module says what each one is. A member's composition is not the belt's one mix. It is read at the member's own
semi-major axis from the taxon kernel (`taxon_kernel.kernel_at_distance`, GEN.473): the disk temperature there sets the water and ammonia ice lines, so the shares of the dry, hydrated, primitive and ice-rich
families change across the belt, with the inner edge drier and the outer edge icier, and the same star gives each member its place on that gradient.

1. **The taxon** (`pick_taxon`) is drawn from the kernel's weights at the member's distance by one uniform number. The differentiated taxa (M the core, A the mantle, V the crust) have weight only when the
   star's aluminium-26 ratio is above the melting ratio (`taxon_kernel.melting_ratio`, about 1.1e-5) and only in the dry family; a star born with less has none, and a belt round it is undifferentiated
   (`is_differentiated` is the member's flag).
2. **The grains** are the taxon's row of the shared table (`small_body_composition.TAXON_DATA`): its rock, metal, water and ice shares and the grain density they mix to.
3. **The density** is the grain density times the filling that the compaction taper gives a body of this mass (`compaction`, GEN.442): the tapered brittle law of Bierson and Nimmo closes the pores between
   `compaction.TAPER_START_KG` and `TAPER_END_KG`, so the filling of a body of the mass is the compacted density over the compaction's own grain density at the same rock fraction. A member above the taper's start is
   a stored body (the ledger's cut is 1e19 kg), so the taper decides it; a lighter one (drawn on demand) takes the taxon's family porosity (`small_body_composition.POROSITY_RANGES`) from a second uniform number. A
   member over `compaction.TAPER_END_KG` is for the structure solver (`compaction.OutOfRange`).
4. **The radius** follows from the mass and the density.

Pure functions over the ledger's draw: the random numbers (one for the taxon, one for the porosity of a light body) come from the caller. Nothing calls them yet (the belt-wide wiring is GEN.473). Stand-ins are in
`docs/design/research/gap-filling-plan.md` section 6.115.
"""

import math
from dataclasses import dataclass

from planetgen.physics import compaction, small_body_composition as sbc, taxon_kernel
from planetgen.util.checks import finite_domain


class BeltMemberError(ValueError):
    """An input the model cannot use."""


@dataclass(frozen=True)
class Member:
    """
    One belt member's composition.

    Attributes:
        taxon (str): The class letter.
        differentiated (bool): Whether it is a melted taxon (M, A or V).
        mass_kg, distance_au (float): As given.
        rock_fraction (float): Rock and metal over the grains, as `compaction` reads it.
        filling (float): One minus the porosity: the compaction taper's for a stored member, the family's drawn one for a light one.
        density_kg_m3, radius_km (float): The bulk density and the radius of the sphere of that mass.
        albedo (float): The taxon's mean albedo.
        makeup (sbc.SmallBodyMakeup): The grains and their groups.
    """

    taxon: str
    differentiated: bool
    mass_kg: float
    distance_au: float
    rock_fraction: float
    filling: float
    density_kg_m3: float
    radius_km: float
    albedo: float
    makeup: object


def shares(distance_au, luminosity_w, aluminium_26_ratio):
    """The taxon weights of a member at `distance_au` from a star of luminosity `luminosity_w` born with this aluminium-26 over aluminium-27 ratio (`taxon_kernel.kernel_at_distance`)."""
    return taxon_kernel.kernel_at_distance(distance_au, luminosity_w, aluminium_26_ratio)


def is_differentiated(taxon):
    """Whether a taxon is one of the melted ones (`taxon_kernel.IGNEOUS_SPLIT`: M, A, V)."""
    return taxon in taxon_kernel.IGNEOUS_SPLIT


def pick_taxon(u, weights):
    """The taxon a uniform number `u` (0 to 1) picks from the weights, in the order of `taxon_kernel.TAXA`."""
    if not 0.0 <= u <= 1.0:
        raise BeltMemberError("u must be from 0 to 1")
    total = math.fsum(weights.values())
    if not total > 0.0:
        raise BeltMemberError("the weights add to nothing")
    running, chosen = 0.0, None
    for taxon in taxon_kernel.TAXA:
        weight = weights.get(taxon, 0.0)
        if weight <= 0.0:
            continue
        chosen = taxon
        running += weight / total
        if u < running:
            return taxon
    return chosen


def _filling(makeup, mass_kg, taxon, u_porosity):
    low, high = compaction.ROCK_FRACTION_RANGE
    rock = min(max(makeup.rock + makeup.metal, low), high)
    if mass_kg < compaction.TAPER_START_KG:
        if not 0.0 <= u_porosity <= 1.0:
            raise BeltMemberError("u_porosity must be from 0 to 1")
        porosity_low, porosity_high = sbc.POROSITY_RANGES[sbc.TAXON_DATA[taxon].family]
        return rock, 1.0 - (porosity_low + u_porosity * (porosity_high - porosity_low))
    radius_km = compaction.compacted_radius_km(mass_kg, rock)
    compacted = compaction.compacted_density_kg_m3(radius_km, rock)
    reference_grain = 1.0 / (rock / compaction.SILICATE_GRAIN_KG_M3 + (1.0 - rock) / compaction.ICE_GRAIN_KG_M3)
    return rock, min(compacted / reference_grain, 1.0)


@finite_domain()
def typical_density_kg_m3(taxon, mass_kg):
    """
    The bulk density of a typical body of `taxon` at `mass_kg`, kg/m^3: the filling of `member` (the compaction taper's for a stored body, the middle of the family's porosity range for a light one) times the taxon's grain
    density. For the size table of a belt, which has no single body to draw.

    Raises:
        compaction.OutOfRange: The mass is over `compaction.TAPER_END_KG`.
    """
    if taxon not in sbc.TAXON_DATA or not mass_kg > 0.0:
        raise BeltMemberError("the taxon must be a known one and the mass positive")
    makeup = sbc.makeup_from_components("belt", sbc.TAXON_DATA[taxon].fractions, 0.0)
    _, filling = _filling(makeup, mass_kg, taxon, 0.5)
    return makeup.grain_density_kg_m3 * filling


@finite_domain()
def member(u_taxon, u_porosity, mass_kg, distance_au, luminosity_w, aluminium_26_ratio):
    """
    The composition, density and radius of one belt member.

    Args:
        u_taxon (float): A uniform number, 0 to 1, for the taxon.
        u_porosity (float): A uniform number, 0 to 1, for the porosity of a body under `compaction.TAPER_START_KG` (not read for a heavier one).
        mass_kg (float): The member's mass (the ledger's).
        distance_au (float): Its semi-major axis.
        luminosity_w (float): The star's luminosity.
        aluminium_26_ratio (float): The star's aluminium-26 over aluminium-27 at the birth of its solids.

    Raises:
        compaction.OutOfRange: The mass is over `compaction.TAPER_END_KG`.
    """
    if not mass_kg > 0.0 or not distance_au > 0.0 or not luminosity_w > 0.0:
        raise BeltMemberError("the mass, distance and luminosity must be positive")
    taxon = pick_taxon(u_taxon, shares(distance_au, luminosity_w, aluminium_26_ratio))
    data = sbc.TAXON_DATA[taxon]
    makeup = sbc.makeup_from_components("belt", data.fractions, 0.0)
    rock, filling = _filling(makeup, mass_kg, taxon, u_porosity)
    density = makeup.grain_density_kg_m3 * filling
    radius_km = (3.0 * mass_kg / (4.0 * math.pi * density)) ** (1.0 / 3.0) / 1.0e3
    return Member(taxon, is_differentiated(taxon), mass_kg, distance_au, rock, filling, density, radius_km, data.albedo, makeup)
