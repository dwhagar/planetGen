# planetgen/physics/captured_composition.py

"""
What a captured small body is made of, and why a family's members differ (GEN.1080; docs/design/systems/captured-asteroid-moons.md section 4.4)
===============================================================================================================================================

A body captured by a planet is a body of the reservoir it came from, so it is made of what that reservoir was made of **at the distance it came from**, not of a mix chosen for moons. This module asks the taxon kernel
(`taxon_kernel.kernel_at_distance`, GEN.473, through `belt_member.member`, GEN.1067) at the source, and then takes away what capture takes away:

1. **The source.** The window that brought the body says where it came from (`source_annulus`). An **encounter** or **exchange** capture takes it from the outer planetesimal disc the instability stirred up, whose limits the
   architecture gives (`irregular_population` draws the count; the disc of the system is the caller's). A **gas-era** capture takes it from the planet's own neighbourhood, the feeding zone `FEEDING_ZONE_HILL_RADII` Hill radii
   either side of its orbit, because the gas that slows a body is only dense near the planet. The body's place in the annulus follows the disc's solids, `a^-SURFACE_DENSITY_SLOPE` (`belt_history`), so a little more of the mass
   is near the inner edge (the mass up to `a` goes as `a^0.5`; `draw_source_distance`).
2. **The taxon and the grains** are `belt_member.member` at that distance and the body's mass: the dry, hydrated, primitive or ice-rich family from the temperature there, the differentiated taxa only if the star's aluminium-26 allowed
   them, the density from the compaction taper or the family's porosity. At 5 AU round a Sun-like star the kernel gives the primitive and ice-rich taxa (P and D); at 2.5 AU a mix with carbonaceous and rocky ones.
3. **What capture takes away.** The existing rule for a captured moon (`moon_composition.CAPTURE_ICE_LOSS`, a parametrisation; [R0269] qualitative) removes a share of the water ice, uniform from 0 to 0.5. The
   grains are renormalised; the porosity stays.
4. **A family shares a parent, but not one composition.** A body broken by a collision after capture leaves fragments that are not alike: the Himalia group's members differ in phyllosilicate content ([R4041], JWST
   abstract, read), and the parent was a heterogeneous rubble. A `Parent` is `PARENT_LUMPS` lumps, each from its own distance in the source annulus (so its own taxon), and the same capture (one ice loss for the whole parent). A
   fragment is one of the lumps (`fragment_lump`), with its own mass and porosity, so a family has the parent's taxa in the parent's proportions and its members differ from each other where the kernel's weights differ across
   the annulus and by the kernel's implanted share. A parent whose annulus is narrow and whose kernel is one-sided (D at 20 to 35 AU is 92 percent) gives a nearly uniform family, which is the physics asking for it.

A body over `compaction.TAPER_END_KG` (a massive exchange capture) is for the structure solver, as for a belt member; `compose` lets `compaction.OutOfRange` through. Pure functions: the random numbers come from the caller.
Nothing calls them yet; the captured-moon wiring (GEN.1083) does. Stand-ins are CP1 to CP4 of `docs/design/research/gap-filling-plan.md` section 6.119.
"""

import math
from dataclasses import dataclass

from planetgen.physics import belt_history, belt_member, irregular_zone, moon_composition, small_body_composition as sbc
from planetgen.util.checks import finite_domain

FEEDING_ZONE_HILL_RADII = 2.0 * math.sqrt(3.0)
"""float: The half-width of the zone a gas-era capture draws from, in the planet's Hill radii (the isolation spacing of planetesimal accretion, recalled, not read). [Stand-in CP1, tier T6, tweak L0.]"""
PARENT_LUMPS = 4
"""int: The lumps a family's parent is made of, each from its own distance (the Himalia group has four known members of different spectra). [Stand-in CP3, tier T6, tweak L3.]"""
WINDOWS = ("gas", "encounter", "exchange")
"""tuple: The windows whose bodies have a composition (`irregular_population`; a temporary capture has no object)."""


class CapturedCompositionError(ValueError):
    """An input the model cannot use."""


@dataclass(frozen=True)
class Captured:
    """
    One captured body's composition.

    Attributes:
        member (belt_member.Member): The body as it was at its source, before capture took its ice.
        window (str): The window that brought it.
        source_distance_au (float): Where it came from.
        ice_loss (float): The share of its water ice capture removed.
        makeup (sbc.SmallBodyMakeup): Its grains after the loss.
        density_kg_m3, radius_km (float): The bulk density (the grain density times the member's filling) and the radius of that mass.
    """

    member: object
    window: str
    source_distance_au: float
    ice_loss: float
    makeup: object
    density_kg_m3: float
    radius_km: float

    @property
    def taxon(self):
        return self.member.taxon


@dataclass(frozen=True)
class Lump:
    """A piece of a parent: the distance it came from and the uniform number that picks its taxon there."""

    distance_au: float
    u_taxon: float


@dataclass(frozen=True)
class Parent:
    """The parent of a family: its lumps and the uniform number of the ice loss its capture gave all of them."""

    lumps: tuple
    u_ice_loss: float
    window: str


@finite_domain()
def source_annulus(window, planet_distance_au, planet_mass_kg, star_mass_kg, disc_inner_au, disc_outer_au):
    """
    The `(inner, outer)` limits, AU, of the reservoir a window drew from: the outer disc's for an encounter or exchange capture, the planet's feeding zone for a gas-era one.

    Args:
        window (str): `gas`, `encounter` or `exchange`.
        planet_distance_au, planet_mass_kg, star_mass_kg (float): The planet's orbit and the masses.
        disc_inner_au, disc_outer_au (float): The planetesimal disc at the window.
    """
    if window not in WINDOWS:
        raise CapturedCompositionError(f"the window {window!r} is none of {WINDOWS}")
    if window == "gas":
        if not planet_distance_au > 0.0 or not planet_mass_kg > 0.0 or not star_mass_kg > 0.0:
            raise CapturedCompositionError("the distance and the masses must be positive")
        half = FEEDING_ZONE_HILL_RADII * irregular_zone.hill_radius_au(planet_distance_au, planet_mass_kg, star_mass_kg)
        return max(planet_distance_au - half, 0.1 * planet_distance_au), planet_distance_au + half
    if not 0.0 < disc_inner_au < disc_outer_au:
        raise CapturedCompositionError("the disc needs 0 < inner < outer")
    return disc_inner_au, disc_outer_au


@finite_domain()
def draw_source_distance(inner_au, outer_au, u):
    """A distance in the annulus drawn by the mass of its solids (`a^-SURFACE_DENSITY_SLOPE`, so the mass up to `a` goes as `a^(2 - slope)`), AU."""
    if not 0.0 < inner_au < outer_au or not 0.0 <= u <= 1.0:
        raise CapturedCompositionError("the annulus needs 0 < inner < outer and u from 0 to 1")
    k = 2.0 - belt_history.SURFACE_DENSITY_SLOPE
    return (inner_au ** k + u * (outer_au ** k - inner_au ** k)) ** (1.0 / k)


@finite_domain()
def ice_loss_share(u):
    """The share of its water ice a capture removes: `moon_composition.CAPTURE_ICE_LOSS` read at `u` (0 to 1)."""
    if not 0.0 <= u <= 1.0:
        raise CapturedCompositionError("u must be from 0 to 1")
    low, high = moon_composition.CAPTURE_ICE_LOSS
    return low + u * (high - low)


@finite_domain()
def compose(window, mass_kg, source_distance_au, luminosity_w, aluminium_26_ratio, u_taxon, u_porosity, u_ice_loss):
    """
    The composition of a body captured from `source_distance_au`.

    Args:
        window (str): `gas`, `encounter` or `exchange` (recorded).
        mass_kg (float): The body's mass.
        source_distance_au (float): Where it came from (`draw_source_distance`).
        luminosity_w, aluminium_26_ratio (float): The star's luminosity and aluminium-26 over aluminium-27 ratio (`belt_member.member`).
        u_taxon, u_porosity, u_ice_loss (float): Uniform numbers, 0 to 1: the taxon, the porosity of a body under `compaction.TAPER_START_KG`, the ice loss.

    Raises:
        compaction.OutOfRange: The mass is over `compaction.TAPER_END_KG`.
    """
    if window not in WINDOWS:
        raise CapturedCompositionError(f"the window {window!r} is none of {WINDOWS}")
    source = belt_member.member(u_taxon, u_porosity, mass_kg, source_distance_au, luminosity_w, aluminium_26_ratio)
    loss = ice_loss_share(u_ice_loss)
    components = dict(source.makeup.components)
    components["water_ice"] = components.get("water_ice", 0.0) * (1.0 - loss)
    makeup = sbc.makeup_from_components("moon", components, 1.0 - source.filling, detail={"origin": "captured", "window": window, "taxon": source.taxon,
                                                                                       "source_au": source_distance_au, "ice_loss": loss})
    density = makeup.grain_density_kg_m3 * source.filling
    radius_km = (3.0 * mass_kg / (4.0 * math.pi * density)) ** (1.0 / 3.0) / 1.0e3
    return Captured(source, window, source_distance_au, loss, makeup, density, radius_km)


@finite_domain()
def draw_parent(window, inner_au, outer_au, u_distances, u_taxa, u_ice_loss):
    """
    The `Parent` of a family: `PARENT_LUMPS` lumps, each at a distance drawn from the annulus (`draw_source_distance`) with its own taxon number, and one ice loss for the whole body.

    Args:
        u_distances, u_taxa (sequence): `PARENT_LUMPS` uniform numbers each.
    """
    if window not in WINDOWS:
        raise CapturedCompositionError(f"the window {window!r} is none of {WINDOWS}")
    if len(u_distances) != PARENT_LUMPS or len(u_taxa) != PARENT_LUMPS:
        raise CapturedCompositionError(f"a parent has {PARENT_LUMPS} lumps: give that many numbers")
    lumps = tuple(Lump(draw_source_distance(inner_au, outer_au, ud), ut) for ud, ut in zip(u_distances, u_taxa))
    return Parent(lumps, u_ice_loss, window)


def fragment_lump(parent, u):
    """The lump a fragment is made of, picked with equal weight by the uniform number `u` (0 to 1)."""
    if not 0.0 <= u <= 1.0:
        raise CapturedCompositionError("u must be from 0 to 1")
    return parent.lumps[min(int(u * len(parent.lumps)), len(parent.lumps) - 1)]


def fragment(parent, u_lump, mass_kg, luminosity_w, aluminium_26_ratio, u_porosity):
    """One fragment of a family: the composition of the lump `u_lump` picks, at the fragment's mass and with the parent's ice loss."""
    lump = fragment_lump(parent, u_lump)
    return compose(parent.window, mass_kg, lump.distance_au, luminosity_w, aluminium_26_ratio, lump.u_taxon, u_porosity, parent.u_ice_loss)
