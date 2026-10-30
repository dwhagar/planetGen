# planetgen/physics/topup_field.py

"""
The light top-up asteroid field: what the Rogues level leaves unmatched becomes thin clouds of small rubble (GEN.1126; docs/design/systems/region-object-store.md section 8.5, last bullets)
===========================================================================================================================================================================

Stage 3 of a sector's generation (`rogue_topup.plan`) makes up the sector's mass shortfall with rogue planets, never more than the Rogues level allows, and hands back a `remainder_kg`.
Boss (5:41 and 5:42 PM EDT 2026-10-30): the remainder goes to interstellar asteroid fields in the same last pass, a distinct, lighter kind than a field with a real parent. This module is that
planning step as a pure function; nothing calls it yet (the region run does).

**What is lighter about it.** The bodies are smaller (`LARGEST_BODY_MASS_KG`, a hundredth of the cut of a field with a history), the rubble is more porous (`BULK_DENSITY_KG_M3`), and the field has
spread over a larger volume (it is placed at a radius of `MIN_RADIUS_LY` to a light year, and its age follows from that radius, `field_dispersal.age_yr_at_radius`). It has no parent ID, no body history and
no collision event: its **origin is the stand-in `unbound remnant`** (`ORIGIN`, registered, GEN.725), and it is given the same *apparent origin* and velocity a top-up rogue has
(`rogue_topup.apparent_history`: where the star that threw it would be now, marked apparent). It retires past a light year like every field (`field_dispersal.RETIREMENT_RADIUS_M`).

**How many and how heavy.** Fields are drawn one at a time, a mass log-uniform between `MIN_FIELD_MASS_KG` and `MAX_FIELD_MASS_KG`, until the remainder is used; a field that would overshoot is cut to
what is left, and the last of `MAX_FIELDS` takes everything that remains however large, so **the placed mass is the remainder exactly**. A remainder under `FLOOR_MASS_KG` (the ledger's mass cut,
below which nothing is stored) is returned as `unmatched_kg` and makes no field.

Random numbers come from the caller's stream (`rng.random()`), so the same stream gives the same plan. Stand-ins are in `docs/design/research/gap-filling-plan.md` section 6.108.
"""

import math
from dataclasses import dataclass

from planetgen import tuning
from planetgen.galaxy.age_law import AgeLaw
from planetgen.physics import belt_ledger, constants, field_dispersal, rogue_topup

ORIGIN = "unbound remnant"
"""str: The origin every top-up field carries: a stand-in, no collision event and no parent claimed. Source tier T9, tweak level L3."""
ORIGIN_APPARENT = rogue_topup.ORIGIN
"""str: The kind of the apparent origin (`apparent`): where the star that threw it would be now, not a system."""
LARGEST_BODY_MASS_KG = 1.0e17
"""float: The mass of the largest body a top-up field keeps (a rubble pile about 50 km across), a hundredth of the ledger's mass cut (1e19 kg). Source tier T9, tweak level L3."""
BULK_DENSITY_KG_M3 = 1500.0
"""float: The bulk density of the rubble, a porous pile, half the density of the silicate pieces of `field_dispersal.DEFAULT_DENSITY_KG_M3`. Source tier T9, tweak level L3."""
SPEED_KM_S_RANGE = (0.3, 3.0)
"""tuple: The rms speed the field spreads at, km/s, log-uniform (slow, so that it stays one cloud for hundreds of thousands of years). Source tier T9, tweak level L3."""
MIN_RADIUS_LY = 0.25
"""float: The least rms radius of a top-up field, light-years: it is drawn between this and the retirement radius of one light year, so it fills more volume than a field with a history. Source tier T9, tweak level L3."""
MIN_FIELD_MASS_KG = 1.0e21
"""float: The lightest field the drawing makes (a sixth of a thousandth of an Earth mass), kg. Source tier T9, tweak level L3."""
MAX_FIELD_MASS_KG = 1.2e24
"""float: The heaviest field drawn before the last one takes the rest (0.2 Earth masses: a collision field keeps 0.27 at most, IF6), kg. Source tier T9, tweak level L3."""
MAX_FIELDS = 12
"""int: The most fields in one sector's top-up; the last takes all that remains. Source tier T9, tweak level L3."""
FLOOR_MASS_KG = belt_ledger.L1_MASS_KG
"""float: The ledger's mass cut (`belt_ledger.L1_MASS_KG`, row BR17): a remainder under this makes no field, kg."""

_LY = constants.LIGHTYEAR_M


class TopUpFieldError(ValueError):
    """An input the planner cannot use."""


@dataclass(frozen=True)
class TopUpField:
    """
    One light field the top-up places.

    Attributes:
        mass_kg (float): Its mass.
        largest_body_mass_kg (float): The mass of its largest body (`LARGEST_BODY_MASS_KG`, or the whole field when that is lighter).
        bulk_density_kg_m3 (float): The density of its rubble.
        rms_speed_ms (float): The speed it spreads at.
        age_yr (float): Its age, from its radius (`field_dispersal.age_yr_at_radius`).
        radius_m (float): Its rms radius, 0.25 to 1 light year.
        position_pc (tuple): Where its centre is in the sector, pc from the sector's corner.
        velocity_km_s (tuple): Its velocity relative to the local standard of rest `(U, V, W)`.
        apparent_origin (rogue_topup.ApparentOrigin): Where it seems to have come from.
        origin (str): Always `ORIGIN`.
        shape (str): Always `ball`.
    """

    mass_kg: float
    largest_body_mass_kg: float
    bulk_density_kg_m3: float
    rms_speed_ms: float
    age_yr: float
    radius_m: float
    position_pc: tuple
    velocity_km_s: tuple
    apparent_origin: rogue_topup.ApparentOrigin
    origin: str = ORIGIN
    shape: str = "ball"


@dataclass(frozen=True)
class FieldTopUp:
    """
    The plan for a sector.

    Attributes:
        fields (tuple): The `TopUpField`s.
        placed_mass_kg (float): Their mass: the remainder, less `unmatched_kg`.
        unmatched_kg (float): The remainder when it is under `FLOOR_MASS_KG` (nothing is stored below the cut), else 0.
    """

    fields: tuple
    placed_mass_kg: float
    unmatched_kg: float


def draw_mass_kg(rng):
    """A field mass, log-uniform between `MIN_FIELD_MASS_KG` and `MAX_FIELD_MASS_KG`."""
    return MIN_FIELD_MASS_KG * (MAX_FIELD_MASS_KG / MIN_FIELD_MASS_KG) ** rng.random()


def draw_speed_ms(rng):
    """A spread speed, m/s, log-uniform over `SPEED_KM_S_RANGE`."""
    low, high = SPEED_KM_S_RANGE
    return 1000.0 * low * (high / low) ** rng.random()


def draw_radius_age(rng, speed_ms):
    """
    `(radius_m, age_yr)` of a field spreading at `speed_ms`: the radius uniform between `MIN_RADIUS_LY` and the retirement radius (one light year), the age the one that gives it
    (`field_dispersal.age_yr_at_radius`), so the step is reversible and the field is never past its retirement.
    """
    low = MIN_RADIUS_LY * _LY
    radius = low + rng.random() * (field_dispersal.RETIREMENT_RADIUS_M - low)
    return radius, field_dispersal.age_yr_at_radius(radius, speed_ms)


def _field(mass_kg, rng, sector_edge_pc, age_law, dispersion, placement, radius_pc):
    position = rogue_topup.place_in_sector(rng, sector_edge_pc, placement)
    speed = draw_speed_ms(rng)
    radius, age = draw_radius_age(rng, speed)
    velocity, origin = rogue_topup.apparent_history(rng, age_law, dispersion, radius_pc)
    return TopUpField(mass_kg, min(LARGEST_BODY_MASS_KG, mass_kg), BULK_DENSITY_KG_M3, speed, age, radius, position, velocity, origin)


def plan(rng, remainder_kg, host_masses_kg, sector_edge_pc=tuning.DEFAULT_SECTOR_EDGE_PC, age_law=rogue_topup.DEFAULT_AGE_LAW, dispersion_km_s=rogue_topup.VELOCITY_DISPERSION_KM_S,
         placement=None, galactocentric_radius_pc=rogue_topup.GALACTOCENTRIC_RADIUS_PC):
    """
    The light fields of one sector's top-up.

    Args:
        rng: A stream with `random()`, 0 up to 1 (the caller's `draw.py` stream).
        remainder_kg (float): The mass the rogues left unmatched (`rogue_topup.TopUp.remainder_kg`).
        host_masses_kg (sequence): The masses of the sector's own host stars, whose systems donate the mass (at least one when anything is to be placed).
        sector_edge_pc (float): The sector's edge.
        age_law (AgeLaw): The star-formation history's age law for the sector's population.
        dispersion_km_s (tuple): The local hosts' velocity dispersion.
        placement (float, optional): Override of `tuning.ROGUE_TOPUP_PLACEMENT`.
        galactocentric_radius_pc (float): The sector's distance from the galactic centre, for the arc the apparent origin lies along.

    Returns:
        FieldTopUp: The fields, the mass placed and the mass under the cut.
    """
    if not remainder_kg >= 0.0 or not math.isfinite(remainder_kg):
        raise TopUpFieldError("the remainder must be a finite mass, not negative")
    if remainder_kg < FLOOR_MASS_KG:
        return FieldTopUp((), 0.0, remainder_kg)
    if not host_masses_kg:
        raise TopUpFieldError("a top-up needs at least one host star to donate the mass")
    fields, remaining = [], remainder_kg
    while remaining > 0.0:
        mass = draw_mass_kg(rng) if len(fields) < MAX_FIELDS - 1 else remaining
        if mass >= remaining or remaining - mass < FLOOR_MASS_KG:
            mass = remaining                                              # a rest under the cut joins the last field
        fields.append(_field(mass, rng, sector_edge_pc, age_law, dispersion_km_s, placement, galactocentric_radius_pc))
        remaining -= mass
    return FieldTopUp(tuple(fields), math.fsum(f.mass_kg for f in fields), 0.0)
