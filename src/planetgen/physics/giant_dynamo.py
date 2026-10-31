"""
The dynamo region of a giant planet from its interior layers (GEN.809, the lock-free half; the wiring step has `magnetism.generate_field` take a giant's
moment from here in place of its three anchors).

A pure module: a `giant_interior.GiantInterior`, the body's radius, the heat it radiates and its spin in; the region that carries the dynamo and the
`dynamo.Field` of that region out. It reads nothing from a body and draws nothing.

1. **Which layer.** A gas giant's field is made in its metallic layer (`metallic_envelope`: hydrogen above the metallization pressure, convecting and
   conducting), from the radius of metallization down to the helium gradient below it. The layers under that, the helium gradient and shell the
   rain leaves and the dilute core (Saturn's stable region, Mankovich and Fuller 2021; Jupiter's, Juno), are stably stratified and carry no
   convection, so they are not part of the region. An ice giant's field is made in the conducting water of its mantle (`ice_mantle`, ionic or
   superionic), a shell over the rock core and under the molecular envelope (Stanley and Bloxham 2004, 2006 for the shell). A brown dwarf is a gas
   giant here.
2. **The region.** Its outer radius is the top of that layer, its inner radius the bottom, its mass the layer's, its density the mass over the shell's
   volume. They are the numbers `dynamo.field` takes (the buoyancy-flux scaling of Olson and Christensen 2006 and the Rossby number that gates the
   dipole).
3. **The power.** `dynamo.giant_power_w` of the heat the giant radiates from its own interior (the stand-in DY2).

The radii of the layers are the structure's (`giant_structure`, GEN.733): an ice giant's layers carry radii only once that is wired, and a body
whose conducting layer has none (the radii it is born with, the whole body) has no region here.

Every number that is not a published value is a stand-in in the gap registry (docs/design/research/gap-filling-plan.md section 6.30).
"""

import math
from dataclasses import dataclass

from planetgen.physics import dynamo, unit_registry as units

_mearth_to_kg = units.converter("Mearth", "kg")

CONVECTING_LAYER = {"gas_giant": "metallic_envelope", "brown_dwarf": "metallic_envelope", "ice_giant": "ice_mantle"}
"""dict: The layer that carries the dynamo of each kind of interior (`giant_interior.INTERIOR_KINDS`)."""


@dataclass(frozen=True)
class Region:
    """
    The shell of a giant that carries its dynamo.

    Attributes:
        layer (str): The layer's kind.
        outer_m, inner_m (float): Its outer and inner radius, m.
        mass_kg (float): Its mass.
    """

    layer: str
    outer_m: float
    inner_m: float
    mass_kg: float

    @property
    def thickness_m(self):
        return self.outer_m - self.inner_m

    @property
    def density_kg_m3(self):
        return self.mass_kg / dynamo.region_volume_m3(self.outer_m, self.inner_m)

    @property
    def inner_over_outer(self):
        """How thin the shell is: 0 for a full sphere, near 1 for a thin shell."""
        return self.inner_m / self.outer_m


def region_of(interior, radius_m):
    """
    The `Region` that carries a giant's dynamo, or `None` when it has none that is placed (no such layer, or one with no radii of its own).

    Args:
        interior: The body's `giant_interior.GiantInterior`.
        radius_m: The body's radius at its 1 bar level, m.
    """
    layer = interior.layer(CONVECTING_LAYER[interior.kind])
    if layer is None or radius_m <= 0.0:
        return None
    top, bottom = layer.top_radius_fraction, layer.bottom_radius_fraction
    if (top == 1.0 and bottom == 0.0) or not 0.0 <= bottom < top <= 1.0:
        return None  # the radii it is born with: no structure has placed it yet
    return Region(layer.kind, top * radius_m, bottom * radius_m, _mearth_to_kg(layer.mass_earth))


@dataclass(frozen=True)
class GiantField:
    """
    A giant's field.

    Attributes:
        region (Region): The shell that makes it.
        power_w (float): The power that drives it.
        field (dynamo.Field): Its moment, Rossby number and dipole share.
    """

    region: Region
    power_w: float
    field: dynamo.Field


def giant_field(interior, radius_m, radiated_w, spin_rad_s):
    """
    The `GiantField` of a giant, or `None` when it has no placed dynamo region.

    Args:
        interior: The body's `giant_interior.GiantInterior`.
        radius_m: Its radius at the 1 bar level, m.
        radiated_w: The heat it radiates from its own interior, W (`dynamo.giant_power_w` takes the share that drives the dynamo).
        spin_rad_s: Its spin rate, rad/s.
    """
    region = region_of(interior, radius_m)
    if region is None:
        return None
    power = dynamo.giant_power_w(radiated_w)
    return GiantField(region, power, dynamo.field(power, region.outer_m, region.inner_m, region.density_kg_m3, spin_rad_s))


def spin_rad_s(rotation_period_hours):
    """The spin rate of a body that turns once in `rotation_period_hours`, rad/s."""
    return 2.0 * math.pi / (rotation_period_hours * 3600.0)
