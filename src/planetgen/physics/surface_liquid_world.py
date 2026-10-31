# planetgen/physics/surface_liquid_world.py

"""
The liquid a world holds on its surface (GEN.1035 wiring; the table and the rules are `surface_liquids`)
=========================================================================================================

Sets `COLUMN_FIELDS` on a planet or moon from its surface temperature, pressure, air, ices, ground and hydrosphere (called by
`surface_composition.generate_surface`, after the hydrosphere and the ground). It draws nothing. A body has at most one surface liquid:

1. **Water.** A world whose hydrosphere is a liquid surface (`surface ocean`, `brine seep`, `hycean`) has `water`, covering `ocean_fraction`
   to the depth the hydrosphere stored. Every other state of water (dry, vapour, ice, ice-covered ocean) leaves the surface to the rules below: a
   frozen or buried water ocean is not a surface liquid.
2. **Lava.** A surface hot enough to be molten (at least the liquidus of the ground's rock: ultramafic ground melts to komatiite, felsic ground to
   rhyolite, the rest to basalt, and a surface above the peridotite liquidus a peridotitic magma ocean) is covered by it all over: the temperature is the
   body's mean, so a mean above the liquidus is a molten surface (`physics.mixture_lava` is the vapour over it). The depth of a magma ocean is not
   known (no rule yet). A lava surface holds in its own vapour, so its boiling window is checked at least at one bar.
3. **Seas of the volatile ices and the air.** Otherwise `surface_liquids.dominant` picks the liquid that is stable at the surface temperature and pressure
   and that the world supplies: a volatile of the air at its saturation pressure condenses, and the ices the body holds (`ice_<species>_mass_fraction`)
   supply their liquid. A methane sea of a cold world with nitrogen air is Titan's. The share of the ice that lies on the surface as a liquid is
   `SURFACE_SHARE_OF_ICE`, spread over the surface like water (`hydrosphere.ocean_fraction`'s hypsometry); the liquid that condenses out of an
   air over saturation lies at the mass the excess partial pressure weighs (`(p - p_sat) / g` per square metre).

Every number that is not a published value is a stand-in in the gap registry (docs/design/research/gap-filling-plan.md section 6.69, rows LQ11 to LQ14).
"""

import math

from planetgen.physics import constants, hydrosphere, ice_species, mixture, saturation, surface_liquids

COLUMN_FIELDS = ("surface_liquid", "surface_liquid_fraction", "surface_liquid_depth_km")
"""tuple: What a body stores of its surface liquid (columns of `planets` and `moons`): the liquid's id in `surface_liquids` (`None` for none), the
share of the surface it covers, and its mean depth over that share, km (`None` when not known)."""

WATER_ID = "water"
WATER_STATES = ("surface ocean", "brine seep", "hycean")
"""tuple: The hydrosphere states whose water is a liquid at the surface."""

SURFACE_SHARE_OF_ICE = 3.0e-4
"""float: STAND-IN LQ11 (L2): the share of a body's ice of a volatile that lies on the surface as its liquid. Tuned so that a Titan-like body with 1 percent of its
mass as methane ice (a guess, the ices are drawn) covers 1.7 percent of its surface, against the 1.6 percent the lakes and seas of Titan cover (Hayes 2016)."""

ICE_LIQUID = {"CH4": "methane", "N2": "nitrogen", "CO2": "carbon_dioxide", "CO": "carbon_monoxide", "HCN": "hydrogen_cyanide",
              "SO2": "sulfur_dioxide", "NH3": "ammonia"}
"""dict: The liquid each stored ice species supplies when it melts."""

LAVA_OF_GROUND = (("ultramafic", "komatiite_melt"), ("felsic", "rhyolite_melt"), ("basaltic", "basalt_melt"))
"""tuple: The lava a ground type of the surface melts to; every other ground (and a ground not drawn) melts to basalt. STAND-IN LQ12."""
LAVA_DEFAULT = "basalt_melt"
MAGMA_OCEAN = "peridotite_magma"
LAVA_IDS = (*(liquid for _ground, liquid in LAVA_OF_GROUND), LAVA_DEFAULT, MAGMA_OCEAN)
"""tuple: Every liquid id that is a lava; the depth of one is not known."""
LAVA_MIN_PRESSURE_PA = 1.0e5
"""float: A lava surface holds in its own vapour, so its window is checked at least at this pressure."""


def _liquid_water(body):
    body.surface_liquid = WATER_ID
    body.surface_liquid_fraction = body.ocean_fraction
    body.surface_liquid_depth_km = body.ocean_depth_km


def _lava(body, temperature_k, pressure_pa):
    """The molten surface of a body, or `None` when it is cooler than the liquidus of its ground's rock (a magma ocean over the peridotite liquidus)."""
    pressure = max(pressure_pa, LAVA_MIN_PRESSURE_PA)
    if surface_liquids.is_stable(MAGMA_OCEAN, temperature_k, pressure):
        return MAGMA_OCEAN
    lava = max(LAVA_OF_GROUND, key=lambda pair: getattr(body, f"ground_{pair[0]}_fraction", None) or 0.0)
    liquid = lava[1] if (getattr(body, f"ground_{lava[0]}_fraction", None) or 0.0) > 0.0 else LAVA_DEFAULT
    return liquid if surface_liquids.is_stable(liquid, temperature_k, pressure) else None


def _reservoir(body):
    """`{liquid id: mass fraction of the body}` of the ices that supply a liquid."""
    out = {}
    for formula, liquid in ICE_LIQUID.items():
        amount = getattr(body, ice_species.FIELD_OF[formula], None) if formula in ice_species.FIELD_OF else None
        if amount and amount > 0.0:
            out[liquid] = amount
    return out


def _sea(body, temperature_k, pressure_pa):
    """`(liquid id, surface mass per square metre in kg)` of the sea a cold or temperate body holds, or `None`."""
    partials = mixture.partial_pressures_pa(body)
    reservoir = _reservoir(body)
    pick = surface_liquids.dominant(temperature_k, pressure_pa, partials, reservoir)
    if pick is None:
        return None
    gravity = (body.gravity or 0.0) * constants.EARTH_GRAVITY
    area_m2 = 4.0 * math.pi * (body.radius * 1000.0) ** 2
    if pick.basis == "reservoir":
        return pick.liquid_id, SURFACE_SHARE_OF_ICE * reservoir[pick.liquid_id] * body.mass / area_m2
    row = surface_liquids.liquid(pick.liquid_id)
    excess = partials[row.saturation] - saturation.pressure_pa(row.saturation, temperature_k)
    return pick.liquid_id, max(excess, 0.0) / gravity if gravity > 0.0 else 0.0


def generate(body):
    """Sets `COLUMN_FIELDS` on a planet or moon (call after its hydrosphere and ground): `None` for a giant, a world with no liquid at its surface and a
    body saved before the columns."""
    for field in COLUMN_FIELDS:
        setattr(body, field, None)
    if body.body_type == "g":
        return
    if getattr(body, "hydrosphere", None) in WATER_STATES:
        _liquid_water(body)
        return
    temperature = body.surface_temperature or 0.0
    pressure = body.atmospheric_pressure or 0.0
    if not temperature > 0.0:
        return
    lava = _lava(body, temperature, pressure)
    if lava is not None:
        body.surface_liquid, body.surface_liquid_fraction = lava, 1.0
        return
    if not pressure > 0.0:
        return                              # in vacuum a liquid boils away
    found = _sea(body, temperature, pressure)
    if found is None:
        return
    liquid_id, mass_per_m2 = found
    gravity = (body.gravity or 0.0) * constants.EARTH_GRAVITY
    layer_km = mass_per_m2 / surface_liquids.density_kg_m3(liquid_id, temperature) / 1000.0
    covered = hydrosphere.ocean_fraction(layer_km, gravity)
    if not covered:
        return
    body.surface_liquid, body.surface_liquid_fraction = liquid_id, covered
    body.surface_liquid_depth_km = layer_km / covered
