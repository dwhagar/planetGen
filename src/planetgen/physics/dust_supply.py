# planetgen/physics/dust_supply.py

"""
The cosmic dust a planet's air takes in, and the carbon it carries, for the V8 species of the sulfur-cloud set (GEN.1542 slice 0)
=================================================================================================================================

V8 of the Venus sulfur-cloud set (`sulfur_cloud_set`, iron compounds of carbon rings in the acid cloud) needs carbon from cosmic dust. Until now the chain was told by its caller whether a body "has"
dust (`carbon_dust`), and no caller knew. A planet's own dust supply is what the belts and comets of its system send it (`influx.planet_influx_kg_per_yr`, GEN.1543: a thin drifting disk, the planet's focused
cross-section, less what a giant in the way ejects). This module turns that mass a year into the one number the gate reads, the carbon the dust lays on each square meter of the planet a year:

    carbon_flux = influx x CARBON_SHARE_OF_DUST / (4 pi R^2)

and says whether it reaches the floor `V8_MIN_CARBON_FLUX_KG_M2_YR` the layer needs to be kept up against its photochemical loss.

Pure functions: no body is read and nothing is stored. The wiring (reading a planet's reservoirs and giants at its fill, and passing `carbon_dust` to `sulfur_cloud_set.chemistry`) is slice 1 of GEN.1542.

What it leaves out (stated): the dust a planet lifts itself (surface dust and sand from composition and wind, volcanic and impact dust from the activity and age of the body, the air's winds). A body has no
volcanic rate, dust cover or wind field stored to read, and the layer sits in the cloud, far above where such dust is lifted on a thick-aired world like Venus; it is the row DUS3 of the gap plan.
"""

import math

CARBON_SHARE_OF_DUST = 0.1
"""float: The share of the infalling dust's mass that is carbon the cloud can take up (organic and carbonaceous grains: the carbon-rich dust of comets and primitive asteroids is of this order, the
mineral part of the dust carries none). A stand-in with no source in the registry. [Stand-in 6.306 DUS1, tier T6, tweak L4.]"""

V8_MIN_CARBON_FLUX_KG_M2_YR = 1.0e-9
"""float: The carbon, kg per square meter a year, the dust must lay on a planet for its cloud to keep the V8 layer against the layer's photochemical loss. Set below the flux a Sun-like system's
main-belt-sized reservoir gives a body at 0.7 to 1 AU (about 3e-9) and above a Mars-like body's (about 9e-10), so a body takes it where the dust is at least that of an inner planet of such a system. No source:
nobody has measured the loss rate of the layer. [Stand-in 6.306 DUS2, tier T6, tweak L4.]"""


def carbon_flux_kg_m2_yr(influx_kg_yr, radius_km):
    """The carbon the dust of `influx_kg_yr` (kg of dust a year, `influx.planet_influx_kg_per_yr`) lays on the sphere of a planet of `radius_km`, kg per square meter a year; 0 for no dust or no radius."""
    if influx_kg_yr <= 0.0 or radius_km <= 0.0:
        return 0.0
    area_m2 = 4.0 * math.pi * (radius_km * 1000.0) ** 2
    return CARBON_SHARE_OF_DUST * influx_kg_yr / area_m2


def v8_supplied(influx_kg_yr, radius_km):
    """Whether the dust of `influx_kg_yr` kg a year gives a planet of `radius_km` the carbon V8 needs: `carbon_flux_kg_m2_yr` at or above `V8_MIN_CARBON_FLUX_KG_M2_YR`."""
    return carbon_flux_kg_m2_yr(influx_kg_yr, radius_km) >= V8_MIN_CARBON_FLUX_KG_M2_YR
