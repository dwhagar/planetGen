# planetgen/physics/crust_porosity.py

"""
The porosity of a world's upper crust, from its gravity and what it is made of (GEN.780)
========================================================================================

A body of `compaction.TAPER_END_KG` or more has no macroporosity (the pores of the whole body are closed, decision S9), but its upper
crust is fractured and porous, and how deep the pores stay open depends on gravity: the pores close at a fixed pressure, so the
depth that holds that pressure is shorter where the gravity is greater. The usual form (Athy 1930; Clifford 1993 for the Martian
megaregolith, Hanna and Phillips 2005) is an exponential:

    phi(z) = phi0 exp(-z / K),    K = K_earth * g_earth / g

with `phi0` the porosity at the surface and `K_earth` the scale depth of the Earth's crust. The value stored (`porosity`) is the mean
over the top `REFERENCE_DEPTH_KM`, which is the megaregolith the Moon and Mars are measured for:

    mean = phi0 K / D (1 - exp(-D / K))

`phi0` follows the composition: `SURFACE_POROSITY_ROCK` for a rocky crust and `SURFACE_POROSITY_ICE` for an icy one, by the share
of the body's mass that is water and ice (a mix in between by that share). With the numbers below the Earth's upper 10 km hold 5.5
percent (observed upper-crust porosity 2 to 10 percent), the Moon 15 (GRAIL: 12 percent of the whole crust, more in the megaregolith)
and Mars 11 (10 to 20 percent). Stand-ins are the registry rows CP1 to CP3 (`docs/design/research/gap-filling-plan.md`).

A body lighter than `compaction.TAPER_END_KG` is a heap of grains and keeps its macroporosity (`bulk_draw`), not this.
"""

import math

from planetgen.physics import constants

SURFACE_POROSITY_ROCK = 0.20
"""float: The porosity of a rocky crust at the surface (Clifford 1993, Hanna and Phillips 2005: 0.2 for basaltic megaregolith)."""

SURFACE_POROSITY_ICE = 0.30
"""float: The porosity of an icy crust at the surface (a stand-in: Bierson and Nimmo's ice law starts at 0.6 and closes to 0.2 by 10 MPa)."""

EARTH_SCALE_DEPTH_KM = 2.82
"""float: The scale depth `K` of the Earth's crust (Clifford 1993, from Athy and Davis)."""

REFERENCE_DEPTH_KM = 10.0
"""float: The depth the stored mean is taken over: the Moon's and Mars' megaregolith."""

EARTH_GRAVITY_M_S2 = 9.80665


def scale_depth_km(gravity_m_s2):
    """`K` of a crust under this gravity: the Earth's, scaled by the gravity."""
    return EARTH_SCALE_DEPTH_KM * EARTH_GRAVITY_M_S2 / gravity_m_s2


def surface_porosity(ice_share):
    """`phi0` for a body whose mass is `ice_share` (0 to 1) water and ice."""
    share = min(max(ice_share, 0.0), 1.0)
    return SURFACE_POROSITY_ROCK + (SURFACE_POROSITY_ICE - SURFACE_POROSITY_ROCK) * share


def mean_porosity(gravity_m_s2, ice_share=0.0):
    """The mean porosity over the top `REFERENCE_DEPTH_KM` of the crust of a body of this gravity and ice share."""
    depth = scale_depth_km(gravity_m_s2)
    return surface_porosity(ice_share) * depth / REFERENCE_DEPTH_KM * (1.0 - math.exp(-REFERENCE_DEPTH_KM / depth))


def gravity_m_s2(mass_kg, radius_km):
    """The surface gravity of a body of this mass and radius."""
    return constants.G * mass_kg / (radius_km * 1000.0) ** 2
