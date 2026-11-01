# planetgen/physics/aquifers.py

"""
Aquifers: liquid water held underground by a world whose surface is dry, frozen or too hot (GEN.654, Boss 2026-10-23)
===================================================================================================================

A world whose surface cannot keep liquid water can still hold it in the pores of its crust, where the heat flow of its
interior puts the temperature in the liquid range. The model is one column of crust, with no draw:

* **Geotherm.** T(z) = T_surface + q z / k, with q the stored internal heat flux (plus an icy moon's tide) and
  k = `tuning.AQUIFER_CONDUCTIVITY_W_M_K`.
* **Where it is liquid.** Below the 273.15 K level, z_cold = (273.15 - T_surface) k / q, and where the water's pressure,
  the surface pressure plus rho g z, exceeds its vapour pressure at T(z) (a hot surface boils the first metres away). Down
  to the hottest water life can use (`tuning.AQUIFER_LIFE_MAX_K`, 122 C) or to where the pores close, whichever is first.
* **Where there are pores.** The crust's porosity falls as phi0 exp(-z/K), K the Martian e-folding depth scaled by the
  inverse of the gravity (`tuning.AQUIFER_POROSITY_SCALE_KM_MARS`); the pores are gone under `tuning.AQUIFER_POROSITY_MIN`.
* **How much water.** The world's surface reservoir (the water it has that is not in the mantle), as a layer of water
  over the whole surface, fills the pores from the bottom up: a small inventory leaves a deep water table and a thin layer, a
  large one fills every pore of the liquid range and the water table is at its top.

The result is the depth to the top of the water-bearing layer and its thickness, in km; a body with none keeps `None` in both.
Every parameter is a DESIGN stand-in (`tuning`), the Martian ones from [R0824], read at second hand.
"""

import math

from planetgen import tuning
from planetgen.physics import constants, saturation

AQUIFER_FIELDS = ("aquifer_depth_km", "aquifer_thickness_km")
"""tuple: What `generate` sets on a planet or moon, as stored in `planets` and `moons`."""

WATER_MELTING_K = 273.15
LIQUID_WATER_DENSITY_KG_M3 = 1000.0
SURFACE_STATES = ("dry", "ice", "vapour")
"""tuple: The hydrosphere states that can have an aquifer: the ones with no liquid water at the surface."""


def _boils(temperature_k, pressure_pa):
    """Whether water at `pressure_pa` is a vapour at `temperature_k` (the vapour pressure is at or over the pressure)."""
    return saturation.pressure_pa("H2O", temperature_k) >= pressure_pa


def layer_km(surface_temperature_k, surface_pressure_pa, gravity_m_s2, flux_w_m2, water_depth_m):
    """
    `(depth_km, thickness_km)` of the water-bearing layer of a world with a surface at `surface_temperature_k` and
    `surface_pressure_pa`, `gravity_m_s2`, a heat flow of `flux_w_m2` and a surface reservoir of water that is
    `water_depth_m` deep over the whole surface; `(None, None)` for none.
    """
    if not (surface_temperature_k and gravity_m_s2 > 0 and flux_w_m2 > 0 and water_depth_m > 0):
        return None, None
    if surface_temperature_k >= tuning.AQUIFER_LIFE_MAX_K:
        return None, None
    gradient = flux_w_m2 / tuning.AQUIFER_CONDUCTIVITY_W_M_K * 1000.0  # K per km
    life_km = (tuning.AQUIFER_LIFE_MAX_K - surface_temperature_k) / gradient
    scale_km = tuning.AQUIFER_POROSITY_SCALE_KM_MARS * tuning.AQUIFER_MARS_GRAVITY_M_S2 / gravity_m_s2
    closure_km = scale_km * math.log(tuning.AQUIFER_POROSITY_SURFACE / tuning.AQUIFER_POROSITY_MIN)
    bottom_km = min(life_km, closure_km)
    top_km = max(0.0, (WATER_MELTING_K - surface_temperature_k) / gradient)

    def boils_at(depth_km):
        pressure = (surface_pressure_pa or 0.0) + LIQUID_WATER_DENSITY_KG_M3 * gravity_m_s2 * depth_km * 1000.0
        return _boils(surface_temperature_k + gradient * depth_km, pressure)

    if top_km == 0.0 and boils_at(0.0):
        if boils_at(bottom_km):
            return None, None
        low, high = 0.0, bottom_km
        for _ in range(60):
            middle = 0.5 * (low + high)
            low, high = (middle, high) if boils_at(middle) else (low, middle)
        top_km = high
    if top_km >= bottom_km:
        return None, None
    # The pores fill from the bottom: the water that fits above `bottom_km` below a depth z is phi0 K (e^-z/K - e^-bottom/K)
    pore_depth_km = tuning.AQUIFER_POROSITY_SURFACE * scale_km
    closed = math.exp(-bottom_km / scale_km)
    needed_km = water_depth_m / 1000.0 / pore_depth_km + closed
    table_km = -scale_km * math.log(needed_km) if needed_km < 1.0 else 0.0
    depth_km = max(top_km, table_km)
    thickness_km = bottom_km - depth_km
    if thickness_km < tuning.AQUIFER_MIN_THICKNESS_KM:
        return None, None
    return depth_km, thickness_km


def generate(body, surface_water_fraction, internal_flux_w_m2):
    """Sets `AQUIFER_FIELDS` on a planet or moon whose hydrosphere has no liquid at the surface (`SURFACE_STATES`) from the
    fraction of its mass that is its surface reservoir (`hydrosphere.generate_hydrosphere`), and clears them otherwise."""
    body.aquifer_depth_km = body.aquifer_thickness_km = None
    if body.hydrosphere not in SURFACE_STATES or not surface_water_fraction:
        return
    gravity = (body.gravity or 0.0) * constants.EARTH_GRAVITY
    area_m2 = 4.0 * math.pi * (body.radius * 1000.0) ** 2
    water_depth_m = surface_water_fraction * body.mass / (LIQUID_WATER_DENSITY_KG_M3 * area_m2)
    flux = internal_flux_w_m2 + (getattr(body, "tidal_flux_w_m2", None) or 0.0)
    body.aquifer_depth_km, body.aquifer_thickness_km = layer_km(
        body.surface_temperature, body.atmospheric_pressure, gravity, flux, water_depth_m)
