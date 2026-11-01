# planetgen/physics/shape_temperature.py

"""
The temperature a body took its shape at (GEN.1201)
===================================================

The roundness engine (`roundness`) gives a body a strength that falls with the temperature the body had when its weight was shaping it: ice creeps at a tenth of its cold strength near
its melting point, and dry olivine creeps below its cold strength from about 1,100 K. Until now no body stored such a temperature, so every strength was the cold strength. This module
works it out from what a body stores, in one place and without drawing from the run's stream:

    T_shape = min(T_surface + dT_accretion + dT_heat, T_ceiling)

* **The surface** (`surface_k`): the body's temperature without any inner heat, the equilibrium temperature the star's light gives it at its distance. It is what a small cold body keeps.
* **Accretion** (`dT_accretion`): the gravitational energy of assembling the body, per kilogram (the binding energy the body stores over its mass), over the heat capacity of its mix. The
  whole of it is kept, as the thermal history of `interior_heat` starts a body at `Ts + E / C` (a magma ocean radiates the rest away: the ceiling). It is a kelvin for a planetesimal, about
  1,500 K for a Moon, and more than the ceiling for a planet.
* **Heat** (`dT_heat`): the heat the body's rock makes while it is still warm enough to flow. Heat leaves a body of radius `R` by conduction in the time `t_c = R^2 / (pi^2 kappa)`
  (the first mode of a conducting sphere), so what the short-lived 26Al and the long-lived isotopes release in `t_c` (`radiogenic.heat_terms`, the stored content of the body decayed by the
  integral `a (1 - exp(-lambda t_c)) / lambda`) stays, and the tides' mean power over the same time (`tidal_w`) with it. A body 10 km across keeps the heat of its first 80 kyr and no more, one 100
  km across that of its first 8 Myr (most of its 26Al). The heat conducting away at the same time is not subtracted: a steady-state balance, `h R^2 / (6 kappa)`, comes out at 1.6 times the
  heat of `t_c`, so the heat of the window is the smaller and stands.
* **The ceiling** (`T_ceiling`): the rock's liquidus at the pressure in the middle of the body (`interior_heat.liquidus_k`): above it the body is a magma ocean, whose temperature the
  radiating surface holds there (the same ceiling `interior_heat` puts on the mantle).

The strength engine caps what it reads: ice creeps no warmer than its melting point (273.15 K) and rock changes its strength only above about 1,000 K, so the number matters for an icy body
that stayed warm and for a rocky one that a short-lived isotope melted, and is the surface temperature for the rest.

What it does not have: the surface is the body's own temperature now, not the young star's at the time the body solidified (the star was brighter before it reached the main sequence); a
body that formed late has less 26Al (the content it stores already carries that); a moon's tides are the mean power of today. The heat of a body's whole history is `interior_heat`'s work,
and this is the one number for the shape. Stand-ins SHT1 and SHT2 of `gap-filling-plan.md` section 6.162. Nothing reads the module until GEN.1201's wiring step lands (it needs GEN.1159's stored
roundness columns).

Use:

    shape_temperature_k(mass_kg, radius_m, ice, metal, silicate, surface_k, binding_energy_j, content, aluminium, tidal_w)
"""

import math

from planetgen.physics import constants, interior_heat, radiogenic
from planetgen.util.checks import finite_domain

ICE_SPECIFIC_HEAT = 1700.0
"""float: J/(kg K), the heat capacity of the ice of a body (ice Ih is 1,500 at 160 K and 2,100 at 273 K). [Stand-in SHT1, tier T7, tweak L0.]"""

WINDOW_SHAPE_FACTOR = math.pi ** 2
"""float: The conduction time of a sphere is `R^2 / (WINDOW_SHAPE_FACTOR kappa)`: the first mode of a conducting sphere decays at `pi^2 kappa / R^2`. [Stand-in SHT2, tier T6, tweak L1.]"""

_GY_S = interior_heat.GY_S


class ShapeTemperatureError(ValueError):
    """A mass, radius, share or temperature argument is out of range."""


def conduction_time_s(radius_m, diffusivity=interior_heat.MANTLE_DIFFUSIVITY):
    """The time heat takes to leave a sphere of this radius, s: `R^2 / (pi^2 kappa)` (the rock's diffusivity of `interior_heat`)."""
    if not (radius_m > 0.0 and diffusivity > 0.0):
        raise ShapeTemperatureError("the radius and the diffusivity must be positive")
    return radius_m ** 2 / (WINDOW_SHAPE_FACTOR * diffusivity)


def heat_capacity_j_kg_k(ice, metal, silicate):
    """The mass-weighted heat capacity of a body's mix, J/(kg K): ice, silicate and metal shares of its mass (their sum is 1 or less; any rest counts as silicate)."""
    rest = max(1.0 - ice - metal - silicate, 0.0)
    return ice * ICE_SPECIFIC_HEAT + (silicate + rest) * interior_heat.MANTLE_SPECIFIC_HEAT + metal * interior_heat.CORE_SPECIFIC_HEAT


def radiogenic_energy_j_kg(content, window_s, silicate_share, aluminium=radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION):
    """The energy, J per kg of the whole body, that its silicate's radioactive decays release in the first `window_s` seconds: `silicate_share * sum a (1 - exp(-lambda t)) / lambda`."""
    if content is None or window_s <= 0.0 or silicate_share <= 0.0:
        return 0.0
    window_gy = window_s / _GY_S
    return silicate_share * _GY_S * math.fsum(a * (-math.expm1(-rate * window_gy)) / rate for a, rate in radiogenic.heat_terms(content, aluminium))


def ceiling_k(mass_kg, radius_m):
    """The hottest a solid body's shape can have been set at, K: the rock's liquidus at the pressure in the middle of the body (`rho g R / 2`, uniform density)."""
    density = mass_kg / (4.0 / 3.0 * math.pi * radius_m ** 3)
    gravity = constants.G * mass_kg / radius_m ** 2
    return interior_heat.liquidus_k(density * gravity * radius_m / 2.0e9)


@finite_domain()
def shape_temperature_k(mass_kg, radius_m, ice, metal, silicate, surface_k, binding_energy_j=None, content=None,
                        aluminium=radiogenic.BULK_SILICATE_EARTH_AL_MASS_FRACTION, tidal_w=0.0):
    """
    The temperature, K, a body of this mass and radius took its shape at.

    Args:
        mass_kg, radius_m (float): The body's mass and mean radius.
        ice, metal, silicate (float): The shares of its mass that are ice of any kind, metal and silicate (any rest counts as silicate).
        surface_k (float): The body's temperature with no inner heat: the star's light at its distance.
        binding_energy_j (float or None): The energy released by assembling it, J (`body_strength`'s stored value); None counts the accretion heat as 0.
        content (radiogenic.Content or None): The body's radioactive content at formation; None for a body that stores none (it makes no radiogenic heat).
        aluminium (float): Aluminium in the bulk silicate, mass fraction (26Al's parent).
        tidal_w (float): The mean power the tides dissipate in the body, W (0 for a body with none).
    """
    if not (mass_kg > 0.0 and radius_m > 0.0 and surface_k > 0.0):
        raise ShapeTemperatureError("the mass, the radius and the surface temperature must be positive")
    for name, share in (("ice", ice), ("metal", metal), ("silicate", silicate)):
        if not 0.0 <= share <= 1.0:
            raise ShapeTemperatureError(f"the {name} share must lie in 0 to 1, got {share!r}")
    if ice + metal + silicate > 1.0 + 1.0e-9:
        raise ShapeTemperatureError("the ice, metal and silicate shares add to more than 1")
    if tidal_w < 0.0:
        raise ShapeTemperatureError("the tidal power must not be negative")
    capacity = heat_capacity_j_kg_k(ice, metal, silicate)
    window = conduction_time_s(radius_m)
    accretion = (binding_energy_j or 0.0) / mass_kg
    heat = radiogenic_energy_j_kg(content, window, silicate, aluminium) + tidal_w * window / mass_kg
    return min(surface_k + (accretion + heat) / capacity, max(ceiling_k(mass_kg, radius_m), surface_k))
