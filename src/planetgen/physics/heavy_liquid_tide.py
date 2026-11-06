# planetgen/physics/heavy_liquid_tide.py

"""
The tide of a deep, global layer of liquid: the fluid Love numbers of a two-layer body (GEN.1052)
==================================================================================================

`ocean_tides` describes a thin sea, and refuses (`OceanTideError`) a liquid that is no thin layer. A global layer of a heavy liquid (a molten iron or mercury layer, a
deep water ocean) over a core is not thin, and the sea's own weight and self-attraction change the tide by an amount that depends on how dense the layer and the core
are. This is the incompressible, hydrostatic limit of that body, solved exactly: a core of density `rho_c` out to the radius fraction `f`, a liquid of density `rho_l` from
there to the surface, both inviscid, each interface an equipotential of the tidal potential plus the potential of the mass that its own displacement moves.

With the tidal potential `-(r / R)^2` (unit equilibrium height of `1 / g`) and `4 pi G` taken as 1, the two interface displacements `zeta_a` (at `f R`) and `zeta_R`
(at the surface) obey

    g_a zeta_a = f^2 + (1/5) (sigma_a f + sigma_R f^2),            g_R zeta_R = 1 + (1/5) (sigma_a f^4 + sigma_R)

with the sheet densities `sigma_a = (rho_c - rho_l) zeta_a` and `sigma_R = rho_l zeta_R`, the gravities `g_a = rho_c f / 3` and `g_R = rho_mean / 3`, and the mean density
`rho_mean = rho_c f^3 + rho_l (1 - f^3)`. Then `h2 = g_R zeta_R` and `k2` is the potential of the two sheets at the surface over the tidal potential.

Three limits are held by the tests. A core of the liquid's own density (a homogeneous fluid body) gives `k2 = 3/2` and `h2 = 5/2`. A core that does not move (`core_rigid`, a rock
mantle against the tide) gives `h2 = 1 / (1 - 3 rho_l / (5 rho_mean))` and `k2 = (3 rho_l / (5 rho_mean)) h2`, the loading factor of the thin sea at **any** depth, so the thin formula of
`ocean_tides` is exact for a rigid floor. A core that yields (a molten core under a deep layer) gives smaller numbers, between a homogeneous fluid's and the rigid floor's.

A layer must be lighter than the core (`rho_l < rho_c`) to stay on top; a liquid denser than the body's mean density cannot be a surface layer at all, because the mean density
is the core's and the layer's together (`layer_is_stable`), and it is a core.

What it does not do: compressibility (a deep layer's density rises with pressure), the rigidity of the core when it yields only partly (`tidal_rheology.response` gives that
for a rocky mantle), rotation, and the dissipation. The numbers are an equilibrium tide's.
"""

import math


class HeavyLiquidError(ValueError):
    """A layer that cannot stay on top, or a body with no layer."""


def mean_density_kg_m3(core_density_kg_m3, core_fraction, liquid_density_kg_m3):
    """The mean density of a core of `core_fraction` of the radius under a layer of liquid."""
    volume = core_fraction ** 3
    return core_density_kg_m3 * volume + liquid_density_kg_m3 * (1.0 - volume)


def layer_is_stable(core_density_kg_m3, liquid_density_kg_m3):
    """Whether a liquid layer stays on a core: the layer is lighter."""
    return liquid_density_kg_m3 < core_density_kg_m3


def two_layer_love(core_density_kg_m3, core_fraction, liquid_density_kg_m3, core_rigid=False):
    """
    `(k2, h2)` of an incompressible body of a core (`core_fraction` of the radius, density `core_density_kg_m3`) under a global liquid layer of `liquid_density_kg_m3`, with a core that
    moves with the tide (`core_rigid` false) or one that does not. `HeavyLiquidError` for a layer heavier than the core or a core fraction outside 0 to 1.
    """
    f = core_fraction
    if not 0.0 < f < 1.0:
        raise HeavyLiquidError(f"a core of {f:g} of the radius is no core under a layer")
    if liquid_density_kg_m3 <= 0.0 or core_density_kg_m3 <= 0.0:
        raise HeavyLiquidError("densities must be positive")
    if not layer_is_stable(core_density_kg_m3, liquid_density_kg_m3) and core_density_kg_m3 != liquid_density_kg_m3:
        raise HeavyLiquidError(f"a layer of {liquid_density_kg_m3:g} kg/m3 over a core of {core_density_kg_m3:g} kg/m3 would sink: it is the core")
    rho_c, rho_l = core_density_kg_m3, liquid_density_kg_m3
    mean = mean_density_kg_m3(rho_c, f, rho_l)
    g_r = mean / 3.0
    if core_rigid:
        # zeta_a = 0: g_R zeta_R = 1 + (1/5) rho_l zeta_R
        zeta_r = 1.0 / (g_r - rho_l / 5.0)
        return (rho_l / 5.0) * zeta_r, g_r * zeta_r
    g_a = rho_c * f / 3.0
    d = rho_c - rho_l
    # g_a za = f^2 + (1/5) (d za f + rho_l zr f^2)      ->  (g_a - d f / 5) za - (rho_l f^2 / 5) zr = f^2
    # g_R zr = 1 + (1/5) (d za f^4 + rho_l zr)          ->  -(d f^4 / 5) za + (g_R - rho_l / 5) zr = 1
    a11, a12, b1 = g_a - d * f / 5.0, -rho_l * f * f / 5.0, f * f
    a21, a22, b2 = -d * f ** 4 / 5.0, g_r - rho_l / 5.0, 1.0
    det = a11 * a22 - a12 * a21
    if det <= 0.0:
        raise HeavyLiquidError("the layer's own weight outweighs its gravity: no equilibrium tide")
    zeta_a = (b1 * a22 - a12 * b2) / det
    zeta_r = (a11 * b2 - a21 * b1) / det
    k2 = (d * zeta_a * f ** 4 + rho_l * zeta_r) / 5.0
    return k2, g_r * zeta_r


def fluid_love(core_density_kg_m3, core_fraction, liquid_density_kg_m3):
    """`(k2, h2)` of the body with a core that yields: `two_layer_love` with `core_rigid` false."""
    return two_layer_love(core_density_kg_m3, core_fraction, liquid_density_kg_m3, core_rigid=False)


def rigid_floor_love(mean_density_kg_m3, liquid_density_kg_m3):
    """`(k2, h2)` of a sea of `liquid_density_kg_m3` over a floor that does not move, on a body of the given mean density: `h2 = 1 / (1 - 3 rho_l / (5 rho_mean))` and `k2 = (3 rho_l / (5 rho_mean)) h2`. The thin-sea loading of
    `ocean_tides`, exact at any depth. `HeavyLiquidError` for a liquid denser than the body's mean density, which cannot be a surface layer."""
    if liquid_density_kg_m3 > mean_density_kg_m3:
        raise HeavyLiquidError(f"a liquid of {liquid_density_kg_m3:g} kg/m3 on a body of mean density {mean_density_kg_m3:g} kg/m3 is denser than the body: it cannot be a surface layer")
    ratio = 3.0 * liquid_density_kg_m3 / (5.0 * mean_density_kg_m3)
    h2 = 1.0 / (1.0 - ratio)
    return ratio * h2, h2
