# planetgen/physics/water_layer.py

"""
The compressed water layer of a rocky or icy body (GEN.461; docs/design/composition/body-composition-parts/C-volatile-and-giants.md
sections 2.3 and 2.4)
===========================================================================================================================

Pure functions, no draws and no state. A body's water lies in a layer on top of its rock and iron. Its depth is not its mass
over 1000 kg/m3: at a water fraction of 0.1 of an Earth mass the layer is 746 km deep, not 1047, because the water at the
bottom is at 11 GPa and nearly twice as dense. The layer is integrated here from the surface down in the mass per area `s`
(kg/m2) above a level, the way the reference solver of `bulk_structure.layered_structure` does it for the whole body:

    dP/ds = g q,    dz/ds = q / rho(P),    g = G (M - f A s) / (R - z)^2,    q = (R / (R - z))^2

`P` is the pressure, `z` the depth, `rho(P)` the 300 K isotherm of water (`bulk_structure.water_density`: IAPWS-95 liquid to
1 GPa, then the ice VI, VII and X reference points of the materials table, GEN.427), `f` the share of the surface the water
covers and `A` the surface area, so the mass of the body below the level is the mass of the body less the water above it, and `q` is the spread of a unit of
surface area onto the smaller sphere at depth (the column mass per area grows as the sphere shrinks).
Four-step Runge-Kutta in `s`.

The liquid ends where the pressure reaches the **liquidus** of the temperature at the bottom (`liquidus_pressure_pa`): the
melting curve of ice III, V and VI recomputed from the [R2204] Gibbs-energy splines (design doc 2.3). Below it the
water is high-pressure ice, ice VI and then VII (`Column.hp_ice_km`), and an ocean with ice under it is sealed off from the
rock.

The density does not jump at the liquidus here: the isotherm is continuous through the freezing, and ice VI is 5 to 12 percent
denser than the liquid at the same pressure, so the high-pressure ice is thinner than this by that much at most. The 300 K
isotherm is also cooler than the bottom of a very deep ocean, which makes the layer slightly too deep by about 1 percent at 10
GPa (both stand-ins until a temperature profile is modelled).
"""

import math
from bisect import bisect_left
from dataclasses import dataclass

from planetgen.physics import bulk_structure, constants

SURFACE_PRESSURE_PA = bulk_structure.P_SURFACE_PA
"""float: The pressure at the top of the water when the caller gives none: the 1 bar level the radius is measured at."""

STEPS = 32
"""int: Runge-Kutta steps over a deep column (the 1 percent of the reference holds from 16)."""

SHALLOW_STEPS = 8
"""int: The steps over a column under `SHALLOW_PRESSURE_PA`, where the density changes by a few percent and most of the
galaxy's bodies are (an Earth's ocean is 40 MPa)."""
SHALLOW_PRESSURE_PA = 1.0e9

LIQUIDUS_K_MPA = (
    (251.2, 210.0), (256.2, 350.0), (271.7, 600.0), (272.1, 620.0), (286.8, 800.0), (300.4, 1000.0), (327.0, 1500.0),
    (347.0, 2000.0), (353.8, 2200.0),
)
"""tuple: Where liquid water stops being stable, `(K, MPa)`: ice III (251.2 K at 210 MPa), ice V (256.2 K at 350 MPa) and ice
VI (272.1 K at 620 MPa to 353.8 K at 2200 MPa) melting, COMPUTED from the [R2204] splines (design doc 2.3). Liquid
is stable below the line. The three published sets of triple points agree within 1 to 3 percent."""

LIQUIDUS_SLOPE_ABOVE_MPA_PER_K = (LIQUIDUS_K_MPA[-1][1] - LIQUIDUS_K_MPA[-2][1]) / (LIQUIDUS_K_MPA[-1][0] - LIQUIDUS_K_MPA[-2][0])
"""float: Beyond the ice VI-VII-liquid triple point (2.2 GPa, 355 K) the line continues with ice VI's last slope: a stand-in for
the ice VII melting curve, which the table does not hold."""


@dataclass(frozen=True)
class Column:
    """A water column from the surface down. Depths are kilometres of the actual (compressed) layer, over the covered area."""
    ice_shell_km: float
    """float: The ice lid at the top (0 with none)."""
    liquid_km: float
    """float: The liquid between the lid and the high-pressure ice (0 when none)."""
    hp_ice_km: float
    """float: The high-pressure ice (VI, VII) under the liquid, or under the lid when the liquid is gone (0 with none)."""
    base_pressure_pa: float
    """float: The pressure at the bottom of the water."""

    @property
    def depth_km(self):
        return self.ice_shell_km + self.liquid_km + self.hp_ice_km


def liquidus_pressure_pa(bottom_temperature_k):
    """The pressure at which water at `bottom_temperature_k` freezes into a high-pressure ice, Pa. Linear between the table's
    points, flat below the ice III triple point (251.2 K, the line through the metastable region is not modelled) and
    continued with ice VI's last slope above 353.8 K."""
    points = LIQUIDUS_K_MPA
    if bottom_temperature_k <= points[0][0]:
        return points[0][1] * 1e6
    if bottom_temperature_k >= points[-1][0]:
        return (points[-1][1] + (bottom_temperature_k - points[-1][0]) * LIQUIDUS_SLOPE_ABOVE_MPA_PER_K) * 1e6
    for (t0, p0), (t1, p1) in zip(points, points[1:]):
        if bottom_temperature_k <= t1:
            return (p0 + (bottom_temperature_k - t0) * (p1 - p0) / (t1 - t0)) * 1e6
    raise AssertionError("unreachable")


MIN_RADIUS_SHARE = 0.05
"""float: The column never reaches below this share of the radius (a body that is nearly all water at a radius that holds
less keeps the integration finite; the caller's radius and water fraction are then inconsistent)."""

ICE_LID_DENSITY_KG_M3 = 917.0
"""float: Ice Ih at the melting point ([R2204]: 916.72 kg/m3), the lid's density."""


def surface_gravity_m_s2(mass_kg, radius_km):
    """GM / R^2 of a body, m/s2."""
    return constants.G * mass_kg / (radius_km * 1000.0) ** 2


def _rates(pressure, depth, s, mass_kg, radius_m, covered, area_m2):
    radius_here = max(radius_m - depth, MIN_RADIUS_SHARE * radius_m)
    gravity = constants.G * max(mass_kg - covered * area_m2 * s, 0.0) / radius_here ** 2
    spread = (radius_m / radius_here) ** 2
    return gravity * spread, spread / bulk_structure.water_density(pressure)


def integrate(mass_kg, radius_km, s_from, s_to, *, covered=1.0, pressure_pa=SURFACE_PRESSURE_PA, depth_m=0.0, steps=STEPS):
    """The column from the mass per area `s_from` to `s_to` (kg/m2 of water above the level) under `pressure_pa` at `depth_m`:
    the nodes `[(s, pressure, depth_m, gravity), ...]`, `steps + 1` of them. `mass_kg` and `radius_km` are the body's (its radius
    at the top of the water); the water above a level is spread over `covered` of the surface."""
    radius_m = radius_km * 1000.0
    area_m2 = 4.0 * math.pi * radius_m ** 2
    h = (s_to - s_from) / steps
    s, pressure, depth = s_from, pressure_pa, depth_m
    nodes = []
    for _ in range(steps):
        k1 = _rates(pressure, depth, s, mass_kg, radius_m, covered, area_m2)
        nodes.append((s, pressure, depth, k1[0] / (radius_m / (radius_m - depth)) ** 2))
        k2 = _rates(pressure + 0.5 * h * k1[0], depth + 0.5 * h * k1[1], s + 0.5 * h, mass_kg, radius_m, covered, area_m2)
        k3 = _rates(pressure + 0.5 * h * k2[0], depth + 0.5 * h * k2[1], s + 0.5 * h, mass_kg, radius_m, covered, area_m2)
        k4 = _rates(pressure + h * k3[0], depth + h * k3[1], s + h, mass_kg, radius_m, covered, area_m2)
        pressure += h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0
        depth += h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0
        s += h
    last = _rates(pressure, depth, s, mass_kg, radius_m, covered, area_m2)[0]
    nodes.append((s, pressure, depth, last / (radius_m / (radius_m - depth)) ** 2))
    return nodes


def _depth_slope(pressure, gravity):
    """dz/dP = 1 / (rho g)."""
    return 1.0 / (bulk_structure.water_density(pressure) * gravity)


def depth_at_pressure_m(nodes, pressure_pa):
    """The depth (m) in `nodes` where the pressure is `pressure_pa`: cubic Hermite in the pressure with the slope 1 / (rho g)
    of the depth at the two nodes around it. The pressure must lie between the first and the last node's."""
    pressures = [node[1] for node in nodes]
    i = min(max(bisect_left(pressures, pressure_pa) - 1, 0), len(nodes) - 2)
    (_, p0, z0, g0), (_, p1, z1, g1) = nodes[i], nodes[i + 1]
    h = p1 - p0
    t = (pressure_pa - p0) / h
    d0, d1 = _depth_slope(p0, g0), _depth_slope(p1, g1)
    return ((2 * t ** 3 - 3 * t ** 2 + 1) * z0 + (t ** 3 - 2 * t ** 2 + t) * h * d0
            + (-2 * t ** 3 + 3 * t ** 2) * z1 + (t ** 3 - t ** 2) * h * d1)


def column(mass_kg, radius_km, water_kg, *, covered=1.0, shell_km=0.0, bottom_temperature_k, pressure_pa=SURFACE_PRESSURE_PA,
           steps=None):
    """The water of a body of `mass_kg` and `radius_km`: `water_kg` of it over `covered` of the surface, under an ice lid
    `shell_km` thick (ice Ih), with `pressure_pa` of air on top and `bottom_temperature_k` where the liquid meets the ice below
    it. Returns the `Column`: lid, liquid, high-pressure ice and the base pressure. A lid heavier than the water is the whole
    column; a liquid column already past the liquidus at the lid's base is high-pressure ice all the way."""
    area_m2 = 4.0 * math.pi * (radius_km * 1000.0) ** 2
    total = water_kg / (covered * area_m2)
    gravity = surface_gravity_m_s2(mass_kg, radius_km)
    lid = ICE_LID_DENSITY_KG_M3 * shell_km * 1000.0
    if steps is None:
        steps = SHALLOW_STEPS if gravity * total < SHALLOW_PRESSURE_PA else STEPS
    if lid >= total:
        return Column(total / ICE_LID_DENSITY_KG_M3 / 1000.0, 0.0, 0.0, pressure_pa + gravity * total)
    start_pressure = pressure_pa + gravity * lid
    nodes = integrate(mass_kg, radius_km, lid, total, covered=covered, pressure_pa=start_pressure, depth_m=shell_km * 1000.0,
                      steps=steps)
    base_pressure, base_depth = nodes[-1][1], nodes[-1][2]
    liquidus = liquidus_pressure_pa(bottom_temperature_k)
    if base_pressure <= liquidus:
        return Column(shell_km, (base_depth - shell_km * 1000.0) / 1000.0, 0.0, base_pressure)
    if start_pressure >= liquidus:
        return Column(shell_km, 0.0, (base_depth - shell_km * 1000.0) / 1000.0, base_pressure)
    freezing_depth = depth_at_pressure_m(nodes, liquidus)
    return Column(shell_km, (freezing_depth - shell_km * 1000.0) / 1000.0, (base_depth - freezing_depth) / 1000.0, base_pressure)
