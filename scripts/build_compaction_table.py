#!/usr/bin/env python3
# scripts/build_compaction_table.py

"""
Builds `src/planetgen/physics/data/compaction_table.json` (GEN.442, docs/design/small-bodies-composition.md section 4):
the radius and bulk density of an ice-and-rock body of a given mass and rock fraction, from the reference integration of
a uniform-composition sphere in hydrostatic balance.

Each layer holds ice grains (920 kg/m3 at zero pressure) and silicate grains (3500 kg/m3), both compressed by a third-order
Birch-Murnaghan equation of state (ice K0 9.2 GPa, K' 5.5; silicate K0 110 GPa, K' 4.0: the design note's stand-ins for the
solver's equation of state, not measurements), and each with the brittle porosity of [R0390] (Icarus 326,
Table 1 and Eqs. 1 and 2, P in MPa): ice phi = max(0.60 exp(-0.1 P), 0.20), silicate phi = min(0.6 P^-0.11, 0.6). The two
porosities are weighted by the grain volumes of the mixture. Between the masses `TAPER_START_KG` (1e19) and `TAPER_END_KG`
(1e21) the porosity is multiplied by a weight that falls from 1 to 0 in log mass (Boss's card of 8:13 PM EDT on 2026-10-21,
decision S9: the endpoints are flagged stand-ins until viscous compaction replaces them); above 1e21 kg the body has no
macroporosity and the structure solver of GEN.429 takes over.

The integration is run once per central pressure, outwards from the centre, and the (mass, radius) pairs are then put on a
common grid of masses, so one body costs one integration and the whole table takes about a minute.

Run it from the repository root with scipy installed: `python scripts/build_compaction_table.py`.
"""

import json
import math
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G = 6.6743e-11
ICE = (920.0, 9.2e9, 5.5)
SILICATE = (3500.0, 110.0e9, 4.0)
TAPER_START_KG = 1e19
TAPER_END_KG = 1e21
ROCK_FRACTIONS = [round(0.30 + 0.05 * i, 2) for i in range(14)]
LOG_MASS_RANGE = (12.0, 21.0)
MASS_STEP_DEX = 0.15
PRESSURE_GRID_PA = np.logspace(-3, 11, 800)
CENTRAL_PRESSURES_PA = np.logspace(-2.5, 9.2, 120)
OUTPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "planetgen", "physics", "data",
                      "compaction_table.json")


def grain_density(grain, pressure_pa):
    """Third-order Birch-Murnaghan density (kg/m3) of a grain `(rho0, K0, K')` at `pressure_pa`."""
    rho0, k0, kp = grain
    if pressure_pa <= 0.0:
        return rho0

    def f(x):
        return 1.5 * k0 * ((x ** (7 / 3) - x ** (5 / 3)) * (1 + 0.75 * (kp - 4) * (x ** (2 / 3) - 1))) - pressure_pa

    return rho0 * brentq(f, 1.0, 3.0)


def ice_porosity(pressure_pa):
    return max(0.60 * math.exp(-0.1 * pressure_pa / 1e6), 0.20)


def silicate_porosity(pressure_pa):
    return min(0.6 * (max(pressure_pa, 1.0) / 1e6) ** -0.11, 0.6)


def taper(mass_kg):
    span = math.log10(TAPER_END_KG) - math.log10(TAPER_START_KG)
    return float(min(max((math.log10(TAPER_END_KG) - math.log10(mass_kg)) / span, 0.0), 1.0))


def density_curve(rock_fraction, weight):
    """Bulk density (kg/m3) on `PRESSURE_GRID_PA` of a layer of this mixture, with the porosity scaled by `weight`."""
    ice_fraction = 1.0 - rock_fraction
    out = []
    for p in PRESSURE_GRID_PA:
        ice_density = grain_density(ICE, p)
        rock_density = grain_density(SILICATE, p)
        ice_pore = weight * ice_porosity(p)
        rock_pore = weight * silicate_porosity(p)
        out.append(1.0 / (ice_fraction / (ice_density * (1.0 - ice_pore)) + rock_fraction / (rock_density * (1.0 - rock_pore))))
    return np.array(out)


def sphere(central_pressure_pa, density_of_pressure):
    """`(mass_kg, radius_m)` of the sphere whose centre is at `central_pressure_pa`."""
    def rhs(r, y):
        p, m = y
        rho = density_of_pressure(p)
        return [-G * m * rho / (r * r), 4.0 * math.pi * r * r * rho]

    def surface(_r, y):
        return y[0] - 1e-6
    surface.terminal = True
    surface.direction = -1
    r0 = 0.01
    rho_c = density_of_pressure(central_pressure_pa)
    solved = solve_ivp(rhs, [r0, 5e8], [central_pressure_pa, 4.0 / 3.0 * math.pi * r0 ** 3 * rho_c], events=surface,
                       rtol=1e-7, atol=[1e-9, 1e-3], method="LSODA")
    return float(solved.y_events[0][0][1]), float(solved.t_events[0][0])


def build():
    log_masses = [LOG_MASS_RANGE[0] + MASS_STEP_DEX * i for i in range(int(round((LOG_MASS_RANGE[1] - LOG_MASS_RANGE[0]) / MASS_STEP_DEX)) + 1)]
    columns = []
    for rock in ROCK_FRACTIONS:
        # The taper depends on the mass of the body, which is what is being solved for: iterate the weight on the mass.
        radii, densities = [], []
        for log_mass in log_masses:
            mass = 10.0 ** log_mass
            weight = taper(mass)
            curve = density_curve(rock, weight)
            log_p = np.log10(PRESSURE_GRID_PA)

            def density_of_pressure(p, curve=curve, log_p=log_p):
                return float(np.interp(math.log10(max(p, 1e-3)), log_p, curve))

            def mass_error(log_pc):
                return math.log(sphere(10.0 ** log_pc, density_of_pressure)[0]) - math.log(mass)

            log_pc = brentq(mass_error, -2.5, 9.6, xtol=1e-5)
            got_mass, radius = sphere(10.0 ** log_pc, density_of_pressure)
            radii.append(radius / 1e3)
            densities.append(got_mass / (4.0 / 3.0 * math.pi * radius ** 3))
        columns.append({"rock_fraction": rock, "radius_km": [round(r, 6) for r in radii],
                        "density_kg_m3": [round(d, 3) for d in densities]})
        print(f"rock {rock}: R {radii[0]:.3f} to {radii[-1]:.1f} km, rho {densities[0]:.0f} to {densities[-1]:.0f}")
    return {
        "title": "Radius and bulk density of an ice-and-rock body by mass and rock fraction (GEN.442)",
        "method": "uniform sphere in hydrostatic balance, third-order Birch-Murnaghan grains, brittle porosity laws mixed by "
                  "volume, porosity tapered to zero between 1e19 and 1e21 kg",
        "sources": {
            "bierson-nimmo-2019": "[R0390] (accepted manuscript Table 1 and Eqs. 1 and 2): ice "
                                  "grains 920 and silicate grains 3500 kg/m3, ice phi = max(0.60 exp(-0.1 P), 0.20), silicate "
                                  "phi = min(0.6 P^-0.11, 0.6), P in MPa",
            "stand-in-eos": "the third-order Birch-Murnaghan constants (ice K0 9.2 GPa K' 5.5, silicate K0 110 GPa K' 4.0) are "
                            "the design note's stand-ins for the structure solver's equation of state (GEN.429)",
            "taper": "the 1e19 to 1e21 kg taper is Boss's decision S9 (8:13 PM EDT on 2026-10-21); its endpoints are two numbers "
                     "read off [R0666]'s statements and are a stand-in until viscous compaction replaces them",
        },
        "taper_start_kg": TAPER_START_KG,
        "taper_end_kg": TAPER_END_KG,
        "log10_mass_kg": [round(m, 6) for m in log_masses],
        "columns": columns,
    }


if __name__ == "__main__":
    table = build()
    with open(OUTPUT, "w", encoding="utf-8") as handle:
        json.dump(table, handle, indent=1)
        handle.write("\n")
    print("wrote", os.path.normpath(OUTPUT))
