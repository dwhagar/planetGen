#!/usr/bin/env python3
# scripts/build_white_dwarf_table.py

"""
Builds `src/planetgen/physics/data/white_dwarf_mass_radius.json` (GEN.439): the zero-temperature white dwarf mass-radius
relation, integrated from the Chandrasekhar equation of state of a relativistic degenerate electron gas.

The gas is cold, fully ionised and has two nucleons per electron (mu_e = 2: a carbon-oxygen or oxygen-neon core). With x the
electron Fermi momentum over m_e c the pressure is P = A f(x) and the density rho = B x^3, with
A = pi m_e^4 c^5 / (3 h^3), B = 8 pi m_u mu_e (m_e c)^3 / (3 h^3) and
f(x) = x (2x^2 - 3) sqrt(1 + x^2) + 3 asinh(x). Hydrostatic equilibrium dP/dr = -G m rho / r^2 and dm/dr = 4 pi r^2 rho are
integrated outward from the centre (central x_c) to the radius where x falls to 1e-4. There is no Coulomb correction, no
temperature and no hydrogen or helium layer, so a real white dwarf is a few percent larger; Sirius B (1.018 solar masses,
0.00803 +- 0.00011 solar radii, [R2208]) comes out within 0.1 percent, and the honest test limit is 5 percent.

Run it from the repository root with scipy installed: `python scripts/build_white_dwarf_table.py`. The constants are
[R4361]; the solar mass and radius are those of `physics/constants.py`.
"""

import json
import math
import os

from scipy.integrate import solve_ivp

G = 6.6743e-11
C = 299792458.0
H = 6.62607015e-34
M_E = 9.1093837015e-31
M_U = 1.66053906660e-27
MU_E = 2.0
SOLAR_MASS_KG = 1.98841e30
SOLAR_RADIUS_M = 695.7e6
CENTRAL_X = (0.3, 12.0)
POINTS = 80
OUTER_X = 1e-4

A = math.pi * M_E ** 4 * C ** 5 / (3.0 * H ** 3)
B = 8.0 * math.pi * M_U * MU_E * (M_E * C) ** 3 / (3.0 * H ** 3)


def star(central_x):
    """`(mass_kg, radius_m)` of the star whose centre has the Fermi momentum `central_x` (in m_e c)."""
    def rhs(r, y):
        x, m = y
        dp_dx = 8.0 * A * x ** 4 / math.sqrt(1.0 + x * x)
        return [-G * m * B * x ** 3 / (r * r) / dp_dx, 4.0 * math.pi * r * r * B * x ** 3]

    def edge(_r, y):
        return y[0] - OUTER_X
    edge.terminal = True
    solved = solve_ivp(rhs, [1.0, 5e7], [central_x, 4.0 / 3.0 * math.pi * B * central_x ** 3], events=edge, rtol=1e-10,
                       atol=[1e-12, 1e10], method="LSODA")
    return solved.y[1][-1], solved.t[-1]


def build():
    rows = []
    for index in range(POINTS):
        central = CENTRAL_X[0] * (CENTRAL_X[1] / CENTRAL_X[0]) ** (index / (POINTS - 1))
        mass, radius = star(central)
        rows.append([round(central, 6), round(mass / SOLAR_MASS_KG, 6), round(radius / SOLAR_RADIUS_M, 8)])
    return rows


def main():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "planetgen", "physics", "data",
                        "white_dwarf_mass_radius.json")
    document = {
        "retrieved": "2026-10-17",
        "title": "The zero-temperature white dwarf mass-radius relation (GEN.439)",
        "sources": [
            "Chandrasekhar's equation of state of a cold relativistic degenerate electron gas with mu_e = 2, integrated by "
            "scripts/build_white_dwarf_table.py ([R4361] constants; scipy LSODA); Sirius B [R2208] is "
            "reproduced to 0.1 percent",
            "no Coulomb or temperature correction and no surface layer, so a real star is a few percent larger"],
        "columns": ["central_x", "mass_solar", "radius_solar"],
        "limiting_mass_solar": round(5.83 / MU_E ** 2, 4),
        "rows": rows_sorted(build())}
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=1)
        handle.write("\n")


def rows_sorted(rows):
    return sorted(rows, key=lambda row: row[1])


if __name__ == "__main__":
    main()
