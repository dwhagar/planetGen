#!/usr/bin/env python3
"""
Builds the fast radius of a small differentiated body (GEN.429): a polynomial fitted to `bulk_structure.layered_structure`.

    ln R/R_earth = sum over i<=3, j<=2, k<=2 of a_ijk u^i c^j w^k    u = ln(M/M_earth), c = metal share of the dry part, w = water share

over 1.5e-4 to 0.3 Earth masses, c from 0 to 0.6 and w from 0 to 1 (mantle Fe/(Fe+Mg) 0.11, core density factor 0.9). Part B's fit
stops at 0.1 Earth masses and part C's at 0.3, so a Mercury, a Mars or a Ganymede has no closed form without it. Run it
to rebuild `SMALL_FIT` in `src/planetgen/physics/bulk_structure.py`; it prints the coefficients and the error on a fresh set.

    python3 scripts/fit_small_body_radius.py
"""
import itertools
import math
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from planetgen.physics import bulk_structure as bs  # noqa: E402

MASS_RANGE = (1.5e-4, 0.3)
CMF_RANGE = (0.0, 0.6)
WATER_RANGE = (0.0, 1.0)
POWERS = [(i, j, k) for i in range(4) for j in range(3) for k in range(3)]


def terms(mass_earth, cmf, water):
    u = math.log(mass_earth)
    return [u ** i * cmf ** j * water ** k for i, j, k in POWERS]


def solver_radius(mass_earth, cmf, water):
    dry = 1.0 - water
    return bs.layered_structure(mass_earth, cmf * dry, (1.0 - cmf) * dry, water, 0.11, 0.9, core_steps=40, mantle_steps=80,
                                water_steps=80)["radius"]


def sample(rng, n):
    rows = []
    for _ in range(n):
        mass = math.exp(rng.uniform(math.log(MASS_RANGE[0]), math.log(MASS_RANGE[1])))
        cmf = rng.uniform(*CMF_RANGE)
        water = rng.choice((0.0, 0.0, rng.uniform(0.0, 1.0)))
        rows.append((mass, cmf, water))
    return rows


def main():
    rng = random.Random(429)
    train = sample(rng, 600)
    x = np.array([terms(*row) for row in train])
    y = np.array([math.log(solver_radius(*row)) for row in train])
    coefficients, *_ = np.linalg.lstsq(x, y, rcond=None)
    test = sample(random.Random(430), 150)
    errors = [math.exp(sum(c * t for c, t in zip(coefficients, terms(*row)))) / solver_radius(*row) - 1.0 for row in test]
    print("rms %.3f%%  max %.3f%%" % (100 * math.sqrt(sum(e * e for e in errors) / len(errors)), 100 * max(abs(e) for e in errors)))
    print("SMALL_FIT = (")
    for (i, j, k), c in zip(POWERS, coefficients):
        print(f"    ({i}, {j}, {k}, {float(c)!r}),")
    print(")")


if __name__ == "__main__":
    main()
