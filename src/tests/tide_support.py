# tests/tide_support.py

"""
A rocky body a test can hand to `interior_heat.tidal_response`: the stored mass and radius, the content of the bulk silicate Earth, a star of an age
and the two draws of its mantle's rheology (GEN.955).
"""

import math
import types

from planetgen.physics import body_spin, constants
from planetgen.physics.mathcheck import _bulk_silicate_earth


def rocky(mass_kg, radius_km, age_gy=4.5, eta0_pa_s=1.0e12, activation_kj_mol=250.0, **fields):
    """A rocky body of the given mass and radius with the bulk silicate Earth's content and a mantle of the given viscosity law."""
    content = _bulk_silicate_earth()
    stored = dict(
        mass=mass_kg, radius=radius_km, body_type="t", water_mass_fraction=0.0, star=types.SimpleNamespace(age=age_gy),
        c_k0_ppm=content.k0_ppm, c_th0_ppm=content.th0_ppm, c_u0_ppm=content.u0_ppm, al26_ratio=content.al26_ratio,
        mantle_eta0_pa_s=eta0_pa_s, mantle_activation_kj_mol=activation_kj_mol)
    stored.update(fields)
    return types.SimpleNamespace(**stored)


def locking_time_s(body, host_mass_kg, semi_major_axis_m, spin_rate, eccentricity=0.0, limit_years=1.0e12):
    """The age at which the torques on `body` have caught a spin born at `spin_rate` (`body_spin.evolve` with a synchronous state, a resonance or a balance), s,
    found by bisection in the log of the age; `math.inf` if it is not caught within `limit_years`."""
    year = constants.SECONDS_PER_YEAR

    def caught(age_s):
        found = body_spin.evolve(body, host_mass_kg, semi_major_axis_m, spin_rate, age_s, eccentricity)
        return found.kind != "free"

    if not caught(limit_years * year):
        return math.inf
    low, high = math.log(1.0e3 * year), math.log(limit_years * year)
    if caught(math.exp(low)):
        return math.exp(low)
    for _ in range(40):
        mid = 0.5 * (low + high)
        low, high = (low, mid) if caught(math.exp(mid)) else (mid, high)
    return math.exp(high)
