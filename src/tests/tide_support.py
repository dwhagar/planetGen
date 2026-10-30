# tests/tide_support.py

"""
A rocky body a test can hand to `interior_heat.tidal_response`: the stored mass and radius, the content of the bulk silicate Earth, a star of an age
and the two draws of its mantle's rheology (GEN.955).
"""

import types

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
