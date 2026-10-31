# planetgen/physics/air_envelope.py

"""
The envelope share of a rocky body with a hydrogen-rich air (GEN.1226; Boss, 1:56 PM EDT 2026-10-31)
======================================================================================================

A rocky body drawn with no envelope stores `envelope_mass_fraction` 0 whatever air it holds, so the body type of `planet_type` (GEN.1219) cannot
see a hydrogen-air "gas dwarf" (class T). This module gives the share to store: the **true air share**, the mass of the air over the mass of the body,
for a body whose air is mostly hydrogen and helium; any other rocky body keeps 0.

The mass of the air is the column the surface pressure holds: `P * 4 pi R^2 / g` (the air's weight over the whole sphere), with g the surface gravity
`G M / R^2`. The bodies this applies to hold 1e-7 to 1e-3 of their mass in air (Venus 1e-4), so by `planet_type` they stay terrestrial; the
number is stored because it is true, not to move a cut.
"""

import math

from planetgen.physics import constants

HYDROGEN_RICH_SHARE = 0.3
"""float: STAND-IN AE1 (L0): an air with this share or more of its pressure in hydrogen and helium is hydrogen-rich (the hydrogen line of the best-fit class table)."""
LIGHT_GASES = ("H2", "He")
"""tuple: The formulas that count as the hydrogen-helium envelope gas."""


def light_share(partials_pa):
    """The share of the total pressure held by hydrogen and helium."""
    total = sum(value for value in partials_pa.values() if value > 0.0)
    if total <= 0.0:
        return 0.0
    return sum(partials_pa.get(formula, 0.0) for formula in LIGHT_GASES) / total


def air_mass_kg(pressure_pa, mass_kg, radius_km):
    """The mass of the air a surface pressure holds on a body, kg (`P 4 pi R^2 / g`)."""
    radius_m = radius_km * constants.KM_TO_M_FACTOR
    gravity = constants.G * mass_kg / radius_m ** 2
    return pressure_pa * 4.0 * math.pi * radius_m ** 2 / gravity


def envelope_share(partials_pa, pressure_pa, mass_kg, radius_km):
    """The envelope share to store for a rocky body: its air's mass over its own when the air is hydrogen-rich, else 0."""
    if not pressure_pa or pressure_pa <= 0.0 or light_share(partials_pa) < HYDROGEN_RICH_SHARE:
        return 0.0
    return min(1.0, air_mass_kg(pressure_pa, mass_kg, radius_km) / mass_kg)
