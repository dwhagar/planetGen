# planetgen/physics/hycean_candidate.py

"""
The Hycean box: which sub-Neptunes could carry a liquid ocean under a hydrogen envelope (GEN.963, registry 6.79)
================================================================================================================

Hycean worlds [R2804] are planets of 1 to 10 Earth masses and 1.4 to 2.6 Earth radii with an equilibrium temperature of 150 to 500 K and a density of
1.5 to 3 g/cc: a rocky core, a water layer that can be liquid at the surface, and a hydrogen envelope. The repository already has the *state*: `hydrosphere.py` calls a
liquid surface ocean under air that is more than half hydrogen (`HYCEAN_MIN_H2_SHARE`) `"hycean"`, and `biogenic_oxygen.py`, `habitability_world.py`, `magnetism.py` and
`bodylook.py` read it. What was missing is the box itself, as a flag a body can be asked about, so that the bodies the box selects can be counted and compared with the
bodies that reach the state. The envelope survival the item names is already `envelope_escape` (GEN.822) and the retention module (GEN.355).

The flag is derived from stored quantities, sits on the existing classes T (gas dwarf), I (ice giant) and V (high-gravity world) and moves no class and no stored value.

Use:

    hycean_candidate(planet_class, mass_earth, radius_earth, teq_k)       # True or False
    in_box(mass_earth, radius_earth, teq_k)                               # the four ranges, whatever the class
"""

import math

from planetgen.physics import constants, unit_registry as units

_mearth_to_kg = units.converter("Mearth", "kg")

HYCEAN_MASS_EARTH = (1.0, 10.0)
"""tuple: The mass of a Hycean candidate, Earth masses (the item, af[R2804]; row HC1)."""

HYCEAN_RADIUS_EARTH = (1.4, 2.6)
"""tuple: The radius, Earth radii (row HC1)."""

HYCEAN_TEMPERATURE_K = (150.0, 500.0)
"""tuple: The equilibrium temperature, K (row HC1)."""

HYCEAN_DENSITY_G_CM3 = (1.5, 3.0)
"""tuple: The bulk density, g/cc (row HC1)."""

HYCEAN_CLASS_LETTERS = ("T", "I", "V")
"""tuple: The classes the flag sits on: the gas dwarf, the ice giant and the high-gravity world (the item)."""


def density_g_cm3(mass_earth, radius_earth):
    """The bulk density, g/cc."""
    volume = 4.0 / 3.0 * math.pi * (radius_earth * constants.EARTH_RADIUS_KM * 1.0e5) ** 3
    return _mearth_to_kg(mass_earth) * 1000.0 / volume


def _inside(value, bounds):
    return bounds[0] <= value <= bounds[1]


def in_box(mass_earth, radius_earth, teq_k):
    """Whether mass, radius, equilibrium temperature and the density they give are all in the Hycean ranges."""
    return (_inside(mass_earth, HYCEAN_MASS_EARTH) and _inside(radius_earth, HYCEAN_RADIUS_EARTH) and _inside(teq_k, HYCEAN_TEMPERATURE_K)
            and _inside(density_g_cm3(mass_earth, radius_earth), HYCEAN_DENSITY_G_CM3))


def hycean_candidate(planet_class, mass_earth, radius_earth, teq_k):
    """Whether a body of class `planet_class` (a letter) is a Hycean candidate: one of `HYCEAN_CLASS_LETTERS` and in the box."""
    return planet_class in HYCEAN_CLASS_LETTERS and in_box(mass_earth, radius_earth, teq_k)
