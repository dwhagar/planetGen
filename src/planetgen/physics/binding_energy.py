# planetgen/physics/binding_energy.py

"""
Gravitational binding energy of a mass object (GEN.1158; /mnt/project-files/reports/round-shape-minimum-mass-2026-10-31.md)
=========================================================================================================================

The energy it takes to pull a body apart to infinity, J, always positive: the energy released when it assembles. It is `structure * G M^2 / R`, and the structure factor depends on how
the mass is spread, so each kind of mass object has its own formula. Boss's rule (8:44 PM and 8:45 PM EDT, 2026-10-30): every mass object that can have it stores it.

| Object | Formula | Structure factor | Source |
|---|---|---|---|
| Rocky or icy body, no core stored (asteroid, comet, small moon, dwarf planet) | uniform sphere | 3/5 | exact for a uniform sphere |
| Differentiated rocky body (planet, big moon) | two uniform shells, core and mantle | `interior_heat.binding_energy_j` | exact for shells |
| Gas giant | polytrope n = 1 | 3/4 | `giant_interior` uses n = 1 for gas giants |
| Brown dwarf, fully convective star (under `FULLY_CONVECTIVE_MASS_SOLAR`) | polytrope n = 3/2 | 6/7 | Lane-Emden |
| Star with a radiative core, giants | polytrope n = 3 (the Eddington standard model) | 3/2 | Lane-Emden |
| White dwarf | polytrope n = 3/2 (non-relativistic degenerate) | 6/7 | [R0730] |
| Neutron star | `BE / (M c^2) = 0.6 beta / (1 - 0.5 beta)`, `beta = G M / (R c^2)` | | [R2573] (ApJ 550, 426, [R2573: Eq. 33], READ) |
| Black hole | `INFINITE` | | Boss: stored as -1, meaning infinite |

For a polytrope of index `n` below 5 the energy is `3 / (5 - n) G M^2 / R` (Lane-Emden; exact for a self-gravitating polytropic sphere).

**-1 means infinite** (`INFINITE`). A real binding energy is always above zero, so a stored -1 can never be a real value. A black hole has no surface to pull a body away from and nothing can
leave it, so it stores -1 (`is_infinite` tells the sentinel apart). Anything that sums, compares or converts energies must skip it.

Limits. The white dwarf's index is 3/2 whatever its mass; near the maximum mass the true value is higher (the index goes to 3 as the electrons become relativistic, a factor 1.75 on the
structure factor, 6/7 against 3/2), so a heavy white dwarf is under by up to 43 percent. The stellar index is 3 for every star over a third of a solar mass, and 3/2 under it; a real solar model has an
effective index near 3.2. A neutron star's fit holds for stars without extreme softening to about 20 percent (Lattimer and Prakash's own words for their older fit; the new one is closer).
Not stored: a galaxy, cluster or nebula.

Use:

    binding_energy_j("rocky", mass_kg, radius_m)
    binding_energy_j("rocky", mass_kg, radius_m, core_mass_kg=m, core_radius_m=r)
    binding_energy_j("star", mass_kg, radius_m)
    binding_energy_j("black-hole", mass_kg)                                    # -1.0
"""

from planetgen.physics import constants, interior_heat
from planetgen.physics import unit_registry as units

_kg_to_msun = units.converter("kg", "Msun")
from planetgen.util.checks import finite_domain

INFINITE = -1.0
"""float: What a black hole stores for its binding energy: -1, meaning infinite (Boss, 8:45 PM EDT 2026-10-30). Never a real energy, which is always above zero."""

FULLY_CONVECTIVE_MASS_SOLAR = 0.35
"""float: A star lighter than this is convective throughout (index 3/2); a heavier one has a radiative core (index 3). [Stand-in BE1, tier T2, tweak L0.]"""

POLYTROPE_INDEX_GAS_GIANT = 1.0
"""float: The polytrope index of a gas giant (the one `giant_interior` and the central pressure formula use). [Stand-in BE2, tier T1, tweak L0.]"""

POLYTROPE_INDEX_CONVECTIVE = 1.5
"""float: The index of a fully convective star, a brown dwarf and a white dwarf (Lane-Emden n = 3/2). [Stand-in BE3, tier T1, tweak L0.]"""

POLYTROPE_INDEX_RADIATIVE = 3.0
"""float: The index of a star with a radiative core (the Eddington standard model). [Stand-in BE4, tier T2, tweak L1.]"""

UNIFORM_FACTOR = 0.6
"""float: The structure factor of a uniform sphere, 3/5."""

KINDS = ("rocky", "gas-giant", "brown-dwarf", "star", "white-dwarf", "neutron-star", "black-hole")
"""tuple: The kinds of mass object `binding_energy_j` knows."""

_MAX_BETA = 2.0


def is_infinite(value):
    """Whether a stored binding energy is the sentinel for infinite (a black hole)."""
    return value == INFINITE


@finite_domain()
def uniform_j(mass_kg, radius_m):
    """The binding energy of a uniform sphere, J: `3/5 G M^2 / R`."""
    if not (mass_kg > 0.0 and radius_m > 0.0):
        raise ValueError("the mass and the radius must be positive")
    return UNIFORM_FACTOR * constants.G * mass_kg ** 2 / radius_m


@finite_domain()
def polytrope_j(mass_kg, radius_m, index):
    """The binding energy of a polytrope of this index (below 5), J: `3 / (5 - n) G M^2 / R`."""
    if not (mass_kg > 0.0 and radius_m > 0.0):
        raise ValueError("the mass and the radius must be positive")
    if not 0.0 <= index < 5.0:
        raise ValueError("a polytrope is bound for an index from 0 up to 5")
    return 3.0 / (5.0 - index) * constants.G * mass_kg ** 2 / radius_m


@finite_domain()
def two_shell_j(mass_kg, radius_m, core_mass_kg, core_radius_m):
    """The binding energy of a core of `core_mass_kg` inside `core_radius_m` under a mantle to `radius_m`, two uniform shells, J (`interior_heat.binding_energy_j`)."""
    if not (mass_kg > 0.0 and radius_m > 0.0):
        raise ValueError("the mass and the radius must be positive")
    if not (0.0 <= core_mass_kg < mass_kg and 0.0 <= core_radius_m < radius_m):
        raise ValueError("the core lies inside the body")
    return interior_heat.binding_energy_j(mass_kg, radius_m, core_mass_kg, core_radius_m)


@finite_domain()
def neutron_star_j(mass_kg, radius_m):
    """The binding energy of a neutron star, J: [R2573: Eq. 33], `BE / (M c^2) = 0.6 beta / (1 - 0.5 beta)` with `beta = G M / (R c^2)`."""
    if not (mass_kg > 0.0 and radius_m > 0.0):
        raise ValueError("the mass and the radius must be positive")
    c2 = constants.SPEED_OF_LIGHT_M_S ** 2
    beta = constants.G * mass_kg / (radius_m * c2)
    if not 0.0 < beta < _MAX_BETA:
        raise ValueError("a neutron star's compactness G M / (R c^2) lies between 0 and 2")
    return mass_kg * c2 * 0.6 * beta / (1.0 - 0.5 * beta)


def star_index(mass_solar):
    """The polytrope index of a star of this mass: 3/2 under `FULLY_CONVECTIVE_MASS_SOLAR`, 3 over it."""
    return POLYTROPE_INDEX_CONVECTIVE if mass_solar < FULLY_CONVECTIVE_MASS_SOLAR else POLYTROPE_INDEX_RADIATIVE


def binding_energy_j(kind, mass_kg, radius_m=None, core_mass_kg=0.0, core_radius_m=0.0):
    """
    The binding energy of a mass object of one of `KINDS`, J (`INFINITE` for a black hole). `radius_m` is the mean radius (not needed for a black hole). A `rocky` body with a core
    (`core_mass_kg` over zero) is two shells, else a uniform sphere.
    """
    if kind not in KINDS:
        raise ValueError(f"unknown kind of mass object {kind!r}")
    if kind == "black-hole":
        return INFINITE
    if radius_m is None:
        raise ValueError(f"a {kind} needs its radius")
    if kind == "rocky":
        return two_shell_j(mass_kg, radius_m, core_mass_kg, core_radius_m) if core_mass_kg > 0.0 else uniform_j(mass_kg, radius_m)
    if kind == "gas-giant":
        return polytrope_j(mass_kg, radius_m, POLYTROPE_INDEX_GAS_GIANT)
    if kind == "brown-dwarf" or kind == "white-dwarf":
        return polytrope_j(mass_kg, radius_m, POLYTROPE_INDEX_CONVECTIVE)
    if kind == "star":
        return polytrope_j(mass_kg, radius_m, star_index(_kg_to_msun(mass_kg)))
    return neutron_star_j(mass_kg, radius_m)
