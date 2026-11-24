# planetgen/physics/fryer_remnants.py

"""
The metal dependence of a stellar black hole's mass (GEN.554; docs/design/systems/black-hole-mass-by-metals.md)
==============================================================================================================

[R1466] fit the baryonic mass a massive single star leaves, from the star's initial mass and its metallicity relative to the
sun `Z`, for the delayed explosion engine (the fits are the paper's Eqs. 5, 7 and 9 for the stars the failure windows hold,
23 solar masses and up; `Z` is capped at 1 where the fits are used, and at the sun the wind-stripped stars of the top of the
mass range fall to the Eq. 9 remnants of 2 to 3 solar masses). The gravitational mass of a black hole is 0.9 of the baryonic mass
(their Eq. 14).

What the generator takes from them is only the trend: how much heavier the holes get as the metals fall, as the
IMF-weighted mean over the stars that fail (`remnant_progenitors.kind_windows`) at the place's failed share
(`black_hole_formation.failed_share`). `black_hole_stretch` is that trend as the factor by which the span of the registered
black hole law grows (1 at the sun and above).
"""

import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.physics import imf, metallicity
from planetgen.physics.units import pc_to_ly

_GRID_STEP_DEX = 0.1
_CELLS_PER_WINDOW = 160
_HERMITE_POINTS = 15


def delayed_baryonic_mass_sol(initial_mass_sol, z_relative):
    """The baryonic remnant mass, solar masses, of a star of `initial_mass_sol` (23 to 150) at metallicity `z_relative`
    (0 to 1, the sun's 1) in the delayed engine: [R1466: Eq. 5] below 30 solar masses, [R1466: Eq. 7] from 30, and from 50 the greater
    of that and [R1466: Eq. 9], the fit for the stars that keep little of their mass."""
    mass = float(initial_mass_sol)
    z = min(max(float(z_relative), 0.0), 1.0)
    if mass < 30.0:
        return 1.1 + 0.2 * math.exp((mass - 11.0) / 4.0) - (2.0 + z) * math.exp(0.4 * (mass - 26.0))
    found = min(33.35 + (4.75 + 1.25 * z) * (mass - 34.0), mass - math.sqrt(z) * (1.3 * mass - 18.35))
    if mass > 50.0:
        found = max(found, 1.8 + 0.04 * (90.0 - mass) if mass < 90.0 else 1.8 + math.log10(mass - 89.0))
    return found


def delayed_black_hole_mass_sol(initial_mass_sol, z_relative):
    """The gravitational mass of the black hole such a star leaves: 0.9 of the baryonic mass ([R1466: Eq. 14])."""
    return 0.9 * delayed_baryonic_mass_sol(initial_mass_sol, z_relative)


def _mean_hole_mass_sol(failed_share, z_relative):
    """The mean mass of the holes the stars of the failure windows leave at `failed_share`, by IMF weight."""
    from planetgen.physics import remnant_progenitors   # a late import: it reads the structure, which imports the physics
    total = weighted = 0.0
    for low, high in remnant_progenitors.kind_windows(remnant_progenitors.BLACK_HOLE, failed_share):
        edges = np.geomspace(low, high, _CELLS_PER_WINDOW + 1)
        for a, b in zip(edges[:-1], edges[1:]):
            weight = imf.number(float(a), float(b))
            total += weight
            weighted += weight * delayed_black_hole_mass_sol(math.sqrt(float(a) * float(b)), z_relative)
    return weighted / total if total > 0.0 else 0.0


@functools.lru_cache(maxsize=1)
def _trend_table():
    """`(fe_h grid, mass ratio to the sun's)`: the mean Fryer hole mass at each [Fe/H] of the grid (the failed share of the
    reference structure there) over its mean at the sun's."""
    from planetgen.galaxy import structure
    from planetgen.physics import black_hole_formation
    reference = structure.REFERENCE_STRUCTURE
    low = tuning.METALLICITY_RANGE_FE_H[0]
    steps = int(math.ceil(-low / _GRID_STEP_DEX))
    grid = np.array([-_GRID_STEP_DEX * (steps - i) for i in range(steps + 1)])
    solar = _mean_hole_mass_sol(black_hole_formation.failed_share(0.0, reference.black_hole_solar_share,
                                                                  reference.black_hole_metal_free_share), 1.0)
    ratios = np.array([_mean_hole_mass_sol(float(black_hole_formation.failed_share(
        fe_h, reference.black_hole_solar_share, reference.black_hole_metal_free_share)), 10.0 ** float(fe_h)) / solar
                       for fe_h in grid])
    # The anchor law stands for the holes of the disk around the sun (the observed ones are there), a spread of [Fe/H] about
    # the disk's mean at the reference radius: the trend is 1 for that mixture, so a place as metal-rich as the sun's is
    # as the anchor is and a richer one lighter.
    from planetgen.galaxy import scales
    mean, sigma = metallicity.mean_fe_h("intermediate", pc_to_ly(scales.reference_radius_pc()))
    nodes, weights = np.polynomial.hermite_e.hermegauss(_HERMITE_POINTS)
    reference = float(np.sum(weights * np.interp(np.clip(mean + sigma * nodes, low, tuning.METALLICITY_RANGE_FE_H[1]),
                                                 grid, ratios)) / np.sum(weights))
    return grid, ratios / reference


def mass_trend(fe_h):
    """The mean mass of the holes a place of [Fe/H] `fe_h` makes over that of the disk around the sun (a number, or an array of
    them): Fryer's delayed engine over the windows that fail, constant from the sun up (its fit stops at the solar metals) and
    7 times higher at a tenth of them than at the sun."""
    grid, ratios = _trend_table()
    found = np.interp(np.minimum(np.asarray(fe_h, dtype=float), 0.0), grid, ratios)
    return float(found) if np.ndim(fe_h) == 0 else found
