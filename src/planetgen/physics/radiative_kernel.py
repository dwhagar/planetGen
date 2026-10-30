# planetgen/physics/radiative_kernel.py

"""
The compiled two-stream sweep of the radiative-convective column (GEN.850)
=========================================================================

`radiative_convective.sweep` for a few profiles at once: the upward and downward infrared flux at every level, layer by layer, band by band and g point by
g point, with the source linear in optical depth between level temperatures. The sweep was most of the 1 to 12 seconds a column took in numpy, and
`scripts/fit_greenhouse_table.py` asks for tens of thousands of columns.

The arithmetic is the numpy function's, line for line (`tests/test_radiative_convective.py` holds the two to a few parts in 1e12).

Like `field_kernel`, the module is **kept out of the web app's import path** (numba's cache needs a writable code tree): only the column solver, which the
web app does not import either, calls it, and `test_the_web_app_never_imports_numba` holds that.
"""

import math

import numpy as np
from numba import njit

PATH_FLOOR = 1.0e-6
"""float: A path under this is optically thin: the slope of the source is `path / 2` (the series of `1 - (1 - exp(-path)) / path`)."""
PATH_TINY = 1.0e-300
"""float: Keeps the division finite where the thin branch is taken instead."""


@njit(cache=True)
def sweep(emission, ground, depth, weights, diffusivity):
    """
    The upward and downward flux, W/m2, at every level: `emission` is (profile, level, band) the blackbody power of each level's temperature, `ground`
    (profile, band) that of the surface, `depth` (layer, band, g) the layers' optical depths shared by the profiles, `weights` (g) the g weights
    (adding to 1). Returns `(flux_up, flux_down)`, each (profile, level); the level index runs from the top (0) to the surface (layers).
    """
    count = emission.shape[0]
    layers = depth.shape[0]
    bands = depth.shape[1]
    points = depth.shape[2]
    flux_up = np.zeros((count, layers + 1))
    flux_down = np.zeros((count, layers + 1))
    transmitted = np.empty((layers, bands, points))
    absorbed = np.empty((layers, bands, points))
    slope = np.empty((layers, bands, points))
    for layer in range(layers):
        for band in range(bands):
            for point in range(points):
                path = diffusivity * depth[layer, band, point]
                transmitted[layer, band, point] = math.exp(-path)
                lost = -math.expm1(-path)
                absorbed[layer, band, point] = lost
                if path > PATH_FLOOR:
                    slope[layer, band, point] = 1.0 - lost / max(path, PATH_TINY)
                else:
                    slope[layer, band, point] = path / 2.0
    stream = np.empty((bands, points))
    for profile in range(count):
        stream[:, :] = 0.0
        for layer in range(layers):
            total = 0.0
            for band in range(bands):
                top = emission[profile, layer, band]
                bottom = emission[profile, layer + 1, band]
                for point in range(points):
                    value = (stream[band, point] * transmitted[layer, band, point] + top * absorbed[layer, band, point]
                             + (bottom - top) * slope[layer, band, point])
                    stream[band, point] = value
                    total += value * weights[point]
            flux_down[profile, layer + 1] = total
        total = 0.0
        for band in range(bands):
            for point in range(points):
                stream[band, point] = ground[profile, band]
                total += ground[profile, band] * weights[point]
        flux_up[profile, layers] = total
        for layer in range(layers - 1, -1, -1):
            total = 0.0
            for band in range(bands):
                top = emission[profile, layer, band]
                bottom = emission[profile, layer + 1, band]
                for point in range(points):
                    value = (stream[band, point] * transmitted[layer, band, point] + bottom * absorbed[layer, band, point]
                             + (top - bottom) * slope[layer, band, point])
                    stream[band, point] = value
                    total += value * weights[point]
            flux_up[profile, layer] = total
    return flux_up, flux_down
