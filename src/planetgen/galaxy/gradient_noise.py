# planetgen/galaxy/gradient_noise.py

"""
Gradient noise (Perlin) for the galaxy's gas and nebula density
===============================================================

3D gradient noise in about [-1, 1], from a permutation table shuffled by its seed, so the same seed always makes the same noise on every
worker. `nebula_density` (MAP.365) pushes the gas about and adds filaments with it. It was part of `nebula_shape` (GEN.75, the metaball shape of a
nebula), which GEN.1267 removed when a nebula's shape became the density contour the map draws.
"""

import numpy as np

from planetgen.util import draw

_GRADIENTS = np.array([
    [1, 1, 0], [-1, 1, 0], [1, -1, 0], [-1, -1, 0],
    [1, 0, 1], [-1, 0, 1], [1, 0, -1], [-1, 0, -1],
    [0, 1, 1], [0, -1, 1], [0, 1, -1], [0, -1, -1],
], dtype=float)
"""ndarray: Classic Perlin gradient directions."""


class GradientNoise:
    """
    3D gradient noise (Perlin) in about [-1, 1], from a permutation table
    shuffled by `seed`, so the same seed always makes the same noise.
    """

    def __init__(self, seed):
        table = list(range(256))
        draw.Stream(seed).shuffle(table)
        self._perm = np.array(table + table, dtype=np.int64)

    def __call__(self, points):
        """Noise at `points`, an array of shape (..., 3); returns shape (...)."""
        p = np.asarray(points, dtype=float)
        floor = np.floor(p)
        frac = p - floor
        cell = floor.astype(np.int64) & 255
        fade = frac * frac * frac * (frac * (frac * 6 - 15) + 10)
        perm = self._perm

        def corner(dx, dy, dz):
            h = perm[perm[perm[cell[..., 0] + dx] + cell[..., 1] + dy] + cell[..., 2] + dz]
            grad = _GRADIENTS[h % 12]
            offset = frac - np.array([dx, dy, dz], dtype=float)
            return np.sum(grad * offset, axis=-1)

        x, y, z = fade[..., 0], fade[..., 1], fade[..., 2]

        def lerp(a, b, t):
            return a + t * (b - a)

        x00 = lerp(corner(0, 0, 0), corner(1, 0, 0), x)
        x10 = lerp(corner(0, 1, 0), corner(1, 1, 0), x)
        x01 = lerp(corner(0, 0, 1), corner(1, 0, 1), x)
        x11 = lerp(corner(0, 1, 1), corner(1, 1, 1), x)
        return lerp(lerp(x00, x10, y), lerp(x01, x11, y), z)
