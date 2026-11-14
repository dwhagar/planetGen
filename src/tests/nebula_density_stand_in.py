"""A stand-in for `planetgen.galaxy.nebula_density.density_at` for the tests of the nebula shapes (MAP.365): smooth lumps and
filaments as a pure function of the absolute position, finer octaves only when `detail` asks for them. The sampler's tests and the
Galaxy Map's browser fixture use it, so they hold the sampler and the drawing, not the model."""

import numpy as np

FORM_LEVEL = 0.35
MAX_DETAIL = 8


def density(points, field, detail, anchors=None):
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    shift = 0.0 if field == "gas" else 1.7
    lumps = np.sin(0.055 * x + 0.7 + shift) * np.cos(0.071 * y - 0.4) * np.sin(0.047 * z + 1.3) + 0.5 * np.sin(0.031 * (x + y) + shift)
    fine = 0.12 * np.sin(0.37 * x + 0.21 * y) * np.sin(0.29 * z) * min(detail, 3) / 3.0
    return np.clip(0.42 + 0.45 * lumps + fine, 0.0, 1.0).astype(np.float32)
