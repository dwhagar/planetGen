# planetgen/galaxy/nebula_shape.py

"""
A nebula's shape (GEN.75)
=========================

Boss (2026-10-03): a nebula on the map should be "a kind of bulbous region
with an irregular shape", not a sphere or an ellipse. The shape is made the
way he suggested:

1. **Metaballs.** Four to eight centres inside an anisotropic ellipsoid
   (stretched along three axes, then turned), each a polynomial blob
   (`(1 - d^2/r^2)^2` inside its radius `r`, nothing outside) whose
   contributions add up to a scalar field `F`.
2. **Domain warp.** Before `F` is read, the point is pushed about by three
   channels of gradient noise, so the blobs' round outline becomes lumpy.
3. **An isovalue.** A point is inside the nebula where `F(warp(point)) >=
   iso`. Marching cubes at that level makes a triangle mesh, at a coarse
   grid for the Galaxy Map and a finer one for a close view.

Everything here is in units of the nebula's own radius: the shape is
scaled so the farthest point of its surface is at distance 1.0 from the
nebula's centre, so `radius_ly` is a bounding sphere and a cheap first
test. Multiply by the radius (and add the centre) to place it.

A shape is plain data (`NebulaShape.to_dict`, JSON-able), drawn from a
seeded `draw.Stream` (`draw_shape`), so the same seed always gives the
same shape on every worker. The noise is built from a permutation table
shuffled by the shape's own seed (no platform random state).

Needs numpy and scikit-image, both already dependencies.
"""

import hashlib
import math

import numpy as np

SHAPE_VERSION = 1
"""int: Written into `to_dict`; a change to how a field is computed bumps it."""

MESH_LODS = {"low": 10, "full": 28}
"""dict: Marching-cubes grid cells per axis for each mesh level of detail."""

_GRID_HALF_WIDTH = 1.8
"""float: The sampling grid reaches this far (natural units) either way
from the centre -- far enough for the largest drawn shape and its warp."""

from planetgen.util import draw

_FIT_GRID = 24
"""int: Grid cells per axis for finding the surface's farthest point."""

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


def _rotation(angles):
    """The 3x3 rotation matrix of three turns (radians) about x, y, then z."""
    ax, ay, az = angles
    cx, sx, cy, sy, cz, sz = math.cos(ax), math.sin(ax), math.cos(ay), math.sin(ay), math.cos(az), math.sin(az)
    rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return rz @ ry @ rx


class NebulaShape:
    """
    One nebula's shape. Build one with `draw_shape` or `from_dict`.

    Attributes:
        centres (list): Metaball centres, `[x, y, z]` in natural units.
        radii (list): Each ball's radius, natural units.
        axes (list): The ellipsoid's stretch along its three axes (the
            field is read at `point / axes` before the balls).
        angles (list): The ellipsoid's turn, radians about x, y, z.
        warp_amplitude (float): How far the noise pushes a point, natural units.
        warp_frequency (float): Noise cells per natural unit.
        warp_octaves (int): Noise layers (each twice the frequency, half the push).
        noise_seed (int): Seeds the noise's permutation table.
        iso (float): The field value at the surface.
        scale (float): Natural units per nebula radius: the surface's
            farthest point sits at `scale` natural units from the centre.
    """

    def __init__(self, centres, radii, axes, angles, warp_amplitude, warp_frequency, warp_octaves,
                 noise_seed, iso, scale=None):
        self.centres = [[float(v) for v in c] for c in centres]
        self.radii = [float(r) for r in radii]
        self.axes = [float(a) for a in axes]
        self.angles = [float(a) for a in angles]
        self.warp_amplitude = float(warp_amplitude)
        self.warp_frequency = float(warp_frequency)
        self.warp_octaves = int(warp_octaves)
        self.noise_seed = int(noise_seed)
        self.iso = float(iso)
        self._noise = GradientNoise(self.noise_seed)
        self._rot = _rotation(self.angles)
        self._centres = np.array(self.centres)
        self._radii = np.array(self.radii)
        self.scale = float(scale) if scale is not None else self._fit()

    # --- the field ---------------------------------------------------------

    def _warped(self, natural):
        """`natural` points pushed about by the warp noise."""
        out = natural.copy()
        amp, freq = self.warp_amplitude, self.warp_frequency
        for octave in range(self.warp_octaves):
            f = freq * (2 ** octave)
            a = amp * (0.5 ** octave)
            for channel in range(3):
                # A different region of the noise for each axis.
                out[..., channel] += a * self._noise(natural * f + 37.1 * (channel + 1) + 5.3 * octave)
        return out

    def field(self, natural):
        """The summed metaball field at `natural` points, shape (..., 3) -> (...)."""
        q = self._warped(np.asarray(natural, dtype=float))
        q = (q @ self._rot) / np.array(self.axes)
        total = np.zeros(q.shape[:-1])
        for centre, radius in zip(self._centres, self._radii):
            d2 = np.sum((q - centre) ** 2, axis=-1) / (radius * radius)
            falloff = np.clip(1.0 - d2, 0.0, None)
            total += falloff * falloff
        return total

    def _fit(self):
        """The distance from the centre of the surface's farthest point, natural units."""
        axis = np.linspace(-_GRID_HALF_WIDTH, _GRID_HALF_WIDTH, _FIT_GRID)
        grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
        inside = self.field(grid) >= self.iso
        if not inside.any():
            return 1.0
        reach = np.linalg.norm(grid[inside], axis=-1).max()
        # One grid step of slack: the grid only samples the surface.
        return float(reach + 2 * _GRID_HALF_WIDTH / (_FIT_GRID - 1))

    # --- in nebula-radius units ---------------------------------------------

    def contains(self, point):
        """Whether `point` (nebula-radius units from its centre) is inside the nebula."""
        p = np.asarray(point, dtype=float)
        if float(np.linalg.norm(p)) > 1.0:
            return False
        return bool(self.field(p * self.scale) >= self.iso)

    def contains_many(self, points):
        """`contains` for an array of points, shape (n, 3) -> (n,) bools."""
        p = np.asarray(points, dtype=float)
        near = np.linalg.norm(p, axis=-1) <= 1.0
        result = np.zeros(p.shape[:-1], dtype=bool)
        if near.any():
            result[near] = self.field(p[near] * self.scale) >= self.iso
        return result

    def mesh(self, lod="full"):
        """
        The surface as triangles by marching cubes, in nebula-radius units.

        Args:
            lod (str): `"low"` (a coarse grid, for the Galaxy Map) or `"full"`.

        Returns:
            tuple: `(vertices, faces)`: float array (n, 3) and int array (m, 3).
        """
        from skimage.measure import marching_cubes

        cells = MESH_LODS[lod]
        axis = np.linspace(-_GRID_HALF_WIDTH, _GRID_HALF_WIDTH, cells + 1)
        grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1)
        values = self.field(grid)
        if values.max() < self.iso or values.min() >= self.iso:
            return np.zeros((0, 3)), np.zeros((0, 3), dtype=int)
        step = axis[1] - axis[0]
        verts, faces, _normals, _values = marching_cubes(values, level=self.iso, spacing=(step, step, step))
        verts = (verts - _GRID_HALF_WIDTH) / self.scale
        return verts, faces

    # --- storage -------------------------------------------------------------

    def to_dict(self):
        """A JSON-able dict that `from_dict` turns back into the same shape."""
        return {
            "version": SHAPE_VERSION,
            "centres": self.centres, "radii": self.radii, "axes": self.axes, "angles": self.angles,
            "warp_amplitude": self.warp_amplitude, "warp_frequency": self.warp_frequency,
            "warp_octaves": self.warp_octaves, "noise_seed": self.noise_seed,
            "iso": self.iso, "scale": self.scale,
        }

    @classmethod
    def from_dict(cls, data):
        if data.get("version") != SHAPE_VERSION:
            raise ValueError(f"unknown nebula shape version {data.get('version')!r}")
        return cls(data["centres"], data["radii"], data["axes"], data["angles"], data["warp_amplitude"],
                   data["warp_frequency"], data["warp_octaves"], data["noise_seed"], data["iso"],
                   scale=data["scale"])


def draw_shape(rng):
    """
    Draws a shape from `rng` (a `draw.Stream`): 4-8 metaball centres in
    an anisotropic ellipsoid, a warp of the Boss-suggested kind, and the
    isovalue, so equal seeds give equal shapes.
    """
    count = rng.randint(4, 8)
    axes = sorted((rng.uniform(0.7, 1.0), rng.uniform(0.85, 1.3), rng.uniform(1.0, 1.7)))
    rng.shuffle(axes)
    centres = []
    for _ in range(count):
        # Uniform in a ball of radius 0.55 (the ellipsoid stretches it).
        while True:
            c = [rng.uniform(-1, 1) for _ in range(3)]
            if sum(v * v for v in c) <= 1:
                break
        centres.append([v * 0.55 for v in c])
    radii = [rng.uniform(0.35, 0.65) for _ in range(count)]
    angles = [rng.uniform(0, math.tau) for _ in range(3)]
    return NebulaShape(
        centres=centres, radii=radii, axes=axes, angles=angles,
        warp_amplitude=rng.uniform(0.12, 0.3), warp_frequency=rng.uniform(1.0, 2.2),
        warp_octaves=rng.randint(2, 3), noise_seed=rng.getrandbits(32),
        iso=rng.uniform(0.25, 0.45),
    )


def shape_for_nebula(nebula_class, radius_ly, density_cm3, temperature_k, extinction_av, dominant_species):
    """
    The shape a nebula with these properties gets (GEN.75): drawn from a
    seed made of the properties themselves, never from the shared random
    state, so asking for it (at generation, or later for a row saved
    before shapes were stored) consumes no random numbers and always
    answers the same.
    """
    text = "|".join((str(nebula_class), repr(float(radius_ly)), repr(float(density_cm3)),
                     repr(float(temperature_k)), repr(float(extinction_av)), str(dominant_species)))
    seed = int.from_bytes(hashlib.sha256(text.encode("utf-8")).digest()[:8], "big")
    return draw_shape(draw.Stream(seed))


SCALAR_COLUMNS = (
    "shape_axis_x", "shape_axis_y", "shape_axis_z", "shape_angle_x", "shape_angle_y", "shape_angle_z",
    "shape_warp_amplitude", "shape_warp_frequency", "shape_warp_octaves", "shape_noise_seed",
    "shape_iso", "shape_scale",
)
"""tuple: The `nebulae` columns holding a shape's single values; its metaballs
are the rows of `nebula_shape_balls` (`to_columns`, `from_columns`)."""


def to_columns(shape):
    """`shape` as `(scalars, balls)`: a dict keyed by `SCALAR_COLUMNS`, and
    `[(ball_index, x, y, z, radius), ...]` for `nebula_shape_balls`."""
    scalars = dict(zip(SCALAR_COLUMNS, (
        *shape.axes, *shape.angles, shape.warp_amplitude, shape.warp_frequency, shape.warp_octaves,
        shape.noise_seed, shape.iso, shape.scale)))
    balls = [(n, c[0], c[1], c[2], r) for n, (c, r) in enumerate(zip(shape.centres, shape.radii))]
    return scalars, balls


def from_columns(scalars, balls):
    """The `NebulaShape` that `to_columns` wrote: `scalars` is any mapping
    with the `SCALAR_COLUMNS` keys, `balls` the ball rows in index order."""
    ordered = sorted(balls, key=lambda b: b[0])
    return NebulaShape(
        centres=[[b[1], b[2], b[3]] for b in ordered], radii=[b[4] for b in ordered],
        axes=[scalars["shape_axis_x"], scalars["shape_axis_y"], scalars["shape_axis_z"]],
        angles=[scalars["shape_angle_x"], scalars["shape_angle_y"], scalars["shape_angle_z"]],
        warp_amplitude=scalars["shape_warp_amplitude"], warp_frequency=scalars["shape_warp_frequency"],
        warp_octaves=scalars["shape_warp_octaves"], noise_seed=scalars["shape_noise_seed"],
        iso=scalars["shape_iso"], scale=scalars["shape_scale"],
    )
