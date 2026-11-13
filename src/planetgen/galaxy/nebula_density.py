# planetgen/galaxy/nebula_density.py

"""
The density a nebula is drawn from on the Galaxy Map (MAP.365)
==============================================================

Boss (12:37 PM EST on 2026-11-04): a nebula's shape comes from the real gas and star-forming density of the galaxy
model, never a random blob per object ID; wispy edges and a dense core that fades unevenly; coarse blobs far out and
more filaments as you zoom in. `density_at` is that density: a pure function of an absolute position, the galaxy seed
and the stored shape, in 0..1. Two callers who ask for the same point agree, so the nodes and the sectors of the map
meet without a seam, and the same galaxy draws the same shape on every worker. The object ID is never read.

A point is inside a nebula where the value is at least `FORM_LEVEL`; the map draws the contour and fades the colour
from the core out with the value itself.

How the value is made, for a field:

* a smooth model density `s` in 0..1.
  `"gas"` (dark and molecular clouds): the seeded cloud field (`nebula_field`, GEN.47) is the anchor, because the
  clouds are where the model placed gas: each cloud adds a polynomial kernel `(1 - d^2 / r^2)^2` out to its radius,
  strengthened where the model gas is dense (`nebula_field.gas_factor`: the arms, the thin layer, so an outline leans
  with them), and a superbubble (`superbubble_field`, GEN.254) clears the gas inside it. A thin background of the
  model gas alone (`GAS_BACKGROUND` of its saturating share) lets filaments reach a little past a cloud's radius.
  `"young"` (emission and reflection nebulae, which form where hot young stars are): the young stars' density
  (`density.population_densities_array`) against its value on the solar circle, saturating, anchored on the nebula's own
  sphere (GEN.1267, MAP.423): the model young density falls under the form level outside about 6 kpc, where a hot star's
  nebula would have no shape at all, so each nebula adds its own polynomial kernel (as a cloud does for the gas) and the model
  density only lifts it where the galaxy really is star-forming.
* octaves of gradient noise (`gradient_noise.GradientNoise`, seeded by the galaxy seed and the octave), the first of
  wavelength `W0_PC` and each next half that. They are added to `s` and weighted by `4 s (1 - s)`, so they bite at the
  edge and in the thin gas, and leave the middle of a dense core alone to fade unevenly. Each octave's amplitude is
  `PERSISTENCE` times the one before, and adding one octave changes the value by at most the weight times that
  octave's amplitude (`octave_bound`), so a finer detail adds filaments and does not move the coarse blobs.

Every constant here shapes the look Boss asked for (soft edge, bulbous outline, how fast filaments come with zoom); none is a
measurement, and MAP.365's drawing part tunes them against the map. The model's own numbers (the gas, the young stars, the
clouds and the superbubbles) are read, never changed.
"""

import functools
import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import density as galaxyDensity, gradient_noise, nebula_field, nebula_volume, seed as galaxySeed, superbubble_field
from planetgen.physics.units import ly_to_pc
from planetgen.util import detmath

FIELDS = ("gas", "young")
"""tuple: The densities a nebula can be drawn from."""

FORM_LEVEL = 0.35
"""float: A point is inside a nebula of a field where `density_at` is at least this."""

W0_PC = 64.0
"""float: The wavelength of the first noise octave, parsecs: about a molecular cloud's diameter (the clouds are 12 to 120 pc across)."""

MAX_DETAIL = 8
"""int: The most octaves: the finest wavelength is `W0_PC / 2**7`, under 1 pc."""

PERSISTENCE = 0.5
"""float: Each octave's amplitude over the one before's."""

NOISE_AMPLITUDE = 0.4
"""float: The first octave's amplitude, in the units of the density (0..1), before the edge weight."""

GAS_BACKGROUND = 0.2
"""float: The most the model gas alone adds to the smooth density (its saturating share is 0..1): under `FORM_LEVEL`
with the noise at its edge weight, so gas alone never fills the sky."""

CLOUD_GAS_FLOOR = 0.8
"""float: The share of a cloud's kernel kept where the model gas is nil; the rest scales with the gas (0.2), so a cloud's
center is always dense enough (0.8) to hold a core through the first octaves."""

YOUNG_SCALE = 3.0
"""float: The young-star density, over its value on the solar circle, at which the smooth `"young"` density reaches half of
`YOUNG_CEILING`."""

YOUNG_CEILING = 0.6
"""float: The most the smooth `"young"` density reaches: under 1 so a star-forming region has a core to be denser than."""

YOUNG_ANCHOR_FLOOR = 0.8
"""float: The share of a young-family nebula's own kernel kept where the model young density is nil (GEN.1267, MAP.423); the rest scales with
the model density (0.2), the same split as `CLOUD_GAS_FLOOR` for a cloud, so a nebula's center is always dense enough to hold a core through
the first octaves and a nebula in a star-forming arm is a little denser than one in the thin outer disk."""

CONTAIN_DETAIL = 3
"""int: The noise octaves the containment test reads (GEN.1267): the finest wavelength is `W0_PC / 4`, 16 pc, so the outline is the coarse
form of the map's, stable and cheap, and a star drawn inside a nebula lies inside the outline the map draws at any closer view (a finer
octave adds filaments at the edge and does not move the coarse blobs, `octave_bound`)."""

CONTAIN_ROUND = 1e-4
"""float: The density is rounded to this before it is compared to `FORM_LEVEL` in the containment test, so a last-bit difference of a platform's
arithmetic cannot put a star on the other side of the outline."""

WARP_PC = 14.0
"""float: How far the smooth `"gas"` density is pushed about by a coarse, smooth warp before the clouds' kernels are read,
parsecs. It makes a cloud's outline lumpy even with no octave added (`detail` 0), as a cloud on the map has to be bulbous
and not a circle (Boss); the warp is a function of position, so the same point is pushed the same way in every node."""

WARP_WAVELENGTH_PC = 96.0
"""float: The wavelength of the warp, parsecs: longer than most clouds, so each is bent as a whole."""

BUBBLE_EDGE = 0.15
"""float: The share of a superbubble's radius over which the gas returns at its wall."""

_NOISE_OFFSET = 1013.0
"""float: How far apart the octaves' lattices are placed, in lattice cells, so that no two share a value at the origin."""


def wavelength_pc(detail):
    """The wavelength of octave `detail` (0 is the first), parsecs: `W0_PC / 2**detail`. A sample at `detail = n` has
    octaves 0 to n - 1, so it resolves features down to about `wavelength_pc(n - 1)`."""
    return W0_PC / 2.0 ** int(detail)


def octave_amplitude(octave):
    """The amplitude of noise octave `octave`: `NOISE_AMPLITUDE * PERSISTENCE ** octave`."""
    return NOISE_AMPLITUDE * PERSISTENCE ** int(octave)


def octave_bound(octave):
    """The most the value can change when octave `octave` is added: the edge weight is at most 1 and a gradient noise is
    inside +-1 (`gradient_noise.GradientNoise` is under 1.0 in size), so its amplitude."""
    return octave_amplitude(octave)


@functools.lru_cache(maxsize=16)
def _noise(galaxy_seed, octave):
    """The gradient noise of octave `octave` for a galaxy, from its own permutation table."""
    return gradient_noise.GradientNoise(galaxySeed.unit_seed(galaxy_seed or nebula_field.FALLBACK_SEED, "nebula-wisp", (octave,)))


@functools.lru_cache(maxsize=256)
def _cell_spheres(galaxy_seed, shape, index):
    """The kept clouds of field cell `index` as an array of rows `(x, y, z, radius_pc)`."""
    rows = [(*center, ly_to_pc(nebula.radius_ly)) for nebula, center in nebula_field.cell_clouds(galaxy_seed, shape, index)]
    return np.array(rows, dtype=float).reshape(-1, 4)


def _spheres_near(galaxy_seed, shape, low, high):
    """Every kept cloud whose sphere can reach the box `low` to `high` (3-tuples, parsecs), as an `(n, 4)` array."""
    edge = tuning.NEBULA_FIELD_CELL_PC
    span = nebula_field.MAX_CLOUD_RADIUS_PC
    first = nebula_field.cell_index(tuple(c - span for c in low))
    last = nebula_field.cell_index(tuple(c + span for c in high))
    parts = []
    for i in range(first[0], last[0] + 1):
        for j in range(first[1], last[1] + 1):
            for k in range(first[2], last[2] + 1):
                spheres = _cell_spheres(galaxy_seed, shape, (i, j, k))
                if len(spheres):
                    parts.append(spheres)
    return np.concatenate(parts) if parts else np.zeros((0, 4))


def _model_gas_share(points, shape):
    """The model's molecular gas at `points` as a saturating share `g / (1 + g)` of the solar circle's, 0..0.75 (the cap is 3)."""
    gas = nebula_field.gas_factor_array(points[..., 0], points[..., 1], points[..., 2], shape)
    return gas / (1.0 + gas)


def _bubble_clearing(galaxy_seed, shape, points, low, high):
    """1 outside every superbubble, falling to 0 over `BUBBLE_EDGE` of the radius inside the wall, at `points`."""
    keep = np.ones(points.shape[:-1], dtype=float)
    middle = tuple((a + b) / 2.0 for a, b in zip(low, high))
    reach = detmath.dist(low, high) / 2.0
    for bubble in superbubble_field.bubbles_reaching(galaxy_seed, shape, middle, reach):
        depth = np.linalg.norm(points - np.array(bubble.center_pc), axis=-1) / bubble.radius_pc
        inside = np.clip((1.0 - depth) / BUBBLE_EDGE, 0.0, 1.0)
        keep *= 1.0 - inside * inside * (3.0 - 2.0 * inside)
    return keep


def _warped(galaxy_seed, points):
    """`points` pushed about by `WARP_PC` of smooth gradient noise on three channels (octaves -1, -2 and -3 of the noise)."""
    lattice = points / WARP_WAVELENGTH_PC
    push = np.stack([_noise(galaxy_seed, -1 - channel)(lattice + _NOISE_OFFSET * (channel + 1)) for channel in range(3)], axis=-1)
    return points + WARP_PC * push


def _smooth_gas(galaxy_seed, shape, points, low, high):
    share = _model_gas_share(points, shape)
    reach = WARP_PC * 1.8
    spheres = _spheres_near(galaxy_seed, shape, tuple(c - reach for c in low), tuple(c + reach for c in high))
    kernel = np.zeros(points.shape[:-1], dtype=float)
    flat = _warped(galaxy_seed, points).reshape(-1, 3)
    if len(spheres):
        flat_kernel = np.zeros(len(flat), dtype=float)
        for sx, sy, sz, radius in spheres:
            offset = flat - np.array([sx, sy, sz])
            near = np.einsum("ij,ij->i", offset, offset) < radius * radius
            if near.any():
                inside = 1.0 - np.einsum("ij,ij->i", offset[near], offset[near]) / (radius * radius)
                flat_kernel[near] = 1.0 - (1.0 - flat_kernel[near]) * (1.0 - inside * inside)     # clouds add as a union
        kernel = flat_kernel.reshape(points.shape[:-1])
    kernel = kernel * (CLOUD_GAS_FLOOR + (1.0 - CLOUD_GAS_FLOOR) * share / 0.75)
    smooth = 1.0 - (1.0 - GAS_BACKGROUND * share / 0.75) * (1.0 - kernel)
    return smooth * _bubble_clearing(galaxy_seed, shape, points, low, high)


def _smooth_young(galaxy_seed, shape, points, anchors=None):
    """
    The smooth `"young"` density at `points`: the model's young-star density, and with `anchors` (`(centers, radii)`: the center `(..., 3)` or
    `(3,)` and the radius, parsecs, of the nebula each point is read for) the nebula's own kernel `(1 - d^2 / r^2)^2` on the warped point, weighted by
    `YOUNG_ANCHOR_FLOOR` and the model density, joined to the model as a union (GEN.1267, MAP.423).
    """
    reference = nebula_field._reference_gas(shape)
    young = np.maximum(galaxyDensity.population_densities_array(points[..., 0], points[..., 1], points[..., 2], shape)["young"], 0.0)
    if reference <= 0.0:
        model = np.zeros(points.shape[:-1])
    else:
        x = young / (YOUNG_SCALE * reference)
        model = YOUNG_CEILING * x / (1.0 + x)
    if anchors is None:
        return model
    centers, radii = anchors
    offset = _warped(galaxy_seed, points) - np.asarray(centers, dtype=float)
    inside = np.clip(1.0 - np.einsum("...i,...i->...", offset, offset) / np.square(np.asarray(radii, dtype=float)), 0.0, 1.0)
    kernel = inside * inside * (YOUNG_ANCHOR_FLOOR + (1.0 - YOUNG_ANCHOR_FLOOR) * model / YOUNG_CEILING)
    return 1.0 - (1.0 - model) * (1.0 - kernel)


def density_at(galaxy_seed, shape, points_pc, field, detail=0, anchors=None):
    """
    The nebula density at `points_pc`, in 0..1 (`float32`).

    Args:
        galaxy_seed (bytes): The galaxy's 16-byte seed.
        shape (GalaxyShape): The stored shape.
        points_pc (array-like): Points `[x, y, z]`, parsecs, galaxy frame, any leading shape with a last axis of 3.
        field (str): `"gas"` or `"young"` (`FIELDS`).
        detail (int): Noise octaves, 0 (the smooth model alone) to `MAX_DETAIL`; larger is clamped, a negative is 0.
        anchors (tuple, optional): For the `"young"` field, `(centers, radii)`: the center `(..., 3)` or `(3,)` and the radius in parsecs of the
            nebula each point is read for, which anchors the density on that nebula (GEN.1267). The `"gas"` field is anchored on the cloud field
            and ignores it. `None` reads the model alone.

    Returns:
        numpy.ndarray: One value per point, the shape of `points_pc` without its last axis. At least `FORM_LEVEL` is inside.

    Raises:
        ValueError: For a field that is not in `FIELDS` or points without a last axis of 3.
    """
    if field not in FIELDS:
        raise ValueError(f"field is one of {', '.join(FIELDS)}, not {field!r}")
    points = np.asarray(points_pc, dtype=float)
    if points.ndim < 1 or points.shape[-1] != 3:
        raise ValueError("points_pc needs a last axis of 3 ([x, y, z] in parsecs)")
    octaves = min(max(int(detail), 0), MAX_DETAIL)
    if points.size == 0:
        return np.zeros(points.shape[:-1], dtype=np.float32)
    flat = points.reshape(-1, 3)
    low, high = tuple(flat.min(axis=0)), tuple(flat.max(axis=0))
    if field == "gas":
        smooth = _smooth_gas(galaxy_seed, shape, points, low, high)
    else:
        smooth = _smooth_young(galaxy_seed, shape, points, anchors)
    smooth = np.clip(smooth, 0.0, 1.0)
    noise = np.zeros_like(smooth)
    for octave in range(octaves):
        lattice = points / wavelength_pc(octave) + _NOISE_OFFSET * octave
        noise += octave_amplitude(octave) * _noise(galaxy_seed, octave)(lattice)
    edge_weight = 4.0 * smooth * (1.0 - smooth)
    return np.clip(smooth + edge_weight * noise, 0.0, 1.0).astype(np.float32)


FIELD_OF_FAMILY = nebula_volume.FIELD_OF_FAMILY
"""dict: Which density a family of nebula is drawn from and contained by: the dark family (molecular gas) from the gas, the rest from the young stars
(`galaxy.nebula_volume`'s own table, so a map and a containment test read one density)."""


def field_of_class(nebula_class):
    """The density (`"gas"` or `"young"`) a nebula of class `nebula_class` (a letter of `tuning.NEBULA_CLASSES`) is drawn from."""
    return FIELD_OF_FAMILY[tuning.NEBULA_CLASSES[nebula_class]["family"]]


def inside(galaxy_seed, shape, points_pc, center_pc, radius_pc, field):
    """
    Which of `points_pc` lie inside the nebula at `center_pc` of radius `radius_pc` (GEN.1267): inside its sphere, and where the density, read for that
    nebula, is at least `FORM_LEVEL` at `CONTAIN_DETAIL` octaves. The map draws the same contour (`galaxy.nebula_volume`), so a star placed inside a
    nebula lies inside the outline a reader sees. Nothing about the nebula but its center, radius and field is read, and the object ID never is.

    Args:
        galaxy_seed (bytes): The galaxy's 16-byte seed.
        shape (GalaxyShape): The stored shape.
        points_pc (array-like): Points `[x, y, z]`, parsecs, galaxy frame, last axis 3.
        center_pc (tuple): The nebula's center, parsecs.
        radius_pc (float): The nebula's radius, parsecs.
        field (str): `"gas"` or `"young"` (`field_of_class`).

    Returns:
        numpy.ndarray: A bool for each point, the shape of `points_pc` without its last axis.
    """
    points = np.asarray(points_pc, dtype=float)
    center = np.asarray(center_pc, dtype=float)
    held = np.einsum("...i,...i->...", points - center, points - center) < radius_pc * radius_pc
    result = np.zeros(points.shape[:-1], dtype=bool)
    if not held.any():
        return result
    near = points[held]
    value = density_at(galaxy_seed, shape, near, field, CONTAIN_DETAIL, anchors=(center, radius_pc))
    result[held] = np.rint(value.astype(float) / CONTAIN_ROUND) >= round(FORM_LEVEL / CONTAIN_ROUND)
    return result
