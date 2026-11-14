# planetgen/galaxy/sky_band.py

"""
The glow band of the unresolved stars in the sky of a system (VIEW.27, 8.2 stream, goal #4907): the light of the stars too
faint to be drawn one by one, as a small surface-brightness image of the whole sky, worked out by marching along each line
of sight through the galaxy's stellar density with the galactic dust in the way, and kept per 50 pc cell of the
observer's position.

`sky-view.md` 3.1: for a line of sight with direction `u` from the observer `P`,

    J(u) = integral of  j0 * density(P + s u) * f_faint(s) * 10^(-0.4 A_V(s))  ds        [L_V,sun / pc^2]
    SB   = -0.19 - 2.5 log10(Omega_arcsec2 * J)                                         [mag / arcsec^2]

* `density` is the galaxy's own stellar density in units of the real local one (`galaxy.density.total_density_array`), so
  the band follows the model galaxy's disk, bulge, arms and halo; `j0` (`LOCAL_V_LUMINOSITY_DENSITY`) is the V light of
  a unit of it, 0.044 L_V,sun per cubic parsec (stand-in SKB1).
* `A_V(s)` is the dust of the smooth disk (`physics.sky_photometry.dust_av_per_pc`, VIEW.20), summed along the same ray
  (the midpoint rule, half a step into the step it dims). Galactic dust only: the dust of nebulae stays with VIEW.13.
* **The subtraction of the resolved stars** is `f_faint(s)`: the share of the V light of a stellar population that comes
  from stars too faint to reach the limiting magnitude `m_lim` at the distance `s`, behind the dust of the ray. A star of
  absolute V magnitude `M` looks like `M + 5 log10(s / 10 pc) + A_V(s)`, so the stars that are NOT drawn are those fainter
  than `M_cut(s) = m_lim - 5 log10(s / 10 pc) - A_V(s)`, and `f_faint` is `faint_share(M_cut)` of the table
  `FAINT_LIGHT_TABLE` (stand-in SKB2). The band therefore adds to the drawn stars without counting them twice: a brighter
  limit (a lower `m_lim`, as on a phone) draws fewer stars and leaves more light in the glow.
* **Steps** are spaced as the square of a uniform fraction (fine near the observer, where the dust and the faint-star
  cut change fastest, and coarse far away where the density is smooth), `STEPS` of them out to `MAX_DISTANCE_PC`.
* **The image** is a cube map, six faces of `FACE_PIXELS` by `FACE_PIXELS` pixels (`face_directions`), each pixel the
  surface brightness in mag/arcsec^2 along the direction of its center (`DARK_MAG` where no light arrives). The shader
  looks it up with a bilinear filter, so a small map is enough for a band this smooth.
* **The cache** keeps one image per cell of `CELL_PC` = 50 pc of the observer's position, worked out from the cell's
  center, so every system in the cell shares it (moving the observer 20 ly changes the band by at most 0.004 mag in the
  Milky Way's disk, the plan's figure). The key (`cache_key`) names the cell, the galaxy (the caller's key: the galaxy's
  id, seed and epoch), the limiting magnitude, the image size, the steps and every constant the image depends on, so
  a change of any of them is a miss and not a stale hit. `BandStore` keeps the images in memory or in `.npy` files.

Nothing in the band is stored in the galaxy; a band is a cache, worked out again when it is missing.
"""

import hashlib
import json
import math
import os
import tempfile

import numpy as np

from planetgen.galaxy import density
from planetgen.physics import sky_photometry
from planetgen.util import detmath

CELL_PC = 50.0
"""float: Side of the observer cells that share a band (parsecs)."""

FACE_PIXELS = 32
"""int: Pixels along an edge of a cube face (6 faces, so 6144 directions): about 3 degrees a pixel."""

STEPS = 400
"""int: Steps of a ray, spaced as the square of a uniform fraction."""

MAX_DISTANCE_PC = 30000.0
"""float: How far a ray goes (parsecs); the light past it, in a galaxy of the Milky Way's size, is negligible."""

DEFAULT_LIMITING_MAGNITUDE = 7.5
"""float: The limiting magnitude of the drawn stars when the caller gives none (VIEW.19's default)."""

LOCAL_V_LUMINOSITY_DENSITY = 0.044
"""float: The V-band luminosity density (L_V,sun per cubic parsec) at a relative density of 1: sky-view.md 3.1, from the 3,081
HYG stars within 25 pc (the literature is about 0.05). Stand-in SKB1."""

SB_ZERO_POINT = -0.19
"""float: The constant of the surface brightness in mag/arcsec^2 for a flux in L_V,sun/pc^2 (sky-view.md 3.1)."""

OMEGA_ARCSEC2 = (math.pi / 648000.0) ** 2
"""float: One square arcsecond in steradians."""

DARK_MAG = 40.0
"""float: The surface brightness stored where no light arrives."""

FAINT_LIGHT_TABLE = (
    (-10.0, 1.0), (-0.19, 0.81), (2.31, 0.35), (4.81, 0.07), (9.0, 0.01), (12.0, 0.0),
)
"""tuple: `(M_cut, share)`: the share of the V light of the stars near the Sun that comes from stars fainter than absolute
V magnitude `M_cut`. The middle three rows are sky-view.md 3.1 (stars dimmer than 1, 10 and 100 L_sun carry 7, 35 and 81
percent of the local V light, from HYG within 25 pc); the ends are invented (stand-in SKB2)."""

MODEL_VERSION = 1
"""int: Raised when the arithmetic of the band changes, so old cached images miss."""

SB_REF = 22.0
"""float: The surface brightness of the page's unit (`html/static/skyband.js` `SB_REF`): a texel of 1.0 is 22 mag/arcsec^2."""

EQUIRECT_WIDTH, EQUIRECT_HEIGHT = 256, 128
"""int: The size of the equirectangular glow image the System View samples (the contract of `skyband.js`)."""

TINT_TEMPERATURE_K = 4500.0
"""float: The temperature of the blackbody whose colour tints the glow (stand-in GLW1): the unresolved light is mostly K and M dwarfs and the cool
giants of the bulge."""


def face_directions(pixels=FACE_PIXELS):
    """
    The unit vector through the center of every pixel of a cube map, shape `(6, pixels, pixels, 3)`. The faces are +x,
    -x, +y, -y, +z, -z; on each, `row` and `column` run through -1..1 across the face (pixel centers).
    """
    if pixels < 1:
        raise ValueError("a face needs at least one pixel")
    centers = (np.arange(pixels) + 0.5) / pixels * 2.0 - 1.0
    column, row = np.meshgrid(centers, centers)
    one = np.ones_like(column)
    faces = np.stack([
        np.stack([one, row, column], axis=-1), np.stack([-one, row, column], axis=-1),
        np.stack([column, one, row], axis=-1), np.stack([column, -one, row], axis=-1),
        np.stack([column, row, one], axis=-1), np.stack([column, row, -one], axis=-1),
    ])
    return faces / np.linalg.norm(faces, axis=-1, keepdims=True)


def _face_coordinates(directions):
    """The cube face (0 to 5, `axis * 2 + negative`) each of `directions` (`(..., 3)`) falls on and its `(column, row)` on that face
    in -1 to 1, the two other axes in the order `face_directions` puts them."""
    d = np.asarray(directions, dtype=float)
    a = np.abs(d)
    axis = np.argmax(a, axis=-1)
    sign = np.take_along_axis(d, axis[..., None], axis=-1)[..., 0] < 0
    major = np.take_along_axis(a, axis[..., None], axis=-1)[..., 0]
    column_axis = np.where(axis == 0, 2, 0)
    row_axis = np.where(axis == 1, 2, 1)
    column = np.take_along_axis(d, column_axis[..., None], axis=-1)[..., 0] / major
    row = np.take_along_axis(d, row_axis[..., None], axis=-1)[..., 0] / major
    return axis * 2 + sign.astype(int), column, row


def face_lookup(band, directions):
    """
    The surface brightness of `band` (the `(6, n, n)` array `band_for` returns) in the directions `directions`, an array
    `(..., 3)`, by the nearest pixel (the shader filters bilinearly; this is for tests and the CPU paths).
    """
    band = np.asarray(band)
    n = band.shape[1]
    face, column, row = _face_coordinates(directions)
    ci = np.clip(((column + 1.0) / 2.0 * n).astype(int), 0, n - 1)
    ri = np.clip(((row + 1.0) / 2.0 * n).astype(int), 0, n - 1)
    return band[face, ri, ci]


def faint_share(absolute_cut):
    """The share of the V light that comes from stars fainter than absolute V magnitude `absolute_cut` (arrays too): the
    table `FAINT_LIGHT_TABLE`, linear between rows, 1 below the first and 0 past the last."""
    cuts = [row[0] for row in FAINT_LIGHT_TABLE]
    shares = [row[1] for row in FAINT_LIGHT_TABLE]
    return np.interp(np.asarray(absolute_cut, dtype=float), cuts, shares)


def step_distances(steps=STEPS, max_distance_pc=MAX_DISTANCE_PC):
    """The midpoint distances `s` and the step widths `ds` of a ray, parsecs: `s = max * f^2` for `f` in (0, 1)."""
    fraction = (np.arange(steps) + 0.5) / steps
    return max_distance_pc * fraction ** 2, 2.0 * max_distance_pc * fraction / steps


def cell_of(observer_pc):
    """The cell `(i, j, k)` of the 50 pc grid that holds `observer_pc` (galactocentric parsecs)."""
    p = np.asarray(observer_pc, dtype=float).reshape(3)
    return tuple(int(v) for v in np.floor(p / CELL_PC))


def cell_center(cell):
    """The center of a cell of `cell_of`, parsecs."""
    return tuple((v + 0.5) * CELL_PC for v in cell)


def cache_key(cell, galaxy_key, limiting_magnitude=DEFAULT_LIMITING_MAGNITUDE, pixels=FACE_PIXELS, steps=STEPS,
              max_distance_pc=MAX_DISTANCE_PC, dust=None):
    """
    The cache key of a band: a hash of the cell, the caller's `galaxy_key` (a string naming the galaxy's id, seed and
    epoch), the limiting magnitude, the image size, the steps, the dust law's overrides (a `physics.dust_model.DustProfile` by its `key()`, the digest of its two tables) and the constants the image
    depends on (`MODEL_VERSION`, the luminosity density, the faint-light table, the cell size).
    """
    payload = {
        "version": MODEL_VERSION, "cell": [int(v) for v in cell], "cell_pc": CELL_PC, "galaxy": str(galaxy_key),
        "m_lim": float(limiting_magnitude), "pixels": int(pixels), "steps": int(steps), "max_pc": float(max_distance_pc),
        "j0": LOCAL_V_LUMINOSITY_DENSITY, "faint": [list(row) for row in FAINT_LIGHT_TABLE],
        "dust": sorted((str(k), v.key() if hasattr(v, "key") else float(v)) for k, v in (dust or {}).items()),
        "dust_default": [sky_photometry.DUST_AV_PER_KPC_AT_SUN, sky_photometry.DUST_SCALE_LENGTH_PC,
                         sky_photometry.DUST_HEIGHT_PC, sky_photometry.DUST_HEIGHT_FLARE_PC_PER_KPC,
                         sky_photometry.DUST_FLARE_START_PC, sky_photometry.DUST_SOLAR_RADIUS_PC],
    }
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()[:32]


def surface_brightness(observer_pc, shape, limiting_magnitude=DEFAULT_LIMITING_MAGNITUDE, pixels=FACE_PIXELS, steps=STEPS,
                       max_distance_pc=MAX_DISTANCE_PC, dust=None, density_at=None, chunk=1024):
    """
    The glow band seen from `observer_pc` (galactocentric parsecs): the ray-march of the module doc, a float32 array
    `(6, pixels, pixels)` of surface brightness in mag/arcsec^2 (`DARK_MAG` where no light arrives).

    Args:
        observer_pc: `(x, y, z)` of the observer.
        shape: The galaxy's `GalaxyShape`.
        limiting_magnitude (float): The faintest star that is drawn; the glow holds the light of the stars past it.
        pixels (int): Pixels along a cube face.
        steps (int): Steps of a ray.
        max_distance_pc (float): How far a ray goes.
        dust (dict): Overrides of the dust law (`sky_photometry.dust_av_per_pc` keywords), for a galaxy of another size.
        density_at: A function `(x, y, z, shape)` of arrays that gives the relative stellar density, in place of
            `density.total_density_array` (for a test with a known density).
        chunk (int): Rays marched at once, to bound the memory.
    """
    dust = dict(dust or {})
    density_at = density_at or density.total_density_array
    observer = np.asarray(observer_pc, dtype=float).reshape(3)
    directions = face_directions(pixels).reshape(-1, 3)
    s, ds = step_distances(steps, max_distance_pc)
    log_distance = 5.0 * np.log10(s / 10.0)
    out = np.empty(len(directions))
    for start in range(0, len(directions), chunk):
        u = directions[start:start + chunk]
        x = observer[0] + u[:, 0:1] * s[None, :]
        y = observer[1] + u[:, 1:2] * s[None, :]
        z = observer[2] + u[:, 2:3] * s[None, :]
        rho = density_at(x, y, z, shape)
        step_av = sky_photometry.dust_av_per_pc(x, y, z, **dust) * ds[None, :]
        a_v = np.cumsum(step_av, axis=1) - step_av / 2.0
        cut = limiting_magnitude - log_distance[None, :] - a_v
        light = LOCAL_V_LUMINOSITY_DENSITY * rho * faint_share(cut) * 10.0 ** (-0.4 * a_v) * ds[None, :]
        out[start:start + chunk] = light.sum(axis=1)
    flux = np.maximum(out, 1e-300)
    mags = SB_ZERO_POINT - 2.5 * np.log10(OMEGA_ARCSEC2 * flux)
    return np.minimum(mags, DARK_MAG).reshape(6, pixels, pixels).astype(np.float32)


class MemoryBandStore:
    """Bands kept in memory, by cache key."""

    def __init__(self):
        self._bands = {}

    def get(self, key):
        """The band stored under `key`, or `None`."""
        return self._bands.get(key)

    def put(self, key, band):
        """Stores `band` under `key`."""
        self._bands[key] = np.array(band, dtype=np.float32)


class FileBandStore:
    """Bands kept as `.npy` files in `directory`, one per cache key; a file is written whole and moved into place, so a
    reader never sees half of one."""

    def __init__(self, directory):
        self.directory = os.fspath(directory)

    def _path(self, key):
        return os.path.join(self.directory, f"{key}.npy")

    def get(self, key):
        """The band stored under `key`, or `None` (also when the file cannot be read)."""
        try:
            return np.load(self._path(key))
        except (OSError, ValueError):
            return None

    def put(self, key, band):
        """Stores `band` under `key`."""
        os.makedirs(self.directory, exist_ok=True)
        handle, temp = tempfile.mkstemp(dir=self.directory, suffix=".tmp")
        try:
            with os.fdopen(handle, "wb") as file:
                np.save(file, np.asarray(band, dtype=np.float32))
            os.replace(temp, self._path(key))
        except BaseException:
            if os.path.exists(temp):
                os.remove(temp)
            raise


def band_for(observer_pc, shape, galaxy_key, store, limiting_magnitude=DEFAULT_LIMITING_MAGNITUDE, pixels=FACE_PIXELS,
             steps=STEPS, max_distance_pc=MAX_DISTANCE_PC, dust=None, density_at=None):
    """
    The band for an observer: the stored image of its 50 pc cell when there is one, else worked out from the cell's center
    (`surface_brightness`), stored and returned. `galaxy_key` names the galaxy for the cache (see `cache_key`).
    """
    cell = cell_of(observer_pc)
    key = cache_key(cell, galaxy_key, limiting_magnitude, pixels, steps, max_distance_pc, dust)
    band = store.get(key)
    if band is not None:
        return band
    band = surface_brightness(cell_center(cell), shape, limiting_magnitude, pixels, steps, max_distance_pc, dust,
                              density_at)
    store.put(key, band)
    return band


def values_of(surface_brightness_mag):
    """The page's unit of each surface brightness (arrays): `10^(-0.4 (SB - SB_REF))`, and 0 where no light arrives (`DARK_MAG` or more)."""
    sb = np.asarray(surface_brightness_mag, dtype=float)
    return np.where(sb >= DARK_MAG, 0.0, 10.0 ** (-0.4 * (sb - SB_REF)))


def tint():
    """The colour of the glow, linear red, green and blue with a luminance (0.2126, 0.7152, 0.0722) of 1: the blackbody of `TINT_TEMPERATURE_K`."""
    from planetgen.galaxy import sky
    rgb = sky.star_colours([TINT_TEMPERATURE_K])[0]
    return rgb / detmath.dot(rgb, (0.2126, 0.7152, 0.0722))


def face_bilinear(values, directions):
    """
    `values` (the `(6, n, n)` array of `values_of` a band) in the directions `directions` (`(..., 3)`), bilinear within the face each direction falls on
    (clamped at the face's edge; a pixel is about 3 degrees, so the seams between faces are not seen).
    """
    values = np.asarray(values, dtype=float)
    n = values.shape[1]
    face, column, row = _face_coordinates(directions)
    x = np.clip((column + 1.0) / 2.0 * n - 0.5, 0.0, n - 1.0)
    y = np.clip((row + 1.0) / 2.0 * n - 0.5, 0.0, n - 1.0)
    x0, y0 = np.floor(x).astype(int), np.floor(y).astype(int)
    x1, y1 = np.minimum(x0 + 1, n - 1), np.minimum(y0 + 1, n - 1)
    fx, fy = x - x0, y - y0
    top = values[face, y0, x0] * (1.0 - fx) + values[face, y0, x1] * fx
    bottom = values[face, y1, x0] * (1.0 - fx) + values[face, y1, x1] * fx
    return top * (1.0 - fy) + bottom * fy


def equirect(band, width=EQUIRECT_WIDTH, height=EQUIRECT_HEIGHT):
    """
    The glow image of the page (`html/static/skyband.js`): a `(height, width, 3)` float32 array of linear light in the page's unit (`values_of`), tinted
    (`tint`), equirectangular in the sky's galactic axes: row 0 at the north pole, column 0 at longitude -180, texel centers half a texel in,
    longitude `atan2(y, x)`, latitude `asin(z)`. `band` is the `(6, n, n)` surface brightness of `band_for`.
    """
    column = (np.arange(width) + 0.5) / width * 2.0 * math.pi - math.pi
    latitude = (0.5 - (np.arange(height) + 0.5) / height) * math.pi
    lon, lat = np.meshgrid(column, latitude)
    directions = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)
    level = face_bilinear(values_of(band), directions)
    return (level[..., None] * tint()[None, None, :]).astype(np.float32)


def to_half_bytes(image):
    """The image as little-endian half floats, red, green and blue a texel, row by row: what `skyband.js` puts in its texture."""
    return np.asarray(image, dtype="<f2").tobytes()
