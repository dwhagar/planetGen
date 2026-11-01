# planetgen/physics/ground_albedo.py

"""
Ground albedo by band (GEN.462, docs/design/body-composition-parts/E-ground-radiogenic-shielding.md sections 2.6 and 2.7)
===========================================================================================================================

The colour of a solid world's ground, as the reflectance in the three visible bands 0.4 to 0.5, 0.5 to 0.6 and 0.6 to 0.7 micrometre
that the materials table holds for 26 rows (laboratory powders, slabs, ices, evaporites, clays, oxides and carbons; [R3425] and the USGS
spectral library, averaged by band). A body's ground is a set of areas, each a row of the table (the ground types and the ice, from GEN.441 and GEN.442); the
band albedo is the area-weighted mean. The caller states the shares of the table's rows; a cover with no band albedo in the table (the CO2, CH4, NH3, N2 and
SO2 frosts: visible albedo not reached) is an error naming the row, never a guess, and a cover that is not a table row (open ocean) comes in as an explicit
`(b1, b2, b3)` in `extra`.

Space weathering darkens an airless surface that has not been refreshed: `A(t) = A_fresh [1 - d (1 - exp(-t / tau_w))]` with `t` the age of the surface
since it was resurfaced, `tau_w` log-uniform in 1 to 2 Gyr and `d` uniform in 0.3 to 0.7 (a parametrisation of the design note, placeholders to fit against lunar and
asteroid reflectance; the exponential decay of the unweathered fraction is [R4713]). A surface under an atmosphere is not darkened here.

A laboratory reflectance factor is neither a geometric nor a Bond albedo: `visible_mean` is the visible-band proxy that the temperature
(GEN.412) and the colour model (GEN.418) read, and the conversion to a Bond albedo over the whole spectrum is an open gap (master section 8 rank 9).

Use:

    albedo = band_albedo({"basalt_powder": 0.7, "water_ice": 0.3})        # BandAlbedo
    weathered = weather(albedo, surface_age_gy=3.0, tau_w_gy=1.5, depth=0.5)
"""

import math
from dataclasses import dataclass

from planetgen.physics import materials
from planetgen.util import draw

BANDS = ("alb_040_050", "alb_050_060", "alb_060_070")
"""tuple: The three property keys of the materials table, in order of wavelength."""
BAND_MICRONS = ((0.4, 0.5), (0.5, 0.6), (0.6, 0.7))

TAU_W_GY = (1.0, 2.0)
"""tuple: The space-weathering time, log-uniform (the published range is 1 Myr to 2 Gyr; this is the design note's placeholder)."""
DEPTH = (0.3, 0.7)
"""tuple: The darkening amplitude `d`, uniform (a placeholder of the design note)."""

SHARE_TOLERANCE = 1e-6


class MissingAlbedo(ValueError):
    """A cover whose row has no visible band albedo in the materials table."""


@dataclass(frozen=True)
class BandAlbedo:
    """The reflectance in the bands `BAND_MICRONS`, each in 0 to 1."""

    b1: float
    b2: float
    b3: float

    @property
    def values(self):
        return (self.b1, self.b2, self.b3)

    @property
    def visible_mean(self):
        """The mean of the three bands: the visible-band proxy for the Bond albedo (the conversion over the whole spectrum is an open gap)."""
        return (self.b1 + self.b2 + self.b3) / 3.0


def row_albedo(name):
    """The `BandAlbedo` of a row of the materials table.

    Raises:
        MissingAlbedo: The row has no value in one of the three bands.
    """
    row = materials.material(name)
    values = []
    for key in BANDS:
        try:
            value = row.prop(key).value
        except KeyError:
            value = None
        if value is None:
            raise MissingAlbedo(f"the materials row {row.id!r} has no visible band albedo ({key}); add it with its source or state the albedo in `extra`")
        values.append(float(value))
    return BandAlbedo(*values)


def band_albedo(cover, extra=None):
    """
    The area-weighted band albedo of a ground cover.

    Args:
        cover (dict): Row of the materials table -> area share (scaled to add up to 1).
        extra (dict or None): Name -> `(b1, b2, b3)` for a cover that is not a table row (the ocean); a name in `extra` is used before the table.

    Raises:
        ValueError: No cover, a negative share, or a share of zero in total.
        MissingAlbedo: A row of the cover has no band albedo.
    """
    extra = extra or {}
    if not cover:
        raise ValueError("no cover")
    total = math.fsum(cover.values())
    if not total > 0.0 or any(share < 0.0 for share in cover.values()):
        raise ValueError(f"the cover shares {cover!r} are negative or add up to zero")
    sums = [0.0, 0.0, 0.0]
    for name, share in cover.items():
        if share == 0.0:
            continue
        albedo = BandAlbedo(*extra[name]) if name in extra else row_albedo(name)
        for i, value in enumerate(albedo.values):
            sums[i] += share / total * value
    return BandAlbedo(*sums)


def weathering_factor(surface_age_gy, tau_w_gy, depth):
    """`1 - d (1 - exp(-t / tau_w))`: the share of the fresh albedo a surface of this age keeps (1 for a fresh surface)."""
    if surface_age_gy < 0.0 or not tau_w_gy > 0.0 or not 0.0 <= depth < 1.0:
        raise ValueError("the age must not be negative, the time positive and the amplitude in [0, 1)")
    return 1.0 - depth * (1.0 - math.exp(-surface_age_gy / tau_w_gy))


def weather(albedo, surface_age_gy, tau_w_gy, depth):
    """The `BandAlbedo` after space weathering (every band darkens by the same factor)."""
    factor = weathering_factor(surface_age_gy, tau_w_gy, depth)
    return BandAlbedo(*(value * factor for value in albedo.values))


def draw_weathering():
    """`(tau_w_gy, depth)` for one body: the time log-uniform in `TAU_W_GY`, then the amplitude uniform in `DEPTH`."""
    tau = math.exp(draw.uniform(math.log(TAU_W_GY[0]), math.log(TAU_W_GY[1])))
    return tau, draw.uniform(*DEPTH)
