# planetgen/generation/object_first.py

"""
The object-first scatter (PERF.58): decide how many objects a shell holds,
then where, instead of walking every band of every shell.

For one shell a certified density majorant `M_b` bounds the density anywhere
in band `b` of that shell. The pass draws one Poisson count `N` with mean
`e * sum(cells_b * M_b)` (`cells_b` the band's wedges weighted by their volume in
edge cubes), gives each of the `N` candidates a band in proportion to that
(a cumulative table), a wedge and a point in the cell, and keeps it with
probability `true density / M_b`. That is exactly a Poisson process with the
true density, so the counts and the spatial distribution match a band-by-band
walk; an empty shell costs only its majorant. Study:
`docs/design/scatter-queue-feasibility.md`.

A sector can take several objects (`sector_capacity`): independent draws, a
count per sector, and an object past the sector's capacity dropped. The
capacity is the smallest of `tuning.OBJECT_FIRST_CAPACITIES` for which a
Poisson count with the sector's expected objects exceeds it with a chance
below `tuning.OBJECT_FIRST_TAIL`.
"""

import itertools
import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import core_gas as coreGas, density, spherical
from planetgen.galaxy import star_formation as starFormation


def _sech_squared(x):
    """`density._sech_squared` on an array (`1 / cosh(x) ** 2`, never overflowing)."""
    e = np.exp(-2.0 * np.minimum(np.abs(x), 350.0))
    return np.minimum(1.0, 4.0 * e / (1.0 + e) ** 2)


def _vertical(z, scale_height_pc):
    """`density._vertical` on an array."""
    return _sech_squared(z / (2.0 * scale_height_pc))


def _disk_share_bound(shape, z_near, z_far, fractions, rho_near, rho_far):
    """
    An upper bound on `sum(share_p * fraction_p)` over the three thin-disk
    populations in a band's cells, wherever in them (`z_near` and `z_far`
    are arrays: the nearest and farthest `|z|` of each band, `rho_near` and `rho_far` the nearest and farthest
    cylindrical radius). A population's weight is its
    share of disk star formation, times its vertical profile (between the
    band's near and far edge: the profile falls with `|z|`), its arm
    contrast (between `1 - A` and `1 + A`) and, for the young and intermediate, the recent star-formation factor
    (GEN.134), whose true smallest and largest values over the radius range are taken
    (`star_formation.recent_factor_range`: it is not monotone, so the ends are not enough). The sum is a linear
    fraction of the weights, so over the box of their ranges it peaks at a corner.
    """
    recent_low, recent_high = starFormation.recent_factor_range(rho_near, rho_far, shape)
    ranges, f = [], []
    for name, ratio in density.population_height_ratios(shape).items():
        height = shape.disk_scale_height_pc * ratio
        share = starFormation.disk_population_share(name, shape)
        arm = abs(tuning.STELLAR_POPULATION_ARM_AMPLITUDE[name])
        low = share * _vertical(z_far, height) / height * max(1.0 - arm, 0.0)
        high = share * _vertical(z_near, height) / height * (1.0 + arm)
        if name in starFormation.RECENT_POPULATIONS:
            low, high = low * recent_low, high * recent_high
        ranges.append((low, high))
        f.append(fractions[name])
    # Far off the plane every profile underflows and the old disk takes it all (`population_densities`).
    best = np.full(np.shape(z_near), float(fractions["old"]))
    for corner in itertools.product(*ranges):
        total = sum(corner)
        mix = sum(weight * fraction for weight, fraction in zip(corner, f))
        with np.errstate(divide="ignore", invalid="ignore"):
            best = np.maximum(best, np.where(total > 0.0, mix / np.where(total > 0.0, total, 1.0), 0.0))
    return best


def band_geometry(shell, high_band, edge_pc):
    """
    For the bands `-high_band .. high_band` of a shell (index `band + high_band`), the extremes of
    each band's cells that bound a density falling with the distance from the axis and with `|z|`:
    `(rho_near, z_near, z_far)` arrays, the smallest cylindrical radius, the smallest `|z|` and the
    largest `|z|` of any point in the band.
    """
    r_in, r_out = spherical.shell_bounds_pc(shell, edge_pc)
    north, south = [], []
    for band in range(-high_band, high_band + 1):
        n, s = spherical.band_colatitude_bounds(shell, band)
        north.append(n)
        south.append(s)
    north, south = np.array(north), np.array(south)
    cos_n, cos_s = np.cos(north), np.cos(south)
    rho_near = r_in * np.minimum(np.sin(north), np.sin(south))
    spans = (cos_n > 0.0) & (cos_s < 0.0)
    z_near = np.where(spans, 0.0, np.where(cos_s >= 0.0, r_in * cos_s, r_in * np.abs(cos_n)))
    z_far = np.where(spans, r_out * np.maximum(cos_n, -cos_s),
                     np.where(cos_s >= 0.0, r_out * cos_n, r_out * np.abs(cos_s)))
    return rho_near, z_near, z_far


def band_majorants(shape, shell, high_band, edge_pc, fractions):
    """
    Per band `-high_band .. high_band` of a shell (index `band + high_band`), an upper bound on
    `sum_p density_p * fractions[p]` (the stars of the pass expected per unit of `e`) anywhere in the
    band's cells. Every component falls with the distance from the axis and with `|z|`, so the band's
    smallest cylindrical radius and its nearest `|z|` bound it; the arm term is bounded by `1 + |A|`,
    the populations' split by `_disk_share_bound`, the halo by its value at the band's nearest point (it falls
    with the spheroidal radius, which both the nearest radius and the nearest height lower).

    Returns:
        numpy.ndarray: One bound per band.
    """
    inner, z_near, z_far = band_geometry(shell, high_band, edge_pc)
    outer = spherical.shell_bounds_pc(shell, edge_pc)[1]
    r_in = spherical.shell_bounds_pc(shell, edge_pc)[0]
    return extent_majorants(shape, inner, z_near, z_far, fractions, outer, r_near=r_in)


def _nuclear_bounds(terms, inner, z_near, z_far, r_near):
    """
    The nuclear cluster's and the nuclear disc's density bounds (unnormalized) over a region whose smallest
    cylindrical radius is `inner`, whose `|z|` runs from `z_near` to `z_far` and whose smallest spherical radius is
    `r_near` (PERF.138). The cluster falls with the spherical radius alone, so it is bounded at the larger of `r_near`
    and the distance of the point `(inner, z_near)`. The disc falls with the cylindrical radius and the height each;
    the nearest point `(inner, z_near)` alone is far inside the region when the two ends come from different cells (a
    polar group has `inner` 0 beside a `z_near` far under its shell's radius), so the heights are cut into
    `NUCLEAR_BOUND_CUTS` slices and each slice is bounded at its lowest height with the smallest radius a point of it
    can have: `rho >= sqrt(r_near^2 - z_top^2)` for a slice that ends at `z_top`. The largest slice bounds the region.
    """
    reach = np.maximum(np.sqrt(inner ** 2 + z_near ** 2), r_near)
    cluster = density._nuclear_cluster_array(reach, terms)
    steps = (np.arange(1, NUCLEAR_BOUND_CUTS + 1) / NUCLEAR_BOUND_CUTS) ** 2
    disc = np.zeros(np.shape(inner))
    low = z_near
    for step in steps:
        top = z_near + (z_far - z_near) * step
        rho_low = np.maximum(inner, np.sqrt(np.maximum(r_near ** 2 - top ** 2, 0.0)))
        disc = np.maximum(disc, density._nuclear_disc_array(rho_low, low, terms))
        low = top
    return cluster, disc


NUCLEAR_BOUND_CUTS = 12
"""int: How many slices of height `_nuclear_bounds` cuts a region into for the nuclear disc."""


def extent_majorants(shape, inner, z_near, z_far, fractions, outer, r_near=0.0):
    """
    `band_majorants` for any extents (arrays): an upper bound on `sum_p density_p * fractions[p]` in a region
    whose smallest cylindrical radius is `inner`, whose largest is at most `outer` (its outer spherical radius is
    one) and whose `|z|` runs from `z_near` to `z_far`. A band is one such region and a group of bands
    (`ReachGroups`) another, bounded at its nearest point. `r_near` (arrays too) is the region's smallest spherical
    radius where it is known (a shell's, a group's `r_in`); the nuclear terms use it (`_nuclear_bounds`).
    """
    terms = density.model_terms(shape)
    thin = np.exp(-inner / shape.disk_scale_length_pc) * _vertical(z_near, shape.disk_scale_height_pc)
    thick = (terms["thick_disk_amplitude"] * np.exp(-inner / terms["thick_disk_scale_length_pc"])
             * _vertical(z_near, terms["thick_disk_scale_height_pc"]))
    along = inner / shape.bulge_scale_radius_pc
    up = z_near / terms["bulge_scale_z_pc"]
    bulge = shape.bulge_amplitude * np.exp(-0.5 * np.sqrt(along ** 4 + up ** 4))
    arm = shape.arm_level * (1.0 + abs(shape.arm_amplitude))
    disk_share = _disk_share_bound(shape, z_near, z_far, fractions, inner, outer)
    k = shape.k_norm
    halo = density._halo(inner, z_near, density.halo_terms(shape))
    nuclear_cluster, nuclear_disc = _nuclear_bounds(terms, inner, z_near, z_far, r_near)
    bound = (np.maximum(k * bulge, 0.0) * fractions["bulge"]
             + np.maximum(k * nuclear_cluster, 0.0) * fractions["nuclear_cluster"]
             + np.maximum(k * nuclear_disc, 0.0) * fractions["nuclear_disc"]
             + np.maximum(k * thick, 0.0) * fractions["old"]
             + np.maximum(k * halo, 0.0) * fractions["halo"]
             + np.maximum(k * thin, 0.0) * arm * disk_share)
    # The central molecular zone's young stars are taken out of the nuclear disc's density into the young population's
    # (GEN.282, `density.population_densities`): the ring's own density bounds what moved.
    bound = bound + k * coreGas.young_bound(inner, outer, z_near, coreGas.core_gas(shape)) * fractions["young"]
    return bound * (1.0 + 1e-12)


class ReachGroups:
    """
    Regions of the sphere beyond the galaxy's outline, each bounded in one piece (GEN.322): a range of colatitude of
    one shell, north or south, so the pass bounds a handful of regions instead of every band of a shell of 100,000.

    `for_shell` cuts the part of one shell beyond the outline into groups of up to
    `tuning.OBJECT_FIRST_GROUP_BANDS` bands: walking from the outline toward the pole, a group takes the next 64
    bands, and the last also takes the polar cap. The first group starts exactly at the outline's last band edge, so
    the groups and the bands inside the outline (`-high_band .. high_band`) tile the shell with no cell in both. A
    shell the outline does not reach (`high_band` below 0) is cut from the equator, the equatorial band split between
    the two halves. `far` cuts many whole shells at once into `pieces` equal steps of colatitude each, for the far
    shells that hold almost nothing, so that a plan does not visit them one by one.

    Attributes (numpy arrays, one entry per group, north groups first):
        shell_low, shell_high (numpy.ndarray): The first and last shell the group spans (one for a group of bands).
        r_in, r_out (numpy.ndarray): The group's inner and outer radius, parsecs.
        inner, z_near, z_far (numpy.ndarray): The group's smallest cylindrical radius and its smallest and largest
            `|z|`, parsecs: `band_geometry`'s for the group.
        cos_low, cos_high (numpy.ndarray): The group's range of `cos(colatitude)`.
        colat_near, colat_far (numpy.ndarray): The group's colatitude nearest the polar axis and nearest the outline (or the
            equator), whichever side: north groups hold `colat_near <= colat < colat_far`, south groups the mirror.
        volume (numpy.ndarray): The group's volume in edge cubes (the sum of its cells').
    """

    def __init__(self, shell_low, shell_high, r_in, r_out, hi, lo, edge_pc):
        """The groups of colatitude `lo .. hi` (north, `hi` on the equator's side; arrays of one entry per group
        with the shells they span and their radii alongside), and their south mirror images."""
        north_inner = r_in * np.sin(lo)
        north_near = r_in * np.cos(hi)
        north_far = r_out * np.cos(lo)
        volume = (r_out ** 3 - r_in ** 3) / 3.0 * (np.cos(lo) - np.cos(hi)) * 2.0 * math.pi / edge_pc ** 3
        self.shell_low = np.concatenate((shell_low, shell_low))
        self.shell_high = np.concatenate((shell_high, shell_high))
        self.r_in = np.concatenate((r_in, r_in))
        self.r_out = np.concatenate((r_out, r_out))
        self.inner = np.concatenate((north_inner, north_inner))
        self.z_near = np.concatenate((north_near, north_near))
        self.z_far = np.concatenate((north_far, north_far))
        self.cos_low = np.concatenate((np.cos(hi), -np.cos(lo)))
        self.cos_high = np.concatenate((np.cos(lo), -np.cos(hi)))
        self.volume = np.concatenate((volume, volume))
        self.colat_near = np.concatenate((lo, np.pi - lo))
        self.colat_far = np.concatenate((hi, np.pi - hi))

    def __len__(self):
        return len(self.volume)

    @classmethod
    def for_shell(cls, shell, high_band, edge_pc, pieces=None):
        """The part of `shell` beyond an outline that reaches `high_band` (see the class), cut in groups of bands, or
        (`pieces`) in that many equal steps of colatitude between the outline and the pole."""
        hi = np.zeros(0)
        if shell > 0:
            layout = spherical.shell_layout(shell)
            if high_band < layout.high_band + 1:
                # `hi` is each group's colatitude on the outline's side: the outline's last band edge toward the pole
                # (the equator for a shell it does not reach), then every 64th band edge on toward the pole.
                first = layout.high_band - high_band if high_band >= 0 else layout.high_band
                indices = np.arange(first, 0, -tuning.OBJECT_FIRST_GROUP_BANDS) if first > 0 else np.zeros(1, dtype=int)
                hi = layout.cap_colatitude + indices * layout.band_step
                if high_band < 0:
                    hi = np.append(math.pi / 2, hi)
                if pieces is not None:
                    hi = hi[0] * (1.0 - np.arange(pieces) / pieces)
        lo = np.append(hi[1:], 0.0) if len(hi) else hi
        r_in, r_out = spherical.shell_bounds_pc(shell, edge_pc)
        count = len(hi)
        return cls(np.full(count, shell), np.full(count, shell), np.full(count, r_in), np.full(count, r_out), hi, lo,
                   edge_pc)

    @staticmethod
    def far_blocks(first_shell, last_shell):
        """The runs of shells `far` cuts `first_shell .. last_shell` into, as `(first, last)` pairs: each as many shells
        as 1/`tuning.OBJECT_FIRST_FAR_STEPS` of its distance, so its density changes by about that much across it."""
        blocks, start = [], first_shell
        while start <= last_shell:
            end = min(start + max(1, start // tuning.OBJECT_FIRST_FAR_STEPS) - 1, last_shell)
            blocks.append((start, end))
            start = end + 1
        return blocks

    @classmethod
    def far(cls, first_shell, last_shell, edge_pc, pieces=tuning.OBJECT_FIRST_FAR_PIECES):
        """Shells `first_shell .. last_shell` (all past the outline and at least 1), cut into runs of shells
        (`far_blocks`) and each run into `pieces` equal steps of colatitude from the equator to the pole: a few
        thousand regions for 60,000 shells, for the far shells that hold almost nothing, so that a plan does not
        visit them one by one."""
        blocks = np.array(cls.far_blocks(first_shell, last_shell), dtype=float).reshape(-1, 2)
        steps = np.pi / 2 * (1.0 - np.arange(pieces) / pieces)
        hi = np.tile(steps, len(blocks))
        lo = np.tile(np.append(steps[1:], 0.0), len(blocks))
        low, high = np.repeat(blocks[:, 0], pieces), np.repeat(blocks[:, 1], pieces)
        return cls(low.astype(int), high.astype(int), (low - 0.5) * edge_pc, (high + 0.5) * edge_pc, hi, lo, edge_pc)

    def probes(self):
        """`(plane, z)` arrays of the two places each group is judged at, mid-shell on its two edges of colatitude:
        `colat_near`'s then `colat_far`'s, one array each."""
        r_mid = (self.r_in + self.r_out) / 2.0
        near = (r_mid * np.sin(self.colat_near), r_mid * np.cos(self.colat_near))
        far = (r_mid * np.sin(self.colat_far), r_mid * np.cos(self.colat_far))
        return near, far


def sector_capacity(expected_objects):
    """
    How many objects a sector may take: the smallest of
    `tuning.OBJECT_FIRST_CAPACITIES` that a Poisson count with mean
    `expected_objects` exceeds with a chance below `tuning.OBJECT_FIRST_TAIL`
    (the largest, if none does).
    """
    capacities = tuning.OBJECT_FIRST_CAPACITIES
    if expected_objects <= 0.0:
        return capacities[0]
    term = math.exp(-expected_objects)
    cumulative = term
    seen = 0
    for capacity in capacities:
        while seen < capacity:
            seen += 1
            term *= expected_objects / seen
            cumulative += term
        if 1.0 - cumulative < tuning.OBJECT_FIRST_TAIL:
            return capacity
    return capacities[-1]

