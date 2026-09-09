# stellarObjects/galaxySkeleton.py

"""
Galaxy-Wide Density Skeleton: Per-Shell Candidate Bands
===========================================================

Finds, for a given `galaxyDensity.GalaxyShape`, the compact structural
facts worth precomputing about the whole galaxy -- not a per-sector
plan. See `docs/design/galaxy-coordinate-system.md` section 9's own
storage-analysis addendum for the reasoning this implements: a sector's
own position, density, and vertices are pure deterministic functions of
`(shell_index, shell_slot_index)` and a handful of shape parameters
(`galaxyDensity`/`sectorGeometry`), so none of that needs precomputing or
storing -- it's cheaper to recompute on demand (`galaxyGen.
ensure_sector_generated`) than to look up. What *does* need precomputing
is where the galaxy has any content at all, which is exactly what
`find_shell_bands` finds, one shell at a time.

**Why a shell's qualifying region is (almost always) one band.** A
sector's `relative_density` is `bulge(r_3d) + disk_radial(r_cyl) *
f_z(z) * arm_factor(r_cyl, theta)` (`galaxyDensity._raw_density`).
Within one shell, `r_3d = r_k` is fixed, so `bulge` is a constant --
either it alone already clears the qualification threshold everywhere in
the shell (`whole_shell_qualifies`, an O(1) check, no scan needed), or it
doesn't and the disk/arm term decides. That term's *exact* maximum over
`theta` at fixed `(r_k, phi)` is `(1 + arm_amplitude)` (`arm_factor`'s own
range), independent of `theta` -- so `bound_relative_density_at_phi` is an
*exact* upper bound (not sampled) on what any slot at that `phi` could
reach, for whatever `theta` it happens to land at. Because `r_cyl = r_k *
sin(phi)` and `z = r_k * cos(phi)` are both symmetric under `phi -> pi -
phi`, this bound is symmetric about the galactic plane (`phi = pi/2`);
for every real-scale parameter choice checked against this module during
its own development, it was also unimodal (strictly falling away from the
plane toward either pole), giving a single contiguous qualifying band
centered on the plane. `find_shell_bands` doesn't *assume* this -- it
scans and reports however many bands it actually finds -- but it explains
why one row per shell is normal, not a coincidence of the parameters
tried so far.

**Why this is a safe superset, not an exact membership list.** The bound
above is exact as an upper bound, but a specific slot's own `theta` might
land it well below that bound (an inter-arm trough near the edge of the
band, say) -- so a slot inside the returned band range isn't guaranteed
to individually qualify. That's deliberate: this module answers "which
region could possibly hold content", cheaply, once, for the whole galaxy;
the *exact* per-slot answer (a single `relative_density` evaluation) is
cheap enough to defer to the moment that slot is actually visited
(`galaxyGen.ensure_sector_generated`), so there is no reason to pay for
it here.

Pure and side-effect-free, like `galaxyGeometry.py`/`galaxyDensity.py` --
no database or I/O; `galaxyPlan.py` owns persistence and parallelizing
the per-shell calls this module exposes.
"""

import math
from collections import namedtuple

from . import program_constants
from .galaxyDensity import _sech_squared
from .galaxyGeometry import shell_radius_pc, shell_sector_count, slot_index_bounds_for_phi_range
from .spaceSector import SpaceSector

PHI_SCAN_SAMPLES = 256
"""int: How many evenly-spaced samples of `phi` in `[0, pi]`
`find_shell_bands` coarse-scans before bisecting each detected sign
change to precision -- cheap insurance against a band too narrow for a
coarser grid to notice at all (the design doc's own note on this exact
tradeoff recommended 500 for a real build after finding 200 adequate;
256 splits the difference, still trivial cost per shell)."""

PHI_BISECT_TOLERANCE = 1e-9
"""float: How closely (radians) `find_shell_bands` refines each band
edge found by the coarse scan above -- far finer than any shell's own
slot-to-slot angular spacing, so it never affects which slot index a
refined edge maps to; only the coarse scan's resolution can miss a band
entirely."""

ShellBand = namedtuple("ShellBand", ["slot_index_min", "slot_index_max"])
"""One contiguous candidate slot-index range within a single shell,
inclusive on both ends -- see the module docstring for what "candidate"
means here."""


def expected_system_count_at_density_1(edge_ly=program_constants.DEFAULT_SECTOR_EDGE_LY):
    """
    The `E` constant every qualification check compares against: how many
    systems a sector this size would hold at `relative_density = 1.0`
    (real local stellar density) -- `SpaceSector.expected_system_count()`,
    evaluated once rather than per call.

    Args:
        edge_ly (float): The sector edge length, in light-years.

    Returns:
        float: Expected system count at density 1.
    """
    return SpaceSector(name="galaxySkeleton-calibration", edge_ly=edge_ly).expected_system_count()


def _bound_raw_density_at_phi(shape, r_k, phi):
    """
    The exact upper bound, over every possible `theta`, of `_raw_density`
    (`galaxyDensity`'s pre-`k_norm` value) at a fixed `(r_k, phi)` -- see
    the module docstring for why `arm_factor`'s own maximum
    (`1 + arm_amplitude`, independent of `theta`) makes this exact rather
    than approximate.

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        r_k (float): This shell's nominal radius, parsecs
                     (`galaxyGeometry.shell_radius_pc`).
        phi (float): Polar angle, radians (`0` = north pole, `pi/2` =
                     galactic plane, `pi` = south pole).

    Returns:
        float: `>= 0`, unnormalized (before `shape.k_norm`).
    """
    r_cyl = r_k * math.sin(phi)
    z = r_k * math.cos(phi)
    bulge = shape.bulge_amplitude * math.exp(-r_k / shape.bulge_scale_radius_pc)
    disk_radial = math.exp(-r_cyl / shape.disk_scale_length_pc)
    f_z = _sech_squared(z / shape.disk_scale_height_pc)
    return bulge + disk_radial * f_z * (1.0 + shape.arm_amplitude)


def bound_relative_density_at_phi(shape, r_k, phi):
    """
    `_bound_raw_density_at_phi`, normalized by `shape.k_norm` -- directly
    comparable to a qualification threshold expressed in
    `relative_density` terms (e.g. `1.0 / expected_system_count_at_density_1(...)`).

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        r_k (float): This shell's nominal radius, parsecs.
        phi (float): Polar angle, radians.

    Returns:
        float: `>= 0`.
    """
    return shape.k_norm * _bound_raw_density_at_phi(shape, r_k, phi)


def whole_shell_qualifies(shape, r_k, threshold_rho):
    """
    Whether every slot in a shell at radius `r_k` trivially qualifies,
    regardless of `phi`/`theta` -- true exactly when the bulge term alone
    (constant across the whole shell, unlike the disk/arm term) already
    clears `threshold_rho` on its own. Since `relative_density = k_norm *
    (bulge + disk_term)` and `disk_term >= 0` always (`arm_factor`'s range
    is `[1 - arm_amplitude, 1 + arm_amplitude]`, both `>= 0` since
    `arm_amplitude` is meant to stay below `1`), `k_norm * bulge >=
    threshold_rho` alone guarantees the true value clears the threshold at
    every point in the shell too -- this is an exact fact, not a bound, so
    no scan is needed when it holds.

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        r_k (float): This shell's nominal radius, parsecs.
        threshold_rho (float): The `relative_density` a sector must reach
            to qualify (typically `1.0 / expected_system_count_at_density_1(...)`).

    Returns:
        bool
    """
    bulge_alone = shape.k_norm * shape.bulge_amplitude * math.exp(-r_k / shape.bulge_scale_radius_pc)
    return bulge_alone >= threshold_rho


def _bisect_crossing_to_qualifying_side(f, phi_neg, phi_pos, tolerance):
    """
    Refines a sign change of `f` to within `tolerance` radians via plain
    bisection, given one point already known to fall short
    (`f(phi_neg) < 0`) and one already known to clear it (`f(phi_pos) >=
    0`) -- `phi_neg`/`phi_pos` may be given in either numeric order.
    Always returns the qualifying-side point, the conservative direction
    for a candidate band boundary (see the module docstring: this band is
    meant to be a safe superset, never a slightly-too-narrow one).

    Args:
        f (callable): `phi -> float`; only its sign at the two starting
                      points is trusted, not monotonicity in between (a
                      real crossing is still found correctly either way,
                      since bisection only relies on the sign disagreeing
                      at the two current bounds).
        phi_neg (float): A point where `f < 0`.
        phi_pos (float): A point where `f >= 0`.
        tolerance (float): Radians; stops once the bracket is this narrow.

    Returns:
        float: `phi_pos`'s refined value -- within `tolerance` of the true
              crossing, on the `f >= 0` side.
    """
    while abs(phi_pos - phi_neg) > tolerance:
        mid = (phi_neg + phi_pos) / 2.0
        if f(mid) < 0:
            phi_neg = mid
        else:
            phi_pos = mid
    return phi_pos


def _merge_adjacent_bands(bands):
    """
    Merges any bands that touch or overlap once expressed as integer slot
    ranges -- two bands found narrowly apart in the continuous `phi` scan
    can end up adjacent or overlapping after `slot_index_bounds_for_phi_range`'s
    own defensive 1-slot buffer, and a caller (`galaxyPlan.py`) should
    never see two rows that describe the same contiguous range.

    Args:
        bands (list): `ShellBand` instances, any order.

    Returns:
        list: `ShellBand` instances, sorted by `slot_index_min`, with any
              touching/overlapping pairs combined into one.
    """
    if not bands:
        return []
    bands = sorted(bands, key=lambda b: b.slot_index_min)
    merged = [bands[0]]
    for band in bands[1:]:
        last = merged[-1]
        if band.slot_index_min <= last.slot_index_max + 1:
            merged[-1] = ShellBand(last.slot_index_min, max(last.slot_index_max, band.slot_index_max))
        else:
            merged.append(band)
    return merged


def find_shell_bands(shape, edge_pc, shell_index, threshold_rho,
                      scan_samples=PHI_SCAN_SAMPLES, tolerance=PHI_BISECT_TOLERANCE):
    """
    This shell's candidate slot-index band(s) -- see the module docstring
    for what "candidate" means and why it's normally exactly one band.

    Coarse-scans `bound_relative_density_at_phi` across `phi` in `[0, pi]`
    at `scan_samples` evenly-spaced points, finds every contiguous run
    that clears `threshold_rho`, refines each run's true edge via
    bisection (`tolerance`), then converts each refined `phi` range to an
    exact slot-index range via `galaxyGeometry.slot_index_bounds_for_phi_range`
    (the same exact closed-form inverse `enumerate_sectors_within_radius`
    already relies on elsewhere in this codebase).

    Args:
        shape (galaxyDensity.GalaxyShape): The galaxy's shape parameters.
        edge_pc (float): The sector edge length, parsecs.
        shell_index (int): The shell index `k`.
        threshold_rho (float): The `relative_density` a sector must reach
            to qualify.
        scan_samples (int): See `PHI_SCAN_SAMPLES`.
        tolerance (float): See `PHI_BISECT_TOLERANCE`.

    Returns:
        list: `ShellBand` instances, in ascending `slot_index_min` order --
              empty if this shell has no qualifying content at all.
    """
    n_k = shell_sector_count(shell_index)
    r_k = shell_radius_pc(shell_index, edge_pc)

    if whole_shell_qualifies(shape, r_k, threshold_rho):
        return [ShellBand(0, n_k - 1)]

    def f(phi):
        return bound_relative_density_at_phi(shape, r_k, phi) - threshold_rho

    phis = [math.pi * j / scan_samples for j in range(scan_samples + 1)]
    samples = [f(phi) for phi in phis]

    bands = []
    i = 0
    n = len(samples)
    while i < n:
        if samples[i] < 0:
            i += 1
            continue
        start = i
        while i < n and samples[i] >= 0:
            i += 1
        end = i - 1

        phi_lo = phis[start] if start == 0 else \
            _bisect_crossing_to_qualifying_side(f, phis[start - 1], phis[start], tolerance)
        phi_hi = phis[end] if end == n - 1 else \
            _bisect_crossing_to_qualifying_side(f, phis[end + 1], phis[end], tolerance)

        i_min, i_max = slot_index_bounds_for_phi_range(phi_lo, phi_hi, n_k)
        bands.append(ShellBand(i_min, i_max))

    return _merge_adjacent_bands(bands)
