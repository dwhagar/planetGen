# stellarObjects/galaxyDensity.py

"""
Galaxy Disk/Spiral Density Model
====================================

Implements the exponential-disk-plus-bulge-plus-spiral-arm density model
from `docs/design/galaxy-disk-density.md` (revision 2) -- the piece of
that design pass that was, until now, only ever prototyped in throwaway
scratch scripts during design discussion, never shipped as real,
tested code.

The model gives a `relative_density` at any galaxy-frame position,
normalized so it equals `1.0` at a chosen calibration point (by
convention, the position on the disk plane, at the inter-arm minimum
azimuth, at `calibration_radius_pc` from the galactic center) -- matching
the same "1.0 = realistic local density" convention `sectorGen.py
--density` already uses elsewhere in this codebase. A sector's own
predicted star count is then `SpaceSector.expected_system_count() *
relative_density(sector_center)`; a sector with a predicted count `< 1` is
one the design doc's own gating rule treats as not worth generating.

Deliberately **not** hardcoded to Milky-Way scale: `GalaxyShape` bundles
every shape parameter as plain fields, so the same formulas describe a
galaxy of any size -- a small toy/demonstration spiral just as validly as
a Milky-Way-scale one. `GalaxyShape.build` computes the one derived value
(`k_norm`, the normalization constant) from the others; nothing here reads
`physical_constants.GALACTIC_CENTER_DISTANCE_LY` or any other Milky-Way-
specific constant, since that would only make sense for a galaxy actually
built to that real scale.

This module is pure and side-effect-free, like `galaxyGeometry.py` --
callers own deciding what to do with a `relative_density` value (whether
that's an occupancy decision, a `--density` multiplier, or just an
analysis figure), and own persistence, if any.
"""

import math
from collections import namedtuple

GalaxyShape = namedtuple(
    "GalaxyShape",
    [
        "disk_scale_length_pc", "disk_scale_height_pc",
        "bulge_scale_radius_pc", "bulge_amplitude",
        "arm_count", "pitch_angle_rad", "arm_amplitude",
        "spiral_reference_radius_pc", "spiral_reference_angle_rad",
        "k_norm",
    ],
)
"""A galaxy's disk/bulge/spiral shape, at whatever physical scale the
caller chooses -- see `GalaxyShape.build` to construct one with `k_norm`
computed automatically rather than guessed."""


def _sech_squared(x):
    """
    `1 / cosh(x)**2`, computed so it never overflows for large `|x|` --
    plain `math.cosh(x) ** 2` raises `OverflowError` once `|x|` exceeds
    ~710 (`cosh` itself grows like `exp(|x|) / 2`), which a small
    `disk_scale_height_pc` relative to a galaxy's own radius reaches
    easily (found directly: a toy-scale shape with
    `disk_scale_height_pc=12` overflows by shell ~2400).

    Rewritten around `exp(-2*|x|)` (always in `(0, 1]`, so it can only
    underflow to a harmless `0.0`, never overflow) rather than `exp(x)`
    directly: `1/cosh(x)**2 == 4*e / (1 + e)**2` where `e = exp(-2*|x|)`
    -- algebraically identical to the direct formula (factor `exp(|x|)`
    out of `cosh(x) = (e^x + e^-x)/2` and simplify), and correctly
    approaches its true limit of `0.0` as `|x| -> infinity` instead of
    crashing partway there.

    Args:
        x (float): Any real number.

    Returns:
        float: `1 / cosh(x)**2`, in `(0, 1]`.
    """
    ax = abs(x)
    if ax > 700.0:  # exp(-1400) underflows to exactly 0.0 anyway
        return 0.0
    e = math.exp(-2.0 * ax)
    return 4.0 * e / (1.0 + e) ** 2


def _raw_density(position_pc, shape):
    """
    `relative_density` before the `k_norm` normalization is applied --
    see `docs/design/galaxy-disk-density.md` revision 2, section 1.

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        shape (GalaxyShape): The galaxy's shape parameters (`k_norm` is
                             ignored here -- this is the pre-normalization
                             value).

    Returns:
        float: `rho_bulge(r_3d) + rho_disk_radial(r_cyl) * f_z(z) *
              arm_factor(r_cyl, theta)`, unnormalized (`>= 0`).
    """
    x, y, z = position_pc
    r_cyl = math.hypot(x, y)
    r_3d = math.sqrt(x * x + y * y + z * z)

    bulge = shape.bulge_amplitude * math.exp(-r_3d / shape.bulge_scale_radius_pc)
    disk_radial = math.exp(-r_cyl / shape.disk_scale_length_pc)
    f_z = _sech_squared(z / shape.disk_scale_height_pc)

    if r_cyl > 1e-9:
        theta = math.atan2(y, x)
        theta_arm = (
            shape.spiral_reference_angle_rad
            + math.log(r_cyl / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad)
        )
        arm_factor = 1 + shape.arm_amplitude * math.cos(shape.arm_count * (theta - theta_arm))
    else:
        arm_factor = 1.0

    return bulge + disk_radial * f_z * arm_factor


def relative_density(position_pc, shape):
    """
    This galaxy's stellar density at `position_pc`, relative to
    `shape`'s own calibration point (`1.0` there, by construction -- see
    `GalaxyShape.build`).

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        shape (GalaxyShape): The galaxy's shape parameters.

    Returns:
        float: `>= 0`. `1.0` at the calibration point; `> 1.0` in denser
              regions (bulge, spiral arm crests near the core); `< 1.0`
              in sparser ones (outer disk, inter-arm, off-plane).
    """
    return shape.k_norm * _raw_density(position_pc, shape)


def build_galaxy_shape(
    disk_scale_length_pc,
    disk_scale_height_pc,
    bulge_scale_radius_pc,
    bulge_amplitude,
    arm_count,
    pitch_angle_rad,
    arm_amplitude,
    calibration_radius_pc=None,
    spiral_reference_angle_rad=0.0,
):
    """
    Builds a `GalaxyShape` with `k_norm` computed automatically, rather
    than left for the caller to guess -- calibrated so `relative_density`
    equals exactly `1.0` at `(calibration_radius_pc, z=0)`, at the
    inter-arm minimum azimuth at that radius (the design doc's own
    calibration choice: real spiral galaxies' analog of "the Sun" sits
    between arms, not on one, so this is the conservative reading of
    "1.0 = typical local density").

    Args:
        disk_scale_length_pc, disk_scale_height_pc, bulge_scale_radius_pc,
        bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude:
            The galaxy's shape parameters -- see `GalaxyShape`.
        calibration_radius_pc (float, optional): The in-plane radius the
            `relative_density = 1.0` calibration point sits at. Defaults
            to `2.82 * disk_scale_length_pc` -- the real Milky Way's own
            solar-radius-to-scale-length ratio, reused here as a
            physically-motivated default at whatever scale this galaxy
            actually is, not a Milky-Way-specific constant.
        spiral_reference_angle_rad (float): See `GalaxyShape` -- arbitrary
            fixed orientation, `0.0` unless there's a reason to prefer
            another.

    Returns:
        GalaxyShape: With `k_norm` filled in.
    """
    if calibration_radius_pc is None:
        calibration_radius_pc = 2.82 * disk_scale_length_pc

    unnormalized = GalaxyShape(
        disk_scale_length_pc=disk_scale_length_pc,
        disk_scale_height_pc=disk_scale_height_pc,
        bulge_scale_radius_pc=bulge_scale_radius_pc,
        bulge_amplitude=bulge_amplitude,
        arm_count=arm_count,
        pitch_angle_rad=pitch_angle_rad,
        arm_amplitude=arm_amplitude,
        spiral_reference_radius_pc=disk_scale_length_pc,
        spiral_reference_angle_rad=spiral_reference_angle_rad,
        k_norm=1.0,
    )

    theta_arm_at_calibration = (
        spiral_reference_angle_rad
        + math.log(calibration_radius_pc / disk_scale_length_pc) / math.tan(pitch_angle_rad)
    )
    theta_interarm = theta_arm_at_calibration + math.pi / arm_count
    calibration_position = (
        calibration_radius_pc * math.cos(theta_interarm),
        calibration_radius_pc * math.sin(theta_interarm),
        0.0,
    )
    k_norm = 1.0 / _raw_density(calibration_position, unnormalized)

    return unnormalized._replace(k_norm=k_norm)


def predicted_star_count(position_pc, shape, expected_system_count_at_density_1):
    """
    This sector's predicted star (system) count -- the quantity the design
    doc's own gating rule compares against `1.0` to decide whether a
    sector is worth generating at all.

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs -- the
                             sector's own center.
        shape (GalaxyShape): The galaxy's shape parameters.
        expected_system_count_at_density_1 (float): A sector's own
            `SpaceSector.expected_system_count()` -- the expected system
            count at `relative_density = 1.0` (real local stellar
            density), before this position's own relative richness is
            applied.

    Returns:
        float: `expected_system_count_at_density_1 *
              relative_density(position_pc, shape)`.
    """
    return expected_system_count_at_density_1 * relative_density(position_pc, shape)
