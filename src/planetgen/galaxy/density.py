# planetgen/galaxy/density.py

"""
Galaxy Disk/Spiral Density Model
====================================

Implements the galaxy density model from
`docs/design/galaxy-disk-density.md`: an exponential thin disk carrying
the spiral arms, an exponential thick disk, and a boxy bar bulge, each
fitted to published Milky Way structure (GEN.118, GEN.119; sources in
`tuning`'s "galaxy density model" block).

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

import functools
import math
from collections import namedtuple

from planetgen import tuning

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
    # min(): at |x| below ~1e-8, rounding can put the quotient one ulp
    # above its true maximum of 1.
    return min(1.0, 4.0 * e / (1.0 + e) ** 2)


def _vertical(z, scale_height_pc):
    """A disk's vertical profile: `sech^2(z / (2 * h))`, the isothermal
    sheet, whose tail away from the plane is `exp(-|z| / h)` -- so `h` is
    the exponential scale height surveys quote (BHG16's 300 pc thin disk,
    900 pc thick disk)."""
    return _sech_squared(z / (2.0 * scale_height_pc))


@functools.lru_cache(maxsize=64)
def model_terms(shape):
    """
    The parts of the model that aren't `GalaxyShape` fields but follow
    from them (`tuning`'s fixed Milky Way ratios): the thick disk's scale
    length, height and amplitude, and the bar's axes and angle. The Galaxy
    Map's prisms read the same dict (`galaxymap3d._density_shape`).

    Returns:
        dict: `thick_disk_amplitude` (raw density at the center, in the
            plane, against the thin disk's 1), `thick_disk_scale_length_pc`,
            `thick_disk_scale_height_pc`, `bulge_scale_y_pc`,
            `bulge_scale_z_pc`, `bar_angle_rad` (galaxy-frame azimuth of
            the bar's near end) with its `bar_cos` and `bar_sin`, and
            `solar_angle_rad` (the Sun's azimuth: the inter-arm minimum
            at the solar radius).
    """
    length = shape.disk_scale_length_pc
    solar_radius = tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH * length
    thick_length = tuning.THICK_DISK_SCALE_LENGTH_RATIO * length
    solar_angle = _interarm_angle(solar_radius, shape)
    bar_angle = solar_angle + math.radians(tuning.BULGE_BAR_ANGLE_DEG)
    return {
        "thick_disk_amplitude": tuning.THICK_DISK_LOCAL_DENSITY_RATIO
        * math.exp(solar_radius / thick_length - solar_radius / length),
        "thick_disk_scale_length_pc": thick_length,
        "thick_disk_scale_height_pc": tuning.THICK_DISK_SCALE_HEIGHT_RATIO * shape.disk_scale_height_pc,
        "bulge_scale_y_pc": tuning.BULGE_AXIS_RATIO_Y * shape.bulge_scale_radius_pc,
        "bulge_scale_z_pc": tuning.BULGE_AXIS_RATIO_Z * shape.bulge_scale_radius_pc,
        "bar_angle_rad": bar_angle,
        "bar_cos": math.cos(bar_angle),
        "bar_sin": math.sin(bar_angle),
        "solar_angle_rad": solar_angle,
    }


def shape_with_terms(shape):
    """`shape`'s fields and its `model_terms` in one plain dict: what the
    Galaxy Map's prisms read (`queryDb.galaxy_density_shape`)."""
    return {**shape._asdict(), **model_terms(shape)}


def _interarm_angle(radius_pc, shape):
    """The azimuth of the inter-arm minimum at `radius_pc`."""
    theta_arm = (
        shape.spiral_reference_angle_rad
        + math.log(radius_pc / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad)
    )
    return theta_arm + math.pi / shape.arm_count


def _bulge(x, y, z, shape, terms):
    """The boxy bar bulge (Dwek et al. 1995 G2; see `tuning.BULGE_*`):
    `bulge_amplitude * exp(-r_s^2 / 2)`, `r_s^4 = ((x'/x0)^2 +
    (y'/y0)^2)^2 + (z/z0)^4` in the bar's own frame."""
    cos_a, sin_a = terms["bar_cos"], terms["bar_sin"]
    along = (x * cos_a + y * sin_a) / shape.bulge_scale_radius_pc
    across = (y * cos_a - x * sin_a) / terms["bulge_scale_y_pc"]
    up = z / terms["bulge_scale_z_pc"]
    in_plane = along * along + across * across
    return shape.bulge_amplitude * math.exp(-0.5 * math.sqrt(in_plane * in_plane + up ** 4))


def bulge_bound(r_cyl, z, shape):
    """The bulge's maximum over azimuth at `(r_cyl, z)` (unnormalized):
    along the bar's long axis. Falls as `r_cyl` and `|z|` grow."""
    terms = model_terms(shape)
    along = r_cyl / shape.bulge_scale_radius_pc
    up = z / terms["bulge_scale_z_pc"]
    return shape.bulge_amplitude * math.exp(-0.5 * math.sqrt(along ** 4 + up ** 4))


def _thick_disk(r_cyl, z, shape, terms):
    """The thick disk (BHG16; see `tuning.THICK_DISK_*`), unnormalized."""
    return (terms["thick_disk_amplitude"] * math.exp(-r_cyl / terms["thick_disk_scale_length_pc"])
            * _vertical(z, terms["thick_disk_scale_height_pc"]))


def _components(position_pc, shape):
    """`(bulge, thin disk without arms, thick disk, arm cosine)` at a
    point, unnormalized."""
    x, y, z = position_pc
    terms = model_terms(shape)
    r_cyl = math.hypot(x, y)
    thin = math.exp(-r_cyl / shape.disk_scale_length_pc) * _vertical(z, shape.disk_scale_height_pc)
    return _bulge(x, y, z, shape, terms), thin, _thick_disk(r_cyl, z, shape, terms), _arm_cosine(x, y, shape)


def _raw_density(position_pc, shape):
    """
    `relative_density` before the `k_norm` normalization is applied --
    see `docs/design/galaxy-disk-density.md`, section 1.

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        shape (GalaxyShape): The galaxy's shape parameters (`k_norm` is
                             ignored here -- this is the pre-normalization
                             value).

    Returns:
        float: `bulge + thin_disk * arm_factor + thick_disk`,
              unnormalized (`>= 0`).
    """
    bulge, thin, thick, arm_cos = _components(position_pc, shape)
    return bulge + thin * (1 + shape.arm_amplitude * arm_cos) + thick


def _arm_cosine(x, y, shape):
    """`cos(arm_count * (theta - theta_arm))` at `(x, y)`: 1 on an arm's
    crest, -1 midway between arms; 0 at the center, where arms are
    undefined."""
    r_cyl = math.hypot(x, y)
    if r_cyl <= 1e-9:
        return 0.0
    theta = math.atan2(y, x)
    theta_arm = (
        shape.spiral_reference_angle_rad
        + math.log(r_cyl / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad)
    )
    return math.cos(shape.arm_count * (theta - theta_arm))


def population_densities(position_pc, shape):
    """
    `relative_density` at `position_pc` split by stellar population
    (young, intermediate, old disk stars and bulge stars; see
    `tuning.STELLAR_POPULATION_*`), so a sector can draw each
    system's age from the mix where it sits.

    The components always sum to `relative_density`, so no sector's total
    changes. The bulge term is the "bulge" population, and the thick disk
    and the halo floor are "old". The thin disk is shared among the three
    disk populations in proportion to how many stars each would put here:
    its share of disk star formation (the length of its age range), times
    its own vertical profile (`sech^2(z/2h)/h`, a thinner disk for younger
    stars) and its own arm contrast (`1 + A cos(arm phase)`, strongest for
    young stars). So young stars crowd the arms near the plane, and far
    off the plane or in the bulge nearly every star is old.

    Returns:
        dict: `{"young", "intermediate", "old", "bulge"}` -> density, each
              `>= 0`, summing to `relative_density(position_pc, shape)`.
    """
    z = position_pc[2]
    bulge, thin, thick, arm_cos = _components(position_pc, shape)
    disk = shape.k_norm * thin * (1 + shape.arm_amplitude * arm_cos)
    model = shape.k_norm * (bulge + thick) + disk
    # The halo floor's share (`tuning.MIN_RELATIVE_DENSITY`, GEN.78): old
    # stars, wherever the model itself falls below it.
    halo = max(model, tuning.MIN_RELATIVE_DENSITY) - model

    ages = tuning.STELLAR_POPULATION_AGE_RANGES_GY
    formation_span = tuning.STAR_FORMATION_AGE_RANGE_GY[1] - tuning.STAR_FORMATION_AGE_RANGE_GY[0]
    weights = {}
    for name, ratio in tuning.STELLAR_POPULATION_SCALE_HEIGHT_RATIO.items():
        height = shape.disk_scale_height_pc * ratio
        share = (ages[name][1] - ages[name][0]) / formation_span
        weights[name] = (share * _vertical(z, height) / height
                         * (1 + tuning.STELLAR_POPULATION_ARM_AMPLITUDE[name] * arm_cos))
    total = sum(weights.values())
    if total <= 0.0:
        # Far enough off the plane that every profile underflows: the
        # thickest (old) disk is all that's left.
        weights, total = {name: 0.0 for name in weights}, 1.0
        weights["old"] = 1.0
    densities = {name: max(disk, 0.0) * weight / total for name, weight in weights.items()}
    densities["old"] += max(shape.k_norm * thick, 0.0) + max(halo, 0.0)
    densities["bulge"] = shape.k_norm * bulge
    return densities


def component_masses(shape, steps=40):
    """
    Each component's share of the model's stars, by integrating its
    density over space (unnormalized units; only the ratios mean
    anything) -- what the math check and the tests compare with the Milky
    Way's measured masses. The arm factor averages to 1 around any
    circle, so the thin disk's mass doesn't depend on it.

    Args:
        shape (GalaxyShape): The galaxy's shape parameters.
        steps (int): Grid steps per axis (more is slower and closer).

    Returns:
        dict: `{"bulge", "thin_disk", "thick_disk"}` -> integrated density.
    """
    terms = model_terms(shape)
    # The bulge on a box in its own frame (the integral doesn't care
    # about the bar's angle), out to 4 scale lengths on each axis.
    axes = (shape.bulge_scale_radius_pc, terms["bulge_scale_y_pc"], terms["bulge_scale_z_pc"])
    cells = [[(i + 0.5) * 4.0 * axis / steps for i in range(steps)] for axis in axes]
    bulge = sum(
        math.exp(-0.5 * math.sqrt((u * u + v * v) ** 2 + w ** 4))
        for u in (x / axes[0] for x in cells[0])
        for v in (y / axes[1] for y in cells[1])
        for w in (z / axes[2] for z in cells[2])
    )
    bulge *= 8.0 * shape.bulge_amplitude * (4.0 / steps) ** 3 * axes[0] * axes[1] * axes[2]

    def disk(amplitude, length, height):
        r_max, z_max = 12.0 * length, 12.0 * height
        dr, dz = r_max / (4 * steps), z_max / (4 * steps)
        radial = sum(2.0 * math.pi * r * math.exp(-r / length) for r in ((i + 0.5) * dr for i in range(4 * steps)))
        vertical = sum(_vertical(z, height) for z in ((i + 0.5) * dz for i in range(4 * steps)))
        return amplitude * radial * dr * 2.0 * vertical * dz

    return {
        "bulge": bulge,
        "thin_disk": disk(1.0, shape.disk_scale_length_pc, shape.disk_scale_height_pc),
        "thick_disk": disk(terms["thick_disk_amplitude"], terms["thick_disk_scale_length_pc"],
                           terms["thick_disk_scale_height_pc"]),
    }


def relative_density(position_pc, shape):
    """
    This galaxy's stellar density at `position_pc`, relative to
    `shape`'s own calibration point (`1.0` there, by construction -- see
    `GalaxyShape.build`).

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        shape (GalaxyShape): The galaxy's shape parameters.

    Returns:
        float: `>= tuning.MIN_RELATIVE_DENSITY` (the halo floor, GEN.78).
              `1.0` at the calibration point; `> 1.0` in denser regions
              (bulge, spiral arm crests near the core); `< 1.0` in sparser
              ones (outer disk, inter-arm, off-plane).
    """
    return max(shape.k_norm * _raw_density(position_pc, shape), tuning.MIN_RELATIVE_DENSITY)


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
            to `tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH *
            disk_scale_length_pc` -- the real Milky Way's own
            solar-radius-to-scale-length ratio (8.2 / 2.6 kpc), reused
            here at whatever scale this galaxy actually is.
        spiral_reference_angle_rad (float): See `GalaxyShape` -- arbitrary
            fixed orientation, `0.0` unless there's a reason to prefer
            another.

    Returns:
        GalaxyShape: With `k_norm` filled in.
    """
    if calibration_radius_pc is None:
        calibration_radius_pc = tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH * disk_scale_length_pc

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

    theta_interarm = _interarm_angle(calibration_radius_pc, unnormalized)
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
