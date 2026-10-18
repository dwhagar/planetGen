# planetgen/galaxy/density.py

"""
Galaxy Disk/Spiral Density Model
====================================

Implements the galaxy density model from
`docs/design/galaxy-disk-density.md`: an exponential thin disk carrying
the spiral arms, an exponential thick disk, a boxy bar bulge and a stellar
halo, each fitted to published Milky Way structure (GEN.118, GEN.119, GEN.234; sources in
`tuning`'s "galaxy density model" block).

The model gives a `relative_density` at any galaxy-frame position,
normalized so it equals `1.0` at a chosen calibration point (by
convention, the position on the disk plane, at the inter-arm minimum
azimuth, at `calibration_radius_pc` from the galactic center) -- matching
the same "1.0 = realistic local density" convention `generation/run_sector.py
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

This module is pure and side-effect-free, like `galaxy/spherical.py` --
callers own deciding what to do with a `relative_density` value (whether
that's an occupancy decision, a `--density` multiplier, or just an
analysis figure), and own persistence, if any.
"""

import functools
import math
from collections import namedtuple

import numpy as np

from planetgen import tuning
from planetgen.galaxy import component_field as componentField
from planetgen.galaxy import star_formation as starFormation
from planetgen.galaxy import structure as galaxyStructure
from planetgen.util import detmath

GalaxyShape = namedtuple(
    "GalaxyShape",
    [
        "disk_scale_length_pc", "disk_scale_height_pc",
        "bulge_scale_radius_pc", "bulge_amplitude",
        "arm_count", "pitch_angle_rad", "arm_amplitude", "arm_level", "nuclear_amplitude",
        "spiral_reference_radius_pc", "spiral_reference_angle_rad",
        "k_norm", "components", "masses", "structure",
    ],
    defaults=(None, None, galaxyStructure.REFERENCE_STRUCTURE),
)
"""A galaxy's disk/bulge/spiral shape, at whatever physical scale the
caller chooses -- see `GalaxyShape.build` to construct one with `k_norm`
computed automatically rather than guessed. `components` is the galaxy's
`galaxy.component_field.ComponentField` (GEN.160, GEN.237: its globular clusters and satellite galaxies), `None` for a
galaxy without; it is not part of the stored shape row. `masses` is `(stellar mass, halo mass)` in solar masses, the
galaxy's own (`galaxy.mass_chain`; GEN.469), which its disk's star formation history is worked out from; `None` takes the
median galaxy of the disk's size (`galaxy.disk_history.population_shares`). `structure` is the galaxy's thick disk, bar,
stellar halo and population motions (`galaxy.structure`, GEN.532, GEN.533). Neither is part of the stored row: both are
drawn again from the galaxy's seed; a shape built with no seed has the Milky Way's published structure."""


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


def _vertical_any(z, scale_height_pc):
    """`_vertical` for a scalar or an array of heights (an array broadcasts, one height per band of a shell)."""
    if np.ndim(z) == 0:
        return _vertical(float(z), scale_height_pc)
    e = np.exp(-2.0 * np.minimum(np.abs(np.asarray(z, dtype=float)) / (2.0 * scale_height_pc), 350.0))
    return np.minimum(1.0, 4.0 * e / (1.0 + e) ** 2)


def thick_disk_amplitude(shape):
    """The thick disk's density at the centre, in the plane, against the thin disk's 1 (raw, before the arm level):
    its mass is `structure.thick_mass_ratio` of the thin disk's (`arm_level` times `2 pi L^2 x 4 h`), so with its own
    length and height (the structure's ratios to the thin disk's) the amplitude is the mass ratio over the product
    of the squared length ratio and the height ratio, times the arm level."""
    structure = shape.structure
    return structure.thick_mass_ratio * shape.arm_level / (structure.thick_length_ratio ** 2 * structure.thick_height_ratio)


def population_height_ratios(shape):
    """The scale height of each thin-disc population over the old one's, for this galaxy (`galaxy.structure`): a dict
    of "young", "intermediate" and "old" (1)."""
    return _height_ratios(shape.structure)


@functools.lru_cache(maxsize=64)
def _height_ratios(structure):
    return galaxyStructure.population_height_ratios(structure)


def bulge_to_total(shape):
    """The bulge's mass over the bulge's and the thin and thick disks', in the model's units (the share
    `shape_draw.bulge_amplitude_for` solves the bulge's amplitude from)."""
    terms = model_terms(shape)
    thin = shape.arm_level * 2.0 * math.pi * shape.disk_scale_length_pc ** 2 * 4.0 * shape.disk_scale_height_pc
    thick = (terms["thick_disk_amplitude"] * 2.0 * math.pi * terms["thick_disk_scale_length_pc"] ** 2
             * 4.0 * terms["thick_disk_scale_height_pc"])
    return terms["bulge_mass"] / (terms["bulge_mass"] + thin + thick)


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
            the bar's near end) with its `bar_cos` and `bar_sin`,
            `solar_angle_rad` (the Sun's azimuth: the inter-arm minimum
            at the solar radius), and the nuclear terms (GEN.280): `bulge_mass` (the bulge's
            integral, in the model's unnormalized units), `nuclear_cluster_mass`, `nuclear_cluster_scale_pc`
            (the Hernquist scale radius), `nuclear_cluster_softening_pc` (the radius inside which the
            cusp is held flat, `tuning.NUCLEAR_RESOLUTION_PC`), `nuclear_disc_mass`, `nuclear_disc_length_pc`,
            `nuclear_disc_height_pc` and `nuclear_disc_density` (the disc's centre).
    """
    length = shape.disk_scale_length_pc
    structure = shape.structure
    solar_radius = tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH * length
    thick_length = structure.thick_length_ratio * length
    thick_height = structure.thick_height_ratio * shape.disk_scale_height_pc
    solar_angle = _interarm_angle(solar_radius, shape)
    bar_angle = solar_angle + math.radians(structure.bar_angle_deg)
    scale_y = structure.bar_axis_ratio_y * shape.bulge_scale_radius_pc
    scale_z = structure.bar_axis_ratio_z * shape.bulge_scale_radius_pc
    bulge_mass = shape.bulge_amplitude * shape.bulge_scale_radius_pc * scale_y * scale_z * BULGE_UNIT_INTEGRAL
    cluster_mass = max(shape.nuclear_amplitude, 0.0) * tuning.NUCLEAR_CLUSTER_BULGE_RATIO * bulge_mass
    cluster_scale = tuning.NUCLEAR_CLUSTER_SCALE_RATIO * shape.bulge_scale_radius_pc
    disc_mass = cluster_mass * tuning.NUCLEAR_DISC_CLUSTER_MASS_RATIO
    disc_length = tuning.NUCLEAR_DISC_LENGTH_RATIO * shape.bulge_scale_radius_pc
    disc_height = tuning.NUCLEAR_DISC_HEIGHT_RATIO * shape.bulge_scale_radius_pc
    return {
        "thick_disk_amplitude": thick_disk_amplitude(shape),
        "thick_disk_scale_length_pc": thick_length,
        "thick_disk_scale_height_pc": thick_height,
        "bulge_scale_y_pc": scale_y,
        "bulge_scale_z_pc": scale_z,
        "bar_angle_rad": bar_angle,
        "bar_cos": math.cos(bar_angle),
        "bar_sin": math.sin(bar_angle),
        "solar_angle_rad": solar_angle,
        "bulge_mass": bulge_mass,
        "nuclear_cluster_mass": cluster_mass,
        "nuclear_cluster_scale_pc": cluster_scale,
        "nuclear_cluster_softening_pc": _cluster_softening(cluster_mass, cluster_scale),
        "nuclear_disc_mass": disc_mass,
        "nuclear_disc_length_pc": disc_length,
        "nuclear_disc_height_pc": disc_height,
        "nuclear_disc_density": disc_mass / (8.0 * math.pi * disc_length ** 2 * disc_height),
    }


BULGE_UNIT_INTEGRAL = 20.648191867705272
"""float: The integral of the bulge's unit profile `exp(-r_s^2 / 2)` over space, in the units of its three scale
lengths (a test compares it with the numerical integral of `component_masses`): the bulge's mass is its amplitude
times its three scale lengths times this."""


def _cluster_density_at(radius_pc, mass, scale_pc):
    """A Hernquist (Dehnen slope 1) cluster of `mass` and scale radius `scale_pc` at `radius_pc`: `M a / (2 pi r
    (r + a)^3)`."""
    return mass * scale_pc / (2.0 * math.pi * radius_pc * (radius_pc + scale_pc) ** 3)


def _cluster_softening(mass, scale_pc):
    """
    The radius below which the nuclear cluster's cusp is held flat (GEN.280): where the Hernquist density equals the
    mean density of the grid's core ball (radius `tuning.NUCLEAR_RESOLUTION_PC`, `3 M / (4 pi R (R + a)^2)`), found
    by bisection (the density falls with the radius). The ball's cell, which a point sample puts at the centre, then
    holds exactly the mass inside it. `0.0` for no cluster.
    """
    if mass <= 0.0:
        return 0.0
    ball = tuning.NUCLEAR_RESOLUTION_PC
    mean = 3.0 * mass / (4.0 * math.pi * ball * (ball + scale_pc) ** 2)
    low, high = 1e-12 * ball, ball
    for _ in range(80):
        middle = math.sqrt(low * high)
        if _cluster_density_at(middle, mass, scale_pc) > mean:
            low = middle
        else:
            high = middle
    return high


def _halo_unit_integral(inner_slope, outer_slope, core_pc, break_pc, axis_ratio):
    """
    The volume integral of the halo's unit-amplitude profile over all space, in pc^3: the flat core, the inner
    power law out to the break and the outer one beyond it, each a closed form in the spheroidal radius `m`
    (a shell at `m` holds `4 pi q m^2 dm`).
    """
    def power(n, low, high):
        # integral of (m / break) ** -n * m ** 2 dm from `low` to `high` (`high` may be infinity)
        k = 3.0 - n
        if abs(k) < 1e-12:
            return break_pc ** n * math.log(high / low)
        top = 0.0 if math.isinf(high) else high ** k
        return break_pc ** n * (top - low ** k) / k

    core = (core_pc / break_pc) ** -inner_slope * core_pc ** 3 / 3.0
    return 4.0 * math.pi * axis_ratio * (core + power(inner_slope, core_pc, break_pc)
                                         + power(outer_slope, break_pc, math.inf))


@functools.lru_cache(maxsize=64)
def halo_terms(shape):
    """
    The stellar halo's parameters for `shape` (GEN.234): a spheroid with a broken power law, `amplitude *
    (max(m, core) / break) ** -slope`, `m = sqrt(R**2 + (z / axis_ratio)**2)`, the slope `halo_inner_slope`
    inside `halo_break_pc` and `halo_outer_slope` beyond. The radii are `tuning`'s ratios to the disk scale
    length. The amplitude is solved, like `k_norm`, so the halo's integral is the galaxy's stellar halo share
    (`galaxy.structure.stellar_halo_system_fraction`, GEN.532) of the visible components' (`component_masses`), in the
    model's unnormalized density.

    Returns:
        dict: `halo_amplitude`, `halo_core_pc`, `halo_break_pc`, `halo_axis_ratio`, `halo_inner_slope`,
            `halo_outer_slope`.
    """
    length = shape.disk_scale_length_pc
    core, break_pc = tuning.HALO_CORE_RATIO * length, tuning.HALO_BREAK_RATIO * length
    unit = _halo_unit_integral(tuning.HALO_INNER_SLOPE, tuning.HALO_OUTER_SLOPE, core, break_pc,
                               tuning.HALO_AXIS_RATIO)
    visible = sum(component_masses(shape, steps=HALO_MASS_STEPS).values())
    return {
        "halo_amplitude": galaxyStructure.stellar_halo_system_fraction(shape.structure) * visible / unit,
        "halo_core_pc": core, "halo_break_pc": break_pc, "halo_axis_ratio": tuning.HALO_AXIS_RATIO,
        "halo_inner_slope": tuning.HALO_INNER_SLOPE, "halo_outer_slope": tuning.HALO_OUTER_SLOPE,
    }


HALO_MASS_STEPS = 60
"""int: `component_masses` steps the halo's amplitude is solved with (finer than the default 40; cached per shape)."""


def _halo_scalar(r_cyl, z, terms):
    """`_halo` at one point, in plain floats (the one-point forms run per cell)."""
    ratio = max(detmath.hypot(r_cyl, z / terms["halo_axis_ratio"]), terms["halo_core_pc"]) / terms["halo_break_pc"]
    return terms["halo_amplitude"] * ratio ** -(terms["halo_inner_slope"] if ratio < 1.0 else terms["halo_outer_slope"])


def _halo(r_cyl, z, terms):
    """The halo at cylindrical radius `r_cyl` and height `z` (arrays, or one of them), unnormalized."""
    m = np.sqrt(np.square(r_cyl) + np.square(z / terms["halo_axis_ratio"]))
    ratio = np.maximum(m, terms["halo_core_pc"]) / terms["halo_break_pc"]
    slope = np.where(ratio < 1.0, terms["halo_inner_slope"], terms["halo_outer_slope"])
    return terms["halo_amplitude"] * np.power(ratio, -slope)


def halo_density(position_pc, shape):
    """The halo's `relative_density` at `position_pc` (`k_norm` applied): what the halo adds to a cell."""
    x, y, z = position_pc
    return shape.k_norm * _halo_scalar(detmath.hypot(x, y), z, halo_terms(shape))


def halo_local_density_scale(shape):
    """
    How much the halo would grow if it were normalised to the local density (`tuning.HALO_LOCAL_DENSITY_FRACTION`
    of the thin disk's density at the Sun) instead of to its mass: the alternative the design note keeps, about 2.7.
    """
    solar = tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH * shape.disk_scale_length_pc
    return tuning.HALO_LOCAL_DENSITY_FRACTION / halo_density((solar, 0.0, 0.0), shape)


def shape_with_terms(shape):
    """`shape`'s fields and its `model_terms` and `halo_terms` in one plain dict: what the
    Galaxy Map's prisms read (`db.query.galaxy_density_shape`)."""
    fields = shape._asdict()
    del fields["components"], fields["masses"], fields["structure"]
    return {**fields, **model_terms(shape), **halo_terms(shape)}


def _interarm_angle(radius_pc, shape):
    """The azimuth of the inter-arm minimum at `radius_pc`."""
    theta_arm = (
        shape.spiral_reference_angle_rad
        + math.log(radius_pc / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad)
    )
    return theta_arm + math.pi / shape.arm_count


def _nuclear_cluster(radius_pc, terms):
    """The nuclear star cluster (GEN.280) at `radius_pc` from the centre, unnormalized: a Hernquist profile held
    flat inside `nuclear_cluster_softening_pc`. `0.0` for a galaxy without one."""
    if not terms["nuclear_cluster_mass"]:
        return 0.0
    return _cluster_density_at(max(radius_pc, terms["nuclear_cluster_softening_pc"]), terms["nuclear_cluster_mass"],
                               terms["nuclear_cluster_scale_pc"])


def _nuclear_cluster_array(radius_pc, terms):
    """`_nuclear_cluster` on an array of radii."""
    radius_pc = np.asarray(radius_pc, dtype=float)
    if not terms["nuclear_cluster_mass"]:
        return np.zeros_like(radius_pc)
    r = np.maximum(radius_pc, terms["nuclear_cluster_softening_pc"])
    scale = terms["nuclear_cluster_scale_pc"]
    return terms["nuclear_cluster_mass"] * scale / (2.0 * math.pi * r * (r + scale) ** 3)


def _nuclear_disc(r_cyl, z, terms):
    """The nuclear stellar disc (GEN.280) at cylindrical radius `r_cyl` and height `z`, unnormalized: an exponential
    disc with a `sech^2` height, like the thin disc's form. `0.0` for a galaxy without one."""
    if not terms["nuclear_disc_mass"]:
        return 0.0
    return (terms["nuclear_disc_density"] * math.exp(-r_cyl / terms["nuclear_disc_length_pc"])
            * _vertical(z, terms["nuclear_disc_height_pc"]))


def _nuclear_disc_array(r_cyl, z, terms):
    """`_nuclear_disc` on arrays (`z` one height or an array that broadcasts against `r_cyl`)."""
    if not terms["nuclear_disc_mass"]:
        return np.zeros_like(np.asarray(r_cyl, dtype=float) + np.asarray(z, dtype=float))
    return (terms["nuclear_disc_density"] * np.exp(-np.asarray(r_cyl, dtype=float) / terms["nuclear_disc_length_pc"])
            * _vertical_any(z, terms["nuclear_disc_height_pc"]))


def nuclear_bound(r_cyl, z, shape):
    """The nuclear terms' sum at `(r_cyl, z)` (unnormalized). Both depend on the radius or the height alone, so
    this is their value at any azimuth and the bound of the sector's outline."""
    terms = model_terms(shape)
    return _nuclear_cluster(detmath.hypot(r_cyl, z), terms) + _nuclear_disc(r_cyl, z, terms)


def _bulge(x, y, z, shape, terms):
    """The bulge, a bar when the galaxy has one (`galaxy.structure`; the form of Dwek et al. 1995 G2):
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
    along the bar's long axis, and the nuclear terms (`nuclear_bound`) on top of it. Falls as `r_cyl` and `|z|`
    grow."""
    terms = model_terms(shape)
    along = r_cyl / shape.bulge_scale_radius_pc
    up = z / terms["bulge_scale_z_pc"]
    return (shape.bulge_amplitude * math.exp(-0.5 * math.sqrt(along ** 4 + up ** 4))
            + nuclear_bound(r_cyl, z, shape))


def _thick_disk(r_cyl, z, shape, terms):
    """The thick disk (`galaxy.structure`), unnormalized."""
    return (terms["thick_disk_amplitude"] * math.exp(-r_cyl / terms["thick_disk_scale_length_pc"])
            * _vertical(z, terms["thick_disk_scale_height_pc"]))


def _components(position_pc, shape):
    """`(bulge, thin disk without arms, thick disk, arm cosine)` at a
    point, unnormalized."""
    x, y, z = position_pc
    terms = model_terms(shape)
    r_cyl = detmath.hypot(x, y)
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
        float: `bulge + thin_disk * arm_factor + thick_disk + halo + the nuclear cluster and disc`,
              unnormalized (`>= 0`).
    """
    bulge, thin, thick, arm_cos = _components(position_pc, shape)
    halo = _halo_scalar(detmath.hypot(position_pc[0], position_pc[1]), position_pc[2], halo_terms(shape))
    nuclear = nuclear_bound(detmath.hypot(position_pc[0], position_pc[1]), position_pc[2], shape)
    return bulge + thin * shape.arm_level * (1 + shape.arm_amplitude * arm_cos) + thick + halo + nuclear


def _arm_cosine(x, y, shape):
    """`cos(arm_count * (theta - theta_arm))` at `(x, y)`: 1 on an arm's
    crest, -1 midway between arms; 0 at the center, where arms are
    undefined."""
    r_cyl = detmath.hypot(x, y)
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
    (young, intermediate, old disk stars, bulge stars and halo stars; see
    `tuning.STELLAR_POPULATION_*`), so a sector can draw each
    system's age from the mix where it sits.

    The components always sum to `relative_density`, so no sector's total
    changes. The bulge term is the "bulge" population, the halo term the
    "halo" population (GEN.234) and the thick disk is "old". The nuclear star cluster and the nuclear stellar disc
    (GEN.280) are the "nuclear_cluster" and "nuclear_disc" populations, on top of the bulge. "cluster" and "satellite"
    are the globular clusters' and the satellite galaxies' systems in the cell's cube over a cell's at density 1
    (`shape.components`, GEN.160 and GEN.237; 0 for a galaxy without). The thin disk is shared among the three
    disk populations in proportion to how many stars each would put here:
    its share of disk star formation (the part of the disk's gas history in its age range, `galaxy.disk_history`), times
    its own vertical profile (`sech^2(z/2h)/h`, a thinner disk for younger
    stars) and its own arm contrast (`1 + A cos(arm phase)`, strongest for
    young stars). So young stars crowd the arms near the plane, and far
    off the plane or in the bulge nearly every star is old. The young and intermediate weights also carry the
    recent star-formation factor of the radius (`galaxy.star_formation`, GEN.134): fewer in the outer disk and in
    the bar's gap, more at its edge; the three still share the same disk density.

    Returns:
        dict: `{"young", "intermediate", "old", "bulge", "halo", "cluster", "satellite", "nuclear_cluster",
              "nuclear_disc"}` -> density, each
              `>= 0`, summing to `relative_density(position_pc, shape)`.
    """
    z = position_pc[2]
    bulge, thin, thick, arm_cos = _components(position_pc, shape)
    disk = shape.k_norm * thin * shape.arm_level * (1 + shape.arm_amplitude * arm_cos)
    halo = shape.k_norm * _halo_scalar(detmath.hypot(position_pc[0], position_pc[1]), z, halo_terms(shape))

    recent = starFormation.recent_factor(detmath.hypot(position_pc[0], position_pc[1]), shape)
    weights = {}
    for name, ratio in population_height_ratios(shape).items():
        height = shape.disk_scale_height_pc * ratio
        share = starFormation.disk_population_share(name, shape)
        weights[name] = (share * _vertical(z, height) / height
                         * (1 + tuning.STELLAR_POPULATION_ARM_AMPLITUDE[name] * arm_cos)
                         * (recent if name in starFormation.RECENT_POPULATIONS else 1.0))
    total = sum(weights.values())
    if total <= 0.0:
        # Far enough off the plane that every profile underflows: the
        # thickest (old) disk is all that's left.
        weights, total = {name: 0.0 for name in weights}, 1.0
        weights["old"] = 1.0
    densities = {name: max(disk, 0.0) * weight / total for name, weight in weights.items()}
    densities["old"] += max(shape.k_norm * thick, 0.0)
    densities["bulge"] = shape.k_norm * bulge
    densities["halo"] = max(halo, 0.0)
    densities.update(_component_densities(position_pc, shape))
    terms = model_terms(shape)
    r_cyl = detmath.hypot(position_pc[0], position_pc[1])
    densities["nuclear_cluster"] = shape.k_norm * _nuclear_cluster(detmath.hypot(r_cyl, z), terms)
    densities["nuclear_disc"] = shape.k_norm * _nuclear_disc(r_cyl, z, terms)
    return densities


def population_densities_array(x, y, z, shape):
    """
    `population_densities` at arrays of points `x`, `y` and a height `z` (one,
    or an array that broadcasts against them: a band's height per row, GEN.278)
    (PERF.89): the same terms in the same order, so a point's densities match
    the one-point form (a test compares them).

    Returns:
        dict: `{"young", "intermediate", "old", "bulge", "halo", "cluster", "satellite", "nuclear_cluster",
            "nuclear_disc"}` -> array of densities.
    """
    terms = model_terms(shape)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r_cyl = np.hypot(x, y)
    thin = np.exp(-r_cyl / shape.disk_scale_length_pc) * _vertical_any(z, shape.disk_scale_height_pc)
    thick = (terms["thick_disk_amplitude"] * np.exp(-r_cyl / terms["thick_disk_scale_length_pc"])
             * _vertical_any(z, terms["thick_disk_scale_height_pc"]))
    along = (x * terms["bar_cos"] + y * terms["bar_sin"]) / shape.bulge_scale_radius_pc
    across = (y * terms["bar_cos"] - x * terms["bar_sin"]) / terms["bulge_scale_y_pc"]
    up = z / terms["bulge_scale_z_pc"]
    in_plane = along * along + across * across
    bulge = shape.bulge_amplitude * np.exp(-0.5 * np.sqrt(in_plane * in_plane + up ** 4))
    with np.errstate(divide="ignore", invalid="ignore"):
        theta_arm = (shape.spiral_reference_angle_rad
                     + np.log(r_cyl / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad))
        arm_cos = np.where(r_cyl <= 1e-9, 0.0, np.cos(shape.arm_count * (np.arctan2(y, x) - theta_arm)))
    disk = shape.k_norm * thin * shape.arm_level * (1 + shape.arm_amplitude * arm_cos)
    halo = shape.k_norm * _halo(r_cyl, z, halo_terms(shape))

    recent = starFormation.recent_factor(r_cyl, shape)
    weights = {}
    for name, ratio in population_height_ratios(shape).items():
        height = shape.disk_scale_height_pc * ratio
        share = starFormation.disk_population_share(name, shape)
        weights[name] = (share * _vertical_any(z, height) / height
                         * (1 + tuning.STELLAR_POPULATION_ARM_AMPLITUDE[name] * arm_cos)
                         * (recent if name in starFormation.RECENT_POPULATIONS else 1.0))
    total = sum(weights.values())
    # Every profile underflowed: the thickest (old) disk is all that's left.
    empty = total <= 0.0
    if np.any(empty):
        weights = {name: np.where(empty, 1.0 if name == "old" else 0.0, weight)
                   for name, weight in weights.items()}
        total = np.where(empty, 1.0, total)
    densities = {name: np.maximum(disk, 0.0) * weight / total for name, weight in weights.items()}
    densities["old"] = densities["old"] + np.maximum(shape.k_norm * thick, 0.0)
    densities["bulge"] = shape.k_norm * bulge
    densities["halo"] = np.maximum(halo, 0.0)
    # The components are not in the array form: its points are band bins far wider than a cluster, and the scatter's
    # per-band weights leave components to their own pass (GEN.161, GEN.240); the one-point form is the cell's.
    densities["cluster"] = np.zeros_like(densities["halo"])
    densities["satellite"] = np.zeros_like(densities["halo"])
    densities["satellite_young"] = np.zeros_like(densities["halo"])
    densities["satellite_intermediate"] = np.zeros_like(densities["halo"])
    densities["nuclear_cluster"] = shape.k_norm * _nuclear_cluster_array(np.sqrt(r_cyl ** 2 + np.square(z)), terms)
    densities["nuclear_disc"] = shape.k_norm * _nuclear_disc_array(r_cyl, z, terms)
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
        dict: `{"bulge", "thin_disk", "thick_disk"}` -> integrated density (the halo's is `halo_mass`,
            kept apart because the halo's amplitude is solved from these three).
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
        "thin_disk": disk(shape.arm_level, shape.disk_scale_length_pc, shape.disk_scale_height_pc),
        "thick_disk": disk(terms["thick_disk_amplitude"], terms["thick_disk_scale_length_pc"],
                           terms["thick_disk_scale_height_pc"]),
    }


def halo_mass(shape):
    """The halo's share of the model's stars (the units of `component_masses`), in closed form: the halo
    term integrated over all space, the stellar halo share of the visible components' by construction."""
    terms = halo_terms(shape)
    return terms["halo_amplitude"] * _halo_unit_integral(
        terms["halo_inner_slope"], terms["halo_outer_slope"], terms["halo_core_pc"], terms["halo_break_pc"],
        terms["halo_axis_ratio"])


def _component_densities(position_pc, shape):
    """The components' part of the density at `position_pc` by kind (GEN.160, GEN.237): the systems the clusters and the
    satellites put in the cell's cube over the systems a cell holds at density 1, `0.0` each for a galaxy without."""
    components = shape.components
    if components is None:
        return dict.fromkeys(componentField.POPULATIONS_IN_THE_GATE, 0.0)
    return components.densities(position_pc)


def _component_density(position_pc, shape):
    """All the components' part of the density at `position_pc`: `_component_densities` summed."""
    components = shape.components
    return 0.0 if components is None else components.relative_density(position_pc)


def relative_density(position_pc, shape):
    """
    This galaxy's stellar density at `position_pc`, relative to
    `shape`'s own calibration point (`1.0` there, by construction -- see
    `GalaxyShape.build`).

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        shape (GalaxyShape): The galaxy's shape parameters.

    Returns:
        float: `1.0` at the calibration point; `> 1.0` in denser regions
              (bulge, spiral arm crests near the core); `< 1.0` in sparser
              ones (outer disk, inter-arm, off-plane).
    """
    return shape.k_norm * _raw_density(position_pc, shape) + _component_density(position_pc, shape)


def total_density_array(x, y, z, shape):
    """
    The sum of `population_densities` (each term counted as none when negative,
    as `bright_stars._densities` does) at arrays of points `x`, `y` and a height `z` (one, or an array that broadcasts against them): the weight a band's angle bin carries in the phenomena scatter,
    for a whole band at once (PERF.63). A shape with no negative term gives
    `relative_density` itself.
    """
    terms = model_terms(shape)
    r_cyl = np.hypot(x, y)
    thin = np.exp(-r_cyl / shape.disk_scale_length_pc) * _vertical_any(z, shape.disk_scale_height_pc)
    thick = (terms["thick_disk_amplitude"] * np.exp(-r_cyl / terms["thick_disk_scale_length_pc"])
             * _vertical_any(z, terms["thick_disk_scale_height_pc"]))
    along = (x * terms["bar_cos"] + y * terms["bar_sin"]) / shape.bulge_scale_radius_pc
    across = (y * terms["bar_cos"] - x * terms["bar_sin"]) / terms["bulge_scale_y_pc"]
    up = z / terms["bulge_scale_z_pc"]
    in_plane = along * along + across * across
    bulge = shape.bulge_amplitude * np.exp(-0.5 * np.sqrt(in_plane * in_plane + up ** 4))
    with np.errstate(divide="ignore", invalid="ignore"):
        theta_arm = (shape.spiral_reference_angle_rad
                     + np.log(r_cyl / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad))
        arm_cos = np.where(r_cyl <= 1e-9, 0.0, np.cos(shape.arm_count * (np.arctan2(y, x) - theta_arm)))
    disk = shape.k_norm * thin * shape.arm_level * (1 + shape.arm_amplitude * arm_cos)
    halo = shape.k_norm * _halo(r_cyl, z, halo_terms(shape))
    nuclear = shape.k_norm * (_nuclear_cluster_array(np.sqrt(r_cyl ** 2 + np.square(z)), terms)
                              + _nuclear_disc_array(r_cyl, z, terms))
    return (np.maximum(disk, 0.0) + np.maximum(shape.k_norm * thick, 0.0)
            + np.maximum(halo, 0.0) + np.maximum(shape.k_norm * bulge, 0.0) + np.maximum(nuclear, 0.0))


def arm_terms(arm_density, interarm_density):
    """
    `(arm_amplitude, arm_level)` for a thin disk that is `arm_density` on an
    arm crest and `interarm_density` midway between arms (ADM.49), in the
    disk's own units (1.4 and 0.6 are the usual `1 +/- 0.4`): the thin disk's
    arm factor is `arm_level * (1 + arm_amplitude * cos)`, so the crest is
    `arm_level * (1 + arm_amplitude)` and the trough `arm_level * (1 - arm_amplitude)`.

    Raises:
        ValueError: Unless `arm_density >= interarm_density >= 0` and the crest is above 0.
    """
    if not (arm_density > 0 and 0 <= interarm_density <= arm_density):
        raise ValueError("the arm density must be above 0 and at least the inter-arm density, which is 0 or more")
    # Rounded so the usual 1.4 and 0.6 give exactly the 0.4 and 1 they stand for.
    return (round((arm_density - interarm_density) / (arm_density + interarm_density), 12),
            round((arm_density + interarm_density) / 2.0, 12))


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
    arm_level=1.0,
    nuclear_amplitude=0.0,
    masses=None,
    structure=galaxyStructure.REFERENCE_STRUCTURE,
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
        arm_level (float): The thin disk's density level, the middle of the
            arm crest and the inter-arm trough (ADM.49; see `arm_terms`).
        nuclear_amplitude (float): How much of the observed nuclear star cluster and nuclear stellar disc
            the centre holds (GEN.280): 1 is the observed pair, scaled to the bulge, 0 is none.
        masses (tuple, optional): The galaxy's `(stellar mass, halo mass)`, solar masses (`GalaxyShape.masses`).
        structure (galaxy.structure.Structure): The galaxy's thick disk, bar, stellar halo and population motions
            (GEN.532, GEN.533); the Milky Way's published one when none is given.

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
        arm_level=arm_level,
        nuclear_amplitude=nuclear_amplitude,
        spiral_reference_radius_pc=disk_scale_length_pc,
        spiral_reference_angle_rad=spiral_reference_angle_rad,
        k_norm=1.0,
        masses=None if masses is None else (float(masses[0]), float(masses[1])),
        structure=structure,
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
