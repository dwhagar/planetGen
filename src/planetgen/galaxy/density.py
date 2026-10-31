# planetgen/galaxy/density.py

"""
Galaxy Disk/Spiral Density Model
====================================

Implements the galaxy density model from
`docs/design/galaxy/galaxy-disk-density.md`: an exponential thin disk carrying
the spiral arms, an exponential thick disk, a boxy bar bulge and a stellar
halo, each fitted to published Milky Way structure (GEN.118, GEN.119, GEN.234; sources in
`tuning`'s "galaxy density model" block).

The model gives a `relative_density` at any galaxy-frame position, in units of the real local stellar
density (`physical_constants.LOCAL_STELLAR_DENSITY_LY3`, the "1.0 = realistic local density" convention
`generation/run_sector.py --density` uses elsewhere in this codebase), normalized so the whole galaxy holds the
systems its stellar mass makes (`structure.stellar_mass_solar` over `mean_system_mass_solar`, GEN.404, GEN.545): no
position of the Sun is read. A sector's own
predicted star count is then `SpaceSector.expected_system_count() *
relative_density(sector_center)`; a sector with a predicted count `< 1` is
one the design doc's own gating rule treats as not worth generating.

Deliberately **not** hardcoded to Milky-Way scale: `GalaxyShape` bundles
every shape parameter as plain fields, so the same formulas describe a
galaxy of any size -- a small toy/demonstration spiral just as validly as
a Milky-Way-scale one. `GalaxyShape.build` computes the one derived value
(`k_norm`, the normalization constant) from the others and the galaxy's stellar mass; nothing here reads a
Sun or home position.

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
from planetgen.galaxy import core_gas as coreGas
from planetgen.galaxy import star_formation as starFormation
from planetgen.galaxy import system_mass as systemMass
from planetgen.galaxy import structure as galaxyStructure
from planetgen.physics import constants as physical_constants
from planetgen.physics.units import ly_to_pc
from planetgen.util import detmath

GalaxyShape = namedtuple(
    "GalaxyShape",
    [
        "disk_scale_length_pc", "disk_scale_height_pc",
        "bulge_scale_radius_pc", "bulge_amplitude",
        "arm_count", "pitch_angle_rad", "arm_amplitude", "arm_level", "nuclear_amplitude",
        "spiral_reference_radius_pc", "spiral_reference_angle_rad",
        "k_norm", "components", "masses", "structure", "system_mass_scale",
        "family", "bar_present", "bar_axis_ratio_b", "bar_axis_ratio_c", "bar_angle_rad", "thick_disk_amplitude",
        "spheroid_amplitude", "spheroid_sersic_n", "spheroid_effective_radius_pc", "spheroid_axis_ratio_b",
        "spheroid_axis_ratio_c", "spheroid_orientation_rad", "body_clump_share",
    ],
    defaults=(None, None, galaxyStructure.REFERENCE_STRUCTURE, 1.0,
              "spiral", None, None, None, None, None, 0.0, 0.0, 0.0, 1.0, 1.0, (0.0, 0.0, 0.0), 0.0),
)
"""A galaxy's disk/bulge/spiral shape, at whatever physical scale the
caller chooses -- see `GalaxyShape.build` to construct one with `k_norm`
computed automatically rather than guessed. `components` is the galaxy's
`galaxy.component_field.ComponentField` (GEN.160, GEN.237: its globular clusters and satellite galaxies), `None` for a
galaxy without; it is not part of the stored shape row. `masses` is `(stellar mass, halo mass)` in solar masses, the
galaxy's own (`galaxy.mass_chain`; GEN.469), which its disk's star formation history is worked out from; `None` takes the
median galaxy of the disk's size (`galaxy.disk_history.population_shares`). `structure` is the galaxy's thick disk, bar,
stellar halo and population motions (`galaxy.structure`, GEN.532, GEN.533). Neither is part of the stored row: both are
drawn again from the galaxy's seed; a shape built with no seed has the Milky Way's published structure.
`system_mass_scale` is how many times the mean system mass the generator derives (GEN.545) is divided by in this galaxy:
1 for the derived richness, more when the Census dial (GEN.607) made the galaxy hold more, lighter systems; it is not
part of the stored row (the plan's pins hold the richness, `with_system_richness` sets it again).

The last thirteen fields are the stored family and component numbers (`galaxy_shape`, DB.66, GEN.379). `family` is the
galaxy's family (`spiral` for every galaxy built so far). The bar's presence, its two axis ratios and angle, and the
thick disk's amplitude are `None` for "as the structure draw gives it from the seed" (the value is then in `structure`),
or a number that pins it (`with_family_columns` puts it in `structure` as well). The Sersic spheroid's fields are the
spheroid family's (GEN.381): an amplitude of 0 has no spheroid, and the model has no spheroid component yet. `body_clump_share`
is the irregular family's young clump share (GEN.385)."""

FAMILY_COLUMNS = ("family", "bar_present", "bar_axis_ratio_b", "bar_axis_ratio_c", "bar_angle_rad", "thick_disk_amplitude",
                  "spheroid_amplitude", "spheroid_sersic_n", "spheroid_effective_radius_pc", "spheroid_axis_ratio_b",
                  "spheroid_axis_ratio_c", "spheroid_orientation_1_rad", "spheroid_orientation_2_rad",
                  "spheroid_orientation_3_rad", "body_clump_share")
"""tuple: The `galaxy_shape` columns DB.66 added, in the order of the `GalaxyShape` fields they fill (the orientation's three
angles are one tuple, `spheroid_orientation_rad`)."""


def with_family_columns(shape, row):
    """
    `shape` with the stored family and component numbers of a `galaxy_shape` row (`FAMILY_COLUMNS`) read in (GEN.379): a
    NULL bar or thick disk number stays `None`, meaning the structure draw's own value, never a fixed number; a number
    is put in `shape.structure` too, so the model and everything that reads the structure agree. The thick disk's
    amplitude sets the structure's mass ratio (`thick_disk_amplitude` is the other way round).
    """
    fields = {name: row[name] for name in FAMILY_COLUMNS if not name.startswith("spheroid_orientation_")}
    fields["bar_present"] = None if fields["bar_present"] is None else bool(fields["bar_present"])
    fields["spheroid_orientation_rad"] = tuple(row[f"spheroid_orientation_{n}_rad"] for n in (1, 2, 3))
    pinned = shape._replace(**fields)
    structure = pinned.structure
    changes = {}
    if pinned.bar_present is not None:
        changes["bar_present"] = pinned.bar_present
    if pinned.bar_axis_ratio_b is not None:
        changes["bar_axis_ratio_y"] = pinned.bar_axis_ratio_b
    if pinned.bar_axis_ratio_c is not None:
        changes["bar_axis_ratio_z"] = pinned.bar_axis_ratio_c
    if pinned.bar_angle_rad is not None:
        changes["bar_angle_deg"] = math.degrees(pinned.bar_angle_rad)
    if pinned.thick_disk_amplitude is not None and pinned.arm_level:
        changes["thick_mass_ratio"] = (pinned.thick_disk_amplitude * structure.thick_length_ratio ** 2
                                       * structure.thick_height_ratio / pinned.arm_level)
    return pinned._replace(structure=structure._replace(**changes)) if changes else pinned


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
    """A disk's vertical profile: `exp(-|z| / h)` (GEN.544), the form surveys fit to star counts, so `h` is the
    exponential scale height they quote (BHG16's 300 pc thin disk, 900 pc thick disk) and the midplane density is the
    one the fits give for the same column (`2 h` times it). It has a cusp in the plane; the stars drawn from it do not
    see one."""
    return math.exp(-min(abs(z) / scale_height_pc, 700.0))


def _vertical_any(z, scale_height_pc):
    """`_vertical` for a scalar or an array of heights (an array broadcasts, one height per band of a shell)."""
    if np.ndim(z) == 0:
        return _vertical(float(z), scale_height_pc)
    return np.exp(-np.minimum(np.abs(np.asarray(z, dtype=float)) / scale_height_pc, 700.0))


def _sech_vertical(z, scale_height_pc):
    """The nuclear stellar disc's vertical profile: `sech^2(z / (2 * h))`, the isothermal sheet, whose tail away from
    the plane is `exp(-|z| / h)` and whose column is `4 h` times its midplane density."""
    return _sech_squared(z / (2.0 * scale_height_pc))


def _sech_vertical_any(z, scale_height_pc):
    """`_sech_vertical` for a scalar or an array of heights."""
    if np.ndim(z) == 0:
        return _sech_vertical(float(z), scale_height_pc)
    e = np.exp(-2.0 * np.minimum(np.abs(np.asarray(z, dtype=float)) / (2.0 * scale_height_pc), 350.0))
    return np.minimum(1.0, 4.0 * e / (1.0 + e) ** 2)


def thick_disk_amplitude(shape):
    """The thick disk's density at the centre, in the plane, against the thin disk's 1 (raw, before the arm level):
    its mass is `structure.thick_mass_ratio` of the thin disk's (`arm_level` times `2 pi L^2 x 4 h`), so with its own
    length and height (the structure's ratios to the thin disk's) the amplitude is the mass ratio over the product
    of the squared length ratio and the height ratio, times the arm level."""
    if shape.thick_disk_amplitude is not None:
        return shape.thick_disk_amplitude
    structure = shape.structure
    return structure.thick_mass_ratio * shape.arm_level / (structure.thick_length_ratio ** 2 * structure.thick_height_ratio)


def population_height_ratios(shape):
    """The scale height of each thin-disc population over the old one's, for this galaxy (`galaxy.structure`): a dict
    of "young", "intermediate" and "old" (1)."""
    return _height_ratios(shape.structure)


@functools.lru_cache(maxsize=64)
def _height_ratios(structure):
    return galaxyStructure.population_height_ratios(structure)


def central_black_hole_mass_solar(shape):
    """The mass, solar masses, of the black hole at the centre of the galaxy `shape` describes (GEN.402): the mass in
    its bulge, or in its stars without a classical bulge, with the galaxy's own draw about the relation."""
    return galaxyStructure.black_hole_mass_solar(shape.structure, bulge_to_total(shape))


def _disk_mass(amplitude, scale_length_pc, scale_height_pc):
    """The mass, in the model's units, of an exponential disk of central density `amplitude`."""
    return amplitude * 2.0 * math.pi * scale_length_pc ** 2 * 2.0 * scale_height_pc


def bulge_to_total(shape):
    """The bulge's mass over the bulge's and the thin and thick disks', in the model's units (the share
    `shape_draw.bulge_amplitude_for` solves the bulge's amplitude from)."""
    terms = model_terms(shape)
    thin = shape.arm_level * 2.0 * math.pi * shape.disk_scale_length_pc ** 2 * 2.0 * shape.disk_scale_height_pc
    thick = (terms["thick_disk_amplitude"] * 2.0 * math.pi * terms["thick_disk_scale_length_pc"] ** 2
             * 2.0 * terms["thick_disk_scale_height_pc"])
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
            `reference_angle_rad` (the reference azimuth the bar's angle and the satellites' frame are
            measured from: the inter-arm minimum at the reference radius), and the nuclear terms (GEN.280): `bulge_mass` (the bulge's
            integral, in the model's unnormalized units), `nuclear_cluster_mass`, `nuclear_cluster_scale_pc`
            (the Hernquist scale radius), `nuclear_cluster_softening_pc` (the radius inside which the
            cusp is held flat, `tuning.NUCLEAR_RESOLUTION_PC`), `nuclear_disc_mass`, `nuclear_disc_length_pc`,
            `nuclear_disc_height_pc` and `nuclear_disc_density` (the disc's centre).
    """
    if shape.spheroid_amplitude:
        raise ValueError("this galaxy stores a Sersic spheroid, which the density model has no component for yet")
    length = shape.disk_scale_length_pc
    structure = shape.structure
    reference_radius = tuning.REFERENCE_RADIUS_SCALE_LENGTHS * length
    thick_length = structure.thick_length_ratio * length
    thick_height = structure.thick_height_ratio * shape.disk_scale_height_pc
    reference_angle = _interarm_angle(reference_radius, shape)
    bar_angle = reference_angle + math.radians(structure.bar_angle_deg)
    scale_y = structure.bar_axis_ratio_y * shape.bulge_scale_radius_pc
    scale_z = structure.bar_axis_ratio_z * shape.bulge_scale_radius_pc
    bulge_mass = shape.bulge_amplitude * shape.bulge_scale_radius_pc * scale_y * scale_z * BULGE_UNIT_INTEGRAL
    thick_amplitude = thick_disk_amplitude(shape)
    main_mass = (bulge_mass + _disk_mass(shape.arm_level, length, shape.disk_scale_height_pc)
                 + _disk_mass(thick_amplitude, thick_length, thick_height))
    # The model's units per solar mass of the visible stars (the galaxy's stars less the stellar halo's): the nuclear
    # cluster (GEN.402) and the nuclear disc (GEN.628) are drawn in solar masses and enter the model through it.
    per_solar = main_mass / (structure.stellar_mass_solar * (1.0 - structure.halo_stellar_ratio))
    cluster_mass = max(shape.nuclear_amplitude, 0.0) * structure.nuclear_cluster_mass_solar * per_solar
    cluster_scale = tuning.NUCLEAR_CLUSTER_SCALE_RATIO * shape.bulge_scale_radius_pc
    disc_mass = max(shape.nuclear_amplitude, 0.0) * structure.nuclear_disc_mass_solar * per_solar
    disc_length = tuning.NUCLEAR_DISC_LENGTH_RATIO * shape.bulge_scale_radius_pc
    disc_height = tuning.NUCLEAR_DISC_HEIGHT_RATIO * shape.bulge_scale_radius_pc
    return {
        "thick_disk_amplitude": thick_amplitude,
        "thick_disk_scale_length_pc": thick_length,
        "thick_disk_scale_height_pc": thick_height,
        "bulge_scale_y_pc": scale_y,
        "bulge_scale_z_pc": scale_z,
        "bar_angle_rad": bar_angle,
        "bar_cos": math.cos(bar_angle),
        "bar_sin": math.sin(bar_angle),
        "reference_angle_rad": reference_angle,
        "bulge_mass": bulge_mass,
        "visible_mass": main_mass,
        "mass_per_solar": per_solar,
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


def _cluster_softened_fraction(scale_pc, softening_pc):
    """The share of the nuclear cluster's nominal (Hernquist) mass that its profile holds once it is flat inside
    `softening_pc` (`_cluster_softening`): the Hernquist mass beyond the softening radius plus the flat core, in closed
    form. `1.0` for no softening."""
    outside = 1.0 - (softening_pc / (softening_pc + scale_pc)) ** 2
    core = 2.0 / 3.0 * scale_pc * softening_pc ** 2 / (softening_pc + scale_pc) ** 3
    return outside + core


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
    (`galaxy.structure.stellar_halo_system_fraction`, GEN.532) of the visible components' (`model_terms`
    `visible_mass`), in the model's unnormalized density.

    Returns:
        dict: `halo_amplitude`, `halo_core_pc`, `halo_break_pc`, `halo_axis_ratio`, `halo_inner_slope`,
            `halo_outer_slope`.
    """
    length = shape.disk_scale_length_pc
    core, break_pc = tuning.HALO_CORE_RATIO * length, tuning.HALO_BREAK_RATIO * length
    unit = _halo_unit_integral(tuning.HALO_INNER_SLOPE, tuning.HALO_OUTER_SLOPE, core, break_pc,
                               tuning.HALO_AXIS_RATIO)
    visible = model_terms(shape)["visible_mass"]
    return {
        "halo_amplitude": galaxyStructure.stellar_halo_system_fraction(shape.structure) * visible / unit,
        "halo_core_pc": core, "halo_break_pc": break_pc, "halo_axis_ratio": tuning.HALO_AXIS_RATIO,
        "halo_inner_slope": tuning.HALO_INNER_SLOPE, "halo_outer_slope": tuning.HALO_OUTER_SLOPE,
    }


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


def shape_with_terms(shape):
    """`shape`'s fields and its `model_terms` and `halo_terms` in one plain dict: what the
    Galaxy Map's prisms read (`db.query.galaxy_density_shape`)."""
    fields = shape._asdict()
    for name in ("components", "masses", "structure", "system_mass_scale", "spheroid_orientation_rad"):
        del fields[name]
    for name in FAMILY_COLUMNS:
        fields.pop(name, None)
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
            * _sech_vertical(z, terms["nuclear_disc_height_pc"]))


def _nuclear_disc_array(r_cyl, z, terms):
    """`_nuclear_disc` on arrays (`z` one height or an array that broadcasts against `r_cyl`)."""
    if not terms["nuclear_disc_mass"]:
        return np.zeros_like(np.asarray(r_cyl, dtype=float) + np.asarray(z, dtype=float))
    return (terms["nuclear_disc_density"] * np.exp(-np.asarray(r_cyl, dtype=float) / terms["nuclear_disc_length_pc"])
            * _sech_vertical_any(z, terms["nuclear_disc_height_pc"]))


def nuclear_bound(r_cyl, z, shape):
    """The nuclear terms' sum at `(r_cyl, z)` (unnormalized). Both depend on the radius or the height alone, so
    this is their value at any azimuth and the bound of the sector's outline."""
    terms = model_terms(shape)
    return _nuclear_cluster(detmath.hypot(r_cyl, z), terms) + _nuclear_disc(r_cyl, z, terms)


def _bulge(x, y, z, shape, terms):
    """The bulge, a bar when the galaxy has one (`galaxy.structure`; the form of Dwek et al. 1995 G2):
    `bulge_amplitude * exp(-r_s^2 / 2)`, `r_s^4 = ((x'/x0)^2 +
    (y'/y0)^2)^2 + (z/z0)^4` in the bar's own frame."""
    if not shape.bulge_amplitude:
        return 0.0
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
    """The thick disk (`galaxy.structure`), unnormalized; exactly `0.0` when its amplitude is."""
    if not terms["thick_disk_amplitude"]:
        return 0.0
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


COMPONENT_NAMES = ("thin_disk", "thick_disk", "bulge", "halo", "nuclear_cluster", "nuclear_disc")
"""tuple: The density model's components, in the order `density_components` lists them (GEN.379)."""

ComponentValues = namedtuple("ComponentValues", (*COMPONENT_NAMES, "thin_profile", "arm_cosine"))
"""One point's value of each component (`COMPONENT_NAMES`, unnormalized: the thin disk with its arms), the thin disk's
profile alone (`thin_profile`, before the arm level and the arms) and the arm cosine the arms were drawn with."""


def component_values(position_pc, shape):
    """
    The value of each component of the density model at a point, unnormalized (`ComponentValues`): what `relative_density`
    and `population_densities` sum and split. A component whose amplitude is 0 is exactly `0.0`.
    """
    terms = model_terms(shape)
    x, y, z = position_pc
    r_cyl = detmath.hypot(x, y)
    profile = math.exp(-r_cyl / shape.disk_scale_length_pc) * _vertical(z, shape.disk_scale_height_pc)
    arm_cos = _arm_cosine(x, y, shape)
    return ComponentValues(
        profile * shape.arm_level * (1 + shape.arm_amplitude * arm_cos),
        _thick_disk(r_cyl, z, shape, terms), _bulge(x, y, z, shape, terms), _halo_scalar(r_cyl, z, halo_terms(shape)),
        _nuclear_cluster(detmath.hypot(r_cyl, z), terms), _nuclear_disc(r_cyl, z, terms), profile, arm_cos)


def _raw_density(position_pc, shape):
    """
    `relative_density` before the `k_norm` normalization is applied --
    see `docs/design/galaxy/galaxy-disk-density.md`, section 1.

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        shape (GalaxyShape): The galaxy's shape parameters (`k_norm` is
                             ignored here -- this is the pre-normalization
                             value).

    Returns:
        float: The sum of `component_values`: `bulge + thin_disk * arm_factor + thick_disk + halo + the nuclear cluster and
              disc`, unnormalized (`>= 0`).
    """
    values = component_values(position_pc, shape)
    return (values.bulge + values.thin_disk + values.thick_disk + values.halo
            + (values.nuclear_cluster + values.nuclear_disc))


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
    its own vertical profile (`exp(-|z|/h)/h`, a thinner disk for younger
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
    values = component_values(position_pc, shape)
    arm_cos = values.arm_cosine
    disk = shape.k_norm * values.thin_profile * shape.arm_level * (1 + shape.arm_amplitude * arm_cos)
    halo = shape.k_norm * values.halo

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
    densities["old"] += max(shape.k_norm * values.thick_disk, 0.0)
    densities["bulge"] = shape.k_norm * values.bulge
    densities["halo"] = max(halo, 0.0)
    densities.update(_component_densities(position_pc, shape))
    r_cyl = detmath.hypot(position_pc[0], position_pc[1])
    densities["nuclear_cluster"] = shape.k_norm * values.nuclear_cluster
    densities["nuclear_disc"] = shape.k_norm * values.nuclear_disc
    core = coreGas.core_gas(shape)
    if core.young_density:      # the ring's young field stars are the nuclear disc's own (GEN.282)
        moved = min(shape.k_norm * coreGas.young_density(r_cyl, z, core), densities["nuclear_disc"])
        densities["young"] += moved
        densities["nuclear_disc"] -= moved
    return densities


ArrayValues = namedtuple("ArrayValues", ("r_cyl", "thin_profile", "thick_disk", "bulge", "arm_cosine", "halo",
                                         "nuclear_cluster", "nuclear_disc"))
"""`component_values` on arrays (`component_arrays`): each component's value at each point, unnormalized, the thin disk
before its arm level and arms (`thin_profile`) and the arm cosine beside it; `r_cyl` is the cylindrical radius."""


def component_arrays(x, y, z, shape):
    """
    `component_values` at arrays of points `x`, `y` and a height `z` (one, or an array that broadcasts against them: a band's
    height per row, GEN.278), in the same terms in the same order (PERF.89). The one place the array forms of the density
    read the components from.
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
    return ArrayValues(
        r_cyl=r_cyl, thin_profile=thin, thick_disk=thick, bulge=bulge, arm_cosine=arm_cos,
        halo=_halo(r_cyl, z, halo_terms(shape)),
        nuclear_cluster=_nuclear_cluster_array(np.sqrt(r_cyl ** 2 + np.square(z)), terms),
        nuclear_disc=_nuclear_disc_array(r_cyl, z, terms))


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
    values = component_arrays(x, y, z, shape)
    r_cyl, arm_cos = values.r_cyl, values.arm_cosine
    disk = shape.k_norm * values.thin_profile * shape.arm_level * (1 + shape.arm_amplitude * arm_cos)
    halo = shape.k_norm * values.halo

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
    densities["old"] = densities["old"] + np.maximum(shape.k_norm * values.thick_disk, 0.0)
    densities["bulge"] = shape.k_norm * values.bulge
    densities["halo"] = np.maximum(halo, 0.0)
    # The components are not in the array form: its points are band bins far wider than a cluster, and the scatter's
    # per-band weights leave components to their own pass (GEN.161, GEN.240); the one-point form is the cell's.
    for name in componentField.POPULATIONS_IN_THE_GATE:
        densities[name] = np.zeros_like(densities["halo"])
    densities["nuclear_cluster"] = shape.k_norm * values.nuclear_cluster
    densities["nuclear_disc"] = shape.k_norm * values.nuclear_disc
    core = coreGas.core_gas(shape)
    if core.young_density:      # the ring's young field stars are the nuclear disc's own (GEN.282)
        moved = np.minimum(shape.k_norm * coreGas.young_density_array(r_cyl, z, core), densities["nuclear_disc"])
        densities["young"] = densities["young"] + moved
        densities["nuclear_disc"] = densities["nuclear_disc"] - moved
    return densities


def integrated_masses(shape, steps=40):
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
        r_max = 12.0 * length
        dr = r_max / (4 * steps)
        radial = sum(2.0 * math.pi * r * math.exp(-r / length) for r in ((i + 0.5) * dr for i in range(4 * steps)))
        return amplitude * radial * dr * 2.0 * height   # the column of `exp(-|z| / h)` is `2 h`

    return {
        "bulge": bulge,
        "thin_disk": disk(shape.arm_level, shape.disk_scale_length_pc, shape.disk_scale_height_pc),
        "thick_disk": disk(terms["thick_disk_amplitude"], terms["thick_disk_scale_length_pc"],
                           terms["thick_disk_scale_height_pc"]),
    }


def thick_disk_units(shape):
    """The thick disk's systems in the model's unnormalized units (the part of `population_units`'s `old` that is not the thin disk's)."""
    terms = model_terms(shape)
    return _disk_mass(terms["thick_disk_amplitude"], terms["thick_disk_scale_length_pc"], terms["thick_disk_scale_height_pc"])


def population_units(shape):
    """
    How many systems each population holds, in the model's unnormalized units (only the ratios mean anything): the thin
    disk shared among `young`, `intermediate` and `old` by its star formation history (`disk_history`), the thick disk
    in `old`, and the bulge, the nuclear cluster and disc and the stellar halo each their own integral. The globular
    clusters and satellite galaxies are not here: their systems come on top of the model (`shape.components`).
    """
    terms = model_terms(shape)
    thin = _disk_mass(shape.arm_level, shape.disk_scale_length_pc, shape.disk_scale_height_pc)
    thick = thick_disk_units(shape)
    return {
        "young": thin * starFormation.disk_population_share("young", shape),
        "intermediate": thin * starFormation.disk_population_share("intermediate", shape),
        "old": thin * starFormation.disk_population_share("old", shape) + thick,
        "bulge": terms["bulge_mass"],
        "nuclear_cluster": terms["nuclear_cluster_mass"],
        "nuclear_disc": terms["nuclear_disc_mass"],
        "halo": halo_mass(shape),
    }


def mean_system_mass_solar(shape):
    """
    The mean mass of the galaxy's systems, solar masses (GEN.545): each population's, worked out from the generator's
    own mass function, remnants and multiplicity over the population's ages (`galaxy.system_mass`), averaged by the
    systems the population holds (`population_units`). The galaxy's stellar mass over this is its number of systems.
    """
    return systemMass.mean_over_solar(population_units(shape), shape)


def total_mass(shape):
    """Everything the model holds, in its unnormalized units: the visible components, the nuclear cluster and disc, and
    the stellar halo (the arm factor averages to 1 around any circle, so the arms add nothing)."""
    terms = model_terms(shape)
    return terms["visible_mass"] + terms["nuclear_cluster_mass"] + terms["nuclear_disc_mass"] + halo_mass(shape)


def halo_mass(shape):
    """The halo's share of the model's stars (the units of `component_masses`), in closed form: the halo
    term integrated over all space, the stellar halo share of the visible components' by construction."""
    terms = halo_terms(shape)
    return terms["halo_amplitude"] * _halo_unit_integral(
        terms["halo_inner_slope"], terms["halo_outer_slope"], terms["halo_core_pc"], terms["halo_break_pc"],
        terms["halo_axis_ratio"])


DensityComponent = namedtuple("DensityComponent", ["name", "population", "amplitude", "mass", "mass_solar", "parameters"])
"""One component of the density model (GEN.379): its `name`, the stellar `population` it fills (`"disk"` is shared by the
young, intermediate and old thin-disk populations), its `amplitude` (the density at its centre in the model's units, or the
lever that scales it; `0` removes the component exactly), its `mass` in the model's units and in solar masses
(`mass_solar`: the units scaled so the components together hold the galaxy's stellar mass), and the `parameters` that give
its shape (scale lengths and heights in parsecs, the bar's angle in radians, the halo's slopes)."""


@functools.lru_cache(maxsize=64)
def density_components(shape):
    """
    The galaxy's density as a list of components (GEN.379), in a fixed order: `thin_disk`, `thick_disk`, `bulge` (the bar
    when the galaxy has one), `halo`, `nuclear_cluster`, `nuclear_disc`. Each component states its amplitude, mass and
    shape (`DensityComponent`); a component with amplitude 0 has mass exactly 0 and adds nothing to the density or to the
    potential, so a galaxy without a thick disk, a bulge or a nucleus is the same list with zeros. The density itself is
    still the sum of the terms above (`relative_density`): this is the list the potential (GEN.338) reads each
    component's mass and shape from, and the one a family with a spheroid (GEN.381) adds its component to.

    The masses are in closed form, so they agree with the numerical integral of each profile (`integrated_masses`, and the
    tests' own quadrature), and they sum to `total_mass(shape)`, which is the galaxy's stellar mass (`build_galaxy_shape`
    solves `k_norm` from it): `mass_solar` is each one's share of `structure.stellar_mass_solar`. The nuclear cluster's mass
    is its nominal (Hernquist) one; its density is held flat inside `softening_pc`, which keeps `softened_mass_fraction` of
    it (about one part in ten thousand of the galaxy's stars, but not the mass the profile's integral gives).

    Returns:
        tuple[DensityComponent]: The components.
    """
    terms, halo = model_terms(shape), halo_terms(shape)
    total = total_mass(shape)
    per_solar = shape.structure.stellar_mass_solar / total
    thin_mass = _disk_mass(shape.arm_level, shape.disk_scale_length_pc, shape.disk_scale_height_pc)
    thick_mass = _disk_mass(terms["thick_disk_amplitude"], terms["thick_disk_scale_length_pc"], terms["thick_disk_scale_height_pc"])
    rows = (
        ("thin_disk", "disk", shape.arm_level, thin_mass,
         {"scale_length_pc": shape.disk_scale_length_pc, "scale_height_pc": shape.disk_scale_height_pc,
          "arm_count": shape.arm_count, "arm_amplitude": shape.arm_amplitude, "pitch_angle_rad": shape.pitch_angle_rad}),
        ("thick_disk", "old", terms["thick_disk_amplitude"], thick_mass,
         {"scale_length_pc": terms["thick_disk_scale_length_pc"], "scale_height_pc": terms["thick_disk_scale_height_pc"]}),
        ("bulge", "bulge", shape.bulge_amplitude, terms["bulge_mass"],
         {"scale_x_pc": shape.bulge_scale_radius_pc, "scale_y_pc": terms["bulge_scale_y_pc"],
          "scale_z_pc": terms["bulge_scale_z_pc"], "bar_angle_rad": terms["bar_angle_rad"]}),
        ("halo", "halo", halo["halo_amplitude"], halo_mass(shape),
         {"core_pc": halo["halo_core_pc"], "break_pc": halo["halo_break_pc"], "axis_ratio": halo["halo_axis_ratio"],
          "inner_slope": halo["halo_inner_slope"], "outer_slope": halo["halo_outer_slope"]}),
        ("nuclear_cluster", "nuclear_cluster", max(shape.nuclear_amplitude, 0.0), terms["nuclear_cluster_mass"],
         {"scale_pc": terms["nuclear_cluster_scale_pc"], "softening_pc": terms["nuclear_cluster_softening_pc"],
          "softened_mass_fraction": _cluster_softened_fraction(terms["nuclear_cluster_scale_pc"],
                                                               terms["nuclear_cluster_softening_pc"])}),
        ("nuclear_disc", "nuclear_disc", max(shape.nuclear_amplitude, 0.0), terms["nuclear_disc_mass"],
         {"scale_length_pc": terms["nuclear_disc_length_pc"], "scale_height_pc": terms["nuclear_disc_height_pc"]}),
    )
    return tuple(DensityComponent(name, population, amplitude, mass, mass * per_solar, parameters)
                 for name, population, amplitude, mass, parameters in rows)


def component_masses(shape):
    """
    Each component's mass in solar masses (GEN.379), `{name: mass}` for the names of `COMPONENT_NAMES`: its profile
    integrated over all space in closed form (`density_components`), scaled so the six together hold the galaxy's stellar
    mass (`shape.structure.stellar_mass_solar`). The potential (GEN.338) reads its bulge, disk and halo masses from here;
    `integrated_masses` is the numerical integral of the same profiles the tests check them against. A component with
    amplitude 0 has mass exactly `0.0`.
    """
    return {component.name: component.mass_solar for component in density_components(shape)}


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
    This galaxy's stellar density at `position_pc`, in units of the real local stellar density
    (`physical_constants.LOCAL_STELLAR_DENSITY_LY3`); the galaxy holds the systems its stellar mass makes
    (`build_galaxy_shape`).

    Args:
        position_pc (tuple): `(x, y, z)`, galaxy-frame parsecs.
        shape (GalaxyShape): The galaxy's shape parameters.

    Returns:
        float: `1.0` where the density is the real local one; `> 1.0` in denser regions
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
    values = component_arrays(x, y, z, shape)
    disk = shape.k_norm * values.thin_profile * shape.arm_level * (1 + shape.arm_amplitude * values.arm_cosine)
    halo = shape.k_norm * values.halo
    nuclear = shape.k_norm * (values.nuclear_cluster + values.nuclear_disc)
    return (np.maximum(disk, 0.0) + np.maximum(shape.k_norm * values.thick_disk, 0.0)
            + np.maximum(halo, 0.0) + np.maximum(shape.k_norm * values.bulge, 0.0) + np.maximum(nuclear, 0.0))


def with_system_richness(shape, systems_per_solar_mass):
    """
    `shape` with the Census dial's richness (GEN.607) set: the galaxy holds `systems_per_solar_mass` systems for each
    solar mass of its stars, so the mean mass of its systems is the inverse, whatever the populations' laws derive.
    `system_mass_scale` is set so that `mean_system_mass_solar` gives it; `k_norm` is not touched (the caller's).
    """
    derived = mean_system_mass_solar(shape._replace(system_mass_scale=1.0))
    return shape._replace(system_mass_scale=derived * systems_per_solar_mass)


def build_galaxy_shape(
    disk_scale_length_pc,
    disk_scale_height_pc,
    bulge_scale_radius_pc,
    bulge_amplitude,
    arm_count,
    pitch_angle_rad,
    arm_amplitude,
    spiral_reference_angle_rad=0.0,
    arm_level=1.0,
    nuclear_amplitude=0.0,
    masses=None,
    structure=galaxyStructure.REFERENCE_STRUCTURE,
    systems_per_solar_mass=None,
):
    """
    Builds a `GalaxyShape` with `k_norm` computed automatically, rather than left for the caller to guess: the
    constant that makes the density, in units of the real local stellar density, hold the galaxy's systems (its
    stellar mass `structure.stellar_mass_solar` over `mean_system_mass_solar`, GEN.404 and GEN.545). No point of the galaxy
    is the calibration point and no Sun position is read.

    Args:
        disk_scale_length_pc, disk_scale_height_pc, bulge_scale_radius_pc,
        bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude:
            The galaxy's shape parameters -- see `GalaxyShape`.
        spiral_reference_angle_rad (float): See `GalaxyShape` -- arbitrary
            fixed orientation, `0.0` unless there's a reason to prefer
            another.
        arm_level (float): The thin disk's density level, the middle of the
            arm crest and the inter-arm trough (ADM.49).
        nuclear_amplitude (float): How much of the observed nuclear star cluster and nuclear stellar disc
            the centre holds (GEN.280): 1 is the observed pair, scaled to the bulge, 0 is none.
        masses (tuple, optional): The galaxy's `(stellar mass, halo mass)`, solar masses (`GalaxyShape.masses`); the halo mass
            is `None` for a galaxy with no halo (`mass_chain.has_halo`).
        structure (galaxy.structure.Structure): The galaxy's thick disk, bar, stellar halo and population motions
            (GEN.532, GEN.533); the Milky Way's published one when none is given.
        systems_per_solar_mass (float, optional): The Census dial's richness (GEN.607): the systems the galaxy holds for each
            solar mass of its stars; `None` takes the mean system mass the generator derives (GEN.545).

    Returns:
        GalaxyShape: With `k_norm` filled in.
    """
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
        masses=None if masses is None else (float(masses[0]), None if masses[1] is None else float(masses[1])),
        structure=structure,
    )

    if systems_per_solar_mass is not None:
        unnormalized = with_system_richness(unnormalized, systems_per_solar_mass)
    systems = structure.stellar_mass_solar / mean_system_mass_solar(unnormalized)
    local_per_pc3 = physical_constants.LOCAL_STELLAR_DENSITY_LY3 / ly_to_pc(1.0) ** 3
    k_norm = systems / (local_per_pc3 * total_mass(unnormalized))

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
