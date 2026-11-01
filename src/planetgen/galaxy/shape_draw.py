# planetgen/galaxy/shape_draw.py

"""
The shape of a spiral galaxy, drawn from its seed (GEN.401; docs/design/galaxy/intergalactic-and-halo-contents.md, section
5.13, rows 4 to 6).

The stellar mass comes from the mass chain (`galaxy.mass_chain`, GEN.400). Everything here follows from it and from
draws about the median laws, each through the one random wrapper (`planetgen.util.draw.Stream`) from
`unit_seed(galaxy_seed, "shape", <link>)`, so one seed gives one shape and two seeds give two:

1. **Disk scale length** from the size-mass relation of late-type galaxies ([R4048], astro-ph/0301527):
   `R50 = gamma M^alpha (1 + M/M0)^(beta - alpha)`, a log-normal draw about it, and `Rd = R50 / 1.678` (the half-light
   radius of a pure exponential disk; Boss chose [R4046]-10-16, 10:09 PM EDT).
2. **Scale height** from the length: `hz = Rd / q`, `q = 7.3` with a 0.13 dex draw (Kregel, van der Kruit and de Grijs
   2002, astro-ph/0204154, Section 5.1).
3. **Stage.** Early (Sa to Sb) or late (Sc and later), the late share falling with mass as the GAMA morphology counts
   of [R2310: Table 1] say. STAND-IN: the table gives cumulative fractions above a mass
   bound, read here as the local fraction, and its classes (S0-Sa, Sab-Scd) are not the two classes of the next link;
   the per-bin fractions are the replacement research (5.13.6 item 6).
4. **Bulge share** `B/T`: `log10(B/D)` drawn about the median of the stage ([R1656: Table 4], K band: -0.51 for T = 1 to 3 and -1.19 for T = 6 to 9) with the table's half width as the scatter.
5. **Bulge and bar scale.** The model has one boxy component, the bulge with its bar. Its major axis is set by the
   bar's length in disk scale lengths, drawn about the stage's mean and deviation ([R1289], astro-ph/0508590,
   Table 5: S0-Sab 1.29 +/- 0.54, Sc-Sd 0.57 +/- 0.44) inside 0.2 to 3 (Erwin's smallest, and [R1486]: no bar reaches 3 scale lengths). The profile `exp(-r_s^2 / 2)` is down to 1.1 percent of its
   centre at three scale radii, taken as the bar's end. The bar's own mass ([R1485]: 0.09 of the total in barred
   galaxies) is not added to `B/T`. Whether a bar is present, the bulge's axis ratios and the bar's angle are
   drawn with the thick disk and the halo (`structure_draw`, GEN.532 and GEN.533) and not stored: DB.66 stores them
   later.
6. **Arms.** The count follows the shares of [R1836] ([R1836: Table 1]: 1, 2, 3, 4 and 5 or more arms
   are 3.9, 64.5, 17.8, 6.5 and 7.3 percent; "5 or more" is taken as 5), independent of mass (their rank correlation
   with total mass is +0.05). The pitch angle falls with mass as the two binned means of Yu and Ho 2019 (ApJ 871, 194)
   ([R4849]: 26.0 degrees at log M* = 9.6, 13.8 at 11.1; no fit is published, so the line through them) with
   their 5 to 6 degree scatter, kept inside Hart's 10 to 40 degrees. The two papers disagree on how strong the mass
   trend is (5.13.6 item 7).

An option the caller sets (`fixed`) replaces its drawn value, and what follows from it is recomputed: the bulge's
amplitude is solved from `B/T` for the disk and bulge actually used. The Structure dial's three options are ratios, so
they keep their meaning when the length moves (GEN.604): `flatness` (the scale length over the scale height, replacing
the absolute height), `bulge_to_total` (replacing the bulge amplitude, which is solved from it) and `bulge_size` (the
bulge's scale radius over the scale length, replacing the absolute radius). A fixed option replaces only its own draw;
the random number of the draw it replaces is left alone, so the same seed with a different length draws the same
residuals.

Only the spiral family exists; GEN.378 adds the others.
"""

import dataclasses
import functools
import math

from planetgen import tuning
from planetgen.galaxy import dials, density, mass_chain, seed as galaxySeed
from planetgen.galaxy import structure as galaxyStructure
from planetgen.galaxy import structure_draw
from planetgen.util import draw

SPIRAL = "spiral"

# --- Row 4: size ([R4048], late-type fit) ---------------------------------------------------------------

SHEN_ALPHA = 0.14
SHEN_BETA = 0.39
SHEN_GAMMA_KPC = 0.10
SHEN_M0_SOLAR = 3.98e10
SHEN_SIGMA_LN_SMALL = 0.47
"""float: The scatter in ln R at low mass (0.20 dex)."""
SHEN_SIGMA_LN_LARGE = 0.34
"""float: The scatter in ln R at high mass (0.15 dex)."""
EXPONENTIAL_HALF_LIGHT_OVER_SCALE = 1.678
"""float: R50 / Rd for a pure exponential disk."""

# --- Row 5: height (Kregel, van [R2484]) -----------------------------------------------------

KREGEL_Q = 7.3
KREGEL_SCATTER_DEX = 0.13

AUTO_SIZE_ATTEMPTS = 32
"""int: The most draws of size and flatness an all-Auto galaxy makes to find a shape the plan represents well."""

# --- Row 6: stage, bulge, bar and arms ---------------------------------------------------------------------------

KELVIN_MASS_BOUNDS = (9.0, 9.5, 10.0, 10.5, 11.0)
"""tuple: The lower bounds, log10 solar masses, of the columns of [R2310: Table 1]."""
KELVIN_S0_SA = (0.18, 0.28, 0.42, 0.48, 0.22)
KELVIN_SAB_SCD = (0.28, 0.36, 0.32, 0.20, 0.06)
"""tuple: Fractions of the galaxies above each bound that are S0-Sa and Sab-Scd (Vmax weighted)."""

LOG_BD_EARLY = -0.51
LOG_BD_EARLY_SIGMA = 0.41
LOG_BD_LATE = -1.19
LOG_BD_LATE_SIGMA = 0.415
"""float: Median `log10(B/D)` and its half width, for T = 1 to 3 (Sa-Sb) and T = 6 to 9, Graham and Worley Table 4."""

BAR_LENGTH_OVER_H_EARLY = (1.29, 0.54)
BAR_LENGTH_OVER_H_LATE = (0.57, 0.44)
BAR_LENGTH_OVER_H_RANGE = (0.2, 3.0)
BAR_LENGTH_PER_SCALE_RADIUS = 3.0
"""float: The bar's semi-major length in scale radii of the model's bulge profile."""

ARM_COUNT_SHARES = ((1, 3.9), (2, 64.5), (3, 17.8), (4, 6.5), (5, 7.3))
"""tuple: `(arms, percent)` from [R1836: Table 1]."""
PITCH_AT_LOW_DEG = 26.0
PITCH_AT_HIGH_DEG = 13.8
PITCH_LOG_MASS_LOW = 9.6
PITCH_LOG_MASS_HIGH = 11.1
PITCH_SCATTER_DEG = 5.5
PITCH_RANGE_DEG = (10.0, 40.0)


@dataclasses.dataclass(frozen=True)
class GalaxyDraw:
    """What one seed draws for a spiral's shape; the values `density.build_galaxy_shape` takes."""

    chain: mass_chain.MassChain
    stage: str
    structure: galaxyStructure.Structure
    disk_scale_length_pc: float
    disk_scale_height_pc: float
    bulge_to_total: float
    bar_length_pc: float
    bulge_scale_radius_pc: float
    bulge_amplitude: float
    arm_count: int
    pitch_angle_deg: float
    compactness: float = 0.0
    """float: The size's distance from the size-mass law at this mass, in standard deviations of its scatter."""


def _stream(galaxy_seed, link):
    return draw.Stream(galaxySeed.unit_seed(galaxy_seed or mass_chain.FALLBACK_SEED, "shape", link))


def size_median_kpc(stellar_mass_solar):
    """[R4048]'s median half-light radius, kpc, of a late-type galaxy of `stellar_mass_solar`."""
    return (SHEN_GAMMA_KPC * stellar_mass_solar ** SHEN_ALPHA
            * (1.0 + stellar_mass_solar / SHEN_M0_SOLAR) ** (SHEN_BETA - SHEN_ALPHA))


def size_sigma_ln(stellar_mass_solar):
    """The scatter in ln R, falling from `SHEN_SIGMA_LN_SMALL` to `SHEN_SIGMA_LN_LARGE` with mass."""
    return SHEN_SIGMA_LN_LARGE + (SHEN_SIGMA_LN_SMALL - SHEN_SIGMA_LN_LARGE) / (1.0 + (stellar_mass_solar / SHEN_M0_SOLAR) ** 2)


def late_share(stellar_mass_solar):
    """The share of the spirals that are late, from [R2310: Table 1]: `Sab-Scd / (S0-Sa + Sab-Scd)` above
    each mass bound, linear in log mass between the bounds and held at the ends."""
    shares = [late / (early + late) for early, late in zip(KELVIN_S0_SA, KELVIN_SAB_SCD)]
    log_mass = math.log10(stellar_mass_solar)
    bounds = KELVIN_MASS_BOUNDS
    if log_mass <= bounds[0]:
        return shares[0]
    for i in range(1, len(bounds)):
        if log_mass <= bounds[i]:
            fraction = (log_mass - bounds[i - 1]) / (bounds[i] - bounds[i - 1])
            return shares[i - 1] + fraction * (shares[i] - shares[i - 1])
    return shares[-1]


def pitch_median_deg(stellar_mass_solar):
    """The pitch angle, degrees, on the line through Yu and Ho's two binned means, kept inside Hart's range."""
    slope = (PITCH_AT_HIGH_DEG - PITCH_AT_LOW_DEG) / (PITCH_LOG_MASS_HIGH - PITCH_LOG_MASS_LOW)
    value = PITCH_AT_LOW_DEG + slope * (math.log10(stellar_mass_solar) - PITCH_LOG_MASS_LOW)
    return min(max(value, PITCH_RANGE_DEG[0]), PITCH_RANGE_DEG[1])


def drawn_disk_scale_length_pc(galaxy_seed, stellar_mass_solar, compactness=None):
    """The thin disc's scale length, parsecs, the seed draws for a galaxy of this stellar mass: the half-light radius of the
    size-mass relation with its scatter, over the exponential disc's half-light-to-scale ratio. A `compactness` (the
    Scale dial's pin, in standard deviations) replaces the draw of the scatter."""
    scatter = (_stream(galaxy_seed, "size").gauss(0.0, size_sigma_ln(stellar_mass_solar)) if compactness is None
               else size_sigma_ln(stellar_mass_solar) * compactness)
    return _length_of_scatter(stellar_mass_solar, scatter)


def _length_of_scatter(stellar_mass_solar, scatter):
    """The thin disc's scale length, parsecs, at `scatter` (a draw of ln R50) from the size-mass relation's median."""
    return (math.exp(math.log(size_median_kpc(stellar_mass_solar)) + scatter)
            / EXPONENTIAL_HALF_LIGHT_OVER_SCALE * 1000.0)


def _drawn_flatness(stream):
    """The thin disc's flatness (scale length over scale height) from the next draw of `stream`: Kregel's median with
    its scatter."""
    return 10.0 ** (math.log10(KREGEL_Q) + stream.gauss(0.0, KREGEL_SCATTER_DEX))


def _auto_size(galaxy_seed, stellar_mass_solar, bar_ratio, crest_of):
    """
    `(disk scale length pc, flatness)` of an all-Auto galaxy: the size and flatness draws of the seed, redrawn from the
    same two streams (the next draw of each) until the shape raises none of the grid's warnings (`dials.shape_warnings`:
    the thin disc's crest under 20 times the one-star density, a scale height under 8 pc, a bulge under 4 pc) and the
    flatness is inside the Structure dial's range, so a default galaxy is one the 4 pc plan represents well (GEN.625).
    The first draw is the one a pinned neighbour would see, so a seed whose first draw is fine keeps it; at most
    `AUTO_SIZE_ATTEMPTS` draws are made, the best crest kept.

    `crest_of(length, flatness)` is the plan's own crest for the galaxy those two make (`galaxy.plan_figures.crest`
    over the shape built at the dials' Auto values), so what is tested is what the plan will print.
    """
    size, height = _stream(galaxy_seed, "size"), _stream(galaxy_seed, "height")
    sigma = size_sigma_ln(stellar_mass_solar)
    flat = dials.lever("flatness")
    best = None
    for _ in range(AUTO_SIZE_ATTEMPTS):
        length = _length_of_scatter(stellar_mass_solar, size.gauss(0.0, sigma))
        q = _drawn_flatness(height)
        if not flat.minimum <= q <= flat.maximum:
            continue
        crest = crest_of(length, q)
        if not dials.shape_warnings(crest, 0.0, length / q, length * bar_ratio / BAR_LENGTH_PER_SCALE_RADIUS):
            return length, q
        if best is None or crest > best[0]:
            best = (crest, length, q)
    if best is None:
        return (_length_of_scatter(stellar_mass_solar, _stream(galaxy_seed, "size").gauss(0.0, sigma)),
                _drawn_flatness(_stream(galaxy_seed, "height")))
    return best[1], best[2]


def _auto_crest(chain, structure, bulge_to_total, bar_ratio, arm_level, arms, pitch_deg, length, q):
    """The plan's crest (`galaxy.plan_figures.crest`) of the galaxy of this chain and structure at disc scale length
    `length` and flatness `q`, every dial at its Auto value and the plan's default sector edge: the figure
    `_auto_size` tests against the grid's warning."""
    from planetgen.galaxy import plan_figures, skeleton
    height = length / q
    bulge_radius = length * bar_ratio / BAR_LENGTH_PER_SCALE_RADIUS
    shape = density.build_galaxy_shape(
        disk_scale_length_pc=length, disk_scale_height_pc=height, bulge_scale_radius_pc=bulge_radius,
        bulge_amplitude=bulge_amplitude_for(bulge_to_total, length, height, bulge_radius, structure, arm_level),
        arm_count=arms, pitch_angle_rad=math.radians(pitch_deg), arm_amplitude=dials.lever("arm_contrast").auto_value,
        arm_level=arm_level, nuclear_amplitude=dials.lever("nuclear_amplitude").auto_value,
        masses=(chain.stellar_mass_solar, chain.halo_mass_solar), structure=structure,
        systems_per_solar_mass=dials.lever("richness").auto_value)
    return plan_figures.crest(shape, skeleton.expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY))


def disk_masses(disk_scale_length_pc, disk_scale_height_pc, structure, arm_level=1.0):
    """`(thin, thick)` disk masses in the model's units (central thin density 1): `2 pi L^2 x 2 h` for the thin
    disk (its `exp(-|z| / h)` profile integrates to `2 h`, GEN.544) times its level, and the thick disk's, whose mass is the
    structure's ratio to the thin disk's (`density.thick_disk_amplitude`)."""
    thin = arm_level * 2.0 * math.pi * disk_scale_length_pc ** 2 * 2.0 * disk_scale_height_pc
    return thin, thin * structure.thick_mass_ratio


def bulge_amplitude_for(bulge_to_total, disk_scale_length_pc, disk_scale_height_pc, bulge_scale_radius_pc, structure,
                        arm_level=1.0):
    """The bulge's central density, against the thin disk's 1, that makes its mass `bulge_to_total` of the bulge plus
    the thin and thick disks (the model's `bulge_mass` formula, `density.model_terms`)."""
    thin, thick = disk_masses(disk_scale_length_pc, disk_scale_height_pc, structure, arm_level)
    bulge_mass = (thin + thick) * bulge_to_total / (1.0 - bulge_to_total)
    unit = (bulge_scale_radius_pc * structure.bar_axis_ratio_y * bulge_scale_radius_pc
            * structure.bar_axis_ratio_z * bulge_scale_radius_pc * density.BULGE_UNIT_INTEGRAL)
    return bulge_mass / unit


def stage_of(galaxy_seed, stellar_mass_solar=None):
    """"late" or "early": the stage the seed's galaxy is (Sc and later, or Sa to Sb), at the mass the seed draws or the
    one pinned."""
    mass = mass_chain.draw_mass_chain(galaxy_seed, stellar_mass_solar).stellar_mass_solar
    return "late" if _stream(galaxy_seed, "stage").random() < late_share(mass) else "early"


def random_nucleus_mass_solar():
    """The mass, solar masses, of the central black hole of a galaxy that no plan names (a black hole or quasar made
    on its own is the nucleus of some other galaxy): a whole galaxy is drawn from the run's random stream and the
    black hole its bulge and stellar mass give is taken (`structure.black_hole_mass_solar`, GEN.402)."""
    drawn = draw_galaxy(draw.getrandbits(128).to_bytes(16, "big"))
    return galaxyStructure.black_hole_mass_solar(drawn.structure, drawn.bulge_to_total)


def structure_for(galaxy_seed, stellar_mass_solar=None, stellar_halo_share=None):
    """The structure of the galaxy a seed makes (`structure_draw.draw_structure`), at the pinned stellar mass and stellar
    halo share when the plan stored them: what a loaded shape is rebuilt with, since no column stores it."""
    return structure_draw.draw_structure(galaxy_seed, stellar_mass_solar, stellar_halo_share)


def _truncated(stream, mean, sigma, low, high):
    for _ in range(50):
        value = stream.gauss(mean, sigma)
        if low <= value <= high:
            return value
    return min(max(mean, low), high)


def draw_galaxy(galaxy_seed, family=SPIRAL, arm_level=1.0, **fixed):
    """
    The shape of the galaxy a seed makes: the same seed gives the same shape.

    Args:
        galaxy_seed (bytes or None): The galaxy's 16-byte seed; `None` uses the mass chain's fallback seed.
        family (str): `SPIRAL` (the only family so far).
        arm_level (float): The thin disk's density level (1 since GEN.611), which the bulge's amplitude is solved
            against.
        **fixed: Any of `stellar_mass_solar`, `compactness` (the size's distance from the size-mass law in sigma),
            `disk_scale_length_pc`, `flatness`, `bulge_to_total`, `bulge_size`, `arm_count`, `pitch_angle_deg`; a value
            given (not `None`) replaces the drawn one. A pinned mass reaches every law that reads it; a pinned
            `disk_scale_length_pc` wins over `compactness`.

    A galaxy with the mass, size, flatness, bulge share and bulge size all on Auto takes its size and flatness from
    `_auto_size`: the seed's draws, redrawn until the shape raises none of the grid's warnings (GEN.625). Pin any of
    those five and the draws are the seed's own, as before.

    Returns:
        GalaxyDraw: The drawn values, with the fixed ones in place of theirs.

    Raises:
        ValueError: For a family other than `SPIRAL`, or an option name that is not one of the six.
    """
    if family != SPIRAL:
        raise ValueError(f"no shape law for the {family!r} family yet")
    unknown = set(fixed) - {"stellar_mass_solar", "compactness", "disk_scale_length_pc", "flatness", "bulge_to_total",
                            "bulge_size", "arm_count", "pitch_angle_deg", "halo_ratio", "stellar_halo_share"}
    if unknown:
        raise ValueError(f"not a shape option: {', '.join(sorted(unknown))}")
    fixed = {name: value for name, value in fixed.items() if value is not None}
    pinned_mass = fixed.get("stellar_mass_solar")
    chain = mass_chain.draw_mass_chain(galaxy_seed, pinned_mass)
    mass = chain.stellar_mass_solar
    if "halo_ratio" in fixed:
        chain = mass_chain.draw_mass_chain(galaxy_seed, mass, halo_mass_of(mass, fixed["halo_ratio"]))

    stage = stage_of(galaxy_seed, pinned_mass)
    late = stage == "late"
    structure = structure_draw.draw_structure(galaxy_seed, pinned_mass, fixed.get("stellar_halo_share"))
    log_bd = _stream(galaxy_seed, "bulge-share").gauss(
        LOG_BD_LATE if late else LOG_BD_EARLY, LOG_BD_LATE_SIGMA if late else LOG_BD_EARLY_SIGMA)
    bulge_to_total = fixed.get("bulge_to_total")
    if bulge_to_total is None:
        bulge_to_total = 1.0 / (1.0 + 10.0 ** -log_bd)
    mean, sigma = BAR_LENGTH_OVER_H_LATE if late else BAR_LENGTH_OVER_H_EARLY
    bar_ratio = _truncated(_stream(galaxy_seed, "bar-length"), mean, sigma, *BAR_LENGTH_OVER_H_RANGE)
    bulge_size = fixed.get("bulge_size")

    arms = fixed.get("arm_count")
    if arms is None:
        point = _stream(galaxy_seed, "arm-count").random() * sum(percent for _, percent in ARM_COUNT_SHARES)
        for arms, percent in ARM_COUNT_SHARES:
            point -= percent
            if point < 0.0:
                break
    pitch = fixed.get("pitch_angle_deg")
    if pitch is None:
        pitch = _truncated(_stream(galaxy_seed, "pitch"), pitch_median_deg(mass), PITCH_SCATTER_DEG, *PITCH_RANGE_DEG)

    length, q = fixed.get("disk_scale_length_pc"), fixed.get("flatness")
    if all(name not in fixed for name in ("stellar_mass_solar", "disk_scale_length_pc", "compactness", "flatness",
                                          "bulge_to_total", "bulge_size")):
        length, q = _auto_size(galaxy_seed, mass, bar_ratio, functools.partial(
            _auto_crest, chain, structure, bulge_to_total, bar_ratio, arm_level, arms, pitch))
    if length is None:
        length = drawn_disk_scale_length_pc(galaxy_seed, mass, fixed.get("compactness"))
    if q is None:
        q = _drawn_flatness(_stream(galaxy_seed, "height"))
    height = length / q

    if bulge_size is None:
        bar_length = length * bar_ratio
        bulge_radius = bar_length / BAR_LENGTH_PER_SCALE_RADIUS
    else:
        bulge_radius = bulge_size * length
        bar_length = bulge_radius * BAR_LENGTH_PER_SCALE_RADIUS
    amplitude = bulge_amplitude_for(bulge_to_total, length, height, bulge_radius, structure, arm_level)

    compactness = (math.log(length * EXPONENTIAL_HALF_LIGHT_OVER_SCALE / 1000.0 / size_median_kpc(mass))
                   / size_sigma_ln(mass))
    return GalaxyDraw(chain, stage, structure, length, height, bulge_to_total, bar_length, bulge_radius, amplitude, arms,
                      pitch, compactness)


# --- The Scale dial (GEN.603): any two of mass, size and density -----------------------------------------------------

@dataclasses.dataclass(frozen=True)
class ScaleResolution:
    """What the Scale dial's pins make of mass, size and density: the three numbers, which were pinned and which
    one the other two gave."""

    stellar_mass_solar: float
    disk_scale_length_pc: float
    compactness: float
    density: float
    pinned: tuple
    derived: str
    """str: `"density"` when mass and size were given, `"size"` when mass and density, `"mass"` when size and density,
    and `"none"` when fewer than two were (the seed draws what was left and the density is a read-out)."""


class ScaleError(ValueError):
    """Pins that cannot make a galaxy: all three of mass, size and density, two ways to set the size, or a solved
    value outside its range."""


def thin_disk_share(structure, bulge_to_total):
    """The thin disk's share of the galaxy's stellar mass: of the stars less the stellar halo's, the part that is
    neither bulge (`bulge_to_total` of bulge, thin and thick disks) nor thick disk (`structure.thick_mass_ratio`
    times the thin disk's)."""
    return (1.0 - structure.halo_stellar_ratio) * (1.0 - bulge_to_total) / (1.0 + structure.thick_mass_ratio)


def central_surface_density(stellar_mass_solar, disk_scale_length_pc, share):
    """The thin disk's central surface density, solar masses per square parsec: `M_thin / (2 pi L^2)` for the
    exponential disk, with `M_thin` the `share` of the stellar mass (`thin_disk_share`)."""
    return stellar_mass_solar * share / (2.0 * math.pi * disk_scale_length_pc ** 2)


def _compactness_of(mass, length_pc):
    return (math.log(length_pc * EXPONENTIAL_HALF_LIGHT_OVER_SCALE / 1000.0 / size_median_kpc(mass))
            / size_sigma_ln(mass))


def _length_of(mass, compactness):
    log_r50 = math.log(size_median_kpc(mass)) + size_sigma_ln(mass) * compactness
    return math.exp(log_r50) / EXPONENTIAL_HALF_LIGHT_OVER_SCALE * 1000.0


def resolve_scale(galaxy_seed, mass=None, compactness=None, length_pc=None, density=None, bulge_to_total=None,
                  mass_range=None, length_range=None):
    """
    Mass, size and density of the galaxy a seed makes under the Scale dial's pins. They are three views of two
    numbers (`central_surface_density`), so any two pinned fix the third and all three are refused.

    Args:
        galaxy_seed (bytes or None): The galaxy's seed.
        mass (float, optional): Pinned stellar mass, solar masses.
        compactness (float, optional): Pinned size, in standard deviations from the size-mass law at the mass.
        length_pc (float, optional): Pinned size as the thin disk scale length, parsecs.
        density (float, optional): Pinned central surface density of the thin disk, solar masses per square parsec.
        bulge_to_total (float, optional): The Structure dial's bulge share, which the thin disk's share of the mass
            depends on; drawn when not given.
        mass_range, length_range (tuple, optional): The `(low, high)` a solved mass or length must lie in.

    Returns:
        ScaleResolution

    Raises:
        ScaleError: For pins that cannot make a galaxy (the message names which pin was dropped and the value the
            other two imply).
    """
    if compactness is not None and length_pc is not None:
        raise ScaleError("the size is set twice, by compactness and by the scale length; give one of them")
    size_pinned = compactness is not None or length_pc is not None
    pinned = tuple(name for name, given in (("mass", mass is not None), ("size", size_pinned),
                                            ("density", density is not None)) if given)

    def share_at(stellar_mass):
        structure = structure_draw.draw_structure(galaxy_seed, stellar_mass)
        if bulge_to_total is not None:
            return thin_disk_share(structure, bulge_to_total)
        late = stage_of(galaxy_seed, stellar_mass) == "late"
        log_bd = _stream(galaxy_seed, "bulge-share").gauss(
            LOG_BD_LATE if late else LOG_BD_EARLY, LOG_BD_LATE_SIGMA if late else LOG_BD_EARLY_SIGMA)
        return thin_disk_share(structure, 1.0 / (1.0 + 10.0 ** -log_bd))

    def length_at(stellar_mass):
        if length_pc is not None:
            return length_pc
        if compactness is not None:
            return _length_of(stellar_mass, compactness)
        return draw_galaxy(galaxy_seed, stellar_mass_solar=stellar_mass).disk_scale_length_pc

    def made(stellar_mass, length, derived):
        return ScaleResolution(
            stellar_mass, length, _compactness_of(stellar_mass, length),
            central_surface_density(stellar_mass, length, share_at(stellar_mass)), pinned, derived)

    if len(pinned) == 3:
        implied = made(mass, length_at(mass), "density")
        raise ScaleError(
            f"mass, size and density are three views of two numbers: give any two. The mass and size given imply a "
            f"central surface density of {implied.density:.4g} solar masses per square parsec; the density given "
            f"({density:.4g}) was dropped")
    if density is None:
        drawn_mass = mass if mass is not None else mass_chain.draw_mass_chain(galaxy_seed).stellar_mass_solar
        return made(drawn_mass, length_at(drawn_mass), "density" if mass is not None and size_pinned else "none")
    if mass is not None or not size_pinned:
        stellar_mass = mass if mass is not None else mass_chain.draw_mass_chain(galaxy_seed).stellar_mass_solar
        length = math.sqrt(stellar_mass * share_at(stellar_mass) / (2.0 * math.pi * density))
        if length_range is not None and not length_range[0] <= length <= length_range[1]:
            raise ScaleError(f"a density of {density:.4g} solar masses per square parsec at {stellar_mass:.4g} solar "
                             f"masses needs a scale length of {length:.4g} pc, outside {length_range[0]:g} to "
                             f"{length_range[1]:g} pc")
        return made(stellar_mass, length, "size")
    # Size and density pinned: the mass is the one that holds both (the density rises with mass at a fixed size, and
    # with mass along the size-mass law, whose exponent is under a half).
    low, high = (mass_range or (1.0e8, 2.0e11))
    log_low, log_high = math.log(low), math.log(high)

    def error(log_mass):
        stellar_mass = math.exp(log_mass)
        return central_surface_density(stellar_mass, length_at(stellar_mass), share_at(stellar_mass)) - density

    if error(log_low) > 0.0 or error(log_high) < 0.0:
        raise ScaleError(f"a density of {density:.4g} solar masses per square parsec at that size needs a stellar "
                         f"mass outside {low:.4g} to {high:.4g} solar masses")
    for _ in range(80):
        middle = 0.5 * (log_low + log_high)
        if error(middle) < 0.0:
            log_low = middle
        else:
            log_high = middle
    stellar_mass = math.exp(0.5 * (log_low + log_high))
    return made(stellar_mass, length_at(stellar_mass), "mass")


DISC_HALF_LIGHT_OVER_R200 = 0.015
"""float: The median ratio of a galaxy's half-light radius to its halo's virial radius ([R2271, R2479],
about 0.015; recalled, not read in the repository, so a placeholder), which the size-halo warning compares a plan with."""

DISC_HALF_LIGHT_OVER_R200_DEX = 1.0
"""float: How far from `DISC_HALF_LIGHT_OVER_R200`, in dex, the plan's ratio may lie before it is warned about; the
relation's scatter is 0.2 to 0.3 dex (recalled), so this is three or four of it (a placeholder)."""


def halo_report(drawn):
    """
    What a plan says about its halo (the Halo dial, GEN.606): `(line, warning)`. The line gives the halo mass over the
    stellar mass, the virial radius and the concentration, or says there is no halo; the warning is a sentence when the
    disc's half-light radius is more than `DISC_HALF_LIGHT_OVER_R200_DEX` from the usual share of the virial radius (a
    pinned size or halo can make a disc the halo could not have held), else `None`.
    """
    chain = drawn.chain
    if chain.halo_mass_solar is None:
        return "No halo: the galaxy is under the halo threshold, so no halo is drawn and no halo law runs", None
    line = (f"Halo: {chain.halo_mass_solar:.3g} solar masses ({chain.halo_mass_solar / chain.stellar_mass_solar:.0f} times "
            f"the stars), virial radius {chain.virial_radius_kpc:.0f} kpc, concentration {chain.concentration:.1f}")
    ratio = drawn.disk_scale_length_pc * EXPONENTIAL_HALF_LIGHT_OVER_SCALE / 1000.0 / chain.virial_radius_kpc
    offset = math.log10(ratio / DISC_HALF_LIGHT_OVER_R200)
    if abs(offset) <= DISC_HALF_LIGHT_OVER_R200_DEX:
        return line, None
    return line, (f"the disc's half-light radius is {ratio:.4f} of the halo's virial radius, {abs(offset):.1f} dex "
                  f"{'above' if offset > 0 else 'below'} the usual {DISC_HALF_LIGHT_OVER_R200} "
                  f"(a pinned size or halo makes a disc this halo would not hold; placeholder limit)")


def halo_mass_of(stellar_mass_solar, halo_ratio):
    """The dark halo mass, solar masses, a pinned `halo_ratio` (halo over stellar mass) gives at a stellar mass: the one
    product every place that reads the pin takes, so a loaded shape and a generated sector agree to the last bit."""
    return halo_ratio * stellar_mass_solar


def edge_limit_pc(edge, disk_scale_length_pc):
    """The distance, parsecs, the Edge dial stops the outline at, or `None` for Auto (the outline's own edge): a length in
    parsecs as it is, one in scale lengths (`12L`) times the thin disc's.

    Raises:
        ValueError: For a relative length of a unit other than `L`.
    """
    if edge is None:
        return None
    if isinstance(edge, dials.RelativeLength):
        if edge.unit != "L":
            raise ValueError(f"{edge.unit!r} is not a unit relative to the galaxy")
        return edge.value * disk_scale_length_pc
    return float(edge)


def pins_of(args, disk_scale_length_pc=None):
    """The shape options the user gave, `{lever name: value}` in the units of the table (the lengths in parsecs): what
    a plan stores beside the resolved shape (DB.81), so the galaxy page can say what was pinned and what was drawn.
    A value that is the row's fixed Auto value counts as pinned when it was given. The Shape selector is a choice, not
    a number (it is stored with the family, DB.66). A length in scale lengths (the Edge) is stored in parsecs, resolved
    against `disk_scale_length_pc`."""
    given = {row.name: getattr(args, row.dest) for row in dials.LEVERS
             if row.kind in ("float", "int", "length") and getattr(args, row.dest, None) is not None}
    return {name: edge_limit_pc(value, disk_scale_length_pc) if isinstance(value, dials.RelativeLength) else value
            for name, value in given.items()}


def pinned_stellar_mass(args, drawn):
    """The stellar mass a plan stores for its galaxy (DB.81): the mass of `drawn` when the Scale levers fixed it, a
    pinned mass or the one a pinned size and density hold, else `None` (the seed's own draw, which a loaded shape draws
    again by itself)."""
    mass_given = dials.resolved(args, "mass") is not None
    size_and_density = dials.resolved(args, "density") is not None and (
        dials.resolved(args, "compactness") is not None or dials.resolved(args, "disk_scale_length") is not None)
    return drawn.chain.stellar_mass_solar if mass_given or size_and_density else None


def default_shape(galaxy_seed):
    """
    The `GalaxyShape` of the galaxy a seed makes with every lever on Auto (what `plan` makes when it is given no shape
    option): what work with no stored shape (a stand-alone command, a test) takes for the galaxy.

    Args:
        galaxy_seed (bytes or None): The galaxy's seed; `None` uses the mass chain's fallback seed.

    Returns:
        GalaxyShape: With `k_norm` and `masses` set.
    """
    return _default_shape(bytes(galaxy_seed or mass_chain.FALLBACK_SEED))


@functools.lru_cache(maxsize=64)
def _default_shape(galaxy_seed):
    import types
    return plan_shape(types.SimpleNamespace(**{row.dest: None for row in dials.LEVERS}), galaxy_seed)[0]


def plan_shape(args, galaxy_seed):
    """
    The `GalaxyShape` a `plan` run builds: the shape options the user gave, the rest drawn from the galaxy seed.

    Args:
        args (argparse.Namespace): The plan's parsed arguments (`add_plan_arguments`); a shape option left at `None`
            is drawn.
        galaxy_seed (bytes or None): The galaxy's seed.

    Returns:
        tuple: `(shape, drawn)`: the `GalaxyShape` with `k_norm` set, and the `GalaxyDraw` it was built from.

    Raises:
        ScaleError: For Scale pins that cannot make a galaxy (`resolve_scale`).
        ValueError: For a pinned or solved value outside a derived limit (`dials.DERIVED_LIMITS`).
    """
    arm_amplitude, arm_level = dials.resolved(args, "arm_contrast"), 1.0
    fixed = {}
    scale_pins = [dials.resolved(args, name) for name in ("mass", "compactness", "disk_scale_length", "density")]
    if any(pin is not None for pin in scale_pins):
        mass_row, length_row = dials.lever("mass"), dials.lever("disk_scale_length")
        scale = resolve_scale(
            galaxy_seed, mass=dials.resolved(args, "mass"), compactness=dials.resolved(args, "compactness"),
            length_pc=dials.resolved(args, "disk_scale_length"), density=dials.resolved(args, "density"),
            bulge_to_total=dials.resolved(args, "bulge_share"),
            mass_range=(mass_row.minimum, mass_row.maximum), length_range=(length_row.minimum, length_row.maximum))
        for row, value in ((mass_row, scale.stellar_mass_solar), (length_row, scale.disk_scale_length_pc)):
            message = row.problem(value)
            if message:
                raise ScaleError(f"the {row.label.lower()} these pins give: {message}")
        fixed = {"stellar_mass_solar": scale.stellar_mass_solar, "disk_scale_length_pc": scale.disk_scale_length_pc}
    drawn = draw_galaxy(
        galaxy_seed, family=dials.resolved(args, "shape") or SPIRAL, arm_level=arm_level, **fixed,
        flatness=dials.resolved(args, "flatness"),
        bulge_to_total=dials.resolved(args, "bulge_share"), bulge_size=dials.resolved(args, "bulge_size"),
        arm_count=dials.resolved(args, "arm_count"), pitch_angle_deg=dials.resolved(args, "pitch_angle_deg"),
        halo_ratio=dials.resolved(args, "halo_ratio"), stellar_halo_share=dials.resolved(args, "stellar_halo_share"))
    if fixed or dials.resolved(args, "flatness") is not None or dials.resolved(args, "bulge_size") is not None:
        for name, value in (("disk_scale_height_pc", drawn.disk_scale_height_pc),
                            ("bulge_scale_radius_pc", drawn.bulge_scale_radius_pc)):
            message = dials.check_derived(name, value)
            if message:
                raise ValueError(message)
    shape = density.build_galaxy_shape(
        disk_scale_length_pc=drawn.disk_scale_length_pc,
        disk_scale_height_pc=drawn.disk_scale_height_pc,
        bulge_scale_radius_pc=drawn.bulge_scale_radius_pc,
        bulge_amplitude=drawn.bulge_amplitude,
        arm_count=drawn.arm_count,
        pitch_angle_rad=math.radians(drawn.pitch_angle_deg),
        arm_amplitude=arm_amplitude,
        arm_level=arm_level,
        nuclear_amplitude=dials.resolved(args, "nuclear_amplitude"),
        masses=(drawn.chain.stellar_mass_solar, drawn.chain.halo_mass_solar),
        structure=drawn.structure,
        systems_per_solar_mass=dials.resolved(args, "richness"),
    )
    return shape, drawn
