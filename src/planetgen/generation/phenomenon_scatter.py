"""
Phenomenon pre-placement (GEN.100): the black holes, neutron stars,
planetary nebulae and supernova remnants of the whole galaxy, its
hypervelocity stars and the nucleus at its center are drawn and placed
right after `planetgen plan`, just before the bright-star scatter, into
`phenomenon_scatter`, while every sector stays unfilled. A sector's fill
later builds each of its rows into the object it stands for (from the
row's own `seed`) at the stored point, and no longer rolls those kinds
itself: what was placed stays where it was put, whichever sectors are
filled and in whatever order.

The walk is the bright-star scatter's (`bright_stars`): shell by shell,
band by band, each band's density averaged over longitude bins, a Poisson
count per kind, each object landing in a wedge of a bin picked in proportion
to the density there. A kind's mean is its rate per star
(`tuning.phenomenon_rate_per_star`) times the stars the same walk expects.
A shell's bands are drawn as numpy arrays (PERF.63, `_band_rows`), about
seven times faster than walking the objects one by one; the random numbers
still come from the shell's `draw.Stream`, in bulk.

The mass cut (GEN.167, docs/design/phenomenon-scatter-mass-cut.md): only
the neutron stars and black holes at or above `--phenomenon-min-mass`
(`tuning.PHENOMENON_MIN_MASS_SOLAR`, 14 solar masses) are scattered, each
kind's mean times its share above the cut (`share_above`), and each built
with its mass drawn above the cut. A sector draws the share below it when
it is filled (`below_cut_draws`, GEN.168), from its own stream, so the two
add up to the whole population. Stellar-mass and intermediate-mass black
holes are separate classes (GEN.166) with `1 - BLACK_HOLE_INTERMEDIATE_MASS_CHANCE`
and `BLACK_HOLE_INTERMEDIATE_MASS_CHANCE` of the black-hole rate.

Not scattered: rogue planets, interstellar comets, asteroid fields and
brown dwarfs (a sector rolls them as before), runaway stars (a flag set at
fill time, since they would be about 1.3 billion rows), emission and
reflection nebulae (they grow around a sector's own hot stars), and
molecular clouds (the galaxy's seeded cloud field already places those
galaxy-wide, `nebula_field`).
"""

import math

import numpy as np

from planetgen import tuning
from planetgen.galaxy import density, nebula_field, profile, remnant_distribution
from planetgen.galaxy import spherical
from planetgen.galaxy.spherical import sector_address_at
from planetgen.galaxy.sector import _sample_poisson_count
from planetgen.generation import bright_stars, multistar
from planetgen.generation.bright_stars import MPC_PER_PC
from planetgen.generation.phenomena import registry
from planetgen.physics import stellar_evolution as stellarEvolution
from planetgen.util import draw
from planetgen.util.poisson import poisson_count
from planetgen.util.random import log_uniform

SCATTERED_KINDS = registry.kinds(scatter_pass="band-density")
"""tuple: The kinds drawn per band from the density (the phenomena registry's
`scatter_pass == "band-density"` rows, GEN.250)."""

MAGNETAR_SUBTYPE = "magnetar"
"""str: The `subtype` of a scattered magnetar (GEN.278), a class of `neutron-star`."""

SCATTER_CLASSES = (
    ("black-hole", "stellar"), ("neutron-star", None), ("planetary-nebula", None),
    ("supernova-remnant", None), ("black-hole", "intermediate"), ("neutron-star", MAGNETAR_SUBTYPE),
)
"""tuple: `(kind, subtype)` in the order each shell's stream draws them (a
new class goes last, so earlier ones keep their draws). A black hole's
subtype is its mass class (GEN.166); a neutron star's, when it has one, is
`magnetar` (GEN.278)."""

FILL_STREAM = "phenomena-fill"
"""str: The per-sector stream name of the below-cut draw (GEN.168)."""

HYPERVELOCITY_KIND = "hypervelocity-star"
"""str: A star ejected by the central black hole, placed with its speed."""

NUCLEUS_KINDS = ("quasar", "black-hole")
"""tuple: What the one nucleus row can be: an active quasar, else the
quiescent supermassive black hole every galaxy has."""

NUCLEUS_ADDRESS = (0, 0, 0)
"""tuple: The `(shell, band, wedge)` the nucleus row sits in (the one cell
per galaxy that rolls for an active nucleus)."""

NUCLEUS_SUBTYPE = "supermassive"
"""str: The `subtype` of a nucleus black hole."""

PHENOMENON_SCATTER_COLUMNS = (
    "shell_index", "band_index", "wedge_index", "kind", "subtype",
    "position_x_mpc", "position_y_mpc", "position_z_mpc",
    "velocity_x_kms", "velocity_y_kms", "velocity_z_kms", "seed", "mass_solar", "map_luminosity_sol",
)
"""tuple: The `phenomenon_scatter` columns a scatter writes, in the order
every row below carries them."""

HYPERVELOCITY_STREAM = "hypervelocity"
NUCLEUS_STREAM = "nucleus"


def _row(shell_index, band_index, wedge, kind, point, rng, subtype=None, velocity=(None, None, None),
         mass_solar=None):
    """One object in `PHENOMENON_SCATTER_COLUMNS` order."""
    return (shell_index, band_index, wedge, kind, subtype, point[0], point[1], point[2],
            velocity[0], velocity[1], velocity[2], rng.getrandbits(63), mass_solar, map_luminosity(mass_solar))


def map_luminosity(mass_solar):
    """The map luminosity (DB.23, GEN.200) stored with a row of `mass_solar`: the luminosity of a main-sequence star
    of that mass by the map's relation; `None` for a row with no mass."""
    return None if mass_solar is None else stellarEvolution.map_luminosity_sol(mass_solar)


def draw_mass(kind, subtype, min_mass_solar, rng):
    """The solar mass (MAP.165) a scattered black hole or neutron star is
    drawn with, from `rng`, at or above the scatter's cut; `None` for a
    kind not drawn by mass. The sector build gives the object this mass, so
    the map and the sector page agree."""
    bounds = mass_range(kind, subtype, min_mass_solar, True)
    if bounds is None:
        return None
    low, high = bounds
    if low >= high:
        return low
    law = mass_law(kind, subtype)
    return log_uniform(low, high, rng=rng) if law[1] else rng.uniform(low, high)


def class_label(kind, subtype):
    """What a log line calls a class: a black hole's mass class in front ("stellar black-hole"), else its kind."""
    return f"{subtype} {kind}" if subtype else kind


def split_label(label):
    """`(kind, subtype)` of a `class_label`; the subtype is `""` for a class without one."""
    subtype, _space, kind = label.rpartition(" ")
    return kind, subtype


EXPECTED_LABELS = tuple(class_label(kind, subtype) for kind, subtype in SCATTER_CLASSES) + (
    HYPERVELOCITY_KIND, class_label("black-hole", NUCLEUS_SUBTYPE), "quasar")
"""tuple: Every class a phenomena scatter can place, so a summary lists the ones that drew none too."""


def mass_law(kind, subtype):
    """`((low, high), logarithmic)`: the solar-mass range a class draws
    from and whether the draw is log-uniform; `None` for a kind not drawn
    by mass (planetary nebulae, supernova remnants)."""
    if kind == "neutron-star":
        return tuning.NEUTRON_STAR_MASS_RANGE_SOLAR, False
    if kind == "black-hole" and subtype == "intermediate":
        return tuning.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR, True
    if kind == "black-hole" and subtype == "stellar":
        return tuning.BLACK_HOLE_MASS_RANGE_SOLAR, False
    return None


def class_rate_per_star(kind, subtype):
    """A class's expected count per star: the kind's rate, split between
    the two black-hole classes by `BLACK_HOLE_INTERMEDIATE_MASS_CHANCE`. The
    magnetars' (GEN.278) is `MAGNETAR_PER_NEUTRON_STAR` of the neutron stars',
    per star of the young population they follow (`MAGNETAR_YOUNG_STAR_SHARE`)."""
    rate = tuning.phenomenon_rate_per_star(kind)
    if subtype == MAGNETAR_SUBTYPE:
        return rate * tuning.MAGNETAR_PER_NEUTRON_STAR / tuning.MAGNETAR_YOUNG_STAR_SHARE
    if kind == "black-hole":
        chance = tuning.BLACK_HOLE_INTERMEDIATE_MASS_CHANCE
        rate *= chance if subtype == "intermediate" else 1.0 - chance
    return rate


def share_above(kind, subtype, min_mass_solar):
    """S_k(c): the share of a class's mass law at or above `min_mass_solar`
    (1 for a kind not drawn by mass). Neutron stars `(2.2 - c) / 1.1`,
    stellar-mass black holes `(20 - c) / 15`, intermediate-mass
    `ln(1e5 / c) / ln(1e3)`, clamped to 0 to 1. The magnetars are all
    scattered whatever the cut (1): the fill never draws one (GEN.278)."""
    law = mass_law(kind, subtype)
    if law is None or subtype == MAGNETAR_SUBTYPE:
        return 1.0
    (low, high), logarithmic = law
    if min_mass_solar <= low:
        return 1.0
    if min_mass_solar >= high:
        return 0.0
    if logarithmic:
        return math.log(high / min_mass_solar) / math.log(high / low)
    return (high - min_mass_solar) / (high - low)


def mass_range(kind, subtype, min_mass_solar, above):
    """The `(low, high)` a class's mass is drawn from on one side of the
    cut: `above` for a scattered object, below for a sector's own. `None`
    for a kind not drawn by mass, or with no cut."""
    law = mass_law(kind, subtype)
    if law is None or min_mass_solar is None:
        return None
    low, high = law[0]
    if subtype == MAGNETAR_SUBTYPE:
        return (low, high) if above else None      # the cut does not touch them (share_above)
    return (max(min_mass_solar, low), high) if above else (low, min(min_mass_solar, high))


PLACED_KINDS = {
    ("black-hole", "stellar"): "black-hole", ("black-hole", "intermediate"): "intermediate-mass-black-hole",
    ("neutron-star", None): "neutron-star", ("planetary-nebula", None): "white-dwarf",
}
"""dict: `(kind, subtype)` of a scattered class that carries a member tree -> the `multistar.Placed` kind of its
placed object (the central star of a planetary nebula is a white dwarf). A supernova remnant, a quasar, a
hypervelocity star and the nucleus carry none (multi-star-systems.md 4.3)."""

NEUTRON_STAR_RECYCLED_SHARE = tuning.PULSAR_RECYCLED_CHANCE * tuning.PULSAR_RATE_SCALE
"""float: The share of scattered neutron stars that are recycled millisecond pulsars (pulsing, then recycled)."""


def row_density_divisor(kind, subtype):
    """
    `1 + p_placed`, the divisor of a class's object rate for the scatter (multi-star-systems.md 6.4): a row is a
    system keyed by its placed object, and the same-class companions of the placed objects are objects of the class
    too, so the rows' rate is the object rate over one plus the expected same-class companions per row. 1 for a
    class with no tree or no same-class companion (the whole of the intermediate-mass black holes).
    """
    placed = PLACED_KINDS.get((kind, subtype))
    if placed is None or kind == "planetary-nebula":      # a white dwarf companion is not a planetary nebula
        return 1.0
    return 1.0 + multistar.same_class_companions(placed, NEUTRON_STAR_RECYCLED_SHARE if placed == "neutron-star" else 0.0)


def _scattered_rates(min_mass_solar):
    """`(kind, subtype, rate per star above the cut)` for each class with
    anything above it, in `SCATTER_CLASSES` order. The rate is of rows, that is of systems keyed by their
    placed object (`row_density_divisor`)."""
    rates = []
    for kind, subtype in SCATTER_CLASSES:
        rate = (class_rate_per_star(kind, subtype) * share_above(kind, subtype, min_mass_solar)
                / row_density_divisor(kind, subtype))
        if rate > 0.0:
            rates.append((kind, subtype, rate))
    return rates


def scatter_shell(shape, shell_index, high_band, edge_pc, expected_at_density_1, seed, skip_addresses=None,
                  min_mass_solar=None):
    """
    One shell of the scatter: every black hole and neutron star at or
    above `min_mass_solar` (default `tuning.PHENOMENON_MIN_MASS_SOLAR`),
    and every planetary nebula and supernova remnant, in the bands
    `-high_band .. high_band` of shell `shell_index`. Each shell draws from its
    own stream (the seed and the shell index), so shells can be drawn in any
    order, or in worker processes, and give the same objects.

    The shell is drawn as arrays (PERF.63, `_band_rows`): one Poisson count per
    band and class, then the bins, wedges, points and masses of all the
    shell's objects at once, from the shell's stream in bulk.

    Yields:
        tuple: One row per object, in `PHENOMENON_SCATTER_COLUMNS` order.
    """
    if min_mass_solar is None:
        min_mass_solar = tuning.PHENOMENON_MIN_MASS_SOLAR
    rates = _scattered_rates(min_mass_solar)
    if not rates or high_band < 0:
        return
    rng = draw.Stream(f"{seed}:phenomena:shell:{shell_index}")
    yield from _band_rows(shape, shell_index, high_band, edge_pc, expected_at_density_1, rng, rates,
                          skip_addresses or set(), min_mass_solar)


GAS_WEIGHTED_KINDS = ("supernova-remnant",)
"""tuple: The scatter kinds whose bin weight also carries the gas fraction (GEN.236): a remnant shines only where its
shell meets gas, so none lies in the halo or a cluster. The sector profile's `kind_factors` holds the same factor for
the fill."""


OLD_STAR_KIND = "planetary-nebula"
"""str: The scatter kind whose bin weight also carries the old stars' share against the solar circle's (GEN.528):
a planetary nebula comes from a star that has had time to leave the asymptotic giant branch, so a place of mostly young
stars makes fewer. The sector profile's `kind_factors` holds the same factor for the fill."""


def _gas_weight(x, y, z, shape):
    """The gas fraction (`nebula_field.gas_factor_array` capped at 1, as `SectorProfile.kind_factors` has it) times the
    share of the stars that are not a cluster's, at arrays of points: the weight a gas-born kind carries on top of the
    stars'."""
    gas = np.minimum(1.0, np.maximum(nebula_field.gas_factor_array(x, y, z, shape), 0.0))
    mix = density.population_densities_array(x, y, z, shape)
    total = sum(np.maximum(value, 0.0) for value in mix.values())
    cluster = np.where(total > 0.0, np.maximum(mix["cluster"], 0.0) / np.where(total > 0.0, total, 1.0), 0.0)
    return gas * (1.0 - cluster)


def _band_rows(shape, shell_index, high_band, edge_pc, expected_at_density_1, rng, rates, skip_addresses,
               min_mass_solar):
    """
    Bands `-high_band .. high_band` of a shell, all at once. Per class each band's Poisson count
    has mean `rate * e * wedges * volume / bins * sum(bin weights)` (the volume in edge cubes); the shell's
    objects then get a longitude bin in proportion to the band's weights (one `searchsorted` over the bands'
    cumulative sums), a wedge under the bin's angle, a uniform point in the wedge's cell, and a mass; a point
    that rounding carried out of its cell is dropped.
    """
    bins = bright_stars.ANGLE_BINS
    bands = np.arange(-high_band, high_band + 1)
    band_count = len(bands)
    wedges = np.array([spherical.wedge_count(shell_index, int(band)) for band in bands], dtype=np.int64)
    volume = np.array([bright_stars.volume_ratio(shell_index, int(band), edge_pc) for band in bands])
    centres = [spherical.band_center_pc(shell_index, int(band), edge_pc) for band in bands]
    plane, z_mid = np.array([c[0] for c in centres]), np.array([c[1] for c in centres])
    colat = [spherical.band_colatitude_bounds(shell_index, int(band)) for band in bands]
    cos_north, cos_south = np.cos([c[0] for c in colat]), np.cos([c[1] for c in colat])
    bin_angle = (np.arange(bins) + 0.5) * (2.0 * math.pi / bins)
    x_bins, y_bins = plane[:, None] * np.cos(bin_angle)[None, :], plane[:, None] * np.sin(bin_angle)[None, :]
    base = density.total_density_array(x_bins, y_bins, z_mid[:, None], shape)
    r_low, r_high = spherical.shell_bounds_pc(shell_index, edge_pc)
    row_of_band = np.arange(band_count)
    tables = {}
    for kind, subtype, rate in rates:
        young = registry.scatter_class_gate(kind, subtype) == "young stars"
        table = (kind, subtype) if young else kind
        if table not in tables:
            if young:
                # The young population and its arm contrast, none in the halo or a cluster (GEN.278); a magnetar is
                # too young to have drifted, so no remnant thickening.
                weights = density.population_densities_array(x_bins, y_bins, z_mid[:, None], shape)["young"]
            else:
                factor = np.full(band_count, 1.0)
                if kind in tuning.REMNANT_SCALE_HEIGHT_RATIO:
                    factor = (remnant_distribution.vertical_factor_array(kind, z_mid, shape.disk_scale_height_pc)
                              * remnant_distribution.radial_factor_array(kind, plane))
                weights = base * factor[:, None]
                if kind in GAS_WEIGHTED_KINDS:
                    weights = weights * _gas_weight(x_bins, y_bins, z_mid[:, None], shape)
                if kind == OLD_STAR_KIND:
                    weights = weights * profile.planetary_nebula_factor_array(x_bins, y_bins, z_mid[:, None], shape)
            cumulative = np.cumsum(weights, axis=1)
            top = cumulative[:, -1]
            # Rows are offset by 2 each, so one sorted array holds every band's cumulative shares.
            shares = cumulative / np.where(top > 0.0, top, 1.0)[:, None] + 2.0 * row_of_band[:, None]
            tables[table] = (top, shares.ravel())
        top, shares = tables[table]
        means = rate * expected_at_density_1 * (wedges * volume / bins) * top
        counts = np.array([poisson_count(float(mean), rng) if mean > 0.0 else 0 for mean in means], dtype=np.int64)
        total = int(counts.sum())
        if not total:
            continue
        band_row = np.repeat(row_of_band, counts)
        uniform = rng.random_array(6 * total).reshape(6, total)
        angle_bin = np.minimum(np.searchsorted(shares, 2.0 * band_row + uniform[0], side="right")
                               - band_row * bins, bins - 1)
        wedge_count = wedges[band_row]
        wedge_width = 2.0 * math.pi / wedge_count
        wedge = np.minimum(((angle_bin + uniform[1]) * (2.0 * math.pi / bins) / wedge_width).astype(np.int64),
                           wedge_count - 1)
        lon = (wedge + uniform[2]) * wedge_width
        # Uniform over the cell's volume: the cube of the radius and the cosine of the colatitude are uniform.
        radius = np.cbrt(r_low ** 3 + uniform[3] * (r_high ** 3 - r_low ** 3))
        cos_colat = cos_south[band_row] + uniform[4] * (cos_north[band_row] - cos_south[band_row])
        in_plane = radius * np.sqrt(np.maximum(0.0, 1.0 - cos_colat * cos_colat))
        mpc = np.stack([np.rint(in_plane * np.cos(lon) * MPC_PER_PC), np.rint(in_plane * np.sin(lon) * MPC_PER_PC),
                        np.rint(radius * cos_colat * MPC_PER_PC)]).astype(np.int64)
        band_index = bands[band_row]
        kept = _stored_in_cells(mpc, shell_index, band_index, wedge, wedge_count, edge_pc)
        masses = _masses(kind, subtype, min_mass_solar, uniform[5])
        seeds = rng.randbits_list(63, total)
        columns = [band_index.tolist(), wedge.tolist(), kept.tolist(), seeds, *mpc.tolist(),
                   masses.tolist() if masses is not None else [None] * total]
        for band, wedge_here, in_cell, seed, x, y, z, mass in zip(*columns):
            if in_cell and (shell_index, band, wedge_here) not in skip_addresses:
                yield (shell_index, band, wedge_here, kind, subtype, x, y, z, None, None, None, seed, mass,
                       map_luminosity(mass))


def _masses(kind, subtype, min_mass_solar, uniform):
    """The solar masses `draw_mass` gives from the uniform draws `uniform` (one per object); `None` for a
    kind not drawn by mass."""
    bounds = mass_range(kind, subtype, min_mass_solar, True)
    if bounds is None:
        return None
    low, high = bounds
    if low >= high:
        return np.full(len(uniform), low)
    if mass_law(kind, subtype)[1]:
        return np.exp(math.log(low) + uniform * (math.log(high) - math.log(low)))
    return low + uniform * (high - low)


def _stored_in_cells(mpc, shell_index, band_index, wedge, wedge_count, edge_pc):
    """Whether each stored point (whole milliparsecs, one column each of `mpc`) is still in its cell after
    rounding: a float test with a margin of a few milliparsecs, and the exact `sector_address_at` for the
    few points within the margin of a face."""
    x, y, z = mpc / MPC_PER_PC
    margin = 4.0 / MPC_PER_PC
    in_plane = np.sqrt(x * x + y * y)
    radius = np.sqrt(in_plane * in_plane + z * z)
    r_low, r_high = spherical.shell_bounds_pc(shell_index, edge_pc)
    colat = np.arctan2(in_plane, z)
    bounds = {int(band): spherical.band_colatitude_bounds(shell_index, int(band)) for band in np.unique(band_index)}
    north = np.array([bounds[int(band)][0] for band in band_index])
    south = np.array([bounds[int(band)][1] for band in band_index])
    width = 2.0 * math.pi / wedge_count
    lon = np.arctan2(y, x) % (2.0 * math.pi)
    reach = np.maximum(radius, 1e-9)
    colat_margin = margin / reach
    lon_margin = margin / np.maximum(in_plane, 1e-9)
    inside = ((radius > r_low + margin) & (radius < r_high - margin)
              & (colat > north + colat_margin) & (colat < south - colat_margin)
              & ((wedge_count == 1) | ((lon > wedge * width + lon_margin) & (lon < (wedge + 1) * width - lon_margin))))
    # A cap (or the core) has one wedge and its margin at the pole is the whole circle; its colatitude test above
    # still holds: only the poles' own edge can fall to the exact test.
    for i in np.nonzero(~inside)[0].tolist():
        inside[i] = (spherical.sector_address_at((x[i], y[i], z[i]), edge_pc)
                     == (shell_index, int(band_index[i]), int(wedge[i])))
    return inside


def below_cut_draws(address, center_pc, shape, expected_stars, min_mass_solar, seed):
    """
    The neutron stars and black holes a sector draws itself when the
    scatter placed only those at or above `min_mass_solar` (GEN.168): per
    class, a Poisson count with mean the class's rate times the sector's
    expected stars (`expected_stars`, the density the scatter used, not
    the sector's rolled star count) times its regional factor at the
    sector's center, times the share below the cut. Drawn from the
    sector's own stream (`FILL_STREAM`), so a sector's objects do not
    depend on which other sectors are filled, or in what order.

    Returns:
        tuple: `(draws, rng)`: `draws` lists `(kind, subtype, mass_range,
            count)` for each class with a count; `rng` is the stream, for
            building them.
    """
    rng = draw.Stream(f"{seed}:{FILL_STREAM}:{address[0]}:{address[1]}:{address[2]}")
    draws = []
    for kind, subtype in SCATTER_CLASSES:
        below = 1.0 - share_above(kind, subtype, min_mass_solar)
        if below <= 0.0:
            continue
        mean = class_rate_per_star(kind, subtype) * expected_stars * below
        if kind in tuning.REMNANT_SCALE_HEIGHT_RATIO:
            mean *= remnant_distribution.placement_factor(kind, center_pc, shape.disk_scale_height_pc)
        count = _sample_poisson_count(mean, rng=rng)
        if count:
            draws.append((kind, subtype, mass_range(kind, subtype, min_mass_solar, above=False), count))
    return draws, rng


def shell_expected(shape, shell_index, high_band, edge_pc, expected_at_density_1, min_mass_solar=None):
    """About how many objects `scatter_shell` places in one shell (the
    progress bar's weight): the same per-band means on evenly spaced
    sample bands, as `bright_stars.shell_expected_stars` does for stars."""
    if min_mass_solar is None:
        min_mass_solar = tuning.PHENOMENON_MIN_MASS_SOLAR
    fractions = {population: 1.0 for population in bright_stars.POPULATIONS}
    rates = _scattered_rates(min_mass_solar)
    total = 0.0
    for _band, plane, z, stars in bright_stars.sampled_band_expectations(
            shape, shell_index, high_band, edge_pc, expected_at_density_1, fractions):
        total += stars * sum(rate * _expectation_factor(shape, kind, subtype, plane, z)
                             for kind, subtype, rate in rates)
    return total


def _expectation_factor(shape, kind, subtype, plane, z):
    """What a class's rate per star is multiplied by at a band centre `plane` and `z` parsecs out (the progress bar's
    estimate of `_band_rows`): the remnant factors of a neutron star or black hole, the gas fraction of a supernova
    remnant (GEN.236), the old stars' share of a planetary nebula (GEN.528), the young population's share of the stars
    for a class that follows it (GEN.278)."""
    if registry.scatter_class_gate(kind, subtype) == "young stars":
        shares = density.population_densities((plane, 0.0, z), shape)
        total = sum(shares.values())
        return shares["young"] / total if total > 0.0 else 0.0
    factor = (remnant_distribution.vertical_factor(kind, z, shape.disk_scale_height_pc)
              * remnant_distribution.radial_factor(kind, plane))
    if kind in GAS_WEIGHTED_KINDS:
        factor *= min(1.0, max(nebula_field.gas_factor((plane, 0.0, z), shape), 0.0))
    if kind == OLD_STAR_KIND:
        factor *= float(profile.planetary_nebula_factor_array(plane, 0.0, z, shape))
    return factor


def nucleus_row(seed):
    """
    The galaxy's nucleus: an active quasar with
    `tuning.QUASAR_ACTIVE_NUCLEUS_CHANCE`, else the quiescent supermassive
    black hole every galaxy has, at the galactic origin in `NUCLEUS_ADDRESS`'s
    cell. Drawn once per scatter, from its own stream.
    """
    rng = draw.Stream(f"{seed}:phenomena:{NUCLEUS_STREAM}")
    if rng.random() < tuning.QUASAR_ACTIVE_NUCLEUS_CHANCE:
        return _row(*NUCLEUS_ADDRESS, "quasar", (0, 0, 0), rng)
    return _row(*NUCLEUS_ADDRESS, "black-hole", (0, 0, 0), rng, subtype=NUCLEUS_SUBTYPE,
                mass_solar=log_uniform(*tuning.BLACK_HOLE_SUPERMASSIVE_MASS_RANGE_SOLAR, rng=rng))


def hypervelocity_rows(extents, edge_pc, seed):
    """
    The galaxy's hypervelocity stars: `tuning.HYPERVELOCITY_STARS_PER_GALAXY`
    (log-uniform between its two ends) stars, each ejected from the central
    black hole in a random direction at `HYPERVELOCITY_STAR_SPEED_RANGE_KMS`
    and placed along that ray, moving outward. One that would land outside
    the outline is redrawn nearer in.

    Args:
        extents (list): `(shell_index, high_band)` per shell.
        edge_pc (float): The sector edge.
        seed (int): The scatter's seed.

    Yields:
        tuple: One row per star, in `PHENOMENON_SCATTER_COLUMNS` order.
    """
    high_bands = dict(extents)
    if not high_bands:
        return
    rng = draw.Stream(f"{seed}:phenomena:{HYPERVELOCITY_STREAM}")
    low, high = tuning.HYPERVELOCITY_STARS_PER_GALAXY
    wanted = int(round(log_uniform(low, high, rng=rng)))
    reach_pc = (max(high_bands) + 0.5) * edge_pc
    placed = 0
    attempts = 0
    while placed < wanted and attempts < wanted * HYPERVELOCITY_ATTEMPTS:
        attempts += 1
        z = rng.uniform(-1.0, 1.0)
        phi = rng.uniform(0.0, 2.0 * math.pi)
        planar = math.sqrt(1.0 - z * z)
        direction = (planar * math.cos(phi), planar * math.sin(phi), z)
        distance = reach_pc * rng.random()
        point = tuple(round(component * distance * MPC_PER_PC) for component in direction)
        address = sector_address_at(tuple(value / MPC_PER_PC for value in point), edge_pc)
        if abs(address[1]) > high_bands.get(address[0], -1):
            continue
        speed = log_uniform(*tuning.HYPERVELOCITY_STAR_SPEED_RANGE_KMS, rng=rng)
        yield _row(*address, HYPERVELOCITY_KIND, point, rng,
                   velocity=tuple(speed * component for component in direction))
        placed += 1


def special_rows(extents, edge_pc, seed, filled=()):
    """The nucleus and the hypervelocity stars, minus any in a cell
    already `filled`: the rows no shell draws."""
    rows = [nucleus_row(seed)] + list(hypervelocity_rows(extents, edge_pc, seed))
    return [row for row in rows if tuple(row[:3]) not in filled]


HYPERVELOCITY_ATTEMPTS = 50
"""int: Draws per wanted star before giving up on one (most land inside
the outline; a thin galaxy loses the ones that fly out of its disk)."""
