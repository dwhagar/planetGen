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

The walk is the bright-star scatter's (`bright_stars`): ring by ring,
each ring's density averaged over angle bins, a Poisson count per kind,
each object landing in a slot of a bin picked in proportion to the density
there. A kind's mean is its rate per star
(`tuning.phenomenon_rate_per_star`) times the stars the same walk expects.

The mass cut (GEN.167, docs/design/phenomenon-scatter-mass-cut.md): only
the neutron stars and black holes at or above `--phenomenon-min-mass`
(`tuning.PHENOMENON_MIN_MASS_SOLAR`, 8 solar masses) are scattered, each
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


from planetgen import tuning
from planetgen.galaxy import remnant_distribution
from planetgen.galaxy.geometry import sector_address_at
from planetgen.galaxy.sector import _sample_poisson_count
from planetgen.generation import bright_stars
from planetgen.generation.bright_stars import MPC_PER_PC, _place_one
from planetgen.galaxy.geometry import layer_center_z_pc, ring_radius_pc
from planetgen.util import draw
from planetgen.util.random import log_uniform

SCATTERED_KINDS = ("black-hole", "neutron-star", "planetary-nebula", "supernova-remnant")
"""tuple: The kinds drawn per ring from the density."""

SCATTER_CLASSES = (
    ("black-hole", "stellar"), ("neutron-star", None), ("planetary-nebula", None),
    ("supernova-remnant", None), ("black-hole", "intermediate"),
)
"""tuple: `(kind, subtype)` in the order each ring's stream draws them (a
new class goes last, so earlier ones keep their draws). A black hole's
subtype is its mass class (GEN.166)."""

FILL_STREAM = "phenomena-fill"
"""str: The per-sector stream name of the below-cut draw (GEN.168)."""

HYPERVELOCITY_KIND = "hypervelocity-star"
"""str: A star ejected by the central black hole, placed with its speed."""

NUCLEUS_KINDS = ("quasar", "black-hole")
"""tuple: What the one nucleus row can be: an active quasar, else the
quiescent supermassive black hole every galaxy has."""

NUCLEUS_ADDRESS = (0, 0, 0)
"""tuple: The `(ring, layer, slot)` the nucleus row sits in (the one cell
per galaxy that rolls for an active nucleus)."""

NUCLEUS_SUBTYPE = "supermassive"
"""str: The `subtype` of a nucleus black hole."""

PHENOMENON_SCATTER_COLUMNS = (
    "ring_index", "layer_index", "ring_slot_index", "kind", "subtype",
    "position_x_mpc", "position_y_mpc", "position_z_mpc",
    "velocity_x_kms", "velocity_y_kms", "velocity_z_kms", "seed",
)
"""tuple: The `phenomenon_scatter` columns a scatter writes, in the order
every row below carries them."""

HYPERVELOCITY_STREAM = "hypervelocity"
NUCLEUS_STREAM = "nucleus"


def _row(ring_index, layer_index, slot, kind, point, rng, subtype=None, velocity=(None, None, None)):
    """One object in `PHENOMENON_SCATTER_COLUMNS` order."""
    return (ring_index, layer_index, slot, kind, subtype, point[0], point[1], point[2],
            velocity[0], velocity[1], velocity[2], rng.getrandbits(63))


def class_label(kind, subtype):
    """What a log line calls a class: a black hole's mass class in front ("stellar black-hole"), else its kind."""
    return f"{subtype} {kind}" if subtype else kind


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
    the two black-hole classes by `BLACK_HOLE_INTERMEDIATE_MASS_CHANCE`."""
    rate = tuning.phenomenon_rate_per_star(kind)
    if kind == "black-hole":
        chance = tuning.BLACK_HOLE_INTERMEDIATE_MASS_CHANCE
        rate *= chance if subtype == "intermediate" else 1.0 - chance
    return rate


def share_above(kind, subtype, min_mass_solar):
    """S_k(c): the share of a class's mass law at or above `min_mass_solar`
    (1 for a kind not drawn by mass). Neutron stars `(2.2 - c) / 1.1`,
    stellar-mass black holes `(20 - c) / 15`, intermediate-mass
    `ln(1e5 / c) / ln(1e3)`, clamped to 0 to 1."""
    law = mass_law(kind, subtype)
    if law is None:
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
    return (max(min_mass_solar, low), high) if above else (low, min(min_mass_solar, high))


def _scattered_rates(min_mass_solar):
    """`(kind, subtype, rate per star above the cut)` for each class with
    anything above it, in `SCATTER_CLASSES` order."""
    rates = []
    for kind, subtype in SCATTER_CLASSES:
        rate = class_rate_per_star(kind, subtype) * share_above(kind, subtype, min_mass_solar)
        if rate > 0.0:
            rates.append((kind, subtype, rate))
    return rates


def scatter_layer(shape, layer_index, outer_ring, edge_pc, expected_at_density_1, seed, skip_addresses=None,
                  min_mass_solar=None):
    """
    One layer of the scatter: every black hole and neutron star at or
    above `min_mass_solar` (default `tuning.PHENOMENON_MIN_MASS_SOLAR`),
    and every planetary nebula and supernova remnant, from ring 0 out to
    `outer_ring` at `layer_index`. Each layer draws from its own stream
    (the seed and the layer index), so layers can be drawn in any order,
    or in worker processes, and give the same objects.

    Yields:
        tuple: One row per object, in `PHENOMENON_SCATTER_COLUMNS` order.
    """
    if min_mass_solar is None:
        min_mass_solar = tuning.PHENOMENON_MIN_MASS_SOLAR
    rates = _scattered_rates(min_mass_solar)
    rng = draw.Stream(f"{seed}:phenomena:{layer_index}")
    skip_addresses = skip_addresses or set()
    for ring_index in range(outer_ring + 1):
        slots, bins = bright_stars._ring_bins(ring_index, layer_index, shape, expected_at_density_1, edge_pc)
        base = [sum(densities.values()) for densities in bins]
        if sum(base) <= 0.0:
            continue
        slots_per_bin = slots / len(bins)
        weights_by_kind = {}
        for kind, subtype, rate in rates:
            if kind not in weights_by_kind:
                weights = _kind_weights(kind, base, ring_index, layer_index, edge_pc, shape.disk_scale_height_pc)
                weights_by_kind[kind] = (weights, sum(weights))
            weights, weight_sum = weights_by_kind[kind]
            count = _sample_poisson_count(rate * expected_at_density_1 * slots_per_bin * weight_sum, rng=rng)
            for _ in range(count):
                spot = _place_one(rng, weights, ring_index, layer_index, slots, edge_pc)
                if spot is None or (ring_index, layer_index, spot[0]) in skip_addresses:
                    continue
                yield _row(ring_index, layer_index, spot[0], kind, spot[1], rng, subtype=subtype)


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


def _kind_weights(kind, base, ring_index, layer_index, edge_pc, thin_height_pc):
    """The per-bin weights for `kind`: the stellar density of each bin times
    the kind's regional factor at the bin's centerline point (GEN.132; a
    kind with no research stays flat)."""
    if kind not in tuning.REMNANT_SCALE_HEIGHT_RATIO:
        return base
    radius = ring_radius_pc(ring_index, edge_pc)
    z = layer_center_z_pc(layer_index, edge_pc)
    count = len(base)
    weights = []
    for k, weight in enumerate(base):
        theta = (k + 0.5) * 2 * math.pi / count
        point = (radius * math.cos(theta), radius * math.sin(theta), z)
        weights.append(weight * remnant_distribution.placement_factor(kind, point, thin_height_pc))
    return weights


def layer_expected(shape, layer_index, outer_ring, edge_pc, expected_at_density_1, min_mass_solar=None):
    """About how many objects `scatter_layer` places in one layer (the
    progress bar's weight): the same per-ring means on evenly spaced
    sample rings, as `bright_stars.layer_expected_stars` does for stars."""
    if min_mass_solar is None:
        min_mass_solar = tuning.PHENOMENON_MIN_MASS_SOLAR
    fractions = {population: 1.0 for population in bright_stars.POPULATIONS}
    stars = bright_stars.layer_expected_stars(shape, layer_index, outer_ring, edge_pc, expected_at_density_1,
                                              fractions)
    # The layer-centre vertical factor stands in for each ring's regional one.
    z = layer_center_z_pc(layer_index, edge_pc)
    return stars * sum(
        rate * remnant_distribution.vertical_factor(kind, z, shape.disk_scale_height_pc)
        for kind, _subtype, rate in _scattered_rates(min_mass_solar))


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
    return _row(*NUCLEUS_ADDRESS, "black-hole", (0, 0, 0), rng, subtype=NUCLEUS_SUBTYPE)


def hypervelocity_rows(extents, edge_pc, seed):
    """
    The galaxy's hypervelocity stars: `tuning.HYPERVELOCITY_STARS_PER_GALAXY`
    (log-uniform between its two ends) stars, each ejected from the central
    black hole in a random direction at `HYPERVELOCITY_STAR_SPEED_RANGE_KMS`
    and placed along that ray, moving outward. One that would land outside
    the outline is redrawn nearer in.

    Args:
        extents (list): `(layer_index, outer_ring_index)` per layer.
        edge_pc (float): The sector edge.
        seed (int): The scatter's seed.

    Yields:
        tuple: One row per star, in `PHENOMENON_SCATTER_COLUMNS` order.
    """
    outer_rings = dict(extents)
    if not outer_rings:
        return
    rng = draw.Stream(f"{seed}:phenomena:{HYPERVELOCITY_STREAM}")
    low, high = tuning.HYPERVELOCITY_STARS_PER_GALAXY
    wanted = int(round(log_uniform(low, high, rng=rng)))
    reach_pc = max(outer_ring + 1 for outer_ring in outer_rings.values()) * edge_pc
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
        if address[0] > outer_rings.get(address[1], -1):
            continue
        speed = log_uniform(*tuning.HYPERVELOCITY_STAR_SPEED_RANGE_KMS, rng=rng)
        yield _row(*address, HYPERVELOCITY_KIND, point, rng,
                   velocity=tuple(speed * component for component in direction))
        placed += 1


def special_rows(extents, edge_pc, seed, filled=()):
    """The nucleus and the hypervelocity stars, minus any in a cell
    already `filled`: the rows no layer draws."""
    rows = [nucleus_row(seed)] + list(hypervelocity_rows(extents, edge_pc, seed))
    return [row for row in rows if tuple(row[:3]) not in filled]


HYPERVELOCITY_ATTEMPTS = 50
"""int: Draws per wanted star before giving up on one (most land inside
the outline; a thin galaxy loses the ones that fly out of its disk)."""
