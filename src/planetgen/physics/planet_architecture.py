# planetgen/physics/planet_architecture.py

"""
Planet architecture (GEN.728, GEN.729, GEN.730, GEN.828, GEN.829, GEN.920, docs/design/composition/planet-architecture.md)
==========================================================================================================================

How many planets a star has, where they sit and how heavy they are, from the star's disc alone. A pure law over the star and its disc;
it reads no class, writes none and draws only from the stream it is given.

1. **The disc.** A disc's solids scatter about the median of its star's mass (`DISC_SOLIDS_SCATTER_DEX`, Andrews et al. 2013). A light
   disc builds nothing the size of a planet; a heavy one builds a chain. No coin flip decides whether a star has planets.
2. **The chain.** From an inner edge (a period near 5 days, Mulders et al. 2015) outwards oligarchs form one feeding zone after another, each
   at the isolation mass (`formation.isolation_mass_kg`) and `separation` mutual Hill radii from the last (`OLIGARCH_SEPARATION_MEDIAN`,
   drawn per system). The outer disc has built less (`PLANET_OUTER_GROWTH_AU`).
3. **The thinning.** Pu and Wu 2015 (arXiv:1502.05449) integrated tight chains and fitted the chance that one survives: a pair of spacing
   `K` mutual Hill radii lasts `tau` orbits of its inner planet with probability `1 - exp(-x^2)`, `x = 0.5 K - 0.6 - 0.35 log10 tau`. Each
   pair of the chain gets an instability time drawn from that law; the earliest to go collides (the two merge, keeping their angular
   momentum), the new neighbours get new times, and the chain is thinned until no pair fails inside `SCULPT_AGE_YR`. What is left is
   the observed pattern: spacings clustered near 12 to 20 mutual Hill radii, fewer and heavier planets than the oligarchs.
4. **The gas.** A core of `GIANT_CORE_CRITICAL_MASS_EARTH` or more beyond the snow line may find gas (`planet_mass.gas_outcome`).
   A giant may then migrate in to a hot Jupiter orbit (`GIANT_MIGRATION_SHARE`, rising with [Fe/H] faster than giants do); the chain inside
   its new orbit is gone, which is why hot Jupiters rarely have companions (Zhu and Dong 2021).

The number of planets is the length of the list this returns. It replaces the uniform count between 0 and a ceiling, the uniform ratio
between 1.3 and 1.8 and the constant 0.4 survival share. Every constant is a row of the stand-in registry (section 6.36 of the gap
filling plan).

Use:

    architecture = build(star_mass_kg, snow_line_au, metal_ratio, fe_h, stream)
    for seed in architecture.planets: seed.distance_au, seed.mass_kg, seed.runaway, seed.merged
"""

import collections
import math

from planetgen import tuning
from planetgen.physics import constants, formation, metallicity, planet_mass, unit_registry as units
from planetgen.util import draw

_kg_to_msun = units.converter("kg", "Msun")

_EARTH_KG = constants.EARTH_MASS_TO_KG
_G = constants.G
_DAY_S = 86400.0
_SOLAR_BLACKBODY_1AU_K = 278.3      # a blackbody at 1 AU from the Sun (Lodders 2003 Table 10)
_YEAR_S = 365.25 * _DAY_S

SURVIVAL_K_OFFSET = 0.6
SURVIVAL_K_SLOPE = 0.5
SURVIVAL_LOG_TIME_SLOPE = 0.35
"""float: The three coefficients of Pu and Wu's survival law, `x = SLOPE K - OFFSET - LOG_TIME_SLOPE log10(tau)` (arXiv:1502.05449 eq. 11, READ).
They fitted 7-planet systems of 3 to 9 Earth masses round a Sun-like star, so the law is used well outside its range here (registry row)."""

SURVIVAL_K_MINIMUM = 1.2
SURVIVAL_K_MINIMUM_SLOPE = 0.7
"""float: Below `SURVIVAL_K_MINIMUM + SURVIVAL_K_MINIMUM_SLOPE log10(tau)` no system survives (arXiv:1502.05449 eq. 11)."""

MAX_LOG_TAU = 40.0
"""float: An instability time beyond this many decades of orbits is never (the draw is capped so the power cannot overflow)."""

Seed = collections.namedtuple("Seed", ("distance_au", "mass_kg", "core_kg", "runaway", "merged", "migrated", "formation_au"))
"""One planet of the architecture: where it ends (`distance_au`), its mass, the solid core it grew from, whether it is a runaway giant,
how many oligarchs merged into it, whether it migrated in to a hot Jupiter orbit and where it formed (`formation_au`: the origin the
composition reads, the snow line it was born past)."""

Architecture = collections.namedtuple("Architecture", ("planets", "disc_scatter", "separation", "oligarchs"))
"""The planets in order of distance, the factor the disc's solids were scattered by, the spacing the oligarchs formed at (mutual Hill
radii) and the number of oligarchs before the thinning."""


def mutual_hill_radius_au(a1_au, m1_kg, a2_au, m2_kg, star_mass_kg):
    """The mutual Hill radius of two planets (Chambers et al. 1996), in AU."""
    return 0.5 * (a1_au + a2_au) * ((m1_kg + m2_kg) / (3.0 * star_mass_kg)) ** (1.0 / 3.0)


def spacing_k(a1_au, m1_kg, a2_au, m2_kg, star_mass_kg):
    """The spacing of two planets in mutual Hill radii."""
    return (a2_au - a1_au) / mutual_hill_radius_au(a1_au, m1_kg, a2_au, m2_kg, star_mass_kg)


def orbital_period_s(a_au, star_mass_kg):
    """The period of a circular orbit of `a_au` round `star_mass_kg`, in seconds."""
    a_m = a_au * constants.AU_TO_M
    return 2.0 * math.pi * math.sqrt(a_m ** 3 / (_G * star_mass_kg))


def survival_probability(k, log_tau):
    """Pu and Wu's chance that a pair of spacing `k` mutual Hill radii is still stable after `10 ** log_tau` orbits of its inner planet."""
    x = SURVIVAL_K_SLOPE * k - SURVIVAL_K_OFFSET - SURVIVAL_LOG_TIME_SLOPE * log_tau
    if k < SURVIVAL_K_MINIMUM + SURVIVAL_K_MINIMUM_SLOPE * log_tau:
        return 0.0
    return 1.0 - math.exp(-x * x)


def instability_log_tau(k, u):
    """
    The decades of orbits a pair of spacing `k` lasts, from `u` a uniform draw in (0, 1): the inverse of `survival_probability` in
    time. A pair that is unstable at once (the law gives it no orbit) returns 0; none lasts beyond `MAX_LOG_TAU`.
    """
    x = math.sqrt(-math.log(u))
    log_tau = (SURVIVAL_K_SLOPE * k - SURVIVAL_K_OFFSET - x) / SURVIVAL_LOG_TIME_SLOPE
    return min(max(log_tau, 0.0), MAX_LOG_TAU)


def period_to_au(period_d, star_mass_kg):
    """The semi-major axis of the orbit of `period_d` days round `star_mass_kg`, in AU."""
    a_m = (_G * star_mass_kg * (period_d * _DAY_S) ** 2 / (4.0 * math.pi ** 2)) ** (1.0 / 3.0)
    return a_m / constants.AU_TO_M


def inner_system_fraction(star_mass_kg):
    """
    The share of stars of this mass with a Kepler-like inner system (a planet of an Earth radius or more with a period under 400 days):
    log-linear in the star's mass between `INNER_SYSTEM_FRACTION_MASSES` and `INNER_SYSTEM_FRACTION_SHARES`, held at the ends. Zhu et al. 2018
    find 30 percent of Sun-like stars and Yang et al. 2020 about 75 percent of stars under 5000 K and 35 percent over 6500 K.
    """
    mass = _kg_to_msun(star_mass_kg)
    anchors = list(zip(tuning.INNER_SYSTEM_FRACTION_MASSES, tuning.INNER_SYSTEM_FRACTION_SHARES))
    if mass <= anchors[0][0]:
        return anchors[0][1]
    for (m0, f0), (m1, f1) in zip(anchors, anchors[1:]):
        if mass <= m1:
            w = (math.log(mass) - math.log(m0)) / (math.log(m1) - math.log(m0))
            return f0 + w * (f1 - f0)
    return anchors[-1][1]


def _poisson(mean, source):
    """A Poisson draw by inversion (one uniform)."""
    u = source.random()
    k, term = 0, math.exp(-mean)
    total = term
    while u > total and k < 100:
        k += 1
        term *= mean / k
        total += term
    return k


def inner_chain(star_mass_kg, scatter, fe_h, stream=None, force=False, extra_sigmas=0.0, floor_au=0.0):
    """
    The oligarchs of a Kepler-like inner system as `[distance_au, mass_kg]` lists, or `[]` for a star without one (module docstring).
    Draws, in this order: whether the star has one (`inner_system_fraction`), the number of planets before the thinning
    (`1 + Poisson(INNER_INITIAL_COUNT_EXTRA * exp(INNER_COUNT_SIGMA_LN * N(0, 1)))`), the system's mass scale, its spacing, the period of the first planet, then per planet its
    mass about the scale and its spacing from the last. The first planet is no closer than `floor_au`. `force` gives the star an inner system whatever the share (it then draws no
    uniform for it), and `extra_sigmas` raises the mean of the starting count by that many widths of its spread (the most a star can
    hold, `build`'s `maximum`). Neighbours are alike in mass (`INNER_NEIGHBOUR_SIGMA_DEX`, peas in a pod,
    Weiss et al. 2018) and the outer one is on average a little heavier (`INNER_OUTWARD_MASS_DEX_PER_PLANET`).
    """
    source = stream or draw
    if not force and source.random() >= inner_system_fraction(star_mass_kg):
        return []
    count = 1 + _poisson(tuning.INNER_INITIAL_COUNT_EXTRA
                         * math.exp(tuning.INNER_COUNT_SIGMA_LN * (source.gauss(0.0, 1.0) + extra_sigmas)), source)
    solar_masses = _kg_to_msun(star_mass_kg)
    scale_kg = (tuning.INNER_PLANET_MASS_EARTH * solar_masses ** tuning.INNER_PLANET_MASS_STELLAR_EXPONENT
                * scatter ** tuning.INNER_PLANET_MASS_DISC_EXPONENT * 10.0 ** (tuning.INNER_PLANET_MASS_METALLICITY_EXPONENT * fe_h)) * _EARTH_KG
    separation = planet_mass.draw_separation(tuning.INNER_SEPARATION_MEDIAN, source)
    distance = max(period_to_au(tuning.OLIGARCH_INNER_EDGE_PERIOD_DAYS
                                * 10.0 ** (tuning.OLIGARCH_INNER_EDGE_SIGMA_DEX * source.gauss(0.0, 1.0)), star_mass_kg), floor_au)
    chain = []
    for i in range(count):
        offset = (i - (count - 1) / 2.0) * tuning.INNER_OUTWARD_MASS_DEX_PER_PLANET
        mass = scale_kg * 10.0 ** (offset + tuning.INNER_NEIGHBOUR_SIGMA_DEX * source.gauss(0.0, 1.0))
        chain.append([distance, mass])
        next_mass = scale_kg * 10.0 ** (offset + tuning.INNER_OUTWARD_MASS_DEX_PER_PLANET)
        spacing = separation * math.exp(tuning.INNER_SPACING_SIGMA_LN * source.gauss(0.0, 1.0))
        kappa = min(spacing * ((mass + next_mass) / (3.0 * star_mass_kg)) ** (1.0 / 3.0), 1.8)
        distance *= (1.0 + kappa / 2.0) / (1.0 - kappa / 2.0)
    return chain


def oligarch_chain(star_mass_kg, snow_line_au, metal_ratio, scatter, separation, inner_au, outer_au):
    """
    The oligarchs an outer disc builds between `inner_au` and `outer_au`, as `[distance_au, mass_kg]` lists in order of distance: each at
    the isolation mass of its feeding zone `separation` mutual Hill radii wide, the next one a spacing further out. The outer disc has
    built less (the taper `PLANET_OUTER_GROWTH_AU` and `PLANET_OUTER_GROWTH_EXPONENT`). Bodies lighter than `PLANET_MIN_MASS_EARTH` are kept in the chain (they
    take part in the thinning) and leave it afterwards (`build`).
    """
    chain = []
    distance = inner_au
    while distance < outer_au and len(chain) < tuning.ABSOLUTE_MAX_SYSTEM_OBJECTS:
        mass = planet_mass.oligarch_core_kg(distance, star_mass_kg, snow_line_au, metal_ratio, scatter, separation)
        chain.append([distance, mass])
        kappa = min(separation * ((2.0 * mass) / (3.0 * star_mass_kg)) ** (1.0 / 3.0), 1.8)
        distance *= (1.0 + kappa / 2.0) / (1.0 - kappa / 2.0)
    return chain


def _merge(first, second):
    """Two colliding bodies `[distance_au, mass_kg, merged]` as one: mass kept to `COLLISION_MASS_RETENTION`, angular momentum conserved."""
    m1, m2 = first[1], second[1]
    mass = (m1 + m2) * tuning.COLLISION_MASS_RETENTION
    sqrt_a = (m1 * math.sqrt(first[0]) + m2 * math.sqrt(second[0])) / (m1 + m2)
    return [sqrt_a * sqrt_a, mass, first[2] + second[2]]


def sculpt(chain, star_mass_kg, age_yr, stream=None):
    """
    Thins `chain` (`[distance_au, mass_kg]` lists) by dynamical instability to `age_yr` (module docstring, step 3): the pair with the
    earliest instability time collides and merges, its two new neighbouring pairs draw new times, until every pair outlasts the age.

    Returns:
        list: `[distance_au, mass_kg, merged]` lists in order of distance, `merged` the number of oligarchs in each.
    """
    source = stream or draw
    bodies = [[a, m, 1] for a, m in chain]
    if len(bodies) < 2:
        return bodies

    def pair_time_s(i, u):
        a1, m1 = bodies[i][0], bodies[i][1]
        a2, m2 = bodies[i + 1][0], bodies[i + 1][1]
        log_tau = instability_log_tau(spacing_k(a1, m1, a2, m2, star_mass_kg), u)
        return orbital_period_s(a1, star_mass_kg) * 10.0 ** log_tau

    uniforms = [1.0 - source.random() for _ in range(len(bodies) - 1)]
    times = [pair_time_s(i, uniforms[i]) for i in range(len(bodies) - 1)]
    age_s = age_yr * _YEAR_S
    while times:
        i = min(range(len(times)), key=times.__getitem__)
        if times[i] > age_s:
            break
        bodies[i:i + 2] = [_merge(bodies[i], bodies[i + 1])]
        del times[i], uniforms[i]
        if i > 0:
            uniforms[i - 1] = 1.0 - source.random()
            times[i - 1] = pair_time_s(i - 1, uniforms[i - 1])
        if i < len(times):
            uniforms[i] = 1.0 - source.random()
            times[i] = pair_time_s(i, uniforms[i])
    return bodies


def migration_probability(fe_h):
    """The chance a giant migrates in to a hot Jupiter orbit: `GIANT_MIGRATION_SHARE` at solar metallicity, rising with [Fe/H] by the extra
    exponent hot Jupiters have over giants (Petigura et al. 2018 against Johnson et al. 2010), never above 1."""
    return min(1.0, tuning.GIANT_MIGRATION_SHARE * 10.0 ** (tuning.GIANT_MIGRATION_METALLICITY_EXTRA_EXPONENT * fe_h))


def _clear_around_giants(seeds, star_mass_kg):
    """
    `seeds` without the planets a giant scattered: those inside `GIANT_CLEARING_SPACING` mutual Hill radii of a giant that moved in, and
    every planet inside a hot Jupiter's orbit (a migrating giant sweeps the inner chain: hot Jupiters rarely have companions, Zhu and Dong
    2021). A giant that did not move clears nothing.
    """
    moved = [s for s in seeds if s.runaway and s.migrated]
    hot = [s for s in moved if s.distance_au < tuning.GIANT_HOT_ORBIT_AU]
    kept = []
    for s in seeds:
        if s in moved:
            kept.append(s)
            continue
        if any(s.distance_au < h.distance_au for h in hot):
            continue
        if any(abs(spacing_k(min(s.distance_au, g.distance_au), s.mass_kg, max(s.distance_au, g.distance_au), g.mass_kg, star_mass_kg))
               < tuning.GIANT_CLEARING_SPACING for g in moved):
            continue
        kept.append(s)
    return _scatter_crowded_giants(kept, star_mass_kg)


def _scatter_crowded_giants(seeds, star_mass_kg):
    """
    `seeds` (in order of distance) without the lighter of any neighbouring pair closer than `GIANT_CLEARING_SPACING` mutual Hill radii that
    holds a giant: each giant moves in on its own, so two can land on one another, and the heavier throws the lighter out (the same spacing
    inside which a moved giant clears the planets around it). Pairs the thinning left apart are not touched.
    """
    kept = list(seeds)
    i = 0
    while i < len(kept) - 1:
        inner, outer = kept[i], kept[i + 1]
        k = spacing_k(inner.distance_au, inner.mass_kg, outer.distance_au, outer.mass_kg, star_mass_kg)
        if (inner.runaway or outer.runaway) and k < tuning.GIANT_CLEARING_SPACING:
            kept.pop(i if inner.mass_kg < outer.mass_kg else i + 1)
            i = max(i - 1, 0)
        else:
            i += 1
    return kept


def sublimation_au(luminosity_w):
    """The distance in AU inside which a blackbody would be hotter than `DUST_SUBLIMATION_K`, from the star's luminosity in watts (0 for a
    dark star)."""
    return (_SOLAR_BLACKBODY_1AU_K / tuning.DUST_SUBLIMATION_K) ** 2 * math.sqrt(max(luminosity_w, 0.0) / constants.SOLAR_LUMINOSITY)


def first_distance_au(star_mass_kg, stream=None, floor_au=0.0):
    """Where a chain that has to start from nothing begins: the inner edge of the disc, a period near `OLIGARCH_INNER_EDGE_PERIOD_DAYS`
    (log-normal in dex by `OLIGARCH_INNER_EDGE_SIGMA_DEX`, Mulders et al. 2015), no closer than `floor_au`. One normal from the stream."""
    source = stream or draw
    return max(period_to_au(tuning.OLIGARCH_INNER_EDGE_PERIOD_DAYS * 10.0 ** (tuning.OLIGARCH_INNER_EDGE_SIGMA_DEX * source.gauss(0.0, 1.0)),
                            star_mass_kg), floor_au)


def next_distance_au(distance_au, mass_kg, star_mass_kg, separation=None):
    """The distance of the next oligarch out from one of `mass_kg` at `distance_au`: `separation` mutual Hill radii further (the
    architecture's own median when not given), as the chains are built."""
    b = tuning.OLIGARCH_SEPARATION_MEDIAN if separation is None else separation
    kappa = min(b * ((2.0 * mass_kg) / (3.0 * star_mass_kg)) ** (1.0 / 3.0), 1.8)
    return distance_au * (1.0 + kappa / 2.0) / (1.0 - kappa / 2.0)


def extra_seed(previous, architecture, star_mass_kg, snow_line_au, metal_ratio, fe_h, floor_au=0.0, stream=None):
    """
    One more planet in an architecture's disc, for a system asked to hold more bodies than the thinning left: the next oligarch out from
    `previous` (a `Seed`, or `None` to start at the disc's inner edge, `first_distance_au`), at the isolation mass of its feeding zone and
    through the gas outcome, in the disc the architecture was built in (`architecture.disc_scatter`, `architecture.separation`).
    """
    source = stream or draw
    if previous is None:
        distance = first_distance_au(star_mass_kg, source, floor_au)
    else:
        distance = next_distance_au(previous.distance_au, previous.mass_kg, star_mass_kg, architecture.separation)
    core = planet_mass.oligarch_core_kg(distance, star_mass_kg, snow_line_au, metal_ratio, architecture.disc_scatter, architecture.separation)
    outcome = planet_mass.gas_outcome(core, distance, star_mass_kg, snow_line_au, fe_h, source)
    return Seed(distance, outcome.mass_kg, outcome.core_kg, outcome.runaway, 1, False, distance)


class Plan:
    """
    The planets a system is to hold, slot by slot: the seeds of an `Architecture` (cut to `limit` when the system is asked for fewer), and
    past the last one the next oligarchs out (`extra_seed`) for a system asked for more. `seed(i)` is the planet of slot `i`.
    """

    def __init__(self, architecture, star_mass_kg, snow_line_au, metal_ratio, fe_h, limit=None, floor_au=0.0):
        self.architecture = architecture
        self.seeds = list(architecture.planets if limit is None else architecture.planets[:limit])
        self._context = (star_mass_kg, snow_line_au, metal_ratio, fe_h, floor_au)

    def __len__(self):
        return len(self.seeds)

    def seed(self, i):
        """The seed of slot `i`, adding oligarchs beyond the last until there is one."""
        while len(self.seeds) <= i:
            self.seeds.append(extra_seed(self.seeds[-1] if self.seeds else None, self.architecture, *self._context))
        return self.seeds[i]


def build(star_mass_kg, snow_line_au, metal_ratio, fe_h, outer_au, stream=None, inner=None, maximum=False, floor_au=0.0):
    """
    The planets of a star (module docstring). Draws, in this order: the disc scatter, the inner system, the outer spacing, the thinning,
    then, for each core past the critical mass beyond the snow line, the gas and migration draws.

    Args:
        star_mass_kg (float): The star's mass, in kilograms.
        snow_line_au (float): The star's own snow line, in AU.
        metal_ratio (float): The star's Z / Z0 (`bulk_draw.metal_ratio_to_solar`).
        fe_h (float): The star's [Fe/H], dex.
        outer_au (float): The disc's outer edge, in AU (the oligarchs stop there).
        stream (draw.Stream or None): The star's own stream; `None` draws from the stream bound now.
        inner (bool or None): `True` gives the star an inner system, `False` none, `None` the observed share (`inner_system_fraction`).
        maximum (bool): The most the disc holds: an inner system, and a starting count two widths of its spread above the mean.
        floor_au (float): The closest a body forms or migrates to, in AU (`sublimation_au`).

    Returns:
        Architecture: The planets by distance, and what built them.
    """
    source = stream or draw
    scatter = planet_mass.draw_disc_scatter(source)
    if inner is False:
        chain_inner = []
    else:
        chain_inner = inner_chain(star_mass_kg, scatter, fe_h, source, force=bool(inner) or maximum, extra_sigmas=2.0 if maximum else 0.0,
                                  floor_au=floor_au)
    separation = planet_mass.draw_separation(tuning.OLIGARCH_SEPARATION_MEDIAN, source)
    start_au = max(period_to_au(tuning.INNER_SYSTEM_PERIOD_LIMIT_DAYS * tuning.OUTER_CHAIN_START_RATIO, star_mass_kg), floor_au)
    if chain_inner:
        start_au = max(start_au, chain_inner[-1][0] * tuning.OUTER_CHAIN_START_RATIO)
    chain = chain_inner + oligarch_chain(star_mass_kg, snow_line_au, metal_ratio, scatter, separation, start_au, outer_au)
    thinned = [body for body in sculpt(chain, star_mass_kg, tuning.SCULPT_AGE_YR, source)
               if body[1] >= (tuning.PLANET_MIN_MASS_EARTH if body[0] < snow_line_au else tuning.PLANET_MIN_MASS_OUTER_EARTH) * _EARTH_KG]
    seeds = []
    hot_au = None
    for distance, core, merged in thinned:
        formed_au = distance
        outcome = planet_mass.gas_outcome(core, distance, star_mass_kg, snow_line_au, fe_h, source)
        migrated = False
        if outcome.runaway:
            if hot_au is None and source.random() < migration_probability(fe_h):
                low, high = tuning.GIANT_MIGRATION_PERIOD_RANGE_DAYS
                distance = max(period_to_au(math.exp(source.uniform(math.log(low), math.log(high))), star_mass_kg), floor_au)
                hot_au, migrated = distance, True
            else:
                distance *= math.exp(-abs(source.gauss(0.0, tuning.GIANT_INWARD_MIGRATION_SIGMA_LN)))
        seeds.append(Seed(distance, outcome.mass_kg, outcome.core_kg, outcome.runaway, merged, migrated, formed_au))
    seeds.sort(key=lambda s: s.distance_au)
    seeds = _clear_around_giants(seeds, star_mass_kg)
    return Architecture(tuple(seeds), scatter, separation, len(chain))
