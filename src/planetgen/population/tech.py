# planetgen/population/tech.py

"""
Tech levels of technological species (POP.7 to POP.9)
======================================================

Six domain indices, each an integer 0 to 7 (Energy, Materials, Information, Medical, Propulsion, Defense),
drawn from a species' age, and the tech level as their weighted sum. See "Tech levels" in
`docs/design/population-and-politics.md` for the model and its evidence, and Boss's Technological Assessment
for the six scales and the weights.

The functions take no database. `draw_indices` draws from a stream of its own (callers use
`draw.Stream(f"tech:{planet_id}")`), never the stream `scan_life_worlds` draws a species' traits and age
from, so the stored ages never change.
"""

import math

from planetgen import tuning

DOMAINS = tuning.TECH_DOMAINS
MAX_INDEX = 7


def _poisson(mean, rng):
    """A Poisson count of mean `mean` (Knuth's product)."""
    limit, product, count = math.exp(-mean), rng.random(), 0
    while product > limit:
        count += 1
        product *= rng.random()
    return count


def _student_t4(rng):
    """A Student t draw with 4 degrees of freedom: a normal over the root of a chi-square / 4."""
    return rng.gauss() / math.sqrt(-0.5 * math.log((1.0 - rng.random()) * (1.0 - rng.random())))


def _soft_minimum(a, b, sharpness):
    """The soft minimum of `a` and `b`: close to the smaller when they are far apart, never above it."""
    low = min(a, b)
    return low - math.log(math.exp(-sharpness * (a - low)) + math.exp(-sharpness * (b - low))) / sharpness


def effective_age(age_years, rng):
    """
    The age a species' technology has had to grow: `age_years` less what collapses took and held to a plateau.

    Collapses are a Poisson process of `TECH_COLLAPSE_RATE_PER_DECADE` per decade of log10 age past
    `CIVILIZATION_MIN_AGE_YEARS`, each at a log-uniform age; a collapse multiplies the age so far by a draw of
    `TECH_COLLAPSE_KEEP`, and growth goes on. With `TECH_STAGNATION_CHANCE` progress also freezes at a plateau.
    """
    low = tuning.CIVILIZATION_MIN_AGE_YEARS
    decades = max(0.0, math.log10(max(age_years, low) / low))
    events = sorted(low * 10.0 ** rng.uniform(0.0, decades)
                    for _ in range(_poisson(tuning.TECH_COLLAPSE_RATE_PER_DECADE * decades, rng)))
    effective, previous = 0.0, 0.0
    for when in events:
        effective = (effective + when - previous) * rng.uniform(*tuning.TECH_COLLAPSE_KEEP)
        previous = when
    effective += max(age_years, low) - previous
    if rng.random() < tuning.TECH_STAGNATION_CHANCE:
        effective = min(effective, 10.0 ** rng.uniform(*tuning.TECH_STAGNATION_LOG10_YEARS))
    return effective


def draw_indices(age_years, rng):
    """
    The six domain indices of a species `age_years` past its industrial transition.

    Each is a logistic of the log10 effective age, shifted by the domain's offset, the species' shared pace
    and the domain's own t(4) noise, softly held under the species' ceiling, clamped to its prerequisite
    domains, and, for a lopsided species, one raised; then rounded to an integer 0 to 7.

    Args:
        age_years (float): `civilization_age_years`, years since the industrial transition.
        rng (draw.Stream): The species' own tech stream.

    Returns:
        dict: Each domain name (`tuning.TECH_DOMAINS`) to its integer index.
    """
    log_age = math.log10(max(effective_age(age_years, rng), 1.0))
    pace = rng.gauss(0.0, tuning.TECH_PACE_SIGMA)
    species_ceiling = min(7.0, max(4.0, rng.gauss(tuning.TECH_CEILING_MEAN, tuning.TECH_CEILING_SD)))
    values = {}
    for domain, offset in zip(DOMAINS, tuning.TECH_DOMAIN_OFFSETS):
        shift = offset + pace + _student_t4(rng) * tuning.TECH_DOMAIN_SIGMA
        raw = MAX_INDEX / (1.0 + math.exp(-tuning.TECH_CURVE_STEEPNESS * (log_age - tuning.TECH_CURVE_CENTRE - shift)))
        ceiling = min(7.0, species_ceiling + rng.gauss(0.0, tuning.TECH_CEILING_DOMAIN_SD))
        values[domain] = _soft_minimum(raw, ceiling, tuning.TECH_CEILING_SOFTNESS)
    for domain, needs, margin in tuning.TECH_PREREQUISITES:
        values[domain] = min(values[domain], min(values[need] for need in needs) + margin)
    if rng.random() < tuning.TECH_SPIKE_CHANCE:
        domain = DOMAINS[rng.randrange(len(DOMAINS))]
        values[domain] = min(float(MAX_INDEX), values[domain] + rng.uniform(*tuning.TECH_SPIKE_SIZE))
    return {domain: min(MAX_INDEX, max(0, int(math.floor(value + 0.5)))) for domain, value in values.items()}


def tech_level(indices):
    """The tech level, 0 to 7 to two places: the weighted sum of the six indices (`tuning.TECH_WEIGHTS`)."""
    return round(sum(weight * indices[domain] for domain, weight in zip(DOMAINS, tuning.TECH_WEIGHTS)), 2)


def band_name(level):
    """The tech band a tech level is in (`tuning.TECH_BANDS`)."""
    name = tuning.TECH_BANDS[0][1]
    for start, band in tuning.TECH_BANDS:
        if level >= start:
            name = band
    return name


def limiting_domain(indices):
    """The lowest of Energy, Materials and Information, the foundation domains; a tie goes in that order."""
    return min(DOMAINS[:3], key=lambda domain: indices[domain])


def is_lopsided(indices):
    """Whether the highest index is at least 3 above the lowest."""
    return max(indices.values()) - min(indices.values()) >= 3
