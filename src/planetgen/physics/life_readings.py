# planetgen/physics/life_readings.py

"""
Life readings (GEN.1210): how much, which way and how safe, for each kind of life on a world
==============================================================================================

A planet or moon with life stores four 0 to 9 scores, one for each kind of life (`life_microbial`, `life_flora`, `life_fauna`,
`life_sentient`, GEN.649). This module reads three things about each kind from the stored values alone, so the page and the API can
say, in three words, "Microbes: abundant, stable, secure." It is a pure function: no database, no random draw, no new column, no epoch.

**Abundance** is the score in the words of the ecologists' DAFOR scale: Absent (0), Rare (1 to 2), Occasional (3 to 4), Frequent (5 to
6), Abundant (7 to 8), Dominant (9). The bands are even steps of the 0 to 9 score.

**Trend** is where the kind is heading, from three things the world stores. The conditions: the living class is scored again with the
star as it will be `TREND_SPAN_GYR` later (the main sequence brightens, `physics/main_sequence.py`); the class's likelihood rising or
falling by more than `TREND_BAND` of itself is Growing or Declining. The stage timeline: a kind whose milestone (`KIND_MILESTONES`) was
reached within `SPREAD_GYR` of the world's age is still Spreading, which reads as Growing unless the conditions are declining. Otherwise
Stable. A world with no star to look at (or no isochrone download) is judged by the timeline alone.

**Security** is how long the living class's conditions last: the class is scored again at `STEP_GYR` steps as the star brightens, and the
first step where its score falls to 0 (or the star leaves the main sequence) is the time left. Dying out: the class already scores 0, or
the time left is `DYING_GYR` or less. Precarious: `PRECARIOUS_GYR` or less. Marginal: `MARGINAL_GYR` or less, or a score of `MARGINAL_SCORE` or less
now. Otherwise Secure. With no star, only the score now is read. The kinds of one world share its living class and so its conditions, so a
world's kinds differ in abundance and trend and, today, not in security (a kind's own limits are not stored).

What the scaling leaves out: the star's brightening warms the surface by a quarter of its fractional change (the temperature of a body that
radiates what it absorbs, `T ~ L^(1/4)`), with the greenhouse, the albedo and the ice cover held where they are, and the water held in
its stored state, as the oxygen rise and the air history already do. Every number is a stand-in, registered with its tier and tweak
level in docs/design/research/stand-in-values.md (the `LR` rows).

Nothing here changes a value. A world with no living class reads Absent for every kind, with no trend and no security.
"""

import math
from collections import namedtuple

from planetgen import tuning
from planetgen.datasets.fetch import DataMissing
from planetgen.physics import (atmosphere, constants, habitability_world, life_classes, life_supply, main_sequence, orbits,
                               stellar_evolution)

KINDS = (("microbial", "Microbes"), ("flora", "Flora"), ("fauna", "Fauna"), ("sentient", "Sentient life"))
"""tuple: Each kind of life, `(stored name, label)`, in the order a page lists them; the score is the body's `life_<name>`."""

ABUNDANCE_BANDS = ((0, "Absent"), (1, "Rare"), (3, "Occasional"), (5, "Frequent"), (7, "Abundant"), (9, "Dominant"))
"""tuple: `(lowest score, word)` for each band of the 0 to 9 score."""

TRENDS = ("Growing", "Stable", "Declining")
SECURITIES = ("Secure", "Marginal", "Precarious", "Dying out")

KIND_MILESTONES = {"microbial": "abiogenesis", "flora": "multicellularity", "fauna": "multicellularity",
                   "sentient": "technological_civilization"}
"""dict: The milestone of the schedule (`generation.evolution.MILESTONE_KEYS`) after which each kind of life spreads."""

TREND_SPAN_GYR = 0.5
"""float: How far ahead, Gyr, the conditions are looked at for the trend."""

TREND_BAND = 0.10
"""float: The share of its own value a class's likelihood must gain or lose over the span to be Growing or Declining."""

SPREAD_GYR = 0.3
"""float: How long, Gyr, after its milestone a kind of life is still spreading."""

STEP_GYR = 0.1
"""float: The step, Gyr, in which the time left is searched."""

DYING_GYR = 0.1
PRECARIOUS_GYR = 0.5
MARGINAL_GYR = 2.0
"""float: Time left, Gyr, within which the security is Dying out (one step), Precarious and Marginal."""

MARGINAL_SCORE = 2
"""int: A living class that scores this or less on its world's conditions today is Marginal at best."""

TEMPERATURE_EXPONENT = 0.25
"""float: How a surface temperature follows the star's luminosity: `T ~ L^0.25` for a body that radiates what it absorbs."""

OPEN_WATER = ("surface ocean", "hycean")
"""tuple: The stored states of water that the runaway greenhouse takes away."""

LUMINOSITY_KEYS = ("energy_flux_w_m2", "star_luminosity_w")
"""tuple: The stored values that follow the star's luminosity directly; `surface_temperature_k` follows its fourth root."""

Reading = namedtuple("Reading", "kind label score abundance trend security years_left_gyr")
"""One kind's reading: its stored score, the abundance word, the trend and security words (`None` for a kind that is Absent, or whose
conditions cannot be read) and the time left in Gyr (`None` for 'longer than the search' or not read)."""


def abundance_word(score):
    """The DAFOR word of a 0 to 9 score; `None` for a score that is not stored."""
    if score is None:
        return None
    word = ABUNDANCE_BANDS[0][1]
    for lowest, name in ABUNDANCE_BANDS:
        if int(score) >= lowest:
            word = name
    return word


def sentence(reading):
    """A reading as the page says it: "Microbes: abundant, stable, secure." ("Fauna: absent." for a kind that is Absent)."""
    words = [word for word in (reading.abundance, reading.trend, reading.security) if word]
    return f"{reading.label}: " + ", ".join(word.lower() for word in words) + "."


def star_of(star):
    """The `star` argument of `readings` from a star object (`mass` in kg, `age` in Gyr, `metallicity_fe_h`) or a mapping that already has
    `mass_sol`, `age_gyr` and `fe_h`; `None` when the star has no mass or age."""
    if star is None:
        return None
    if isinstance(star, dict):
        return star if star.get("mass_sol") and star.get("age_gyr") else None
    mass, age = getattr(star, "mass", None), getattr(star, "age", None)
    if not mass or not age:
        return None
    return {"mass_sol": mass / constants.SOLAR_MASS_TO_KG, "age_gyr": age, "fe_h": getattr(star, "metallicity_fe_h", None) or 0.0}


def _luminosity_ratio(star, years_gyr):
    """The star's luminosity `years_gyr` from now over its luminosity now, or `None` if the isochrones are not downloaded; `math.inf` once
    the star has left the main sequence."""
    mass, age, fe_h = star["mass_sol"], star["age_gyr"], star.get("fe_h") or 0.0
    try:
        if age + years_gyr >= stellar_evolution.main_sequence_lifetime_gy(mass, fe_h):
            return math.inf
        now = main_sequence.state(mass, age, fe_h).luminosity_sol
        later = main_sequence.state(mass, age + years_gyr, fe_h).luminosity_sol
    except DataMissing:
        return None
    return later / now if now > 0 else None


def _ocean_lost(values):
    """Whether the starlight in `values` is over the runaway-greenhouse limit of its host (`orbits.runaway_greenhouse_flux`, Kopparapu et al.
    2014) for a world with an open surface ocean: past it the sea boils off and the water classes that need it score nothing."""
    teff, flux, mass = values.get("star_teff_k"), values.get("energy_flux_w_m2"), values.get("mass_kg")
    if values.get("hydrosphere") not in OPEN_WATER or not teff or not flux or not mass:
        return False
    earth_light = constants.SOLAR_LUMINOSITY / (4.0 * math.pi * constants.AU_TO_M ** 2)
    insolation = flux / habitability_world.LIGHT_SHARE / earth_light
    return insolation > orbits.runaway_greenhouse_flux(teff, mass / constants.EARTH_MASS_TO_KG)


def _likelihood(values, row, ratio):
    """The likelihood of class `row` on the world with its star `ratio` times as luminous, or `None` if the world cannot be scored."""
    scaled = dict(values)
    if ratio != 1.0:
        if scaled.get("surface_temperature_k"):
            scaled["surface_temperature_k"] = scaled["surface_temperature_k"] * ratio ** TEMPERATURE_EXPONENT
        for key in LUMINOSITY_KEYS:
            if scaled.get(key):
                scaled[key] = scaled[key] * ratio
    if scaled.get("surface_temperature_k") is None or atmosphere.GAS_FIELD not in scaled:
        return None
    if ratio != 1.0 and _ocean_lost(scaled):
        scaled["hydrosphere"] = "vapour"
    for quantity, supplied in life_supply.supplies(scaled).items():
        scaled.setdefault(quantity, supplied)
    scaled.update(life_classes.derived_quantities(scaled))
    found = life_classes.quantities(scaled)
    return life_classes.score_class(row, found, habitability_world.world_from(scaled), scaled, ).likelihood


UNKNOWN = math.nan
"""float: The time left of a world whose star or whose isochrones are not known."""


def _is_known(years_left):
    return years_left is None or not math.isnan(years_left)


def _years_left(values, row, star):
    """Gyr until the class stops scoring on this world as its star brightens: `0.0` if it scores nothing now, `None` if it lasts through
    `MARGINAL_GYR`, `UNKNOWN` if the star or the isochrones are not known."""
    steps = int(round(MARGINAL_GYR / STEP_GYR))
    for n in range(1, steps + 1):
        years = round(n * STEP_GYR, 6)
        ratio = _luminosity_ratio(star, years)
        if ratio is None:
            return UNKNOWN
        if ratio == math.inf:
            return years
        likelihood = _likelihood(values, row, ratio)
        if likelihood is None:
            return UNKNOWN
        if life_classes.score_of(likelihood) == 0:
            return years
    return None


def _conditions(values, living, star):
    """`(score now, change over the span, years left)` for the living class on this world: `None` for what cannot be read (`UNKNOWN` for the
    time left)."""
    row = life_classes.BY_ID[living]
    now = _likelihood(values, row, 1.0)
    if now is None:
        return None, None, UNKNOWN
    score_now = life_classes.score_of(now)
    if star is None:
        return score_now, None, UNKNOWN
    ratio = _luminosity_ratio(star, TREND_SPAN_GYR)
    change = None
    if ratio is not None:
        later = 0.0 if ratio == math.inf else _likelihood(values, row, ratio)
        if later is not None and now > 0:
            change = (later - now) / now
        elif later is not None:
            change = 0.0 if later == 0 else 1.0
    return score_now, change, (0.0 if score_now == 0 else _years_left(values, row, star))


def _spreading(values, kind):
    """Whether the kind's milestone was reached within `SPREAD_GYR` of the world's age on its stored schedule."""
    key = KIND_MILESTONES[kind]
    reached = values.get(f"life_age_{key}_gyr")
    age = values.get("life_system_age_gyr")
    return reached is not None and age is not None and 0 <= age - reached < SPREAD_GYR


def _trend(values, kind, change):
    if change is not None and change < -TREND_BAND:
        return "Declining"
    if _spreading(values, kind) or (change is not None and change > TREND_BAND):
        return "Growing"
    return "Stable"


def _security(score_now, years_left):
    if score_now is None:
        return None
    known = years_left is not None and _is_known(years_left)
    if score_now == 0 or (known and years_left <= DYING_GYR):
        return "Dying out"
    if known and years_left <= PRECARIOUS_GYR:
        return "Precarious"
    if (known and years_left <= MARGINAL_GYR) or score_now <= MARGINAL_SCORE:
        return "Marginal"
    return "Secure"


def readings(values, star=None):
    """
    The four `Reading`s of a world (kinds in `KINDS` order), from its stored `values` (a column-named mapping that holds the air under
    `atmosphere.GAS_FIELD`, as `db.store.with_life_inputs` and `life_classes.values_from_body` make it) and, if known, its star
    (`star_of`). `[]` for a gas giant and for a body whose kind scores are not stored.

    Args:
        values (Mapping): The body's stored values: `life_microbial`, `life_flora`, `life_fauna`, `life_sentient`, `life_class`,
            `life_class_variant`, `life_stage`, the milestone ages `life_age_<key>_gyr`, `life_system_age_gyr` and what
            `life_classes.score_classes` reads.
        star (dict, optional): `mass_sol`, `age_gyr` and `fe_h` of the host star; without it the trend comes from the timeline and the
            security from today's score.

    Returns:
        list: One `Reading` for each kind.
    """
    scores = {kind: values.get(f"life_{kind}") for kind, _label in KINDS}
    if values.get("body_type") == "g" or all(score is None for score in scores.values()):
        return []
    living = life_classes.living_class_id(values.get("life_class"), values.get("life_class_variant"), values.get("life_stage"))
    outlook = (None, None, UNKNOWN)
    if living in life_classes.BY_ID and any(scores.values()):
        outlook = _conditions(values, living, star_of(star))
    score_now, change, years_left = outlook
    found = []
    for kind, label in KINDS:
        score = scores[kind]
        if living not in life_classes.BY_ID or not score:
            found.append(Reading(kind, label, score, None if score is None else "Absent", None, None, None))
            continue
        left = years_left if years_left is not None and _is_known(years_left) else None
        found.append(Reading(kind, label, int(score), abundance_word(score), _trend(values, kind, change), _security(score_now, years_left), left))
    return found
