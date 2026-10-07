# planetgen/physics/stellar_evolution.py

"""
Star Population Model
=====================

Draws a random star as physics rather than from a spectral-letter table:

1. a mass from the Kroupa (2001) initial mass function (`sample_imf_mass_sol`),
2. an age from the star-formation history (`sample_star_age_gy`),
3. the star's present state from how far that age is into its own
   main-sequence lifetime (`evolve_star`): main sequence (V), subgiant
   (IV), giant (III/II) or supergiant (IB/IAB/IA/0), then a white dwarf
   (VII) or, for a massive star, a neutron star or black hole.

The spectral letter, subclass and Yerkes class then follow from the
resulting temperature and luminosity (`spectral_letter_and_subclass`), so
a hot star is always luminous for its size, a supergiant is always
massive, and white dwarfs and red giants turn up at their real rates.

Every constant lives in `program_constants`' "Star population model"
section.
"""

import math
import random

from planetgen.physics import constants
from planetgen import tuning


YERKES_CLASS_NAMES = {
    "0": "Hypergiant", "IA+": "Luminous Supergiant", "IA": "Supergiant",
    "IAB": "Intermediate-size Luminous Supergiant", "IB": "Less Luminous Supergiant",
    "II": "Bright Giant", "III": "Giant", "IV": "Subgiant", "V": "Main Sequence",
    "VI": "Subdwarf", "VII": "White Dwarf", "D": "White Dwarf",  # D is an alias for VII
}
"""Display name of each Yerkes luminosity class, used in `Star.type`."""


def main_sequence_lifetime_gy(mass_sol):
    """t_MS = 10 Gy * M^-2.5, the same law the rest of the generator uses."""
    return tuning.SOLAR_MS_LIFESPAN_GY * mass_sol ** tuning.MS_LIFESPAN_MASS_EXPONENT


def sample_imf_mass_sol(min_mass_sol=None, max_mass_sol=None, rng=random):
    """
    Draws an initial mass (Msun) from the Kroupa broken power law
    (`IMF_BREAKS_SOL`/`IMF_SLOPES`), optionally truncated to
    `[min_mass_sol, max_mass_sol]`.

    Each segment's number weight is its power-law integral (with the
    segments joined continuously at the breaks); a segment is picked by
    weight and the mass drawn from it by inverse transform.
    """
    breaks = tuning.IMF_BREAKS_SOL
    slopes = tuning.IMF_SLOPES
    lo = breaks[0] if min_mass_sol is None else max(min_mass_sol, breaks[0])
    hi = breaks[-1] if max_mass_sol is None else min(max_mass_sol, breaks[-1])
    if not lo < hi:
        raise ValueError(f"empty IMF range [{lo}, {hi}] Msun")

    segments = []
    norm = 1.0  # continuity factor k_i, with k_0 = 1
    for i, alpha in enumerate(slopes):
        a, b = breaks[i], breaks[i + 1]
        if i > 0:
            norm *= breaks[i] ** (alpha - slopes[i - 1])
        a, b = max(a, lo), min(b, hi)
        if a < b:
            segments.append((a, b, alpha, norm * _power_law_integral(a, b, alpha)))

    total = sum(weight for *_, weight in segments)
    pick = rng.random() * total
    for a, b, alpha, weight in segments:
        if pick < weight:
            break
        pick -= weight
    return _sample_power_law(a, b, alpha, rng)


def _power_law_integral(a, b, alpha):
    if alpha == 1.0:
        return math.log(b / a)
    return (b ** (1 - alpha) - a ** (1 - alpha)) / (1 - alpha)


def _sample_power_law(a, b, alpha, rng):
    u = rng.random()
    if alpha == 1.0:
        return a * (b / a) ** u
    p = 1 - alpha
    return (a ** p + u * (b ** p - a ** p)) ** (1 / p)


def age_window_gy(low, high, age_bias=None):
    """
    Narrows `[low, high]` for `SystemConfig.AGE`: "young" keeps the first
    `YOUNG_STAR_AGE_LIFESPAN_RATIO` of the window, "old" drops the first
    `OLD_STAR_AGE_LIFESPAN_RATIO`, anything else keeps it whole.
    """
    span = high - low
    if age_bias == "young":
        high = low + span * tuning.YOUNG_STAR_AGE_LIFESPAN_RATIO
    elif age_bias == "old":
        low = low + span * tuning.OLD_STAR_AGE_LIFESPAN_RATIO
    return low, high


def population_age_range_gy(population=None):
    """
    The age range (Gy) stars of `population` are drawn from: one of
    `STELLAR_POPULATION_AGE_RANGES_GY`'s keys ("young", "intermediate",
    "old", "bulge"), or None for the whole disk's star-formation history
    (`STAR_FORMATION_AGE_RANGE_GY`).
    """
    if population is None:
        return tuning.STAR_FORMATION_AGE_RANGE_GY
    try:
        return tuning.STELLAR_POPULATION_AGE_RANGES_GY[population]
    except KeyError:
        raise ValueError(f"unknown stellar population {population!r}; expected one of "
                         f"{sorted(tuning.STELLAR_POPULATION_AGE_RANGES_GY)} or None") from None


def sample_star_age_gy(age_bias=None, rng=random, population=None):
    """A star's age, uniform over its population's age range
    (`population_age_range_gy`), biased by `age_bias`."""
    low, high = age_window_gy(*population_age_range_gy(population), age_bias)
    return rng.uniform(low, high)


def main_sequence_luminosity_sol(mass_sol):
    """The piecewise mass-luminosity relation (`MS_MASS_LUMINOSITY_PIECES`)."""
    for piece in tuning.MS_MASS_LUMINOSITY_PIECES:
        if piece["max_mass_sol"] is None or mass_sol < piece["max_mass_sol"]:
            return piece["coeff"] * mass_sol ** piece["exponent"]


def main_sequence_radius_sol(mass_sol):
    """R = M^0.8 below 1 Msun, M^0.57 above."""
    exponent = (tuning.MS_RADIUS_EXPONENT_BELOW_1_SOL if mass_sol < 1.0
                else tuning.MS_RADIUS_EXPONENT_ABOVE_1_SOL)
    return mass_sol ** exponent


def effective_temperature_k(luminosity_sol, radius_sol):
    """Stefan-Boltzmann in solar units: T = T_sun * (L / R^2)^(1/4)."""
    return tuning.SUN_EFFECTIVE_TEMPERATURE_K * (luminosity_sol / radius_sol ** 2) ** 0.25


def white_dwarf_radius_km(mass_sol):
    """The generator's existing white dwarf mass-radius relation."""
    return constants.WHITE_DWARF_BASE_RADIUS_KM * mass_sol ** constants.WHITE_DWARF_MASS_RADIUS_EXPONENT


def spectral_letter_and_subclass(temperature_k):
    """
    The spectral letter whose `TEMP_RANGES` holds `temperature_k` (clamped
    to O above the table and M below it), and the 0-9 subclass within it
    (0 hottest), the same subclass rule the generator already uses.
    """
    ranges = constants.TEMP_RANGES
    letter = None
    for candidate, (low, high) in ranges.items():
        if low <= temperature_k <= high:
            letter = candidate
            break
    if letter is None:
        letter = "O" if temperature_k > ranges["O"][1] else "M"
    low, high = ranges[letter]
    clamped = min(max(temperature_k, low), high)
    subclass = constants.SUBCLASS_MAX_VALUE - round(
        (clamped - low) / (high - low) * constants.SUBCLASS_MAX_VALUE)
    return letter, subclass


def _log_uniform(low, high, rng):
    return math.exp(rng.uniform(math.log(low), math.log(high)))


def giant_luminosity_bands(mass_sol):
    """
    The luminosity bands (Lsun) a giant of initial mass `mass_sol` is drawn
    from, as `(low, high, share)`, dimmest first, each log-uniform within
    itself: a bright giant's one band, or a giant's main band and its short
    bright tip (`GIANT_TIP_FRACTION`, GEN.79).
    """
    pc = tuning
    if mass_sol >= pc.BRIGHT_GIANT_MIN_MASS_SOL:
        return ((*pc.BRIGHT_GIANT_LUMINOSITY_RANGE_SOL, 1.0),)
    return ((*pc.GIANT_LUMINOSITY_RANGE_SOL, 1.0 - pc.GIANT_TIP_FRACTION),
            (*pc.GIANT_TIP_LUMINOSITY_RANGE_SOL, pc.GIANT_TIP_FRACTION))


def _band_share_above(low, high, share, min_luminosity_sol):
    """The part of one log-uniform band's `share` at or above the value."""
    if min_luminosity_sol is None or min_luminosity_sol <= low:
        return share
    if min_luminosity_sol >= high:
        return 0.0
    return share * math.log(high / min_luminosity_sol) / math.log(high / low)


def giant_bright_chance(mass_sol, min_luminosity_sol):
    """The chance a giant's luminosity (`giant_luminosity_bands`) is at
    least `min_luminosity_sol`."""
    return sum(_band_share_above(low, high, share, min_luminosity_sol)
               for low, high, share in giant_luminosity_bands(mass_sol))


def sample_giant_luminosity_sol(mass_sol, rng=random, min_luminosity_sol=None):
    """
    A giant's luminosity from `giant_luminosity_bands`, or from their part
    at or above `min_luminosity_sol` (clamped to the top when none is).
    One `rng.random()` draw, whichever band it lands in.
    """
    bands = giant_luminosity_bands(mass_sol)
    weights = [_band_share_above(low, high, share, min_luminosity_sol) for low, high, share in bands]
    total = sum(weights)
    if total <= 0.0:
        return bands[-1][1]
    pick = rng.random() * total
    # The last band with any weight, unless the pick lands before it.
    index = max(i for i, weight in enumerate(weights) if weight > 0.0)
    for i, weight in enumerate(weights):
        if weight > 0.0 and pick < weight:
            index = i
            break
        pick -= weight
    low, high, _share = bands[index]
    if min_luminosity_sol is not None:
        low = min(max(low, min_luminosity_sol), high)
    fraction = min(max(pick / weights[index], 0.0), 1.0)
    return math.exp(math.log(low) + (math.log(high) - math.log(low)) * fraction)


def _supergiant_yerkes_class(luminosity_sol):
    yerkes = "IB"
    for name, threshold in tuning.SUPERGIANT_YERKES_THRESHOLDS_SOL.items():
        if luminosity_sol >= threshold:
            yerkes = name
    return yerkes


def evolve_star(mass_sol, age_gy, rng=random, min_luminosity_sol=None):
    """
    The present state of a star of initial mass `mass_sol` at `age_gy`.

    `min_luminosity_sol` draws a giant's or supergiant's luminosity (the
    only random ones) from the part of its range at or above that value;
    `stellarPopulation.sample_bright_stars` uses it once it has picked a
    mass and age that can reach it.

    Returns:
        dict or None: None when the star has already collapsed to a neutron
        star or black hole (initial mass at or above
        `SUPERGIANT_MIN_MASS_SOL`, past its giant phase). Otherwise:
        `yerkes_class`, `mass_sol` (present mass: a white dwarf's final
        mass), `luminosity_sol`, `temperature_k`, `radius_km` (a white
        dwarf's from its mass-radius relation, None for every other class,
        whose radius follows from luminosity and temperature),
        `lifespan_gy` (the age at which it becomes a remnant; infinite for a
        white dwarf) and `phase_end_gy` (the age at which it leaves its
        present phase; infinite for a white dwarf).
    """
    pc = tuning
    t_ms = main_sequence_lifetime_gy(mass_sol)
    t_sub = t_ms * pc.SUBGIANT_PHASE_END_MS_FRACTION
    t_giant = t_ms * pc.GIANT_PHASE_END_MS_FRACTION
    massive = mass_sol >= pc.SUPERGIANT_MIN_MASS_SOL
    l_ms = main_sequence_luminosity_sol(mass_sol)
    t_eff_ms = effective_temperature_k(l_ms, main_sequence_radius_sol(mass_sol))
    state = {"mass_sol": mass_sol, "radius_km": None, "lifespan_gy": t_giant}

    if age_gy < t_ms:
        state.update(yerkes_class="V", luminosity_sol=l_ms, temperature_k=t_eff_ms, phase_end_gy=t_ms)
    elif massive and age_gy < t_giant:
        growth_low, growth_high = pc.SUPERGIANT_LUMINOSITY_GROWTH_RANGE
        if min_luminosity_sol is not None:
            growth_low = min(max(growth_low, min_luminosity_sol / l_ms), growth_high)
        luminosity = l_ms * rng.uniform(growth_low, growth_high)
        if rng.random() < pc.RED_SUPERGIANT_FRACTION:
            temperature = rng.uniform(*pc.RED_SUPERGIANT_TEMPERATURE_RANGE_K)
        else:
            temperature = rng.uniform(*pc.BLUE_SUPERGIANT_TEMPERATURE_RANGE_K)
        state.update(yerkes_class=_supergiant_yerkes_class(luminosity), luminosity_sol=luminosity,
                     temperature_k=temperature, phase_end_gy=t_giant)
    elif age_gy < t_sub:
        progress = (age_gy - t_ms) / (t_sub - t_ms)
        end_temperature = min(t_eff_ms, pc.SUBGIANT_END_TEMPERATURE_K)
        state.update(yerkes_class="IV",
                     luminosity_sol=l_ms * (1 + (pc.SUBGIANT_MAX_BRIGHTENING - 1) * progress),
                     temperature_k=t_eff_ms + (end_temperature - t_eff_ms) * progress,
                     phase_end_gy=t_sub)
    elif age_gy < t_giant:
        bright = mass_sol >= pc.BRIGHT_GIANT_MIN_MASS_SOL
        state.update(yerkes_class="II" if bright else "III",
                     luminosity_sol=sample_giant_luminosity_sol(mass_sol, rng, min_luminosity_sol),
                     temperature_k=rng.uniform(*pc.GIANT_TEMPERATURE_RANGE_K),
                     phase_end_gy=t_giant)
    elif massive:
        return None
    else:
        final_mass = min(max(pc.WD_IFMR_SLOPE * mass_sol + pc.WD_IFMR_INTERCEPT_SOL, pc.WD_MASS_RANGE_SOL[0]),
                         pc.WD_MASS_RANGE_SOL[1])
        cooling_age = max(age_gy - t_giant, pc.WD_MIN_COOLING_AGE_GY)
        luminosity = pc.WD_COOLING_L0_SOL * (final_mass / 0.6) * cooling_age ** pc.WD_COOLING_EXPONENT
        luminosity = min(max(luminosity, pc.WD_LUMINOSITY_RANGE_SOL[0]), pc.WD_LUMINOSITY_RANGE_SOL[1])
        radius_km = white_dwarf_radius_km(final_mass)
        radius_sol = radius_km * constants.KM_TO_M_FACTOR / constants.SOLAR_RADIUS_M
        state.update(yerkes_class="VII", mass_sol=final_mass, luminosity_sol=luminosity,
                     temperature_k=effective_temperature_k(luminosity, radius_sol), radius_km=radius_km,
                     lifespan_gy=float("inf"), phase_end_gy=float("inf"))
    return state


def sample_living_star(age_bias=None, large_star=False, rng=random, habitable_host=False,
                       population=None, max_luminosity_sol=None):
    """
    Draws `(initial_mass_sol, age_gy, state)` for a random star that hasn't
    collapsed (see `evolve_star`), redrawing mass and age together until it
    hasn't. With `large_star`, the IMF is truncated at `LARGE_STAR_MIN_MASS_SOL`
    and the age drawn inside the star's own lifetime (and the default
    star-formation window), so a large star is always still shining.
    With `habitable_host` (a habitable world is required), the star must
    be at least `LIFE_MIN_STAR_AGE_GY` old and not a white dwarf, whose
    progenitor would have engulfed its habitable zone.

    `population` draws the age from that population's range instead of the
    whole disk's (see `population_age_range_gy`), and `max_luminosity_sol`
    redraws any star but a white dwarf at or above that luminosity (a
    sector's dim stars, once its bright ones were pre-placed).
    """
    pc = tuning
    if large_star and population is not None:
        oldest_large_star_gy = main_sequence_lifetime_gy(pc.LARGE_STAR_MIN_MASS_SOL) * pc.GIANT_PHASE_END_MS_FRACTION
        if population_age_range_gy(population)[0] >= oldest_large_star_gy:
            # Every star this massive in so old a population has died; a
            # large star is asked for anyway (+large_star), so it keeps the
            # whole disk's ages rather than fail.
            population = None
    if habitable_host and population is not None and population_age_range_gy(population)[1] <= pc.LIFE_MIN_STAR_AGE_GY:
        population = None  # likewise, no star this young can host life yet
    for _ in range(pc.STAR_MODEL_MAX_REDRAWS):
        if large_star:
            mass = sample_imf_mass_sol(min_mass_sol=pc.LARGE_STAR_MIN_MASS_SOL, rng=rng)
            low, high = population_age_range_gy(population)
            high = min(high, main_sequence_lifetime_gy(mass) * pc.GIANT_PHASE_END_MS_FRACTION)
            if high <= low:
                continue  # an old population has no living star this massive
            age = rng.uniform(*age_window_gy(low, high, age_bias))
        else:
            mass = sample_imf_mass_sol(rng=rng)
            age = sample_star_age_gy(age_bias, rng, population)
        if habitable_host and age < pc.LIFE_MIN_STAR_AGE_GY:
            continue
        state = evolve_star(mass, age, rng)
        if state is None or (habitable_host and state["yerkes_class"] == "VII"):
            continue
        # A white dwarf is never pre-placed, so the cap never removes one
        # (the hottest are clamped to exactly the lowest allowed cap).
        if (max_luminosity_sol is not None and state["yerkes_class"] != "VII"
                and state["luminosity_sol"] >= max_luminosity_sol):
            continue
        return mass, age, state
    raise ValueError(f"no living star in {pc.STAR_MODEL_MAX_REDRAWS} draws")


def star_params(initial_mass_sol, age_gy, state):
    """
    The stored form of a population-model star: everything `Star` needs to
    rebuild it without re-rolling (`Star.from_params`), in the units the
    database stores. The temperature is rounded to the nearest hundred
    kelvin, as for every generated star; the radius follows from
    Stefan-Boltzmann at that temperature (a white dwarf's from its
    mass-radius relation).

    Returns:
        dict: `type` (e.g. "G2V Yellow Main Sequence Star"),
        `yerkes_class`, `mass_kg`, `radius_km`, `temperature_k`,
        `luminosity_w`, `age_gy`, `lifespan_gy`, `initial_mass_sol` and
        `phase_end_age_gy`.
    """
    temperature = int(round(state["temperature_k"], tuning.ROUND_TEMPERATURE_NEAREST_HUNDRED))
    spectral_class, subclass = spectral_letter_and_subclass(temperature)
    luminosity_w = state["luminosity_sol"] * constants.SOLAR_LUMINOSITY
    if state["radius_km"] is not None:
        radius_km = state["radius_km"]
    else:
        radius_km = math.sqrt(luminosity_w / (constants.FOUR_PI * constants.STEFAN_BOLTZMANN_CONSTANT
                                              * temperature ** 4)) / constants.KM_TO_M_FACTOR
    yerkes = state["yerkes_class"]
    color = constants.SPECTRAL_CLASS_COLORS[spectral_class]
    return {
        "type": f"{spectral_class}{subclass}{yerkes} {color} {YERKES_CLASS_NAMES[yerkes]} Star",
        "yerkes_class": yerkes,
        "mass_kg": state["mass_sol"] * constants.SOLAR_MASS_TO_KG,
        "radius_km": radius_km,
        "temperature_k": temperature,
        "luminosity_w": luminosity_w,
        "age_gy": age_gy,
        "lifespan_gy": state["lifespan_gy"],
        "initial_mass_sol": initial_mass_sol,
        "phase_end_age_gy": state["phase_end_gy"],
    }
