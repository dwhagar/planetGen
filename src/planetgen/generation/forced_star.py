"""
Stars of a forced type (`SystemConfig.STAR_TYPE`, GEN.941), drawn from the star model.

A forced type such as `K2III` is not a recipe of its own: it asks for a star the population model
(`stellar_evolution.evolve_star`) can make, so the star returned is one of the model's, a mass and an age whose state has that
letter, subclass and Yerkes class. This module scans the model once per [Fe/H] (a grid of masses and ages, weighted by the
mass function and the time a star spends in each cell), then picks a cell of the type, jitters its mass and age inside the cell
and evolves the star for real, so the result has the model's mass, luminosity, age, lifespan and phase end.

A type the model never makes (a B giant, an M hypergiant) is still a type the person asked for: the star is the model's nearest
one of that Yerkes class by temperature, with the temperature set to the request (`ForcedStar.exact` is False and the caller
logs it). Two classes are aliases (`IA+` is the hypergiant `0`, `D` the white dwarf `VII`). A subdwarf (`VI`) is a metal-poor
main-sequence star (a halo star at `SUBDWARF_FE_H_RANGE`) or, for the hot letters, a blue horizontal-branch star
(`SUBDWARF_HOT_MIN_TEMPERATURE_K`); GEN.943 decides what the class means beyond that.
"""
import functools
import math
from collections import defaultdict, namedtuple

import numpy as np

from planetgen import tuning
from planetgen.physics import constants, imf, post_main_sequence
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw

ForcedStar = namedtuple("ForcedStar", "initial_mass_sol age_gy state exact fe_h")
_Cell = namedtuple("_Cell", "weight mass_sol age_gy mass_ratio age_step_gy temperature_k present_mass_sol luminosity_sol")
Summary = namedtuple("Summary", "mass_sol luminosity_sol share")

ALIASES = {"IA+": "0", "D": "VII"}
"""A forced Yerkes class the model names otherwise."""

_FE_H_BUCKETS = (-2.0, -1.5, -1.0, -0.5, 0.0, 0.5)
_MODEL_CLASSES = ("0", "IA", "IAB", "IB", "II", "III", "IV", "V", "VII")


def _bucket(fe_h):
    return min(_FE_H_BUCKETS, key=lambda b: abs(b - fe_h))


def _key(state):
    """The (letter, subclass) of a model state at the temperature the star is stored with."""
    temperature = int(round(state["temperature_k"], tuning.ROUND_TEMPERATURE_NEAREST_HUNDRED))
    return evolution.spectral_letter_and_subclass(temperature)


def _mass_grid():
    low, high = imf.lowest_mass_sol(), imf.highest_mass_sol()
    steps = tuning.FORCED_STAR_MASS_STEPS
    return [low * (high / low) ** (i / (steps - 1)) for i in range(steps)]


def _ages(line, mass_sol, fe_h):
    """`(age, step)` pairs over one mass's life up to the age of the universe: the main sequence, the stretch after it and, for
    a star that leaves a white dwarf, the white dwarf's cooling (log-spaced)."""
    universe = tuning.UNIVERSE_AGE_GY
    t_ms, t_end = line.main_sequence_gy, line.end_gy
    n_ms, n_post, n_wd = tuning.FORCED_STAR_AGE_STEPS
    out = []
    top = min(t_ms, universe)
    out += [((k + 0.5) * top / n_ms, top / n_ms) for k in range(n_ms)]
    if t_ms < universe:
        top = min(t_end, universe)
        step = (top - t_ms) / n_post
        out += [(t_ms + (k + 0.5) * step, step) for k in range(n_post)]
    if t_end < universe and mass_sol < tuning.SUPERGIANT_MIN_MASS_SOL:
        low, high = tuning.WD_MIN_COOLING_AGE_GY, universe - t_end
        ratio = (high / low) ** (1.0 / n_wd)
        for k in range(n_wd):
            cooling = low * ratio ** (k + 0.5)
            out.append((t_end + cooling, cooling * (ratio - 1.0 / ratio) / 2.0))
    return out


@functools.lru_cache(maxsize=None)
def _table(bucket):
    """
    The model's cells at one [Fe/H]: `{(letter, subclass, yerkes): [_Cell, ...]}`, each cell a mass and an age with the weight
    a random star of the population would have there (the mass function per log mass times the time the cell spans), the
    cell's size and its temperature. Built once per process from a stream of its own.
    """
    rng = draw.Stream(f"forced-star-table:{bucket}")
    masses = _mass_grid()
    ratio = (masses[-1] / masses[0]) ** (1.0 / (len(masses) - 1))
    table = defaultdict(list)
    for mass in masses:
        try:
            line = post_main_sequence.timeline(mass, bucket)
        except ValueError:
            continue  # no track after the main sequence at this mass: the main sequence alone is skipped too
        per_log_mass = float(imf.density(np.array([mass]))[0]) * mass
        for age, step in _ages(line, mass, bucket):
            state = evolution.evolve_star(mass, age, rng, fe_h=bucket)
            if state is None or state["yerkes_class"] not in _MODEL_CLASSES:
                continue
            letter, subclass = _key(state)
            cell = _Cell(per_log_mass * step, mass, age, ratio, step, state["temperature_k"], state["mass_sol"],
                         state["luminosity_sol"])
            table[(letter, subclass, state["yerkes_class"])].append(cell)
            if state.get("phase") == "horizontal-branch" and state["temperature_k"] >= tuning.SUBDWARF_HOT_MIN_TEMPERATURE_K:
                table[(letter, subclass, "VI")].append(cell)
            if state["yerkes_class"] == "V" and age >= tuning.SUBDWARF_MIN_AGE_GY and bucket <= tuning.SUBDWARF_FE_H_RANGE[1]:
                table[(letter, subclass, "VI")].append(cell)
    return dict(table)


def _pick(cells, age_bias, rng):
    """One of `cells` by weight; with an age bias, from the older (or younger) half of the weight only."""
    if age_bias in ("old", "young") and len(cells) > 1:
        ordered = sorted(cells, key=lambda cell: cell.age_gy, reverse=age_bias == "old")
        half = sum(cell.weight for cell in ordered) / 2.0
        kept, total = [], 0.0
        for cell in ordered:
            kept.append(cell)
            total += cell.weight
            if total >= half:
                break
        cells = kept
    pick = rng.random() * sum(cell.weight for cell in cells)
    for cell in cells:
        pick -= cell.weight
        if pick <= 0.0:
            return cell
    return cells[-1]


def target_temperature_k(letter, subclass):
    """The temperature a letter and subclass stand for: subclass 0 at the hot end of the letter's `TEMP_RANGES`, 9 at the cool
    end, to the hundred kelvin stars are stored at."""
    low, high = constants.TEMP_RANGES[letter]
    kelvin = low + (constants.SUBCLASS_MAX_VALUE - subclass) * (high - low) / constants.SUBCLASS_MAX_VALUE
    return int(round(kelvin, tuning.ROUND_TEMPERATURE_NEAREST_HUNDRED))


def _nearest_cells(yerkes, target_k, bucket):
    """The cells of `yerkes` closest in temperature to `target_k` (those within `FORCED_STAR_NEAREST_SLACK` of the best)."""
    cells = [(abs(math.log(cell.temperature_k / target_k)), cell)
             for (_, _, found), group in _table(bucket).items() if found == yerkes for cell in group]
    if not cells:
        raise ValueError(f"the star model makes no star of Yerkes class {yerkes}")
    best = min(distance for distance, _ in cells)
    return [cell for distance, cell in cells if distance <= best + tuning.FORCED_STAR_NEAREST_SLACK]


def sample(letter, subclass, yerkes, age_bias=None, rng=draw, population=None, fe_h_peaks=None,
           galactic_center_dist_ly=None, fe_h=None):
    """
    A star of the forced type `letter` `subclass` `yerkes` from the model: a `ForcedStar` with the initial mass, the age, the
    `evolve_star` state and the [Fe/H]. `age_bias` ("old" or "young") picks the older or younger half of the type's
    cells. The birth [Fe/H] is `fe_h` or drawn as for any star (`stellar_evolution.birth_fe_h`), except for a subdwarf's, uniform inside
    `SUBDWARF_FE_H_RANGE`.

    Raises:
        ValueError: when the model makes no star of the Yerkes class at all.
    """
    model_class = ALIASES.get(yerkes, yerkes)
    if fe_h is not None:
        pass
    elif model_class == "VI":
        fe_h = rng.uniform(*tuning.SUBDWARF_FE_H_RANGE)
    else:
        fe_h = evolution.birth_fe_h(rng, population, fe_h_peaks, galactic_center_dist_ly)
    bucket = _bucket(fe_h)
    cells = _table(bucket).get((letter, subclass, model_class))
    exact = bool(cells)
    if not exact:
        cells = _nearest_cells(model_class, target_temperature_k(letter, subclass), bucket)
    cell = _pick(cells, age_bias, rng)

    for _ in range(tuning.FORCED_STAR_JITTER_TRIES):
        mass = cell.mass_sol * cell.mass_ratio ** (rng.random() - 0.5)
        age = min(max(cell.age_gy + (rng.random() - 0.5) * cell.age_step_gy, 0.0), tuning.UNIVERSE_AGE_GY)
        state = evolution.evolve_star(mass, age, rng, fe_h=fe_h)
        if state is None or state["yerkes_class"] not in _MODEL_CLASSES:
            continue
        if model_class == "VI":
            # A subdwarf is the main sequence of a halo star or a hot horizontal-branch star, read as class VI.
            hot = state.get("phase") == "horizontal-branch" and state["temperature_k"] >= tuning.SUBDWARF_HOT_MIN_TEMPERATURE_K
            if not (hot or (state["yerkes_class"] == "V" and age >= tuning.SUBDWARF_MIN_AGE_GY)):
                continue
        elif state["yerkes_class"] != model_class:
            continue
        if exact and _key(state) != (letter, subclass):
            continue
        return ForcedStar(mass, age, state, exact, fe_h)
    # The jittered stars fell out of the type (a state with a draw in it, a horizontal-branch colour or a white dwarf's
    # atmosphere, moves with the draw); the cell's own star is in it, at the temperature the table recorded.
    state = evolution.evolve_star(cell.mass_sol, cell.age_gy, rng, fe_h=bucket)
    state["temperature_k"] = cell.temperature_k
    return ForcedStar(cell.mass_sol, cell.age_gy, state, exact, bucket)


def _summary(cells, total_weight=None):
    """The ranges of `cells`, and their share of `total_weight` (None without one)."""
    return Summary((min(c.present_mass_sol for c in cells), max(c.present_mass_sol for c in cells)),
                   (min(c.luminosity_sol for c in cells), max(c.luminosity_sol for c in cells)),
                   None if total_weight is None else sum(c.weight for c in cells) / total_weight)


def _population_weight():
    return sum(cell.weight for (_, _, yerkes), group in _table(0.0).items() if yerkes != "VI" for cell in group)


def letter_summary(letter):
    """The present mass and luminosity ranges (solar) of the model's main-sequence stars of one spectral letter at solar
    [Fe/H], and the share of all the model's stars (every class, remnants included) that are of the letter (the class reference
    pages, GEN.941)."""
    cells = [cell for (found, _, yerkes), group in _table(0.0).items() if found == letter and yerkes == "V" for cell in group]
    everything = [cell for (found, _, yerkes), group in _table(0.0).items() if found == letter and yerkes != "VI"
                  for cell in group]
    ranges = _summary(cells)
    return Summary(ranges.mass_sol, ranges.luminosity_sol, sum(c.weight for c in everything) / _population_weight())


def class_summary(yerkes):
    """The present mass and luminosity ranges (solar) of the model's stars of one Yerkes class at solar [Fe/H] (a subdwarf's at
    the halo's lower edge of `SUBDWARF_FE_H_RANGE`; a Wolf-Rayet star's from the end of the tracks of the stars that become one),
    and the class's share of all the model's stars (the class reference pages, GEN.941). Raises ValueError for a class it never
    makes."""
    if yerkes == "W":
        masses = [tuning.WR_MIN_INITIAL_MASS_SOL, 40.0, 60.0, 100.0, imf.highest_mass_sol()]
        present = [post_main_sequence.timeline(m, 0.0).state(post_main_sequence.timeline(m, 0.0).end_gy - 1.0e-9).mass_sol
                   for m in masses]
        return Summary((min(present), max(present)), tuning.WR_LUMINOSITY_RANGE_SOL, None)
    if yerkes == "VI":
        bucket = _bucket(tuning.SUBDWARF_FE_H_RANGE[0])
        cells = {id(cell): cell for (_, _, found), group in _table(bucket).items() if found == "VI" for cell in group}.values()
        return _summary(list(cells))
    cells = [cell for (_, _, found), group in _table(0.0).items() if found == ALIASES.get(yerkes, yerkes) for cell in group]
    if not cells:
        raise ValueError(f"the star model makes no star of Yerkes class {yerkes}")
    return _summary(cells, _population_weight())
