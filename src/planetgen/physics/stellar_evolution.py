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
   (VII) or, for a massive star, a neutron star or black hole. A star lighter
   than `HELIUM_FLASH_MASS_SOL` burns helium in its core between its red giant
   branch and its white dwarf: a horizontal-branch star (red clump, RR Lyrae,
   blue or extreme, by [Fe/H] and age) and then an asymptotic giant (GEN.286).

The spectral letter, subclass and Yerkes class then follow from the
resulting temperature and luminosity (`spectral_letter_and_subclass`), so
a hot star is always luminous for its size, a supergiant is always
massive, and white dwarfs and red giants turn up at their real rates.

Every constant lives in `program_constants`' "Star population model"
section.
"""

import math

from planetgen.galaxy import scales
from planetgen.physics import constants, giant_branch, metallicity, ms_lifetime, remnant_mass, wolf_rayet
from planetgen import tuning
from planetgen.util import draw
from planetgen.physics import unit_registry as units

_km_to_m = units.converter("km", "m")
_m_to_km = units.converter("m", "km")


YERKES_CLASS_NAMES = {
    "0": "Hypergiant", "IA+": "Luminous Supergiant", "IA": "Supergiant",
    "IAB": "Intermediate-size Luminous Supergiant", "IB": "Less Luminous Supergiant",
    "II": "Bright Giant", "III": "Giant", "IV": "Subgiant", "V": "Main Sequence",
    "VI": "Subdwarf", "VII": "White Dwarf", "D": "White Dwarf",  # D is an alias for VII
    "W": "Wolf-Rayet",
}
"""Display name of each Yerkes luminosity class, used in `Star.type`."""

PHASE_NAMES = {"horizontal-branch": "Horizontal Branch Giant", "asymptotic-giant": "Asymptotic Giant"}
"""Display name of a giant's phase, in place of its class's name, in `Star.type`: a horizontal-branch star and an
asymptotic giant are both class III."""

HORIZONTAL_BRANCH_COLOURS = ("blue", "rr-lyrae", "red", "extreme")
"""The horizontal branch's colours, in the order of `tuning.HORIZONTAL_BRANCH_COLOURS`' shares."""


def main_sequence_lifetime_gy(mass_sol, fe_h=0.0):
    """The main-sequence lifetime, Gy, of a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` (default 0, the
    reference composition), from `ms_lifetime` (Hurley, Pols and Tout 2000, GEN.467), the one the rest of the
    generator uses."""
    return ms_lifetime.lifetime_gy(mass_sol, fe_h)


def sample_imf_mass_sol(min_mass_sol=None, max_mass_sol=None, rng=draw):
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


def imf_number(low_sol, high_sol):
    """The number of stars in `[low_sol, high_sol]` solar masses per unit of the IMF (`IMF_BREAKS_SOL` and `IMF_SLOPES`, the
    segments joined continuously at the breaks); only the ratio of two means anything. 0 for an empty range."""
    breaks, slopes = tuning.IMF_BREAKS_SOL, tuning.IMF_SLOPES
    norm, total = 1.0, 0.0
    for i, alpha in enumerate(slopes):
        if i > 0:
            norm *= breaks[i] ** (alpha - slopes[i - 1])
        a, b = max(low_sol, breaks[i]), min(high_sol, breaks[i + 1])
        if a < b:
            total += norm * _power_law_integral(a, b, alpha)
    return total


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


def sample_star_age_gy(age_bias=None, rng=draw, population=None, age_law=None):
    """A star's age over its population's age range (`population_age_range_gy`), biased by `age_bias`: uniform, or with
    the density of `age_law` (a `galaxy.age_law.AgeLaw` over that range, GEN.546) when given."""
    low, high = age_window_gy(*population_age_range_gy(population), age_bias)
    return rng.uniform(low, high) if age_law is None else age_law.sample(low, high, rng)


def _piecewise_luminosity_sol(pieces, mass_sol):
    for piece in pieces:
        if piece["max_mass_sol"] is None or mass_sol < piece["max_mass_sol"]:
            return piece["coeff"] * mass_sol ** piece["exponent"]


def main_sequence_luminosity_sol(mass_sol):
    """The piecewise mass-luminosity relation (`MS_MASS_LUMINOSITY_PIECES`)."""
    return _piecewise_luminosity_sol(tuning.MS_MASS_LUMINOSITY_PIECES, mass_sol)


def map_luminosity_sol(mass_sol):
    """
    The luminosity (Lsun) a main-sequence star of `mass_sol` would have by the
    map's relation (`MAP_MASS_LUMINOSITY_PIECES`, GEN.200): the one function
    that gives a faint but massive object (a neutron star, a black hole) the
    visibility on the Galaxy Map of the star of its mass, so major
    gravitational places are always on the map. It is not the star
    generators' relation (`main_sequence_luminosity_sol`): the two differ
    above 2 Msun (1.5 M^3.5 against 1.4 M^3.5, and a linear law from 20 Msun
    against 55), by about 7% from 2 to 20 Msun and by a factor of 12 at 4 million Msun.

    Args:
        mass_sol (float): The mass, solar masses (positive).

    Returns:
        float: The luminosity, solar luminosities; increasing with mass.
    """
    return _piecewise_luminosity_sol(tuning.MAP_MASS_LUMINOSITY_PIECES, mass_sol)


def phenomenon_map_luminosity_sol(map_type, mass_sol=None):
    """
    The luminosity (Lsun) the maps draw a point phenomenon at (MAP.208): that of
    the main-sequence star of its mass (`map_luminosity_sol`), not its own
    physical light, so a 10 Msun black hole is about as bright as a 10 Msun
    star. A black hole is `BLACK_HOLE_MAP_BRIGHTNESS` of that; an object with no
    mass takes `MASSLESS_MARKER_MASS_SOL`. The one function the three maps share:
    the Galaxy Map's tiles and the Sector Map's lights both start from it.

    Args:
        map_type (str): "black_hole", "neutron_star", "quasar" or anything else.
        mass_sol (float, optional): The mass, solar masses.

    Returns:
        float: The luminosity, solar luminosities.
    """
    mass = mass_sol if mass_sol and mass_sol > 0 else tuning.MASSLESS_MARKER_MASS_SOL
    share = tuning.BLACK_HOLE_MAP_BRIGHTNESS if map_type == "black_hole" else 1.0
    return map_luminosity_sol(mass) * share


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



def giant_luminosity_law(mass_sol):
    """
    The luminosity law of an intermediate-mass or heavier giant of initial mass `mass_sol` as `(low, high, slope)`
    (Lsun): the time it spends per unit of ln L falls as L^-slope between `low` and `high` (GEN.287). Not used below
    `HELIUM_FLASH_MASS_SOL`, whose red giant branch is `giant_branch.branch` (GEN.502). A giant lighter than
    `BRIGHT_GIANT_MIN_MASS_SOL` is on the red giant branch: from the base (`SUBGIANT_MAX_BRIGHTENING` times its
    main-sequence luminosity, where its subgiant phase ended) to the branch tip, slope `GIANT_BRANCH_TIME_SLOPE`. A
    bright giant is spread log-uniform (slope 0) over its own range.
    """
    pc = tuning
    if mass_sol >= pc.BRIGHT_GIANT_MIN_MASS_SOL:
        return (*pc.BRIGHT_GIANT_LUMINOSITY_RANGE_SOL, 0.0)
    tip = pc.GIANT_BRANCH_TIP_LUMINOSITY_SOL
    base = pc.SUBGIANT_MAX_BRIGHTENING * main_sequence_luminosity_sol(mass_sol)
    return (min(base, tip / 2.0), tip, pc.GIANT_BRANCH_TIME_SLOPE)


def luminosity_law_share_above(low, high, slope, min_luminosity_sol):
    """The part of the time spent between `low` and `high` (the law of
    `giant_luminosity_law`) spent at or above the value; 1.0 for `None`."""
    if min_luminosity_sol is None or min_luminosity_sol <= low:
        return 1.0
    if min_luminosity_sol >= high:
        return 0.0
    if slope == 0.0:
        return math.log(high / min_luminosity_sol) / math.log(high / low)
    return ((min_luminosity_sol ** -slope - high ** -slope) / (low ** -slope - high ** -slope))


def luminosity_law_quantile(low, high, slope, fraction, min_luminosity_sol=None):
    """The luminosity below which `fraction` (0 to 1) of the time spent at
    or above `min_luminosity_sol` (`None`: from `low`) is spent."""
    if min_luminosity_sol is not None:
        low = min(max(low, min_luminosity_sol), high)
    fraction = min(max(fraction, 0.0), 1.0)
    if slope == 0.0 or low >= high:
        return math.exp(math.log(low) + (math.log(high) - math.log(low)) * fraction)
    return (low ** -slope - fraction * (low ** -slope - high ** -slope)) ** (-1.0 / slope)


def giant_bright_chance(mass_sol, min_luminosity_sol, fe_h=0.0):
    """The chance a giant's luminosity is at least `min_luminosity_sol`: the share of the red giant branch above it
    for a star lighter than `HELIUM_FLASH_MASS_SOL` (`giant_branch`, GEN.502), else of `giant_luminosity_law`."""
    if has_core_helium_phases(mass_sol):
        return giant_branch.branch(mass_sol, fe_h).share_above(min_luminosity_sol)
    return luminosity_law_share_above(*giant_luminosity_law(mass_sol), min_luminosity_sol)


def sample_giant_luminosity_sol(mass_sol, rng=draw, min_luminosity_sol=None, fe_h=0.0):
    """
    A giant's luminosity from its red giant branch (`giant_branch`) below `HELIUM_FLASH_MASS_SOL`, else from
    `giant_luminosity_law`, or from its part at or above `min_luminosity_sol` (clamped to the top when none is).
    One `rng.random()` draw.
    """
    if has_core_helium_phases(mass_sol):
        return giant_branch.branch(mass_sol, fe_h).quantile(rng.random(), min_luminosity_sol)
    return luminosity_law_quantile(*giant_luminosity_law(mass_sol), rng.random(), min_luminosity_sol)


def subgiant_end_luminosity_sol(mass_sol, fe_h=0.0):
    """The luminosity a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` reaches as its subgiant phase ends:
    the base of its red giant branch below `HELIUM_FLASH_MASS_SOL` (`giant_branch`, GEN.502), else
    `SUBGIANT_MAX_BRIGHTENING` times its main-sequence luminosity."""
    l_ms = main_sequence_luminosity_sol(mass_sol)
    if has_core_helium_phases(mass_sol):
        return max(giant_branch.base_luminosity_sol(mass_sol, fe_h), 1.01 * l_ms)
    return tuning.SUBGIANT_MAX_BRIGHTENING * l_ms


def subgiant_end_fraction(mass_sol, fe_h=0.0):
    """The age a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` leaves the Hertzsprung gap, as a multiple of
    its main-sequence lifetime: the paper's time to the base of the giant branch below `HELIUM_FLASH_MASS_SOL`
    (GEN.502), else `SUBGIANT_PHASE_END_MS_FRACTION`."""
    if has_core_helium_phases(mass_sol):
        return giant_branch.subgiant_end_fraction(mass_sol, fe_h)
    return tuning.SUBGIANT_PHASE_END_MS_FRACTION


def giant_end_fraction(mass_sol, fe_h=0.0):
    """The age the giant phase ends (helium ignition below `HELIUM_FLASH_MASS_SOL`, GEN.502) of a star of birth
    [Fe/H] `fe_h`, as a multiple of the main-sequence lifetime; `GIANT_PHASE_END_MS_FRACTION` for the heavier stars."""
    if has_core_helium_phases(mass_sol):
        return giant_branch.giant_end_fraction(mass_sol, fe_h)
    return tuning.GIANT_PHASE_END_MS_FRACTION


def horizontal_branch_bright_chance(min_luminosity_sol):
    """The chance a horizontal-branch star (one luminosity range, whatever
    its colour) is at least `min_luminosity_sol` bright."""
    return luminosity_law_share_above(*tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL, 0.0, min_luminosity_sol)


def agb_bright_chance(min_luminosity_sol):
    """The chance an asymptotic giant is at least `min_luminosity_sol`
    bright: the time-at-luminosity law of the red giant branch between
    `AGB_LUMINOSITY_RANGE_SOL`."""
    return luminosity_law_share_above(*tuning.AGB_LUMINOSITY_RANGE_SOL, tuning.GIANT_BRANCH_TIME_SLOPE,
                                      min_luminosity_sol)


def has_core_helium_phases(mass_sol):
    """Whether a star of initial mass `mass_sol` has a horizontal branch and an asymptotic giant branch: below
    `HELIUM_FLASH_MASS_SOL` (and not a supergiant)."""
    return mass_sol < tuning.HELIUM_FLASH_MASS_SOL


def remnant_age_fraction(mass_sol, fe_h=0.0):
    """The age at which a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` ends its last
    giant phase, as a multiple of its main-sequence lifetime: the giant
    phase's end, or for a star with core helium burning the end of its
    asymptotic giant branch."""
    pc = tuning
    if has_core_helium_phases(mass_sol):
        return (giant_end_fraction(mass_sol, fe_h)
                + pc.HORIZONTAL_BRANCH_DURATION_GY / main_sequence_lifetime_gy(mass_sol, fe_h) + pc.AGB_MS_FRACTION)
    return pc.GIANT_PHASE_END_MS_FRACTION


def remnant_formation_time_gy(mass_sol, fe_h=0.0):
    """The age at which a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` has become its remnant, Gy: the end
    of its last giant phase."""
    return main_sequence_lifetime_gy(mass_sol, fe_h) * remnant_age_fraction(mass_sol, fe_h)


def youngest_progenitor_sol(age_yr, fe_h=0.0):
    """The lightest star of birth [Fe/H] `fe_h` whose giant phases have ended by this age (a star lighter than
    `HELIUM_FLASH_MASS_SOL` also burns helium on the horizontal branch and climbs the asymptotic giant branch, and its
    giant branch has the length of `giant_branch`, GEN.502)."""
    age_gy = age_yr / 1.0e9
    flash = tuning.HELIUM_FLASH_MASS_SOL

    if remnant_formation_time_gy(flash * (1.0 - 1.0e-9), fe_h) <= age_gy:
        low, high = ms_lifetime.MIN_MASS_SOL, flash
        if remnant_formation_time_gy(low, fe_h) <= age_gy:
            return low
        for _ in range(50):
            middle = math.sqrt(low * high)
            if remnant_formation_time_gy(middle, fe_h) > age_gy:
                low = middle
            else:
                high = middle
        return high
    return max(ms_lifetime.mass_for_lifetime_gy(age_gy / tuning.GIANT_PHASE_END_MS_FRACTION, fe_h), flash)


def horizontal_branch_shares(fe_h, age_gy):
    """
    The shares `(blue, rr_lyrae, red, extreme)` of the horizontal branch at
    a star's [Fe/H] and age: `tuning.HORIZONTAL_BRANCH_COLOURS` linear
    between its anchors, from the effective [Fe/H] (higher for a younger
    star: the second parameter). They sum to 1.
    """
    pc = tuning
    effective = fe_h + pc.HORIZONTAL_BRANCH_AGE_SHIFT_DEX_PER_GY * (pc.HORIZONTAL_BRANCH_REFERENCE_AGE_GY - age_gy)
    table = pc.HORIZONTAL_BRANCH_COLOURS
    if effective <= table[0][0]:
        shares = table[0][1]
    elif effective >= table[-1][0]:
        shares = table[-1][1]
    else:
        for (lo_fe, lo_shares), (hi_fe, hi_shares) in zip(table, table[1:]):
            if effective <= hi_fe:
                weight = (effective - lo_fe) / (hi_fe - lo_fe)
                shares = tuple(a + (b - a) * weight for a, b in zip(lo_shares, hi_shares))
                break
    total = sum(shares)
    return tuple(share / total for share in shares)


def sample_horizontal_branch_colour(fe_h, age_gy, rng=draw):
    """A horizontal-branch star's colour, one of `HORIZONTAL_BRANCH_COLOURS` (one `rng.random()` draw)."""
    pick = rng.random()
    shares = horizontal_branch_shares(fe_h, age_gy)
    for colour, share in zip(HORIZONTAL_BRANCH_COLOURS, shares):
        if pick < share:
            return colour
        pick -= share
    return HORIZONTAL_BRANCH_COLOURS[-1]


def _supergiant_yerkes_class(luminosity_sol):
    yerkes = "IB"
    for name, threshold in tuning.SUPERGIANT_YERKES_THRESHOLDS_SOL.items():
        if luminosity_sol >= threshold:
            yerkes = name
    return yerkes


def birth_fe_h(rng=draw, population=None, fe_h_peaks=None, galactic_center_dist_ly=None):
    """
    A star's birth [Fe/H] (GEN.651, GEN.553's order: composition first): drawn from its `population` (or
    `fe_h_peaks`, a cluster's or satellite's own distribution) at `galactic_center_dist_ly` from the centre (the Sun's
    distance without one) by `metallicity.sample_fe_h`. A cluster or satellite star with no distribution of its own is
    drawn as the halo's.
    """
    if population in ("cluster", "satellite", "satellite_young", "satellite_intermediate") and not fe_h_peaks:
        population = "halo"
    distance = scales.default_galactocentric_ly() if galactic_center_dist_ly is None else galactic_center_dist_ly
    return metallicity.sample_fe_h(population, distance, rng, peaks=fe_h_peaks)


def evolve_star(mass_sol, age_gy, rng=draw, min_luminosity_sol=None, fe_h=None, population=None, fe_h_peaks=None,
                galactic_center_dist_ly=None):
    """
    The present state of a star of initial mass `mass_sol` at `age_gy`.

    `min_luminosity_sol` draws a giant's or supergiant's luminosity (the
    only random ones) from the part of its range at or above that value;
    `star_population.sample_bright_stars` uses it once it has picked a
    mass and age that can reach it.

    The star's birth [Fe/H] sets its lifetimes and its giant branch (`ms_lifetime`, `giant_branch`, GEN.651) and a
    horizontal-branch star's colour. `fe_h` is the star's own; without it one is drawn first, before any other draw
    here (`birth_fe_h`: the star's `population` or `fe_h_peaks`, at `galactic_center_dist_ly` from the centre, the Sun's
    distance without one), and either way returned as `metallicity_fe_h` for the star to keep.

    Returns:
        dict or None: None when the star has already collapsed to a neutron
        star or black hole (initial mass at or above
        `SUPERGIANT_MIN_MASS_SOL`, past its giant phase). Otherwise:
        `yerkes_class`, `mass_sol` (present mass: a white dwarf's final
        mass), `luminosity_sol`, `temperature_k`, `radius_km` (a white
        dwarf's from its mass-radius relation, None for every other class,
        whose radius follows from luminosity and temperature),
        `metallicity_fe_h` (the birth [Fe/H]), `lifespan_gy` (the age at which it becomes a remnant; infinite for a
        white dwarf) and `phase_end_gy` (the age at which it leaves its
        present phase; infinite for a white dwarf). A horizontal-branch star
        or an asymptotic giant (class III) also has `phase`
        (`"horizontal-branch"` or `"asymptotic-giant"`, a key of
        `PHASE_NAMES`). A star younger than its contraction time
        (`pre_main_sequence.state`, GEN.242) is a `"protostar"` (under 0.5
        Myr) or a `"PMS"` (pre-main-sequence) class, with `disk_class` ("0",
        "I", "II" or "III") in place of a Yerkes numeral, and `phase_end_gy`
        the age it reaches the main sequence.
    """
    pc = tuning
    if fe_h is None:
        fe_h = birth_fe_h(rng, population, fe_h_peaks, galactic_center_dist_ly)
    t_ms = main_sequence_lifetime_gy(mass_sol, fe_h)
    t_sub = t_ms * subgiant_end_fraction(mass_sol, fe_h)
    t_giant = t_ms * giant_end_fraction(mass_sol, fe_h)
    t_horizontal = t_agb = t_giant
    if has_core_helium_phases(mass_sol):
        t_horizontal = t_giant + pc.HORIZONTAL_BRANCH_DURATION_GY
        t_agb = t_horizontal + t_ms * pc.AGB_MS_FRACTION
    massive = mass_sol >= pc.SUPERGIANT_MIN_MASS_SOL
    l_ms = main_sequence_luminosity_sol(mass_sol)
    t_eff_ms = effective_temperature_k(l_ms, main_sequence_radius_sol(mass_sol))
    state = {"mass_sol": mass_sol, "radius_km": None, "lifespan_gy": t_agb, "metallicity_fe_h": fe_h}

    wr_start = wolf_rayet.start_gy(mass_sol, t_ms, t_giant, fe_h) if massive else None
    young = pre_main_sequence.state(mass_sol, age_gy * 1.0e9)
    if young is not None:
        # GEN.242: still contracting, a protostar (Class 0 and I) or a pre-main-sequence star.
        state.update(young)
    elif age_gy < t_ms:
        state.update(yerkes_class="V", luminosity_sol=l_ms, temperature_k=t_eff_ms, phase_end_gy=t_ms)
    elif wr_start is not None and wr_start <= age_gy < t_giant:
        # GEN.243: the last stretch of a very massive star's life, in place of the end of its supergiant window.
        state.update(wolf_rayet.state(mass_sol, age_gy, wr_start, t_giant, rng, min_luminosity_sol))
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
                     luminosity_sol=l_ms + (subgiant_end_luminosity_sol(mass_sol, fe_h) - l_ms) * progress,
                     temperature_k=t_eff_ms + (end_temperature - t_eff_ms) * progress,
                     phase_end_gy=t_sub)
    elif age_gy < t_giant:
        bright = mass_sol >= pc.BRIGHT_GIANT_MIN_MASS_SOL
        state.update(yerkes_class="II" if bright else "III",
                     luminosity_sol=sample_giant_luminosity_sol(mass_sol, rng, min_luminosity_sol, fe_h),
                     temperature_k=rng.uniform(*pc.GIANT_TEMPERATURE_RANGE_K),
                     phase_end_gy=t_giant)
    elif age_gy < t_horizontal:
        colour = sample_horizontal_branch_colour(fe_h, age_gy, rng)
        state.update(yerkes_class="III", phase="horizontal-branch",
                     luminosity_sol=luminosity_law_quantile(*pc.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL, 0.0,
                                                            rng.random(), min_luminosity_sol),
                     temperature_k=rng.uniform(*pc.HORIZONTAL_BRANCH_COLOUR_TEMPERATURE_RANGE_K[colour]),
                     phase_end_gy=t_horizontal)
    elif age_gy < t_agb:
        state.update(yerkes_class="III", phase="asymptotic-giant",
                     luminosity_sol=luminosity_law_quantile(*pc.AGB_LUMINOSITY_RANGE_SOL, pc.GIANT_BRANCH_TIME_SLOPE,
                                                            rng.random(), min_luminosity_sol),
                     temperature_k=rng.uniform(*pc.GIANT_TEMPERATURE_RANGE_K),
                     phase_end_gy=t_agb)
    elif massive:
        return None
    else:
        final_mass = remnant_mass.white_dwarf_mass_sol(mass_sol)
        cooling_age = max(age_gy - t_agb, pc.WD_MIN_COOLING_AGE_GY)
        luminosity = pc.WD_COOLING_L0_SOL * (final_mass / 0.6) * cooling_age ** pc.WD_COOLING_EXPONENT
        luminosity = min(max(luminosity, pc.WD_LUMINOSITY_RANGE_SOL[0]), pc.WD_LUMINOSITY_RANGE_SOL[1])
        radius_km = white_dwarf_radius_km(final_mass)
        radius_sol = _km_to_m(radius_km) / constants.SOLAR_RADIUS_M
        state.update(yerkes_class="VII", mass_sol=final_mass, luminosity_sol=luminosity,
                     temperature_k=effective_temperature_k(luminosity, radius_sol), radius_km=radius_km,
                     lifespan_gy=float("inf"), phase_end_gy=float("inf"))
    return state


def sample_living_star(age_bias=None, large_star=None, rng=draw, habitable_host=False,
                       population=None, max_luminosity_sol=None, max_mass_sol=None, fe_h_peaks=None, age_law=None,
                       galactic_center_dist_ly=None, fe_h=None):
    """
    Draws `(initial_mass_sol, age_gy, state)` for a random star that hasn't
    collapsed (see `evolve_star`), redrawing mass and age together until it
    hasn't. With `large_star`, the IMF is truncated at `LARGE_STAR_MIN_MASS_SOL`
    and the age drawn inside the star's own lifetime (and the default
    star-formation window), so a large star is always still shining;
    `large_star=False` (`-large_star`) redraws any star born that massive.
    With `habitable_host` (a habitable world is required), the star must
    be at least `LIFE_MIN_STAR_AGE_GY` old and not a white dwarf, whose
    progenitor would have engulfed its habitable zone.

    `population` draws the age from that population's range instead of the
    whole disk's (see `population_age_range_gy`), and `max_luminosity_sol`
    redraws any star but a white dwarf at or above that luminosity (a
    sector's dim stars, once its bright ones were pre-placed); likewise
    `max_mass_sol` redraws any star born with at least that mass (the stars
    the scatter's mass pass placed). `fe_h_peaks` is the place's [Fe/H]
    distribution and `galactic_center_dist_ly` the star's distance from the centre, the place its birth [Fe/H] is
    drawn for (`birth_fe_h`, GEN.651: each try draws it first, before the mass and the age, and the star's lifetimes
    and giant branch follow it; `fe_h` pins it instead, a companion's, born with its primary). `age_law` (a `galaxy.age_law.AgeLaw` over the
    population's range, GEN.546) spreads the age inside that range by the disk's star formation rate instead of
    evenly; it goes with the population, so a fallback to the whole disk's ages drops it.
    """
    pc = tuning
    if large_star and population is not None:
        oldest_large_star_gy = main_sequence_lifetime_gy(pc.LARGE_STAR_MIN_MASS_SOL) * remnant_age_fraction(
            pc.LARGE_STAR_MIN_MASS_SOL)
        if population_age_range_gy(population)[0] >= oldest_large_star_gy:
            # Every star this massive in so old a population has died; a
            # large star is asked for anyway (+large_star), so it keeps the
            # whole disk's ages rather than fail.
            population, age_law = None, None
    if habitable_host and population is not None and population_age_range_gy(population)[1] <= pc.LIFE_MIN_STAR_AGE_GY:
        population, age_law = None, None  # likewise, no star this young can host life yet
    for _ in range(pc.STAR_MODEL_MAX_REDRAWS):
        birth = fe_h if fe_h is not None else birth_fe_h(rng, population, fe_h_peaks, galactic_center_dist_ly)
        if large_star:
            mass = sample_imf_mass_sol(min_mass_sol=pc.LARGE_STAR_MIN_MASS_SOL, rng=rng)
            if max_mass_sol is not None and mass >= max_mass_sol:
                continue
            low, high = population_age_range_gy(population)
            high = min(high, remnant_formation_time_gy(mass, birth))
            if high <= low:
                continue  # an old population has no living star this massive
            low, high = age_window_gy(low, high, age_bias)
            age = rng.uniform(low, high) if age_law is None else age_law.sample(low, high, rng)
        else:
            mass = sample_imf_mass_sol(rng=rng)
            if large_star is False and mass >= pc.LARGE_STAR_MIN_MASS_SOL:
                continue
            if max_mass_sol is not None and mass >= max_mass_sol:
                continue
            age = sample_star_age_gy(age_bias, rng, population, age_law)
        if habitable_host and age < pc.LIFE_MIN_STAR_AGE_GY:
            continue
        state = evolve_star(mass, age, rng, fe_h=birth, population=population, fe_h_peaks=fe_h_peaks)
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
        `phase_end_age_gy`; and `metallicity_fe_h`, the birth [Fe/H] the star keeps (GEN.651).
    """
    temperature = int(round(state["temperature_k"], tuning.ROUND_TEMPERATURE_NEAREST_HUNDRED))
    spectral_class, subclass = spectral_letter_and_subclass(temperature)
    luminosity_w = state["luminosity_sol"] * constants.SOLAR_LUMINOSITY
    if state["radius_km"] is not None:
        radius_km = state["radius_km"]
    else:
        radius_km = math.sqrt(luminosity_w / (constants.FOUR_PI * constants.STEFAN_BOLTZMANN_CONSTANT
                                              * temperature ** 4))
        radius_km = _m_to_km(radius_km)
    yerkes = state["yerkes_class"]
    color = constants.SPECTRAL_CLASS_COLORS[spectral_class]
    if yerkes == wolf_rayet.CLASS:
        # GEN.243: the spectral text is the subtype and subclass (WN6), in place of a letter, and the class says Wolf-Rayet.
        subtype = state["wr_subtype"]
        star_type = f"{subtype}{wolf_rayet.subclass(subtype, temperature)} {color} {wolf_rayet.CLASS_NAME} Star"
    elif yerkes in pre_main_sequence.STAGE_CLASSES:
        # No Yerkes numeral: a contracting star is named for its stage, and a protostar for its Class.
        name = (f"Class {state['disk_class']} Protostar" if yerkes == pre_main_sequence.PROTOSTAR
                else pre_main_sequence.CLASS_NAMES[yerkes])
        star_type = f"{spectral_class}{subclass} {color} {name}"
    else:
        star_type = (f"{spectral_class}{subclass}{yerkes} {color} "
                     f"{PHASE_NAMES.get(state.get('phase'), YERKES_CLASS_NAMES[yerkes])} Star")
    return {
        "type": star_type,
        "yerkes_class": yerkes,
        "mass_kg": state["mass_sol"] * constants.SOLAR_MASS_TO_KG,
        "radius_km": radius_km,
        "temperature_k": temperature,
        "luminosity_w": luminosity_w,
        "age_gy": age_gy,
        "lifespan_gy": state["lifespan_gy"],
        "initial_mass_sol": initial_mass_sol,
        "phase_end_age_gy": state["phase_end_gy"],
        "metallicity_fe_h": state["metallicity_fe_h"],
    }


# Last, because `pre_main_sequence` reads the main-sequence relations above.
from planetgen.physics import pre_main_sequence  # noqa: E402
