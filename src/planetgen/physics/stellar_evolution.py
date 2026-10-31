# planetgen/physics/stellar_evolution.py

"""
Star Population Model
=====================

Draws a random star as physics rather than from a spectral-letter table:

1. a mass from the primary mass function, the single-star initial mass function less the companions of the heavier
   primaries (`imf.sample`, `imf.PRIMARY`, GEN.928),
2. an age from the star-formation history (`sample_star_age_gy`),
3. the star's present state from how far that age is into its own
   main-sequence lifetime (`evolve_star`): main sequence (V), subgiant
   (IV), giant (III/II) or supergiant (IB/IAB/IA/0), then a white dwarf
   (VII) or, for a massive star, a neutron star or black hole. After the main sequence the star's
   luminosity, temperature and radius are read from the MIST tracks (`post_main_sequence`, GEN.935). A star lighter
   than `HELIUM_FLASH_MASS_SOL` burns helium in its core between its red giant
   branch and its white dwarf: a horizontal-branch star (red clump, RR Lyrae,
   blue or extreme, by [Fe/H] and age) and then an asymptotic giant (GEN.286); a heavier one burns it
   as a core-helium-burning giant (or, from `SUPERGIANT_MIN_MASS_SOL`, a supergiant) and then climbs the asymptotic
   giant branch.

The spectral letter, subclass and Yerkes class then follow from the
resulting temperature and luminosity (`spectral_letter_and_subclass`), so
a hot star is always luminous for its size, a supergiant is always
massive, and white dwarfs and red giants turn up at their real rates.

Every constant lives in `program_constants`' "Star population model"
section.
"""

import functools
import math

from planetgen.datasets import fetch, mamajek_sequence, white_dwarf_models
from planetgen.galaxy import scales
from planetgen.physics import (
    age_metallicity, constants, imf, main_sequence, metallicity, ms_lifetime, post_main_sequence, remnant_mass, stellar_abundances,
    wolf_rayet)
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

PHASE_NAMES = {"horizontal-branch": "Horizontal Branch Giant", "asymptotic-giant": "Asymptotic Giant",
               "core-helium-burning": "Core Helium Burning Giant", "classical-cepheid": "Classical Cepheid"}
"""Display name of a giant's phase, in place of its class's name, in `Star.type`: a horizontal-branch star, a
core-helium-burning giant and an asymptotic giant are all class III (II from `BRIGHT_GIANT_MIN_MASS_SOL`); a star on the
track inside the instability strip is a classical Cepheid (`classical_cepheid`)."""

HORIZONTAL_BRANCH_COLOURS = ("blue", "rr-lyrae", "red", "extreme")
"""The horizontal branch's colours, in the order of `tuning.HORIZONTAL_BRANCH_COLOURS`' shares."""


def main_sequence_lifetime_gy(mass_sol, fe_h=0.0):
    """The main-sequence lifetime, Gy, of a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` (default 0, the
    reference composition), from `ms_lifetime` (Hurley, Pols and Tout 2000, GEN.467), the one the rest of the
    generator uses."""
    return ms_lifetime.lifetime_gy(mass_sol, fe_h)


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


@functools.lru_cache(maxsize=1 << 19)
def main_sequence_bright_from_gy(mass_sol, min_luminosity_sol, t_ms_gy, fe_h=0.0):
    """
    The age (Gy) a main-sequence star of initial mass `mass_sol` and birth [Fe/H] `fe_h` first shines at least
    `min_luminosity_sol` bright, a main sequence that brightens with age (GEN.652): 0.0 when it is that bright on
    arrival, `None` when it is not even as it leaves (`t_ms_gy`, the lifetime).
    """
    # The two edges are read from the tables (`main_sequence.zams_state`, `end_state`), which settles all but the few
    # masses whose main sequence crosses the brightness; only those are searched.
    if main_sequence.end_state(mass_sol, fe_h).luminosity_sol < min_luminosity_sol:
        return None
    if main_sequence.zams_state(mass_sol, fe_h).luminosity_sol >= min_luminosity_sol:
        return 0.0
    # Regula falsi with the Illinois step on log luminosity, which is smooth and rising in age: a few rows of the
    # table, not sixteen halvings of it.
    target = math.log(min_luminosity_sol)
    low, high = 0.0, t_ms_gy
    f_high = math.log(main_sequence.state(mass_sol, high, fe_h).luminosity_sol) - target
    if f_high < 0.0:
        return None
    f_low = math.log(main_sequence.state(mass_sol, low, fe_h).luminosity_sol) - target
    if f_low >= 0.0:
        return 0.0
    side = 0
    for _ in range(MAIN_SEQUENCE_BRIGHT_ITERATIONS):
        mid = (low * f_high - high * f_low) / (f_high - f_low)
        f_mid = math.log(main_sequence.state(mass_sol, mid, fe_h).luminosity_sol) - target
        if abs(f_mid) < MAIN_SEQUENCE_BRIGHT_TOLERANCE:
            return mid
        if f_mid < 0.0:
            low, f_low = mid, f_mid
            if side == -1:
                f_high *= 0.5
            side = -1
        else:
            high, f_high = mid, f_mid
            if side == 1:
                f_low *= 0.5
            side = 1
    return high


MAIN_SEQUENCE_BRIGHT_ITERATIONS = 12
"""int: Steps at most when the main sequence is searched for the age a star reaches a brightness."""

MAIN_SEQUENCE_BRIGHT_TOLERANCE = 3.0e-4
"""float: How close (in the natural log of the luminosity) the age found is to the brightness asked for."""


def map_luminosity_sol(mass_sol):
    """
    The luminosity (Lsun) a main-sequence star of `mass_sol` would have by the
    map's relation (`MAP_MASS_LUMINOSITY_PIECES`, GEN.200): the one function
    that gives a faint but massive object (a neutron star, a black hole) the
    visibility on the Galaxy Map of the star of its mass, so major
    gravitational places are always on the map. It is not the star
    generators' relation (`main_sequence.state`, GEN.652): the two differ
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
    the Galaxy Map's nodes and the Sector Map's lights both start from it.

    Args:
        map_type (str): "black_hole", "neutron_star", "quasar" or anything else.
        mass_sol (float, optional): The mass, solar masses.

    Returns:
        float: The luminosity, solar luminosities.
    """
    mass = mass_sol if mass_sol and mass_sol > 0 else tuning.MASSLESS_MARKER_MASS_SOL
    share = tuning.BLACK_HOLE_MAP_BRIGHTNESS if map_type == "black_hole" else 1.0
    return map_luminosity_sol(mass) * share


def effective_temperature_k(luminosity_sol, radius_sol):
    """Stefan-Boltzmann in solar units: T = T_sun * (L / R^2)^(1/4)."""
    return tuning.SUN_EFFECTIVE_TEMPERATURE_K * (luminosity_sol / radius_sol ** 2) ** 0.25


def white_dwarf_radius_km(mass_sol):
    """The generator's existing white dwarf mass-radius relation."""
    return constants.WHITE_DWARF_BASE_RADIUS_KM * mass_sol ** constants.WHITE_DWARF_MASS_RADIUS_EXPONENT


def spectral_letter_and_subclass(temperature_k):
    """
    The spectral letter and the 0-9 subclass (0 hottest) of the type whose temperature is nearest `temperature_k` in the Pecaut and
    Mamajek dwarf sequence (`datasets.mamajek_sequence`, GEN.943): the Sun's 5,770 K is G2, 3,900 K is the K/M boundary, and a
    temperature above O0 or below M9 is O0 or M9.
    """
    return mamajek_sequence.spectral_letter_and_subclass(temperature_k)


def is_subdwarf(luminosity_sol, temperature_k):
    """
    Whether a main-sequence star is a subdwarf (Yerkes VI, GEN.943): it lies `SUBDWARF_LUMINOSITY_DEFICIT_DEX` or more below the
    luminosity of a solar-metallicity dwarf of its temperature (the Pecaut and Mamajek sequence), as a metal-poor halo star does.
    """
    return math.log10(luminosity_sol) < mamajek_sequence.dwarf_log_luminosity(temperature_k) - tuning.SUBDWARF_LUMINOSITY_DEFICIT_DEX


def yerkes_class_of(state):
    """The Yerkes class a model state is stored under: its own, except that a main-sequence star under-luminous for its temperature
    is a subdwarf (`VI`)."""
    if state["yerkes_class"] == "V" and is_subdwarf(state["luminosity_sol"], state["temperature_k"]):
        return "VI"
    return state["yerkes_class"]


def luminosity_law_share_above(low, high, slope, min_luminosity_sol):
    """The part of the time spent between `low` and `high` (a time at each luminosity that falls as L^-slope per unit of
    ln L; log-uniform for slope 0) spent at or above the value; 1.0 for `None`."""
    if min_luminosity_sol is None or min_luminosity_sol <= low:
        return 1.0
    if min_luminosity_sol >= high:
        return 0.0
    if slope == 0.0:
        return math.log(high / min_luminosity_sol) / math.log(high / low)
    return ((min_luminosity_sol ** -slope - high ** -slope) / (low ** -slope - high ** -slope))


def luminosity_law_quantile(low, high, slope, fraction, min_luminosity_sol=None):
    """The luminosity below which `fraction` (0 to 1) of the time spent at or above `min_luminosity_sol` (`None`:
    from `low`) is spent."""
    if min_luminosity_sol is not None:
        low = min(max(low, min_luminosity_sol), high)
    fraction = min(max(fraction, 0.0), 1.0)
    if slope == 0.0 or low >= high:
        return math.exp(math.log(low) + (math.log(high) - math.log(low)) * fraction)
    return (low ** -slope - fraction * (low ** -slope - high ** -slope)) ** (-1.0 / slope)


def horizontal_branch_bright_chance(min_luminosity_sol):
    """The chance a horizontal-branch star (one luminosity range, whatever
    its colour) is at least `min_luminosity_sol` bright."""
    return luminosity_law_share_above(*tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL, 0.0, min_luminosity_sol)


def has_core_helium_phases(mass_sol):
    """Whether a star of initial mass `mass_sol` has a horizontal branch and an asymptotic giant branch of the lighter
    stars' kind (its helium ignites in a flash): below `HELIUM_FLASH_MASS_SOL`."""
    return post_main_sequence.has_core_helium_phases(mass_sol)


def remnant_age_fraction(mass_sol, fe_h=0.0):
    """The age at which a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` ends its last
    giant phase, as a multiple of its main-sequence lifetime: the end of its asymptotic giant branch, or for a star
    that does not reach one the end of its core helium burning (`post_main_sequence.timeline`)."""
    line = post_main_sequence.timeline(mass_sol, fe_h)
    return line.end_gy / line.main_sequence_gy


def remnant_formation_time_gy(mass_sol, fe_h=0.0):
    """The age at which a star of initial mass `mass_sol` and birth [Fe/H] `fe_h` has become its remnant, Gy: the end
    of its last giant phase."""
    return post_main_sequence.timeline(mass_sol, fe_h).end_gy


_LONGEST_POST_MAIN_SEQUENCE = 1.5
"""float: No star's life after its main sequence is longer than this share of the main-sequence lifetime (the intermediate-mass
stars' is 1.2 to 1.3), which brackets the search of `youngest_progenitor_sol`."""


def youngest_progenitor_sol(age_yr, fe_h=0.0, settled=None):
    """The lightest star of birth [Fe/H] `fe_h` whose giant phases have ended by this age (a star lighter than
    `HELIUM_FLASH_MASS_SOL` also burns helium on the horizontal branch and climbs the asymptotic giant branch, and its
    giant branch has the length of `giant_branch`, GEN.502).

    `settled(low, high)`, when given, is asked before each halving of the search whether the answer is already as good as
    it needs to be, and the search then ends with `high`: a reader that only sorts the answer among a fixed list of masses
    stops as soon as no mass of the list lies between the bounds (the answer is somewhere in them), after a few steps
    instead of fifty. Without it the search runs on to the end."""
    age_gy = age_yr / 1.0e9
    low = max(ms_lifetime.mass_for_lifetime_gy(age_gy, fe_h), ms_lifetime.MIN_MASS_SOL)
    high = max(ms_lifetime.mass_for_lifetime_gy(age_gy / _LONGEST_POST_MAIN_SEQUENCE, fe_h), low)
    if remnant_formation_time_gy(low, fe_h) <= age_gy:
        return low
    for _ in range(50):
        if settled is not None and settled(low, high):
            return high
        middle = math.sqrt(low * high)
        if remnant_formation_time_gy(middle, fe_h) > age_gy:
            low = middle
        else:
            high = middle
    return high


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
    `fe_h_peaks`, a cluster's or satellite's own distribution, GEN.662) at `galactic_center_dist_ly` from the centre (the
    Sun's distance without one) by `metallicity.sample_fe_h`. A cluster or satellite star with no distribution of its own
    takes the generic cluster's or satellite's (`metallicity.generic_peaks`), never the halo's.
    """
    distance = scales.default_galactocentric_ly() if galactic_center_dist_ly is None else galactic_center_dist_ly
    return metallicity.sample_fe_h(population, distance, rng, peaks=fe_h_peaks or metallicity.generic_peaks(population))


def classical_cepheid(state, initial_mass_sol):
    """Whether a star on its track (a `post_main_sequence.State`) born with `initial_mass_sol` is a classical Cepheid: a
    giant or supergiant of at most `CEPHEID_MAX_INITIAL_MASS_SOL` that is luminous enough (`CEPHEID_MIN_LUMINOSITY_SOL`)
    and inside the instability strip (`CEPHEID_TEMPERATURE_RANGE_K`)."""
    low, high = tuning.CEPHEID_TEMPERATURE_RANGE_K
    return (state.window not in (post_main_sequence.HORIZONTAL_BRANCH, post_main_sequence.SUBGIANT)
            and initial_mass_sol <= tuning.CEPHEID_MAX_INITIAL_MASS_SOL
            and state.luminosity_sol >= tuning.CEPHEID_MIN_LUMINOSITY_SOL and low <= state.temperature_k <= high)


def _evolved_state(line, age_gy, mass_sol, fe_h, rng):
    """
    The class, luminosity, temperature and phase of a star past its main sequence at `age_gy`, from its
    `post_main_sequence.Timeline`: the track's own values, except on the horizontal branch of a star lighter than
    `HELIUM_FLASH_MASS_SOL`, whose colour is drawn from its [Fe/H] and age and its luminosity from the branch's range.
    """
    pc = tuning
    place = line.state(age_gy)
    window = place.window
    if window == post_main_sequence.HORIZONTAL_BRANCH:
        colour = sample_horizontal_branch_colour(fe_h, age_gy, rng)
        return dict(yerkes_class="III", phase="horizontal-branch", mass_sol=place.mass_sol,
                    luminosity_sol=luminosity_law_quantile(*pc.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL, 0.0, rng.random()),
                    temperature_k=rng.uniform(*pc.HORIZONTAL_BRANCH_COLOUR_TEMPERATURE_RANGE_K[colour]))
    found = dict(luminosity_sol=place.luminosity_sol, temperature_k=place.temperature_k, mass_sol=place.mass_sol)
    if window == post_main_sequence.SUPERGIANT:
        found["yerkes_class"] = _supergiant_yerkes_class(place.luminosity_sol)
    elif window == post_main_sequence.SUBGIANT:
        found["yerkes_class"] = "IV"
    else:
        found["yerkes_class"] = "II" if mass_sol >= pc.BRIGHT_GIANT_MIN_MASS_SOL else "III"
    if window == post_main_sequence.CORE_HELIUM_BURNING:
        found["phase"] = "core-helium-burning"
    elif window == post_main_sequence.ASYMPTOTIC_GIANT:
        found["phase"] = "asymptotic-giant"
    if classical_cepheid(place, mass_sol):
        found["phase"] = "classical-cepheid"
    return found


def evolve_star(mass_sol, age_gy, rng=draw, min_luminosity_sol=None, fe_h=None, population=None, fe_h_peaks=None,
                galactic_center_dist_ly=None, white_dwarf_mass_sol=None):
    """
    The present state of a star of initial mass `mass_sol` at `age_gy`.

    After its main sequence the star's luminosity and temperature are those of its MIST track at that age
    (`post_main_sequence`, GEN.935): one point for one age, not a draw. `min_luminosity_sol` is for the one luminosity that
    is drawn (a Wolf-Rayet star's, from the part of its range at or above that value; `star_population.sample_bright_stars`
    uses it once it has picked a mass and age that can reach it).

    The star's birth [Fe/H] sets its lifetimes and its tracks (`ms_lifetime`, `giant_branch`, `post_main_sequence`, GEN.651)
    and a horizontal-branch star's colour. `fe_h` is the star's own; without it one is drawn first, before any other draw
    here (`birth_fe_h`: the star's `population` or `fe_h_peaks`, at `galactic_center_dist_ly` from the centre, the Sun's
    distance without one), and either way returned as `metallicity_fe_h` for the star to keep.

    Returns:
        dict or None: None when the star has already collapsed to a neutron
        star or black hole (initial mass at or above
        `SUPERGIANT_MIN_MASS_SOL`, past its giant phase). Otherwise:
        `yerkes_class`, `mass_sol` (present mass: MIST's on the main sequence and, with the red giant
        branch's loss of `tuning.RED_GIANT_BRANCH_REIMERS_ETA`, after it, GEN.940; a white dwarf's final mass), `luminosity_sol`, `temperature_k`, `radius_km` (a white
        dwarf's from its mass-radius relation, None for every other class,
        whose radius follows from luminosity and temperature),
        `metallicity_fe_h` (the birth [Fe/H]), `lifespan_gy` (the age at which it becomes a remnant; infinite for a
        white dwarf) and `phase_end_gy` (the age at which it leaves its
        present phase; infinite for a white dwarf). A horizontal-branch star, a core-helium-burning giant,
        an asymptotic giant or a Cepheid (class III or II) also has `phase`
        (`"horizontal-branch"`, `"core-helium-burning"`, `"asymptotic-giant"` or `"classical-cepheid"`, a key of
        `PHASE_NAMES`). A star younger than its contraction time
        (`pre_main_sequence.state`, GEN.242) is a `"protostar"` (under 0.5
        Myr) or a `"PMS"` (pre-main-sequence) class, with `disk_class` ("0",
        "I", "II" or "III") in place of a Yerkes numeral, and `phase_end_gy`
        the age it reaches the main sequence.

    `white_dwarf_mass_sol` is the mass of the white dwarf the star leaves when a close companion stripped its giant of its envelope
    on the red giant branch (a helium core, GEN.996, `helium_core_white_dwarf`), in place of the initial-final relation's; it
    applies only to a star that has become a white dwarf.
    """
    pc = tuning
    if fe_h is None:
        fe_h = birth_fe_h(rng, population, fe_h_peaks, galactic_center_dist_ly)
    line = post_main_sequence.timeline(mass_sol, fe_h)
    t_ms, t_end = line.main_sequence_gy, line.end_gy
    massive = mass_sol >= pc.SUPERGIANT_MIN_MASS_SOL
    state = {"mass_sol": mass_sol, "radius_km": None, "lifespan_gy": t_end, "metallicity_fe_h": fe_h}

    wr_start = wolf_rayet.start_gy(mass_sol, t_ms, t_end, fe_h) if massive else None
    young = pre_main_sequence.state(mass_sol, age_gy * 1.0e9, fe_h)
    if young is not None:
        # GEN.242: still contracting, a protostar (Class 0 and I) or a pre-main-sequence star.
        state.update(young)
    elif age_gy < t_ms:
        on_sequence = main_sequence.state(mass_sol, age_gy, fe_h)
        state.update(yerkes_class="V", luminosity_sol=on_sequence.luminosity_sol, mass_sol=on_sequence.mass_sol,
                     temperature_k=on_sequence.temperature_k, phase_end_gy=t_ms)
    elif wr_start is not None and wr_start <= age_gy < t_end:
        # GEN.243: the last stretch of a very massive star's life, in place of the end of its supergiant window.
        present = line.state(age_gy).mass_sol
        state.update(wolf_rayet.state(mass_sol, present, age_gy, wr_start, t_end, rng, min_luminosity_sol))
    elif age_gy < t_end:
        state.update(_evolved_state(line, age_gy, mass_sol, fe_h, rng))
        state["phase_end_gy"] = line.window_at(age_gy).end_gy
    elif massive:
        return None
    else:
        final_mass = remnant_mass.white_dwarf_mass_sol(mass_sol) if white_dwarf_mass_sol is None else white_dwarf_mass_sol
        cooling_age = cooling_age_gy(mass_sol, age_gy, fe_h)
        state.update(_white_dwarf_state(final_mass, cooling_age, rng))
    return state


def cooling_age_gy(initial_mass_sol, age_gy, fe_h):
    """How long, Gy, a white dwarf of this progenitor has been cooling at `age_gy`: the age less the end of the progenitor's life, at least
    `tuning.WD_MIN_COOLING_AGE_GY`."""
    return max(age_gy - post_main_sequence.timeline(initial_mass_sol, fe_h).end_gy, tuning.WD_MIN_COOLING_AGE_GY)


def helium_core_state(initial_mass_sol, age_gy, fe_h, core_mass_sol, rng):
    """The state fields of the white dwarf a star of `initial_mass_sol` solar masses has become at `age_gy` when it is a helium core of
    `core_mass_sol` solar masses (GEN.996): `evolve_star`'s white dwarf fields (mass, luminosity, temperature, radius, atmosphere), to
    update the state of the same star with. It cools from the end of its progenitor's life like any other white dwarf of it."""
    return _white_dwarf_state(core_mass_sol, cooling_age_gy(initial_mass_sol, age_gy, fe_h), rng)


def _white_dwarf_state(final_mass_sol, cooling_age_gy, rng):
    """
    A white dwarf of `final_mass_sol` that began cooling `cooling_age_gy` ago (GEN.742): the Bedard et al. 2020 cooling sequences
    (`datasets.white_dwarf_models.cooling_state`) give its effective temperature, luminosity, radius and the crystallised share of its
    core, on a hydrogen atmosphere (a thick hydrogen layer, `tuning.STELLAR_WD_DA_SHARE` of white dwarfs) or a helium one (thin).
    When the sequences are not downloaded and the run disregards missing data (GEN.599) it is the Mestel power law of the cooling age
    and the Chandrasekhar-type radius instead, labelled in `fetch.stand_ins_used`.
    """
    if fetch.locate(white_dwarf_models.DATASET).use_stand_in("the Mestel power law of the white-dwarf cooling"):
        luminosity = tuning.WD_COOLING_L0_SOL * (final_mass_sol / 0.6) * cooling_age_gy ** tuning.WD_COOLING_EXPONENT
        luminosity = min(max(luminosity, tuning.WD_LUMINOSITY_RANGE_SOL[0]), tuning.WD_LUMINOSITY_RANGE_SOL[1])
        radius_km = white_dwarf_radius_km(final_mass_sol)
        radius_sol = _km_to_m(radius_km) / constants.SOLAR_RADIUS_M
        return dict(yerkes_class="VII", mass_sol=final_mass_sol, luminosity_sol=luminosity,
                    temperature_k=effective_temperature_k(luminosity, radius_sol), radius_km=radius_km,
                    lifespan_gy=float("inf"), phase_end_gy=float("inf"))
    hydrogen = rng.random() < tuning.STELLAR_WD_DA_SHARE
    cooled = white_dwarf_models.cooling_state(final_mass_sol, cooling_age_gy, "thick" if hydrogen else "thin")
    return dict(yerkes_class="VII", mass_sol=final_mass_sol, luminosity_sol=cooled.luminosity_sol, temperature_k=cooled.teff_k,
                radius_km=cooled.radius_solar * constants.SOLAR_RADIUS_M / 1000.0, lifespan_gy=float("inf"),
                phase_end_gy=float("inf"), crystallised_fraction=cooled.crystallised_fraction,
                white_dwarf_atmosphere="DA" if hydrogen else "non-DA")


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
            mass = imf.sample(min_mass_sol=pc.LARGE_STAR_MIN_MASS_SOL, rng=rng, curve=imf.PRIMARY)
            if max_mass_sol is not None and mass >= max_mass_sol:
                draw.redrawn("star.model_redraw")
                continue
            low, high = population_age_range_gy(population)
            high = min(high, remnant_formation_time_gy(mass, birth))
            if high <= low:
                draw.redrawn("star.model_redraw")
                continue  # an old population has no living star this massive
            low, high = age_window_gy(low, high, age_bias)
            age = age_metallicity.tilted_age(
                birth, population, lambda: rng.uniform(low, high) if age_law is None else age_law.sample(low, high, rng), rng)
        else:
            mass = imf.sample(rng=rng, curve=imf.PRIMARY)
            if mass < stellar_abundances.hydrogen_burning_limit_for_fe_h(birth):
                continue    # GEN.936: under the limit of its place it is a brown dwarf, not a star a system is built around
            if large_star is False and mass >= pc.LARGE_STAR_MIN_MASS_SOL:
                draw.redrawn("star.model_redraw")
                continue
            if max_mass_sol is not None and mass >= max_mass_sol:
                draw.redrawn("star.model_redraw")
                continue
            age = age_metallicity.tilted_age(birth, population, lambda: sample_star_age_gy(age_bias, rng, population, age_law), rng)
        if habitable_host and age < pc.LIFE_MIN_STAR_AGE_GY:
            draw.redrawn("star.model_redraw")
            continue
        state = evolve_star(mass, age, rng, fe_h=birth, population=population, fe_h_peaks=fe_h_peaks)
        if state is None or (habitable_host and state["yerkes_class"] == "VII"):
            draw.redrawn("star.model_redraw")
            continue
        # A white dwarf is never pre-placed, so the cap never removes one
        # (the hottest are clamped to exactly the lowest allowed cap).
        if (max_luminosity_sol is not None and state["yerkes_class"] != "VII"
                and state["luminosity_sol"] >= max_luminosity_sol):
            draw.redrawn("star.model_redraw")
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
    mass-radius relation, and then the luminosity follows from the radius and
    the stored temperature instead, GEN.365).

    Returns:
        dict: `type` (e.g. "G2V Yellow Main Sequence Star"),
        `yerkes_class`, `mass_kg`, `radius_km`, `temperature_k`,
        `luminosity_w`, `age_gy`, `lifespan_gy`, `initial_mass_sol` and
        `phase_end_age_gy`; and `metallicity_fe_h`, the birth [Fe/H] the star keeps (GEN.651).
    """
    temperature = int(round(state["temperature_k"], tuning.ROUND_TEMPERATURE_NEAREST_HUNDRED))
    spectral_class, subclass = spectral_letter_and_subclass(state["temperature_k"])
    luminosity_w = state["luminosity_sol"] * constants.SOLAR_LUMINOSITY
    if state["radius_km"] is not None:
        radius_km = state["radius_km"]
        # GEN.365: a white dwarf's radius is the degenerate star's (its mass sets it) and its temperature the one stored (rounded, or the
        # forced type's), so its luminosity is what those two give, and L = 4 pi R^2 sigma T^4 holds for the stored numbers.
        luminosity_w = constants.FOUR_PI * (radius_km * 1000.0) ** 2 * constants.STEFAN_BOLTZMANN_CONSTANT * temperature ** 4
    else:
        radius_km = math.sqrt(luminosity_w / (constants.FOUR_PI * constants.STEFAN_BOLTZMANN_CONSTANT
                                              * temperature ** 4))
        radius_km = _m_to_km(radius_km)
    yerkes = yerkes_class_of(state)
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
        "mass_kg": _msun_to_kg(state["mass_sol"]),
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

_msun_to_kg = units.converter("Msun", "kg")
