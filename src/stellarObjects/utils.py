# stellarObjects/utils.py

"""
Utilities
=========

This module contains utility functions for the planetGen package, including
mathematical calculations, string formatting, and name generation. These
functions support various aspects of the celestial body generation process,
from calculating physical properties to creating unique and plausible names.
"""

import functools
import inspect
import math
import random

from planetgen.generation.config import SystemConfig
from planetgen.names.wordlists import (
    BAD_CONSONANTS, COMPANION_SUFFIXES, DICTIONARY_WORDS, DIMINUTIVE_PREFIXES, GREEK_LETTERS,
    NSFW_WORDS, ROMAN_NUMERALS_BY_VALUE, SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES,
    UNIVERSAL_PHONEMES, VOWELS, WORD_SIZE_MEAN,
)
from planetgen.physics import constants as physical_constants
from planetgen import tuning

def _numeric_leaves(value):
    if isinstance(value, (tuple, list)):
        for item in value:
            yield from _numeric_leaves(item)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield value


def finite_domain(*, allow_inf=(), clamped=()):
    """
    Decorator giving a numeric physics helper one explicit contract: it
    returns finite numbers or raises `ValueError` -- never a
    `ZeroDivisionError`/`OverflowError` from an out-of-range input, and
    never a silent NaN/inf that would otherwise surface far away (in a
    stored row, a rendered page).

    * Every float argument must be finite, except those named in
      `allow_inf` (which may be +/-inf, e.g. `vis_viva_speed_kms`'s
      parabolic `semi_major_axis_au = math.inf`) and those named in
      `clamped` (which the helper clamps into a valid range itself, so any
      non-NaN value is fine; NaN still surfaces through the result check).
    * A `ZeroDivisionError`/`OverflowError` raised inside (a zero or
      subnormal divisor, a power overflowing) becomes a `ValueError`.
    * A non-finite number anywhere in the result becomes a `ValueError`.

    Plain unit conversions and formatting helpers deliberately don't use
    this: they pass non-finite values through unchanged.
    """
    allow_inf = frozenset(allow_inf)
    clamped = frozenset(clamped)

    def decorate(func):
        params = list(inspect.signature(func).parameters)
        where = func.__name__

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for name, value in (*zip(params, args), *kwargs.items()):
                if isinstance(value, float) and not math.isfinite(value):
                    if name in clamped or (name in allow_inf and math.isinf(value)):
                        continue
                    raise ValueError(f"{where}: {name} must be a finite number, got {value!r}")
            try:
                result = func(*args, **kwargs)
            except (ZeroDivisionError, OverflowError) as exc:
                raise ValueError(f"{where}: input out of range ({exc})") from exc
            for leaf in _numeric_leaves(result):
                if not math.isfinite(leaf):
                    raise ValueError(f"{where}: result out of range ({leaf!r}) for {args!r}")
            return result

        return wrapper

    return decorate


def get_star_spectral_class(star):
    """
    Returns the uppercase spectral class character (e.g. 'G') for a star.

    Works for both a plain `Star` and a `BinaryStarProxy`, without importing
    `BinaryStarProxy` here (which would create a circular import) — a proxy
    is detected by duck-typing on its `_primary` attribute, and its primary
    star's type is used as the representative spectral class.

    Args:
        star: A `Star` or `BinaryStarProxy` instance.

    Returns:
        str: The uppercase spectral class character.
    """
    reference_star = star._primary if hasattr(star, '_primary') else star
    return reference_star.type[0].upper()


def get_star_evolutionary_profile(star):
    """
    Returns the `STAR_EVOLUTION`-shaped profile describing which life
    chemicals and evolutionary paces are plausible for planets orbiting
    `star`, correctly accounting for its Yerkes luminosity class.

    `STAR_EVOLUTION` is keyed by spectral letter (e.g. 'G'), which for a
    main-sequence (Yerkes 'V') star fully determines both its current
    temperature and its total lifespan -- so that class's fixed entry
    applies directly.

    For any evolved or remnant star (giants, supergiants, subgiants, bright
    giants, hypergiants, subdwarfs, white dwarfs), the spectral letter only
    reflects the star's *current* temperature/color, not a lifespan -- e.g.
    an "F5VII" white dwarf merely glows at F-like temperature; it did not
    live and die as an F-type main-sequence star. Using STAR_EVOLUTION['F']
    directly for such a star would misapply a main-sequence lifespan/scale
    table to an unrelated evolutionary history. `potentially_viable_chemicals`
    (driven by the star's current emission spectrum) is still taken from the
    current-temperature letter, since that part genuinely is about current
    color. `supported_evolutionary_scales`, which is really a proxy for how
    much time was available for a biosphere to develop, is instead derived
    from the star's own already-computed age/lifespan (already correctly
    yerkes-class-aware -- see `Star._calculate_initial_star_age_and_lifespan`),
    by checking which evolutionary paces could reach their technological
    civilization milestone within that time budget.

    Args:
        star: A `Star` or `BinaryStarProxy` instance.

    Returns:
        dict: A `STAR_EVOLUTION`-entry-shaped dict (with at least
              `potentially_viable_chemicals` and `supported_evolutionary_scales`
              keys), or `{}` if the spectral letter has no entry.
    """
    reference_star = star._primary if hasattr(star, '_primary') else star
    spectral_class_char = reference_star.type[0].upper()
    base_info = tuning.STAR_EVOLUTION.get(spectral_class_char, {})
    if not base_info:
        return {}

    if reference_star.yerkes_class == "V":
        return base_info

    # A white dwarf's lifespan is infinite (it just cools forever), so "how
    # much time has been available so far" is better represented by its age.
    # Every other evolved class has a finite lifespan, used as-is.
    time_budget = reference_star.age if reference_star.lifespan == float('inf') else reference_star.lifespan

    reachable_scales = [
        scale for scale in ["fast", "normal", "slow"]
        if tuning.EVOLUTIONARY_TIMELINES[scale]['technological_civilization'] <= time_budget
    ]
    if not reachable_scales:
        # Even the fastest pace doesn't fit -- still return it so callers have
        # something to work with, mirroring the fallback in get_evolutionary_timeline.
        reachable_scales = ["fast"]

    return {**base_info, "supported_evolutionary_scales": reachable_scales}


# The distance ladder, smallest first: (label, meters). `format_distance`
# shows a value in the largest unit it is at least 1 of.
DISTANCE_LADDER = (
    ("km", physical_constants.KM_M),
    ("AU", physical_constants.AU_M),
    ("mpc", physical_constants.MILLIPARSEC_M),
    ("cpc", physical_constants.CENTIPARSEC_M),
    ("ly", physical_constants.LIGHTYEAR_M),
    ("pc", physical_constants.PARSEC_M),
    ("kpc", physical_constants.KILOPARSEC_M),
    ("Mpc", physical_constants.MEGAPARSEC_M),
    ("Gpc", physical_constants.GIGAPARSEC_M),
)
# Units that carry a second, more familiar unit in parentheses.
PARSEC_UNITS = frozenset({"mpc", "cpc", "pc", "kpc", "Mpc", "Gpc"})
# A parsec value's parenthetical is in ly when the distance is at least
# this many ly, else in AU when at least this many AU, else in km (Boss,
# 2026-09-30: "4.2 pc (13.7 ly)", "2.4 mpc (495 AU)").
DISTANCE_PAREN_MIN_LY = 0.01
DISTANCE_PAREN_MIN_AU = 0.01


SCIENTIFIC_MIN_INTEGER_DIGITS = 7
"""int: A number shown with no decimals and this many digits or more is
shown in scientific notation instead (UX.36, Boss 2026-10-02: whole
numbers stay plain past the 6th digit)."""

SCIENTIFIC_MIN_DECIMAL_INTEGER_DIGITS = 5
"""int: A number shown with decimals and this many digits before the
decimal point or more is shown in scientific notation instead (UX.20,
kept by UX.36: "the 4th when decimals are used")."""

SCIENTIFIC_SIGNIFICANT_FIGURES = 3
"""int: Significant figures in `scientific_text`'s mantissa."""

_SUPERSCRIPT_DIGITS = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")


def scientific_text(value, significant=SCIENTIFIC_SIGNIFICANT_FIGURES):
    """
    `value` in scientific notation as plain text with Unicode superscripts,
    "1.23 × 10⁶", which reads right in HTML, in a `data-*` attribute read
    back with `.textContent`, in Markdown and in wikitext alike.
    """
    mantissa, exponent = f"{value:.{significant - 1}e}".split("e")
    return f"{mantissa} \u00d7 10{str(int(exponent)).translate(_SUPERSCRIPT_DIGITS)}"


def format_number(value, spec=",.0f"):
    """
    The site's one number formatter (UX.20): `value` formatted with `spec`
    (a `format()` spec, default whole and comma-grouped), unless that shows
    too many digits before the decimal point (`_shows_too_many_digits`),
    when it is `scientific_text(value)` instead. Counts and measurements
    alike go through it; IDs, years in dates, designations, page numbers
    and coordinates don't. `static/numberformat.js` mirrors it.

    Returns:
        str: e.g. "9,999", "1.5", "1.23 × 10⁴". NaN and infinities are
             formatted as-is.
    """
    text = format(value, spec)
    if not math.isfinite(value):
        return text
    if "e" in text or _shows_too_many_digits(text):
        return scientific_text(value)
    return text


def _shows_too_many_digits(text):
    """True when formatted `text` should be scientific instead: 7 or more
    whole digits with no decimals shown, 5 or more with decimals (UX.36)."""
    whole, _point, decimals = text.lstrip("-+").partition(".")
    limit = SCIENTIFIC_MIN_DECIMAL_INTEGER_DIGITS if decimals else SCIENTIFIC_MIN_INTEGER_DIGITS
    return len(whole.replace(",", "")) >= limit


def _three_figures(value):
    """`value` to three significant figures, comma-grouped, trailing zeros
    dropped: 4.2, 13.7, 0.499, 495, 12,300."""
    if value == 0 or not math.isfinite(value):
        return f"{value:g}"
    magnitude = math.floor(math.log10(abs(value)))
    decimals = max(0, 2 - magnitude)
    text = f"{round(value, decimals):,.{decimals}f}"
    if _shows_too_many_digits(text):
        return scientific_text(value)
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


def _format_in_unit(meters, label, unit_m):
    if label == "km" and abs(meters) >= 1e6:
        # Whole kilometers read better than three figures ("384,400 km").
        return f"{format_number(round(meters / unit_m))} km"
    return f"{_three_figures(meters / unit_m)} {label}"


def format_distance_m(meters):
    """
    Formats a distance (or a non-body radius) in the most meaningful unit on
    the ladder km < AU < mpc < cpc < ly < pc < kpc < Mpc < Gpc: the largest
    unit the value is at least 1 of. Values in a parsec unit add one
    parenthetical: ly when the distance is at least `DISTANCE_PAREN_MIN_LY`,
    else AU when at least `DISTANCE_PAREN_MIN_AU`, else km.

    Every page and text output passes its distances through this (or its
    `format_distance_km`/`_au`/`_ly`/`_pc` wrappers) rather than formatting
    them itself; `static/distance.js` mirrors it for the maps. Planet, moon
    and star radii are the exception: see `format_body_radius_km`.

    Args:
        meters (float): The distance, in meters.

    Returns:
        str: e.g. "384,400 km", "1.52 AU", "4.2 pc (13.7 ly)",
             "2.4 mpc (495 AU)". `None` gives an en dash.
    """
    if meters is None:
        return "\u2013"
    meters = float(meters)
    if not math.isfinite(meters):
        return f"{meters:g} km"
    size = abs(meters)
    label, unit_m = DISTANCE_LADDER[0]
    for candidate_label, candidate_m in DISTANCE_LADDER:
        # A hair of tolerance so a value computed as exactly 1 unit (e.g.
        # 1 pc from ly * LY_TO_AU / AU_PER_PARSEC) isn't shown as 1,000 mpc.
        if size >= candidate_m * (1 - 1e-12):
            label, unit_m = candidate_label, candidate_m
    text = _format_in_unit(meters, label, unit_m)
    if label in PARSEC_UNITS:
        text += f" ({distance_parenthetical(meters)})"
    return text


def distance_parenthetical(meters):
    """
    The familiar unit a parsec-family value carries in parentheses: ly when
    the distance is at least `DISTANCE_PAREN_MIN_LY`, else AU when at least
    `DISTANCE_PAREN_MIN_AU`, else km. (No parsec value is under 0.01 AU,
    1 mpc being about 206 AU, but the rule is kept as Boss wrote it.)
    """
    size = abs(meters)
    if size >= DISTANCE_PAREN_MIN_LY * physical_constants.LIGHTYEAR_M:
        return _format_in_unit(meters, "ly", physical_constants.LIGHTYEAR_M)
    if size >= DISTANCE_PAREN_MIN_AU * physical_constants.AU_M:
        return _format_in_unit(meters, "AU", physical_constants.AU_M)
    return _format_in_unit(meters, "km", physical_constants.KM_M)


def format_distance_km(km):
    """`format_distance_m` for a value in kilometers (the schema's unit)."""
    return format_distance_m(None if km is None else km * physical_constants.KM_M)


def format_distance_au(au):
    """`format_distance_m` for a value in AU (generation's orbit unit)."""
    return format_distance_m(None if au is None else au * physical_constants.AU_M)


def format_distance_ly(ly):
    """`format_distance_m` for a value in light-years."""
    return format_distance_m(None if ly is None else ly * physical_constants.LIGHTYEAR_M)


def format_distance_pc(pc):
    """`format_distance_m` for a value in parsecs (galaxy geometry)."""
    return format_distance_m(None if pc is None else pc * physical_constants.PARSEC_M)


# The speed ladder (UX.13), slowest first: (label, km/s per unit, shown
# from this many km/s up). Boss, 2026-10-01: "km/h on the low speed end to
# Mm/s on the high speed end", then multiples of c from a tenth of light
# speed. `static/speed.js` mirrors it.
SPEED_LADDER = (
    ("km/h", 1 / 3600, 0.0),
    ("km/s", 1.0, 1.0),
    ("Mm/s", 1e3, 1e3),
    ("c", physical_constants.SPEED_OF_LIGHT_KMS, 0.1 * physical_constants.SPEED_OF_LIGHT_KMS),
)


def format_speed_kms(kms):
    """
    Formats a speed in the most meaningful unit on the ladder km/h < km/s <
    Mm/s < c: km/h below 1 km/s, km/s up to 1,000 km/s, Mm/s up to a tenth
    of light speed, then multiples of c. One unit, three significant
    figures, in `format_distance_m`'s number style.

    Every page and text output passes its plain speeds through this rather
    than formatting them itself; `static/speed.js` mirrors it. Warp and
    fold factors (`navigation.py`) keep their own "x c" columns.

    Args:
        kms (float): The speed, in km/s.

    Returns:
        str: e.g. "36 km/h", "29.8 km/s", "4.5 Mm/s", "0.25 c". `None`
             gives an en dash.
    """
    if kms is None:
        return "–"
    kms = float(kms)
    if not math.isfinite(kms):
        return f"{kms:g} km/s"
    size = abs(kms)
    label, unit_kms, _ = SPEED_LADDER[0]
    for candidate_label, candidate_kms, threshold_kms in SPEED_LADDER:
        if size >= threshold_kms * (1 - 1e-12):
            label, unit_kms = candidate_label, candidate_kms
    return f"{_three_figures(kms / unit_kms)} {label}"


def format_speed_ms(ms):
    """`format_speed_kms` for a value in m/s."""
    return format_speed_kms(None if ms is None else ms / physical_constants.KM_M)


# The time-period ladder (UX.14), shortest first: (plural label, singular
# label, seconds). A Julian year of 365.25 days, as SECONDS_PER_YEAR.
# `static/period.js` mirrors it.
PERIOD_LADDER = (
    ("µs", "µs", 1e-6),
    ("ms", "ms", 1e-3),
    ("s", "s", 1.0),
    ("minutes", "minute", 60.0),
    ("hours", "hour", 3600.0),
    ("days", "day", 86400.0),
    ("years", "year", physical_constants.SECONDS_PER_YEAR),
    ("ky", "ky", physical_constants.SECONDS_PER_YEAR * 1e3),
    ("My", "My", physical_constants.SECONDS_PER_YEAR * 1e6),
    ("Gy", "Gy", physical_constants.SECONDS_PER_YEAR * 1e9),
)


def format_duration_seconds(seconds):
    """
    Formats a time period or duration in the most meaningful unit on the
    ladder µs < ms < s < minutes < hours < days < years < ky < My < Gy: the
    largest unit the value is at least 1 of, as `format_distance_m` does.
    One unit, three significant figures, singular when the shown value is
    exactly 1 ("1 year", "1 day"). `static/period.js` mirrors it.

    Args:
        seconds (float): The duration, in seconds.

    Returns:
        str: e.g. "12.5 ms", "45 minutes", "27.3 days", "1.88 years",
             "236 My". Zero is "0 s"; `None` gives an en dash.
    """
    if seconds is None:
        return "–"
    seconds = float(seconds)
    if not math.isfinite(seconds):
        return f"{seconds:g} years"
    if seconds == 0:
        return "0 s"
    size = abs(seconds)
    plural, singular, unit_s = PERIOD_LADDER[0]
    for candidate in PERIOD_LADDER:
        if size >= candidate[2] * (1 - 1e-12):
            plural, singular, unit_s = candidate
    number = _three_figures(seconds / unit_s)
    return f"{number} {singular if number == '1' else plural}"


def format_period_years(years):
    """`format_duration_seconds` for a value in (Julian) years -- every
    orbital period the generator stores is in years or Gy."""
    return format_duration_seconds(None if years is None else years * physical_constants.SECONDS_PER_YEAR)


def _whole_or_tenths(value):
    """`value` to whole units, or to one decimal when under 10 in size
    ("15", "-40", "5.9", "0"), in `format_number`'s style."""
    text = format_number(value, ",.1f" if abs(value) < 10 else ",.0f")
    if text.endswith(".0"):
        text = text[:-2]
    return "0" if text in ("-0", "+0") else text


def format_temperature_k(kelvin):
    """
    A surface temperature in kelvin with Celsius and Fahrenheit alongside
    (Boss, 2026-10-01: "temperatures should be reported in K, C, and F for
    surface conditions"). Kelvin to three significant figures in
    `format_distance_m`'s number style; °C and °F to whole degrees, or one
    decimal when under 10 in size. Star effective temperatures stay in K
    alone and don't use this.

    Args:
        kelvin (float): The temperature, in K.

    Returns:
        str: e.g. "288 K (15 °C, 59 °F)", "737 K (464 °C, 867 °F)". `None`
             gives an en dash.
    """
    if kelvin is None:
        return "–"
    kelvin = float(kelvin)
    if not math.isfinite(kelvin):
        return f"{kelvin:g} K"
    celsius = kelvin - physical_constants.CELSIUS_ZERO_K
    fahrenheit = celsius * 9 / 5 + 32
    return f"{_three_figures(kelvin)} K ({_whole_or_tenths(celsius)} °C, {_whole_or_tenths(fahrenheit)} °F)"


def _three_figures_or_tiny(value):
    """`_three_figures`, but scientific below a thousandth, where a
    secondary unit would otherwise run to a row of zeros."""
    if value != 0 and abs(value) < 1e-3:
        return scientific_text(value)
    return _three_figures(value)


# The pressure ladder, smallest first: (label, pascals). A pressure is
# shown in the largest unit it is at least 1 of.
PRESSURE_LADDER = (
    ("Pa", 1.0),
    ("kPa", 1e3),
    ("MPa", 1e6),
    ("GPa", 1e9),
)


def format_pressure_pa(pascals):
    """
    A pressure with a metric primary on the Pa < kPa < MPa < GPa ladder and
    customary atm and psi alongside, all to three significant figures
    (Boss, 2026-10-01: "Atmospheric pressure and surface conditions should
    show customary units as well as a secondary to help contextualize the
    metric values given").

    Args:
        pascals (float): The pressure, in Pa (the schema's
            `atmospheric_pressure_pa`).

    Returns:
        str: e.g. "101 kPa (1 atm, 14.7 psi)", "9.2 MPa (90.8 atm,
             1,334 psi)". `None` gives an en dash.
    """
    if pascals is None:
        return "–"
    pascals = float(pascals)
    if not math.isfinite(pascals):
        return f"{pascals:g} Pa"
    size = abs(pascals)
    label, unit_pa = PRESSURE_LADDER[0]
    for candidate_label, candidate_pa in PRESSURE_LADDER:
        if size >= candidate_pa * (1 - 1e-12):
            label, unit_pa = candidate_label, candidate_pa
    atm = _three_figures_or_tiny(pascals / physical_constants.STANDARD_ATMOSPHERE_PA)
    psi = _three_figures_or_tiny(pascals / physical_constants.PSI_PA)
    return f"{_three_figures(pascals / unit_pa)} {label} ({atm} atm, {psi} psi)"


def format_pressure_atm(atm):
    """`format_pressure_pa` for a value in standard atmospheres."""
    return format_pressure_pa(None if atm is None else atm * physical_constants.STANDARD_ATMOSPHERE_PA)


def format_body_radius_km(system_config: SystemConfig, radius_km, precision=None):
    """
    A planet, moon or star radius: always km in scientific notation (Boss,
    2026-09-30), the one exception to `format_distance_m`'s ladder.

    Args:
        system_config (SystemConfig): Picks HTML or wikitext notation.
        radius_km (float): The radius, in kilometers.
        precision (int, optional): Decimal places of the mantissa; defaults
            to `tuning.SCIENTIFIC_NOTATION_DECIMAL_PLACES`.

    Returns:
        str: e.g. "6.371 × 10^3 km".
    """
    if radius_km is None:
        return "\u2013"
    if precision is None:
        precision = tuning.SCIENTIFIC_NOTATION_DECIMAL_PLACES
    return f"{to_scientific_notation(system_config, radius_km, precision)} km"


def format_length_km(system_config: SystemConfig, value_km, threshold, round_digits, scientific_precision=None):
    """
    Formats a length in kilometers, switching between a comma-grouped plain
    number and scientific notation based on a threshold.

    Args:
        system_config (SystemConfig): The shared SystemConfig object.
        value_km (float): The length to format, in kilometers.
        threshold (float): The value above which scientific notation is used.
        round_digits (int): Decimal places for the plain-number form (and the
                            default precision for the scientific-notation form).
        scientific_precision (int, optional): Decimal places for the
                                              scientific-notation form, if
                                              different from `round_digits`.

    Returns:
        str: The formatted length string, including the " km" unit.
    """
    if value_km <= threshold:
        return f"{format_number(round(value_km, round_digits), ',')} km"
    precision = scientific_precision if scientific_precision is not None else round_digits
    return f"{to_scientific_notation(system_config, value_km, precision)} km"


def format_relative_to_sol(system_config: SystemConfig, value, sol_constant, unit, low_percent_precision=2):
    """
    Formats a physical quantity (mass, luminosity, ...) with its value in
    scientific notation alongside a comparison to the Sol reference value:
    a percentage when the ratio is small, or a "×" multiplier when it's large.

    Args:
        system_config (SystemConfig): The shared SystemConfig object.
        value (float): The quantity's value, in SI units.
        sol_constant (float): The Sol-reference value for the same quantity
                              (e.g. `physical_constants.SOLAR_MASS_TO_KG`).
        unit (str): The unit label for the SI value (e.g. "kg", "W").
        low_percent_precision (int, optional): Decimal places used for the
                                               percentage when the ratio is
                                               below `PERCENT_SOL_THRESHOLD_LOW`.
                                               Defaults to 2.

    Returns:
        str: The formatted string, e.g. "1.23 × 10^30 kg (1.00× Sol)".
    """
    sol_val = value / sol_constant
    sci_notation = to_scientific_notation(system_config, value)
    if sol_val < tuning.PERCENT_SOL_THRESHOLD_LOW:
        return f"{sci_notation} {unit} ({sol_val * tuning.PERCENT_MULTIPLIER:.{low_percent_precision}f}% of Sol)"
    elif sol_val < tuning.PERCENT_SOL_THRESHOLD_HIGH:
        return f"{sci_notation} {unit} ({sol_val * tuning.PERCENT_MULTIPLIER:.1f}% of Sol)"
    else:
        return f"{sci_notation} {unit} ({format_number(sol_val, ',.1f')}× Sol)"


def ly_to_milliparsecs(ly):
    """
    Converts a distance in light-years to milliparsecs (mpc).

    For the database persistence layer only -- sector-scale position/
    geometry columns (star_systems.position_x/y/z_mpc, sectors.edge_mpc;
    see planetgen/db/schema.sql) are stored in milliparsecs specifically,
    distinct from every other distance-shaped column in the schema (which
    is kilometers). Generation/physics code keeps its own native light-year
    units for sector geometry throughout and never calls this.

    Args:
        ly (float): The distance in light-years.

    Returns:
        float: The distance in milliparsecs.
    """
    au = ly * physical_constants.LY_TO_AU
    return au / physical_constants.AU_PER_MILLIPARSEC


def milliparsecs_to_ly(mpc):
    """
    Converts a distance in milliparsecs (mpc) back to light-years -- the
    inverse of `ly_to_milliparsecs`, for reconstructing live objects
    (native light-year sector geometry) from database rows.

    Args:
        mpc (float): The distance in milliparsecs.

    Returns:
        float: The distance in light-years.
    """
    au = mpc * physical_constants.AU_PER_MILLIPARSEC
    return au * physical_constants.AU_TO_LY


def mpc_to_pc(mpc):
    """
    Converts a distance in milliparsecs (mpc) to parsecs (pc). Exact --
    milliparsecs and parsecs are the same unit family, `1/1000` apart, so
    this is a single power-of-ten scaling with no AU round-trip (contrast
    `milliparsecs_to_ly`, which does need one).

    Args:
        mpc (float): The distance in milliparsecs.

    Returns:
        float: The distance in parsecs.
    """
    return mpc / 1000


def pc_to_mpc(pc):
    """
    Converts a distance in parsecs (pc) to milliparsecs (mpc) -- the
    inverse of `mpc_to_pc`. Exact, same reasoning.

    Args:
        pc (float): The distance in parsecs.

    Returns:
        float: The distance in milliparsecs.
    """
    return pc * 1000


def pc_to_ly(pc):
    """
    Converts a distance in parsecs (pc) to light-years (ly), for
    human-readable display (e.g. "~48,923 ly from the galactic core") --
    see `docs/design/galaxy-coordinate-system.md` section 2. Not used by
    the database persistence layer itself (which stores galaxy-scale
    distances in parsecs directly); this exists purely for prose alongside
    the same "display string next to the raw stored value" treatment
    `table_*` columns get elsewhere in this schema.

    Args:
        pc (float): The distance in parsecs.

    Returns:
        float: The distance in light-years.
    """
    return pc * (physical_constants.AU_PER_PARSEC / physical_constants.LY_TO_AU)


def ly_to_pc(ly):
    """
    Converts a distance in light-years (ly) to parsecs (pc) -- the inverse
    of `pc_to_ly`. See that function's docstring.

    Args:
        ly (float): The distance in light-years.

    Returns:
        float: The distance in parsecs.
    """
    return ly * (physical_constants.LY_TO_AU / physical_constants.AU_PER_PARSEC)


def ly_to_au(ly):
    """
    Converts a distance in light-years (ly) to astronomical units (AU),
    for human-readable display at a stellar-phenomenon's own AU scale
    (`planetgen/web/maps/phenomenonmap.py`'s diagram) -- the AU counterpart to
    `pc_to_ly`/`ly_to_pc` above, using the same
    `physical_constants.LY_TO_AU` this package's other ly<->AU
    conversions already share.

    Args:
        ly (float): The distance in light-years.

    Returns:
        float: The distance in astronomical units.
    """
    return ly * physical_constants.LY_TO_AU


def au_to_ly(au):
    """
    Converts a distance in astronomical units (AU) to light-years (ly) --
    the inverse of `ly_to_au`. See that function's docstring.

    Args:
        au (float): The distance in astronomical units.

    Returns:
        float: The distance in light-years.
    """
    return au * physical_constants.AU_TO_LY


def to_scientific_notation(system_config: SystemConfig, number, precision=2):
    """
    Converts a number to scientific notation with the specified precision.

    This function is used for formatting large numbers in a compact and
    standardized way, suitable for data templates and displays. It checks the
    global `config.MARKDOWN` flag to determine the output format.

    Args:
        system_config (SystemConfig): The shared SystemConfig object.
        number (float): The number to convert.
        precision (int, optional): The number of decimal places to show. Defaults to 2.

    Returns:
        str: The number in scientific notation, either as a wikitext template
             or as HTML for markdown.
    """
    if number == 0:
        return "0"
    if not math.isfinite(number):
        # NaN/inf have no exponent; shown as-is rather than raising.
        return str(number)
    # Python's own `e` formatting rounds the mantissa *before* picking the
    # exponent, so 9.999 comes out as 1.00e+01 (not "10.00 x 10^0"), and
    # it's exact for subnormals (5e-324), where dividing by 10**exponent
    # underflowed to a division by zero.
    mantissa, exponent_text = f"{number:.{precision}e}".split("e")
    exponent = int(exponent_text)
    if system_config.MARKDOWN:
        return f"{mantissa} × 10<sup>{exponent}</sup>"
    else:
        output = f"Exp|{mantissa}|{exponent:d}"
        return "{{" + output + "}}"

def format_age_string(age_gy, precision=2):
    """
    Formats an age in billions of years (GY) into a human-readable string,
    dynamically choosing between "Million Years" and "Billion Years".

    Args:
        age_gy (float): The age in billions of years.
        precision (int, optional): The number of decimal places to show. Defaults to 2.

    Returns:
        str: A formatted string (e.g., "12.50 Billion Years" or "500 Million Years").
    """
    if age_gy >= 1.0:
        return f"{age_gy:.{precision}f} Billion Years"
    else:
        return f"{age_gy * 1000:.{precision}f} Million Years"

def calculate_object_mass(object_class, object_radius, planet_classes, planet_density, object_density=None):
    """
    Calculates the mass of a celestial object in kilograms.

    This function computes the mass based on the object's radius and density.
    If the density is not provided, it is randomly determined based on the
    object's class and type.

    Args:
        object_class (str): The class of the object (e.g., 'M', 'N').
        object_radius (float): The radius of the object in kilometers.
        planet_classes (dict): A dictionary defining the properties of planet classes.
        planet_density (dict): A dictionary of density ranges for planet types.
        object_density (float, optional): The density of the object in g/cm³.

    Returns:
        tuple: A tuple containing the volume in km³ and the mass in kg.
    """
    if object_density is None:
        min_density, max_density = planet_density[planet_classes[object_class]['type']]
        p_density = random.uniform(min_density, max_density)
    else:
        p_density = object_density

    volume_km3 = (4 / 3) * math.pi * object_radius ** 3
    volume_m3 = volume_km3 * physical_constants.KM_TO_M_FACTOR ** 3
    mass = volume_m3 * p_density * 1000
    return volume_km3, mass


def power_law_share(low, high, slope):
    """The unnormalized weight of [low, high] under dN/dlogM ~ M^-slope
    (`slope` > 0): `low^-slope - high^-slope`."""
    return low ** -slope - high ** -slope


def sample_power_law(low, high, slope):
    """A value in [low, high] drawn from dN/dlogM ~ M^-slope (`slope` > 0)
    by inverting its cumulative distribution: one `random.random()`."""
    low_term, high_term = low ** -slope, high ** -slope
    return (low_term + random.random() * (high_term - low_term)) ** (-1 / slope)


@finite_domain(clamped=("mode_fraction",))
def sample_bounded_bell(min_val, max_val, mode_fraction, spread_divisor=3.0, max_attempts=1000):
    """
    Draws a random value in [min_val, max_val] from a bounded bell-curve
    (Gaussian) distribution peaking at `min_val + mode_fraction * (max_val -
    min_val)`, instead of a flat uniform draw across the whole range.

    `mode_fraction` (0.0-1.0) is "what fraction through the available range
    is the statistically most common (modal) value" -- e.g. a class whose
    real-world single-body analog sits 27% of the way from its declared
    radius minimum to its maximum uses `mode_fraction=0.27` so generated
    instances cluster around that real value instead of being spread flatly
    across the whole declared range (see `tuning.PLANET_CLASSES`'
    own `size_mode` values and the real-world analogs their docstrings
    cite).

    The standard deviation self-adjusts to whichever bound is nearer the
    mode (reaching it in `spread_divisor` steps, ~3-sigma by default) so a
    mode pinned near one edge of the range still produces a legible bell
    shape instead of wasting most of its probability mass outside the
    range entirely. Uses rejection sampling (redraw until the result falls
    in [min_val, max_val]) rather than clamping a raw Gaussian draw to the
    bounds, which would pile spillover probability mass up at the boundary
    and destroy the bell shape there; `max_attempts` is a termination
    safety net (falls back to the clamped mean if never satisfied, which
    the self-adjusting spread should make effectively unreachable in
    practice for any reasonable `mode_fraction`).

    Args:
        min_val (float): Lower bound of the available range for this draw
                         (not necessarily a class's full declared range --
                         e.g. a moon's available radius may be further
                         capped by its parent's Hill sphere).
        max_val (float): Upper bound of the available range for this draw.
        mode_fraction (float): 0.0-1.0, how far through [min_val, max_val]
                               the distribution's peak sits.
        spread_divisor (float, optional): How many standard deviations
                                          reach the nearer bound from the
                                          mode. Defaults to 3.0.
        max_attempts (int, optional): Rejection-sampling attempt cap.

    Returns:
        float: The sampled value, guaranteed within [min_val, max_val].
    """
    if max_val <= min_val:
        return min_val
    mode_fraction = min(1.0, max(0.0, mode_fraction))
    mean = min_val + mode_fraction * (max_val - min_val)
    span = max_val - min_val
    # Floored so a mode pinned exactly at an edge (fraction 0.0 or 1.0)
    # doesn't collapse the spread to zero.
    nearer_bound_distance = max(min(mean - min_val, max_val - mean), span * 0.02)
    stdev = nearer_bound_distance / spread_divisor
    for _ in range(max_attempts):
        value = random.gauss(mean, stdev)
        if min_val <= value <= max_val:
            return value
    return min(max(mean, min_val), max_val)


@finite_domain()
def calculate_habitable_zone(luminosity):
    """
    Calculates the inner and outer boundaries of the habitable zone for a star.

    The habitable zone is defined as the region around a star where liquid
    water could exist on a planet's surface. This calculation is based on the
    star's luminosity.

    Args:
        luminosity (float): The luminosity of the star in Watts.

    Returns:
        tuple: A tuple containing the inner and outer radii of the habitable
               zone in AU.
    """
    solar_lum = luminosity / physical_constants.SOLAR_LUMINOSITY
    inner_radius = math.sqrt(solar_lum / 1.1)
    outer_radius = math.sqrt(solar_lum / 0.53)
    return (inner_radius, outer_radius)


@finite_domain()
def calculate_hill_sphere(distance_m, body_mass_kg, central_mass_kg):
    """
    Calculates the Hill sphere radius for a celestial body.

    The Hill sphere is the region around a celestial body where its own
    gravity is the dominant force for attracting satellites. This function
    calculates the radius of this sphere.

    Args:
        distance_m (float): The distance (semi-major axis) between the smaller
                            body and the larger central body, in meters.
        body_mass_kg (float): The mass of the smaller body (e.g., a planet) in kilograms.
        central_mass_kg (float): The mass of the larger central body (e.g., a star) in kilograms.

    Returns:
        float: The radius of the Hill sphere in meters.
    """
    return distance_m * (body_mass_kg / (3 * central_mass_kg)) ** (1 / 3)


@finite_domain(clamped=("companion_mass_fraction", "eccentricity"))
def holman_wiegert_critical_semimajor_axis(binary_separation_au, companion_mass_fraction, eccentricity):
    """
    Holman, M. & Wiegert, P. (1999), AJ 117:621, "Long-Term Stability of
    Planets in Binary Systems" -- the empirical fit for an S-type (wide)
    binary's critical semi-major axis: the largest orbit around ONE star of
    the pair that remains long-term stable against the other star's
    periodic gravitational perturbation.

        a_crit / a_bin = 0.464 - 0.380*mu - 0.631*e + 0.586*mu*e
                         + 0.150*e^2 - 0.198*mu*e^2

    `mu` is the *perturbing companion's* mass fraction of the pair's total
    mass, `mu = M_companion / (M_this_star + M_companion)` -- this is
    evaluated once per star, using that star's own companion, so calling
    this twice for one pair (once from each star's perspective) generally
    yields two different `a_crit` values unless the two masses are equal.

    Valid over roughly `mu` in [0.1, 0.9] and `e` in [0.0, 0.8] (Holman &
    Wiegert's own numerical grid doesn't extend meaningfully further) --
    inputs are clamped to `physical_constants.HOLMAN_WIEGERT_MU_RANGE`/
    `HOLMAN_WIEGERT_ECCENTRICITY_RANGE` rather than extrapolated, since a
    saturated-at-the-boundary estimate is more useful than either an
    exception or a silently invalid extrapolation.

    Args:
        binary_separation_au (float): The binary pair's own semi-major
                                      axis (separation), in AU.
        companion_mass_fraction (float): `mu`, as defined above (0-1).
        eccentricity (float): The binary orbit's own eccentricity (0-1).

    Returns:
        float: The critical semi-major axis, in AU (same unit as
              `binary_separation_au`) -- this star's own maximum stable
              planetary orbit given the companion's influence.

    Example (equal-mass, circular pair -- a standard reference case):
        `holman_wiegert_critical_semimajor_axis(1.0, 0.5, 0.0)` gives
        `0.464 - 0.380*0.5 = 0.274` exactly, consistent with the commonly
        cited ~0.27-0.30 * a_bin figure for this case.
    """
    mu_min, mu_max = physical_constants.HOLMAN_WIEGERT_MU_RANGE
    e_min, e_max = physical_constants.HOLMAN_WIEGERT_ECCENTRICITY_RANGE
    mu = min(max(companion_mass_fraction, mu_min), mu_max)
    e = min(max(eccentricity, e_min), e_max)

    ratio = (0.464 - 0.380 * mu - 0.631 * e
             + 0.586 * mu * e + 0.150 * e ** 2 - 0.198 * mu * e ** 2)
    return ratio * binary_separation_au



@finite_domain(clamped=("secondary_mass_fraction", "eccentricity"))
def holman_wiegert_circumbinary_a_crit_au(binary_separation_au, secondary_mass_fraction, eccentricity=0.0):
    """
    Holman, M. & Wiegert, P. (1999), AJ 117:621 -- the empirical fit for a
    P-type (circumbinary) orbit's critical semi-major axis: the SMALLEST
    orbit around both stars of a close pair that stays long-term stable.

        a_crit / a_bin = 1.60 + 5.10*e - 2.22*e^2 + 4.12*mu - 4.27*e*mu
                         - 5.09*mu^2 + 4.61*e^2*mu^2

    `mu` is the lighter star's fraction of the pair's total mass. Inputs
    are clamped to `tuning.HOLMAN_WIEGERT_P_TYPE_MU_RANGE`/
    `HOLMAN_WIEGERT_P_TYPE_ECCENTRICITY_RANGE` (the fit's tested grid)
    rather than extrapolated. A circular equal-mass pair gives about
    2.39 * a_bin.

    Args:
        binary_separation_au (float): The pair's own separation, in AU.
        secondary_mass_fraction (float): `mu`, as defined above.
        eccentricity (float): The pair's orbital eccentricity.

    Returns:
        float: The innermost stable circumbinary orbit, in AU.
    """
    mu_min, mu_max = tuning.HOLMAN_WIEGERT_P_TYPE_MU_RANGE
    e_min, e_max = tuning.HOLMAN_WIEGERT_P_TYPE_ECCENTRICITY_RANGE
    mu = min(max(secondary_mass_fraction, mu_min), mu_max)
    e = min(max(eccentricity, e_min), e_max)

    ratio = (1.60 + 5.10 * e - 2.22 * e ** 2 + 4.12 * mu - 4.27 * e * mu
             - 5.09 * mu ** 2 + 4.61 * e ** 2 * mu ** 2)
    return ratio * binary_separation_au

@finite_domain()
def mutual_hill_radius_au(mass1_kg, mass2_kg, distance1_au, distance2_au, central_mass_kg):
    """
    Gladman (1993), Icarus 106:247, "Dynamical stability of the outer solar
    system and the delivery of comets" -- the mutual Hill radius of two
    orbiting bodies:

        R_H,mutual = ((m1 + m2) / (3 * M_central))^(1/3) * ((a1 + a2) / 2)

    Gladman's own derivation assumes both bodies orbit the SAME central
    mass -- used here (see `systemData.StarSystem._validate_cross_star_clearance`)
    across two planets that orbit *different* stars of a wide binary, this
    is a physically-motivated extension of the criterion's spirit rather
    than a literal textbook application; see that method's docstring for
    the specific choice of `central_mass_kg` and its justification.

    Args:
        mass1_kg (float): First body's own mass, in kg.
        mass2_kg (float): Second body's own mass, in kg.
        distance1_au (float): First body's distance from whatever it
                              orbits, in AU.
        distance2_au (float): Second body's distance from whatever it
                              orbits, in AU.
        central_mass_kg (float): The mass, in kg, both distances above are
                                 measured against (see docstring above for
                                 how this generator chooses it when the two
                                 bodies orbit different stars).

    Returns:
        float: The mutual Hill radius, in AU (same unit as
              `distance1_au`/`distance2_au`).
    """
    return ((mass1_kg + mass2_kg) / (3 * central_mass_kg)) ** (1 / 3) * ((distance1_au + distance2_au) / 2)


def sample_wide_binary_separation_au():
    """
    Draws an S-type (wide) binary's separation (semi-major axis), log-
    uniformly between `tuning.WIDE_BINARY_SEPARATION_MIN_AU`
    and `WIDE_BINARY_SEPARATION_MAX_AU` -- see those constants' own
    docstring for why log-uniform (not linear-uniform) sampling is used.
    Uses the same `math.exp(random.uniform(math.log(...), math.log(...)))`
    idiom already used for moon-distance placement
    (`planetPhysics.generate_moons`).

    Returns:
        float: A separation, in AU.
    """
    low = tuning.WIDE_BINARY_SEPARATION_MIN_AU
    high = tuning.WIDE_BINARY_SEPARATION_MAX_AU
    return math.exp(random.uniform(math.log(low), math.log(high)))


def sample_wide_binary_eccentricity():
    """
    Draws an S-type (wide) binary's orbital eccentricity from a "thermal"
    distribution, `f(e) = 2e`, capped at
    `tuning.WIDE_BINARY_ECCENTRICITY_MAX` -- see that constant's
    own docstring for why wide pairs (unlike the close/P-type pair) keep a
    realistic, generally non-zero eccentricity.

    Derivation: the thermal PDF `f(e) = 2e` restricted to `[0, e_max]` and
    renormalized is still exactly proportional to `e` (just rescaled), so
    its CDF is `F(e) = (e / e_max)^2` and the closed-form inverse-CDF
    sample is `e = e_max * sqrt(u)`, `u ~ Uniform(0, 1)` -- no rejection
    sampling needed.

    Returns:
        float: An eccentricity, in [0, `WIDE_BINARY_ECCENTRICITY_MAX`).
    """
    return tuning.WIDE_BINARY_ECCENTRICITY_MAX * math.sqrt(random.random())


@finite_domain()
def mutual_hill_radius_m(distance1_m, distance2_m, mass1_kg, mass2_kg, central_mass_kg):
    """
    Calculates the *mutual* Hill radius of two bodies that orbit the same
    primary -- the length scale real orbital-dynamics stability criteria
    (Gladman 1993; Chambers, Wetherill & Boslough 1996; Smith & Lissauer
    1999/2009) use for how close two adjacent planets' orbits can safely
    be, as distinct from `calculate_hill_sphere`'s single-body sphere of
    gravitational dominance (the right tool for "how far can a satellite
    orbit *this* body," not for "how close can two planets orbit each
    other").

    R_H,mutual = ((a1 + a2) / 2) * ((m1 + m2) / (3 * M_central)) ** (1/3)

    -- i.e. the single-body formula generalized to use the *pair's*
    combined mass and average distance, rather than either body's own
    mass and distance alone. See
    `tuning.MUTUAL_HILL_RADII_SEPARATION` for how this
    generator turns this length into an actual minimum separation.

    Args:
        distance1_m (float): The first body's distance from the shared
                             primary, in meters.
        distance2_m (float): The second body's distance from the shared
                             primary, in meters.
        mass1_kg (float): The first body's mass, in kilograms.
        mass2_kg (float): The second body's mass, in kilograms.
        central_mass_kg (float): The shared primary's mass, in kilograms.

    Returns:
        float: The mutual Hill radius, in meters.
    """
    avg_distance_m = (distance1_m + distance2_m) / 2
    return avg_distance_m * ((mass1_kg + mass2_kg) / (3 * central_mass_kg)) ** (1 / 3)


@finite_domain()
def snow_line_au(luminosity_w):
    """
    The snow line (ice condensation point, ~170K) for a star of a given
    luminosity, in AU -- real protoplanetary disk temperature falls off
    with stellar flux, i.e. with distance^-2, the same physical reasoning
    `calculate_habitable_zone` already uses for its own sqrt(luminosity)
    boundaries (see `physical_constants.SNOW_LINE_AU_AT_1_LSUN`).

    Args:
        luminosity_w (float): The star's luminosity, in Watts.

    Returns:
        float: The snow line's distance from the star, in AU.
    """
    solar_lum = luminosity_w / physical_constants.SOLAR_LUMINOSITY
    return physical_constants.SNOW_LINE_AU_AT_1_LSUN * math.sqrt(solar_lum)


@finite_domain()
def disk_surface_density_scale(star_mass_kg):
    """
    How much a star's own protoplanetary disk's solid surface density
    should be scaled relative to the Sun's (MMSN) value, based on the
    real, observed disk-mass-vs-stellar-mass relation -- see
    `tuning.DISK_MASS_STELLAR_MASS_EXPONENT`'s docstring for
    the literature basis. Feeds `mmsn_surface_density_gcm2`'s own
    `density_scale` argument.

    Args:
        star_mass_kg (float): The star's mass, in kilograms.

    Returns:
        float: A dimensionless scale factor, 1.0 for a solar-mass star.
    """
    solar_masses = star_mass_kg / physical_constants.SOLAR_MASS_TO_KG
    return solar_masses ** tuning.DISK_MASS_STELLAR_MASS_EXPONENT


@finite_domain()
def mmsn_surface_density_gcm2(distance_au, snow_line_au, density_scale=1.0):
    """
    The Minimum Mass Solar Nebula's solid surface density at a given
    distance from the star (Hayashi 1981) -- see
    `physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2`'s docstring
    for the model and its snow-line ice-boost jump.

    Args:
        distance_au (float): Distance from the star, in AU.
        snow_line_au (float): This star's own snow line (see
                              `snow_line_au`), in AU.
        density_scale (float, optional): A star-dependent scale factor
                                         (see `disk_surface_density_scale`)
                                         applied on top of the Sun's own
                                         MMSN normalization. Defaults to
                                         1.0 (a solar-mass star's disk).

    Returns:
        float: The solid surface density at this distance, in g/cm^2.
    """
    density = (
        physical_constants.MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2
        * density_scale
        * distance_au ** physical_constants.MMSN_SURFACE_DENSITY_EXPONENT
    )
    if distance_au >= snow_line_au:
        density *= physical_constants.SNOW_LINE_ICE_BOOST_FACTOR
    return density


@finite_domain()
def isolation_mass_kg(distance_au, surface_density_gcm2, star_mass_kg):
    """
    The oligarchic-growth isolation mass (Lissauer 1993; Kokubo & Ida
    2000, 2002): the mass a growing embryo reaches once it has cleared its
    own feeding zone of width `b` mutual Hill radii
    (`tuning.MUTUAL_HILL_RADII_SEPARATION`) -- the same `b`
    `StarSystem._mutual_min_distance_au` uses for adjacent-planet spacing,
    so this count estimate and that spacing rule are provably consistent.

    Derivation: M_iso = 2*pi*a*(b*R_H)*Sigma, where R_H = a*(M_iso /
    (3*M_star))^(1/3) is the embryo's *own* Hill radius (it hasn't met a
    neighbor yet, so this is the single-body form, not the mutual one).
    Substituting and solving for M_iso (it appears on both sides) gives
    the closed form below:

        M_iso = (2*pi * b * Sigma * a^2)^(3/2) / (3 * M_star)^(1/2)

    Args:
        distance_au (float): Distance from the star, in AU.
        surface_density_gcm2 (float): Local solid surface density at this
                                      distance (see
                                      `mmsn_surface_density_gcm2`), in
                                      g/cm^2.
        star_mass_kg (float): The star's mass, in kilograms.

    Returns:
        float: The isolation mass, in kilograms.
    """
    distance_m = distance_au * physical_constants.AU_TO_M
    surface_density_kgm2 = surface_density_gcm2 * 10  # 1 g/cm^2 = 10 kg/m^2
    b = tuning.MUTUAL_HILL_RADII_SEPARATION
    base = 2 * math.pi * b * surface_density_kgm2 * distance_m ** 2
    return base ** (3 / 2) / (3 * star_mass_kg) ** 0.5


# -inf is just "not positive" (the documented (0.0, 0.0)); +inf fails the result check.
@finite_domain(allow_inf=("distance_ly",))
def calculate_galactic_orbit(distance_ly):
    """
    Estimates a star system's circular orbital speed and orbital period
    around the galactic center, given its distance from it.

    Uses `physical_constants.GALACTIC_ROTATION_FLAT_VELOCITY_KMS`/
    `GALACTIC_ROTATION_CORE_RADIUS_PC`'s simple rotation-curve model (see
    that module's comment for the physical justification and calibration
    against Sol's own distance) rather than a Keplerian point-mass orbit
    around `physical_constants.MILKY_WAY_MASS` -- the latter would put a
    Sol-distance orbit at nearly 800 km/s, ~4x the real value, since most
    of the galaxy's mass isn't actually enclosed within that radius the way
    a naive point-mass calculation assumes. This function is deliberately
    independent of the orbiting body's own mass (unlike
    `calculate_hill_sphere`), matching real orbital mechanics at galactic
    scale: essentially every star's mass is negligible next to the
    galaxy's, so orbital speed at a given radius is the same for any star
    there, not a function of that star's own mass.

    Args:
        distance_ly (float): The star system's distance from the galactic
                             center, in light-years.

    Returns:
        tuple: `(orbital_speed_kms, orbital_period_gy)` -- circular orbital
              speed in km/s, and orbital period in billions of years (Gy),
              the same unit `Star.age`/`lifespan` already use. Both `0.0`
              for a system placed exactly at the galactic center (r = 0,
              where a circular orbit is degenerate).
    """
    if distance_ly <= 0:
        return 0.0, 0.0

    distance_pc = ly_to_pc(distance_ly)
    core_radius_pc = physical_constants.GALACTIC_ROTATION_CORE_RADIUS_PC
    orbital_speed_kms = (
        physical_constants.GALACTIC_ROTATION_FLAT_VELOCITY_KMS
        * distance_pc / math.sqrt(distance_pc ** 2 + core_radius_pc ** 2)
    )

    circumference_m = 2 * math.pi * distance_ly * physical_constants.LY_TO_M
    orbital_period_s = circumference_m / (orbital_speed_kms * physical_constants.KM_TO_M_FACTOR)
    orbital_period_gy = (orbital_period_s / physical_constants.SECONDS_PER_YEAR) / 1e9

    return orbital_speed_kms, orbital_period_gy


def generate_galactic_orbit_fields(galactic_center_dist_ly=None, galactic_orbital_phase_deg=None):
    """
    Generates the four `galactic_orbital_*` fields every gravitationally-
    bound-to-the-galaxy body this generator tracks needs: a `Star`
    (`Star.__init__`), a compact remnant (`compactRemnant.CompactRemnant.
    _finish_init`), and every standalone exotic phenomenon (`nebulaData.
    Nebula`, `supernovaRemnantData.SupernovaRemnant`, `roguePlanetData.
    RoguePlanet`/`InterstellarComet`, `asteroidFieldData.AsteroidField`) --
    a rogue planet or a comet passing through is unbound from any specific
    star, not from the galaxy itself, so it still orbits the galactic
    center on the same timescale a lone star does, via the exact same
    mass-independent `calculate_galactic_orbit` formula. Extracted here so
    every caller computes this identically instead of re-deriving it.

    Args:
        galactic_center_dist_ly (float, optional): This body's actual
            distance from the galactic center, in light-years -- see
            `Star.calculate_system_perimeter`'s docstring for the same
            parameter/fallback convention. `None` (the default) uses the
            fixed `physical_constants.GALACTIC_CENTER_DISTANCE_LY` constant.
        galactic_orbital_phase_deg (float, optional): This body's current
            angular position around its galactic orbit, in degrees -- see
            `Star.__init__`'s docstring for the same parameter. `None`
            (the default) rolls a fresh random value in `[0, 360)`.

    Returns:
        tuple: `(galactic_orbital_speed_kms, galactic_orbital_period_gy,
              galactic_orbital_phase_deg, galactic_min_update_interval_years)`.
    """
    if galactic_center_dist_ly is None:
        galactic_center_dist_ly = physical_constants.GALACTIC_CENTER_DISTANCE_LY
    speed_kms, period_gy = calculate_galactic_orbit(galactic_center_dist_ly)
    phase_deg = galactic_orbital_phase_deg if galactic_orbital_phase_deg is not None else random.uniform(0, 360)
    min_update_interval_years = minimum_update_interval_years(period_gy * 1e9)
    return speed_kms, period_gy, phase_deg, min_update_interval_years


def format_galactic_orbit(speed_kms, period_gy):
    """
    Formats a `(galactic_orbital_speed_kms, galactic_orbital_period_gy)`
    pair as the display string every "Galactic Orbit" table row uses (e.g.
    `Star.get_table_properties`, `doubleStar.BinaryStarProxy.
    get_table_properties`, `compactRemnant.BlackHole`/`NeutronStar`, and
    every standalone exotic phenomenon) -- extracted so all of them render
    it identically instead of re-deriving the same f-string.

    Args:
        speed_kms (float): Circular orbital speed, km/s.
        period_gy (float): Orbital period, billions of years.

    Returns:
        str: e.g. `"206 km/s (236 My per orbit)"`, both on the shared
             ladders (`format_speed_kms`, `format_period_years`).
    """
    return f"{format_speed_kms(speed_kms)} ({format_period_years(period_gy * 1e9)} per orbit)"


@finite_domain()
def circular_orbital_speed_kms(distance_au, period_years):
    """
    Tangential speed of a circular orbit, given its radius and period:
    `v = 2*pi*r / T`. Exact (constant at every point of the orbit), since
    this generator only ever models circular orbits for planets and moons
    (see `planetPhysics.generate_orbital_motion_properties`) -- unlike
    `calculate_galactic_orbit`, which has to *assume* a rotation-curve
    model to get a speed at all, a planet's/moon's period is already known
    exactly from Kepler's third law (`planetPhysics.
    calculate_orbital_period_years`), so speed here is a direct
    geometric consequence of the two, not a separate physical model.

    Args:
        distance_au (float): Orbital radius (semi-major axis), in AU.
        period_years (float): Orbital period, in years.

    Returns:
        float: Orbital speed, in km/s.
    """
    circumference_km = 2 * math.pi * distance_au * physical_constants.AU_TO_KM
    period_seconds = period_years * physical_constants.SECONDS_PER_YEAR
    return circumference_km / period_seconds


@finite_domain()
def orbital_position_au(distance_au, inclination_deg, ascending_node_deg, phase_deg):
    """
    Converts a circular orbit's elements -- radius, inclination, ascending
    node, and current phase -- into a 3D Cartesian position relative to
    the orbit's primary (the body actually being orbited: a star/binary
    system center for a planet, a planet for a moon).

    Standard orbital-plane-to-reference-frame rotation, specialized for a
    circular orbit: `orbital_phase_deg` already plays the role of the
    argument of latitude `u = omega + true_anomaly` directly (no separate
    argument-of-periapsis term, since a circular orbit has no periapsis to
    measure one from -- see `planetPhysics.generate_orbital_motion_properties`'s
    docstring). `inclination_deg`/`ascending_node_deg` orient the orbital
    plane itself; `phase_deg` is where the body sits within it:

        u = radians(phase_deg), i = radians(inclination_deg), Om = radians(ascending_node_deg)
        x = r * (cos(Om)*cos(u) - sin(Om)*sin(u)*cos(i))
        y = r * (sin(Om)*cos(u) + cos(Om)*sin(u)*cos(i))
        z = r * sin(u)*sin(i)

    At `i = 0` (an uninclined orbit) this correctly collapses to a flat
    circle in the primary's own reference plane (`z = 0` always); `node`
    and `phase` become degenerate there (no inclined plane left for the
    ascending node to describe the crossing of), so only their sum
    matters: `(r*cos(node+u), r*sin(node+u), 0)` -- the standard,
    physically expected behavior at this degenerate case, same as in real
    orbital mechanics, not a bug.

    Args:
        distance_au (float): Orbital radius, in AU.
        inclination_deg (float): Orbital plane tilt, in degrees, relative
                                 to the primary's reference plane.
        ascending_node_deg (float): Longitude of the ascending node, in
                                    degrees.
        phase_deg (float): Current argument of latitude (position angle
                           around the orbit), in degrees.

    Returns:
        tuple: `(x_au, y_au, z_au)`, relative to the primary, in the same
              reference frame `inclination_deg`/`ascending_node_deg` are
              measured against.
    """
    u = math.radians(phase_deg)
    i = math.radians(inclination_deg)
    node = math.radians(ascending_node_deg)

    cos_u, sin_u = math.cos(u), math.sin(u)
    cos_i = math.cos(i)
    cos_node, sin_node = math.cos(node), math.sin(node)

    x = distance_au * (cos_node * cos_u - sin_node * sin_u * cos_i)
    y = distance_au * (sin_node * cos_u + cos_node * sin_u * cos_i)
    z = distance_au * sin_u * math.sin(i)

    return x, y, z


def calculate_reflex_offset(parent_mass_kg, children):
    """
    A parent body's own displacement from its nominal fixed point, caused
    by the combined gravitational pull of everything orbiting it -- the
    "wobble"/reflex-motion half of a proper two-body (barycentric)
    treatment, mirrored on the SQL side by the correlated-subquery
    `UPDATE`s in `_db.advance_orbital_phases`.

    For a single child, this is the exact two-body barycentric formula:
    `offset = -(child_mass / (parent_mass + child_mass)) * relative_vector`,
    where `relative_vector` is the child's own already-stored position
    relative to the parent (e.g. `Planet.position_x/y/z`,
    `BinaryStarProxy.binary_mutual_position_x/y/z`) -- deliberately never
    recomputed or changed by this function, since a large amount of
    existing physics (insolation, Hill sphere, tidal locking) depends on
    that vector remaining the *true* separation, not a barycenter-reduced
    one; this only computes the *parent's* own small displacement.

    For multiple children (a star with several planets, a planet with
    several moons), each child's individual pairwise pull is summed --
    the standard linear-superposition approximation real radial-velocity
    work uses for multi-planet reflex motion, exact to leading order
    whenever the parent is much more massive than any single child (true
    for every relationship in this generator).

    Args:
        parent_mass_kg (float): The parent body's own mass, in kg.
        children (list): `(mass_kg, x_au, y_au, z_au)` tuples, one per
                         body orbiting the parent, each already in the
                         parent's own reference frame.

    Returns:
        tuple: `(x_au, y_au, z_au)`, the parent's own offset from its
              nominal fixed point -- `(0.0, 0.0, 0.0)` if `children` is
              empty.
    """
    offset_x = offset_y = offset_z = 0.0
    for child_mass_kg, x_au, y_au, z_au in children:
        mu_child = child_mass_kg / (parent_mass_kg + child_mass_kg)
        offset_x -= mu_child * x_au
        offset_y -= mu_child * y_au
        offset_z -= mu_child * z_au
    return offset_x, offset_y, offset_z


@finite_domain()
def minimum_update_interval_years(period_years):
    """
    The shortest `elapsed_years` worth advancing a body's
    `orbital_phase_deg` for at all -- below this, the phase delta added is
    smaller than `orbital_phase_deg`'s own floating-point resolution, so
    `MOD(orbital_phase_deg + delta, 360)` is guaranteed to round right
    back to the exact value already stored: a wasted write that changes
    nothing. Exists specifically to guard `_db.advance_orbital_phases`
    against that silent no-op, not as a narrative/display stat -- unlike
    this package's other derived quantities, there is no scale at which a
    human would want to read this number (it lands in the nanosecond
    range for any realistic orbital period, see below).

    Derivation: `orbital_phase_deg` ranges over `[0, 360)`, stored as an
    IEEE 754 double (MySQL `DOUBLE`, Python `float` -- identical
    representation). The coarsest (least precise) representable step
    anywhere in that range is the unit-in-the-last-place at magnitudes
    just under 360 -- `math.ulp(360.0)` -- used here as a single,
    domain-wide conservative bound rather than a per-row value that would
    depend on the body's current phase (tighter near 0, coarser near 360)
    and so would itself need updating every time phase does, for no real
    benefit. A phase delta at or above this many degrees is guaranteed to
    change the stored value, regardless of where in `[0, 360)` the
    current phase happens to sit; anything smaller might not.

    `elapsed_years / period_years * 360 >= ulp_deg`
    `elapsed_years >= period_years * ulp_deg / 360`

    Worked example: a 1-year period gives a floor around 5e-9 seconds --
    roughly 16 orders of magnitude below `planetgen.cli.orbits`'s own "once a
    month or so" real-world cadence (see that module's docstring), so
    this guard exists for correctness/defensiveness (a future caller
    advancing time in much smaller steps, e.g. a fast-forward simulation)
    rather than because today's actual usage pattern ever comes close to
    triggering it.

    Args:
        period_years (float): The body's own orbital period, in years.
                              Always positive and finite for a real
                              generated planet/moon (Kepler's third law on
                              a positive distance and mass).

    Returns:
        float: The minimum `elapsed_years` worth calling
              `_db.advance_orbital_phases` for, for a body with this
              period.
    """
    ulp_deg = math.ulp(360.0)
    return period_years * ulp_deg / 360


def split_into_syllables(name):
    """
    Splits a word into a list of syllables.

    This is a basic syllable splitting function that helps in the process of
    generating new names by rearranging syllables from existing names.

    Args:
        name (str): The word to be split into syllables.

    Returns:
        list: A list of strings, where each string is a syllable.
    """
    syllables = []
    current_syllable = ""
    for i, char in enumerate(name):
        current_syllable += char
        if char in VOWELS and i < len(name) - 1 and name[i+1] not in VOWELS:
            syllables.append(current_syllable)
            current_syllable = ""
    if current_syllable:
        syllables.append(current_syllable)
    return syllables


_NSFW_PLAIN_WORDS = tuple(w for w in NSFW_WORDS if w and "'" not in w and " " not in w)
_NSFW_SPACED_WORDS = tuple(w for w in NSFW_WORDS if w and ("'" in w or " " in w))


def is_name_valid(name):
    """
    Checks if a generated name is valid based on multiple criteria.

    This function ensures that the generated name is not a common English word,
    does not contain any offensive terms, and follows basic phonetic rules.
    The validation checks are as follows:
    - The name should not exist in the NLTK dictionary of words.
    - The name should not contain any substring from the NSFW (Not Safe For Work) word list,
      also checked with its apostrophes and spaces removed.
    - No word of the name may be a `nameUniqueness` decoration word (`_DECORATION_WORDS`).
    - The name should not contain more than two consecutive vowels.
    - The name should not contain more than two consecutive consonants.
    - The name should not contain any of the defined bad consonant clusters.

    These checks help in generating names that are unique, appropriate, and sound plausible.

    Args:
        name (str): The name to validate. The function expects a lowercase string.

    Returns:
        bool: Returns `True` if the name is valid according to all the rules,
              otherwise returns `False`.
    """
    name_lower = name.lower()
    if name_lower in DICTIONARY_WORDS:
        return False
    if any(token in _DECORATION_WORDS for token in name_lower.split()):
        # A word that is also a nameUniqueness decoration ("Liten",
        # "Ohana", ...) would be stripped off an undecorated name by
        # `nameUniqueness.strip_decoration`, grouping it with another base.
        return False
    # Checked with apostrophes/spaces removed too: an apostrophe spliced in
    # by a phoneme chunk ("pak'i", "k'ike") must not hide an offensive word.
    # A word with no apostrophe or space that appears in any variant also
    # appears in the fully squeezed one, so only the few words that carry
    # one of those need every variant (keeps this hot check one pass).
    squeezed = name_lower.replace("'", "")
    fully_squeezed = squeezed.replace(" ", "")
    if any(word in fully_squeezed for word in _NSFW_PLAIN_WORDS):
        return False
    if _NSFW_SPACED_WORDS:
        variants = (name_lower, squeezed, name_lower.replace(" ", ""))
        if any(word in variant for word in _NSFW_SPACED_WORDS for variant in variants):
            return False
    
    vowel_count = 0
    consonant_count = 0
    for char in name_lower:
        if char in VOWELS:
            vowel_count += 1
            consonant_count = 0
        elif char.isalpha():
            consonant_count += 1
            vowel_count = 0
        else:
            vowel_count = 0
            consonant_count = 0
        if vowel_count > 2 or consonant_count > 2:
            return False

    for cluster in BAD_CONSONANTS:
        if cluster in name_lower:
            return False

    return True


_DECORATION_WORDS = frozenset(
    word.lower() for word in (*GREEK_LETTERS, *DIMINUTIVE_PREFIXES, *COMPANION_SUFFIXES,
                              *ROMAN_NUMERALS_BY_VALUE.values())
)
"""Every word `nameUniqueness` decorates a name with (lowercase) -- see
`is_name_valid`/`generate_phoneme_salad_name`, which keep generated words
from ever being one of them."""


def split_long_word(name):
    """
    Splits a long word into two, capitalizing the second word.

    This function improves the readability of long generated names by splitting
    them into two parts, creating a compound name effect.

    A name can contain an embedded apostrophe (from a base name like
    "Hi'iaka" or a spliced-in `UNIVERSAL_PHONEMES` chunk like "ch'"/"b'a").
    If the naive midpoint split landed right next to one of those, one half
    would end up starting or ending with a bare "'" (e.g. "Amech' Snesis") --
    so the split point is nudged past any apostrophe it would otherwise cut
    beside. If nudging would consume the whole rest of the name, the split
    is skipped and the name is returned unsplit rather than produce an
    empty half.

    Args:
        name (str): The long word to be split.

    Returns:
        str: The split and capitalized name, or the original name if not
             long enough (or if avoiding an apostrophe boundary leaves no
             valid split point).
    """
    if len(name) > WORD_SIZE_MEAN:
        if " " in name:
            # Already more than one word (a base name like "El Nath");
            # splitting again could put a second space beside the first.
            return name
        split_point = len(name) // 2
        while split_point < len(name) and (name[split_point - 1] == "'" or name[split_point] == "'"):
            split_point += 1
        if not (0 < split_point < len(name)):
            return name
        return name[:split_point] + " " + name[split_point:].capitalize()
    return name


UNIVERSAL_PHONEME_CHANCE = 0.4
"""
Odds that `generate_phoneme_salad_name` splices an extra chunk from
`names.UNIVERSAL_PHONEMES` into a generated name's syllable pool. Applies
uniformly to stars, planets, moons, and sectors -- every one of those
funnels through this same function -- rather than needing each type's own
prefix/suffix lists to be extended individually.
"""

MAX_NAME_GENERATION_ATTEMPTS = 10_000
"""
Hard cap on `generate_phoneme_salad_name`'s own retry loop -- see that
function's `Raises` doc. A real generation run finds a valid name within a
handful of attempts; this is only a backstop against `is_name_valid`
rejecting every candidate forever (previously an unconditional `while
True` with no way out at all).
"""


def generate_phoneme_salad_name(name_list, prefix_list, suffix_list, allow_split=True, syllable_fraction=1.0, max_length=None):
    """
    Generates a unique, phonetically pleasing name from a list of base names.

    This function creates new names by taking a base name, shuffling its
    syllables, and adding a prefix and suffix. It includes logic to ensure
    the resulting name is phonetically plausible and passes validation checks.
    With `UNIVERSAL_PHONEME_CHANCE` odds, it also splices in one chunk from
    `names.UNIVERSAL_PHONEMES` -- a cross-linguistic phoneme pool shared by
    every name type -- widening the cultural range names are drawn from
    beyond whatever's in `name_list` itself.

    Args:
        name_list (list): A list of base names to choose from.
        prefix_list (list): A list of possible prefixes.
        suffix_list (list): A list of possible suffixes.
        allow_split (bool): Whether a long result may be split into two
                            space-separated words via `split_long_word`
                            (e.g. `"Xyleth Anore"`). Default `True` for
                            stars/planets/moons, where that reads fine as
                            one name. Callers that combine multiple calls
                            into one already-multi-word name (e.g.
                            `sectorGen.generate_sector_name`, which joins
                            two of these into a two-word sector name) must
                            pass `False` here, or a single call splitting
                            internally would silently make the combined
                            result 3-4 words instead of the intended 2.
        syllable_fraction (float): Fraction (0-1] of the shuffled base
                            name's syllables to actually keep before the
                            prefix/suffix are attached, trimming from the
                            end. Default `1.0` keeps every syllable
                            (unchanged behavior for stars/planets/moons).
                            `sectorGen.generate_sector_name` passes `0.5`
                            here so sector names -- built by joining two
                            of these calls into one two-word name -- come
                            out roughly half as long per word; base names
                            long enough to still shrink always keep at
                            least 1 syllable.
        max_length (int): Optional hard cap, in characters, on the fully
                            assembled name (base syllables + prefix +
                            suffix, before capitalization). `None` (the
                            default) leaves names uncapped. Chopping
                            happens after the suffix is attached, so it's
                            the backstop against `syllable_fraction`
                            alone not being enough -- prefixes, suffixes,
                            and the occasional spliced-in
                            `UNIVERSAL_PHONEMES` chunk are fixed-ish
                            overhead that doesn't shrink with
                            `syllable_fraction`, so a long base name can
                            still produce a longer-than-intended result
                            without this. `sectorGen.generate_sector_name`
                            passes `7` here alongside `syllable_fraction=0.5`
                            to reliably keep each half of a sector name
                            short.

    Returns:
        str: A newly generated, unique name.

    Raises:
        RuntimeError: If `MAX_NAME_GENERATION_ATTEMPTS` consecutive draws
            all fail `is_name_valid` -- confirmed (via a hard-timeout test,
            `test_bughunt_name_exhaustion.py`) to otherwise loop forever
            with no way out. In real generation this is never reached (a
            valid name is essentially always found within a handful of
            attempts); this guard only matters if `is_name_valid`'s own
            filters (or a `name_list`/`prefix_list`/`suffix_list` this
            narrow) were ever misconfigured to reject everything, in which
            case a whole generation run should fail loudly and immediately
            rather than hang with no explanation.
    """
    for _attempt in range(MAX_NAME_GENERATION_ATTEMPTS):
        name = random.choice(name_list)

        syllables = split_into_syllables(name)
        if len(syllables) > 1:
            random.shuffle(syllables)

        if syllable_fraction < 1.0 and len(syllables) > 1:
            keep = max(1, round(len(syllables) * syllable_fraction))
            syllables = syllables[:keep]

        if random.random() < UNIVERSAL_PHONEME_CHANCE:
            syllables.insert(random.randint(0, len(syllables)), random.choice(UNIVERSAL_PHONEMES))

        name = "".join(syllables)

        prefix = random.choice(prefix_list)
        if prefix[-1] in VOWELS and name[0].lower() in VOWELS:
            name = prefix + name[1:]
        elif prefix[-1] not in VOWELS and name[0].lower() not in VOWELS:
            if (prefix[-1] + name[0].lower()) in BAD_CONSONANTS:
                name = prefix + "'" + name
            else:
                name = prefix + name
        else:
            name = prefix + name

        suffix = random.choice(suffix_list)
        if name[-1] in VOWELS and suffix[0].lower() in VOWELS:
            name = name + suffix[1:]
        elif name[-1] not in VOWELS and suffix[0].lower() not in VOWELS:
            if (name[-1] + suffix[0].lower()) in BAD_CONSONANTS:
                name = name + "'" + suffix
            else:
                name = name + suffix
        else:
            name = name + suffix

        if max_length is not None and len(name) > max_length:
            # Chopping at a raw character index can land right after an
            # embedded apostrophe (from a base name like "Hi'iaka" or a
            # spliced-in UNIVERSAL_PHONEMES chunk like "ch'"), leaving the
            # truncated name ending in a bare "'" -- so back off past it.
            cutoff = max_length
            while cutoff > 0 and name[cutoff - 1] == "'":
                cutoff -= 1
            name = name[:cutoff] if cutoff > 0 else name[:max_length]

        name = name.lower()

        if is_name_valid(name):
            if allow_split:
                name = split_long_word(name)
                if any(part.lower() in _DECORATION_WORDS for part in name.split()):
                    # A split half that is itself a decoration word
                    # ("Xxxxx Ohana") -- see `_DECORATION_WORDS`.
                    continue
            name = name[0].upper() + name[1:]
            if "'" in name:
                # Two apostrophes can land adjacent here (e.g. a base name like
                # "Hi'iaka" splits into a syllable starting with "'", and a
                # spliced-in UNIVERSAL_PHONEMES chunk like "ch'" ends with "'";
                # shuffling can place them next to each other), producing an
                # empty string between them once split -- guard against
                # indexing that empty part rather than assuming every part is
                # non-empty.
                parts = name.split("'")
                name = "'".join([part[0].upper() + part[1:] if part else part for part in parts])
            return name

    raise RuntimeError(
        f"generate_phoneme_salad_name: no valid name found in {MAX_NAME_GENERATION_ATTEMPTS} attempts "
        f"(name_list={name_list!r}) -- is_name_valid is rejecting every candidate."
    )


def generate_sector_name():
    """
    Generates a random two-word sector name, each word independently drawn
    from the same phoneme-salad name generator used for star/planet/moon
    names -- using the sector-flavored `SECTOR_NAMES`/`SECTOR_PREFIXES`/
    `SECTOR_SUFFIXES` base lists instead, so generated sectors draw on real
    astronomical regions (galactic arms, superclusters, nebulae) and
    science-fiction sector names rather than reusing star names verbatim.
    No literal "Sector" suffix. `generate.py`'s `sector`/`galaxy` subcommands
    override this entirely via `--name`/`-n`, which hard-sets the whole
    name instead; `planetgen.db.store`'s name-uniqueness machinery
    (`nameUniqueness.py`) also calls this directly, to draw an entirely
    fresh sector name on the rare occasion a collision exhausts every
    decoration this project has for one -- both reasons this lives here,
    in `stellarObjects`, rather than in `generate.py` itself, which
    `store.py` can't import (it would be a backwards/circular dependency --
    `generate.py` already imports `planetgen.db.store`).

    Returns:
        str: A newly generated sector name, e.g. "Voranthis Kelmoor" --
        always exactly two words.
    """
    # allow_split=False: generate_phoneme_salad_name can itself split a
    # long result into two words (e.g. "Xyleth Anore"). Since this
    # function already joins two independent calls into one name, leaving
    # splitting on could silently produce 3-4 words instead of 2.
    # syllable_fraction=0.5 trims each word's base syllables by about
    # half before the prefix/suffix are attached -- many SECTOR_NAMES
    # entries (e.g. "Sagittarius", "Metropolis") are long real place
    # names, and two of them joined together made for unwieldy sector
    # names. max_length=7 backstops that: prefixes, suffixes, and the
    # occasional spliced-in universal phoneme are fixed-ish overhead that
    # doesn't shrink with syllable_fraction, so a long base name could
    # still slip through longer than intended without a hard cap too.
    first_word = generate_phoneme_salad_name(SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES, allow_split=False, syllable_fraction=0.5, max_length=7)
    second_word = generate_phoneme_salad_name(SECTOR_NAMES, SECTOR_PREFIXES, SECTOR_SUFFIXES, allow_split=False, syllable_fraction=0.5, max_length=7)
    return f"{first_word} {second_word}"


def to_paragraph(sentences):
    """
    Converts a list of sentences into a single, well-formed paragraph.

    This function is a simple utility to combine a list of strings into a
    single paragraph, which is used for generating descriptive text for
    celestial bodies.

    Args:
        sentences (list): A list of strings, where each string is a sentence.

    Returns:
        str: A single string representing the combined paragraph.
    """
    return " ".join(sentences)


def properties_to_string(system_config: SystemConfig, properties, template_name, markdown_header=None, markdown_key_map=None):
    """
    Converts a dictionary of properties into either a Markdown table or a Wiki template block,
    based on the `system_config.MARKDOWN` flag.

    Args:
        system_config (SystemConfig): The shared SystemConfig object.
        properties (dict): A dictionary where keys are property names (str) in Wikitext format
                           and values are their corresponding values (str or any type convertible to str).
        template_name (str): The name of the template to use for Wiki format (e.g., "Planet Data").
        markdown_header (str, optional): An optional header to prepend to the Markdown table.
                                         Defaults to None.
        markdown_key_map (dict, optional): A dictionary mapping Wikitext property keys to their
                                            desired Markdown table header names. If a key is not
                                            found in this map, it will be converted from
                                            lower_snake_case to Title Case for Markdown.

    Returns:
        str: A string formatted as either a Markdown table or a Wiki template block.
    """
    output_lines = []

    if system_config.MARKDOWN:
        if markdown_header:
            output_lines.append(markdown_header)
        output_lines.extend([
            "| Property | Value |",
            "|---|---|"
        ])
        for prop_key, value in properties.items():
            # Determine the Markdown header name
            if markdown_key_map and prop_key in markdown_key_map:
                markdown_prop_name = markdown_key_map[prop_key]
            else:
                # Default to converting lower_snake_case to Title Case
                markdown_prop_name = prop_key.replace('_', ' ').title()
            output_lines.append(f"| {markdown_prop_name} | {value} |")
        return "\n".join(output_lines)
    else:
        output_lines.append(f"{{{{{template_name}")
        for prop_key, value in properties.items():
            # For Wikitext, the keys are already in the correct format
            output_lines.append(f"|{prop_key}={value}")
        output_lines.append("}}")
        return "\n".join(output_lines)