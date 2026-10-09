# planetgen/util/format.py

"""
Format
======

The text every number, distance, speed, duration, temperature, pressure
and age is shown in, for the generator's wikitext and Markdown and for the
web pages (`planetgen.web.lib.fmt` wraps these with its HTML dash for a
missing value), and the property table every body writes.
"""

import decimal
import math
import re

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.physics import constants as physical_constants


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


_FIXED_SPEC_RE = re.compile(r"^(,?)\.(\d+)f$")


_ROUNDING_CONTEXT = decimal.Context(prec=400)
"""decimal.Context: Precise enough for any float (1e308 with a few decimals)."""


def round_half_up(value, decimals):
    """
    `value` rounded to `decimals` places, ties away from zero, judged on the
    shortest decimal that reads back as `value` (UX.79): 9.995 gives 10.00,
    1.005 gives 1.01, as the browser's own rounding does. Python's own
    `round` and `format` go by the binary value instead (9.995 is really
    9.99499...), so the two copies of a number disagreed.

    Returns:
        decimal.Decimal: Finite `value`, rounded.
    """
    return decimal.Decimal(repr(float(value))).quantize(
        decimal.Decimal(1).scaleb(-decimals), rounding=decimal.ROUND_HALF_UP, context=_ROUNDING_CONTEXT)


def _no_negative_zero(text):
    """`text` without a minus sign when every digit in it is zero ("-0" is "0", UX.80)."""
    if text.startswith("-") and not any(char in "123456789" for char in text):
        return text[1:]
    return text


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
    fixed = _FIXED_SPEC_RE.match(spec)
    if fixed and math.isfinite(value):
        text = _no_negative_zero(format(round_half_up(value, int(fixed.group(2))), spec))
    else:
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
    if value == 0:
        return "0"
    if not math.isfinite(value):
        return f"{value:g}"
    magnitude = math.floor(math.log10(abs(value)))
    decimals = max(0, 2 - magnitude)
    text = f"{round_half_up(value, decimals):,.{decimals}f}"
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
