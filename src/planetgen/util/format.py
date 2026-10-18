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
from planetgen.physics import constants as physical_constants
from planetgen.util import ladders


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


def _in_rung(rung, value):
    """`value` (in base units) shown in `rung`'s unit, with the unit's
    symbol: three significant figures, or a whole number from the rung's
    `whole_from`; the singular form when the number is exactly 1."""
    number = (value - rung.offset) / rung.factor
    if rung.whole_from is not None and abs(value) >= rung.whole_from:
        return f"{format_number(number)} {rung.symbol}"
    text = three_figures_or_tiny(number)
    return f"{text} {rung.singular if rung.singular and text == '1' else rung.symbol}"


def pick_rung(ladder, value):
    """
    The rung `value` (in the ladder's base unit) is shown in: the largest
    one it is at least `start` of, with a hair of tolerance so a value
    computed as exactly 1 unit (e.g. 1 pc from ly * LY_TO_AU / AU_PER_PARSEC)
    isn't shown as 1,000 of the one below it. Zero is in the ladder's `zero`
    rung, the lowest when it names none.
    """
    size = abs(value)
    if size == 0 and ladder.zero is not None:
        return ladder.rung(ladder.zero)
    rung = ladder.rungs[0]
    for candidate in ladder.rungs:
        start = candidate.factor if candidate.start is None else candidate.start
        if size >= start * (1 - 1e-12):
            rung = candidate
    return rung


def format_quantity(ladder, value, companions=True, system=ladders.AUTOMATIC):
    """
    The site's one unit formatter (UX.23): `value`, in the ladder's base
    unit, in the most meaningful unit on `ladder` (`planetgen.util.ladders`),
    to three significant figures in `format_number`'s style. A rung that
    carries an aside adds one more familiar unit in parentheses, and the
    ladder's companion units (°C and °F beside K, atm and psi beside Pa)
    follow in their own parentheses unless `companions` is false. `system`
    is the reader's unit preference (UX.78): "auto" (the ladders as they
    are), "metric" or "customary".
    `static/unitladder.js` mirrors it; tests/test_unit_ladders.py checks the
    two against one table.

    Returns:
        str: e.g. "384,400 km", "4.2 pc (13.7 ly)", "27.3 days". `None` and
             NaN give an en dash; an infinity "inf" or "-inf" in the
             ladder's `nonfinite` unit.
    """
    if value is None:
        return "–"
    value = float(value)
    if math.isnan(value):
        return "–"
    if math.isinf(value):
        return f"{value:g} {ladder.nonfinite}"
    ladder = ladder.for_system(system)
    rung = pick_rung(ladder, value)
    text = _in_rung(rung, value)
    if rung.with_aside:
        text += f" ({_in_rung(ladder.rung(_aside_symbol(ladder, value)), value)})"
    if companions:
        shown = [f"{_companion_number(companion, value)} {companion.symbol}"
                 for companion in ladder.companions if _companion_applies(companion, value)]
        if shown:
            text += f" ({', '.join(shown)})"
    return text


def _companion_applies(companion, value):
    """True when a companion unit is shown beside `value`: its size is in the companion's range."""
    size = abs(value)
    if companion.min_size is not None and size < companion.min_size:
        return False
    return companion.max_size is None or size < companion.max_size


def _companion_number(companion, value):
    """`value` in a companion unit, without the symbol."""
    number = (value - companion.offset) / companion.factor
    return _whole_or_tenths(number) if companion.whole else three_figures_or_tiny(number)


def _aside_symbol(ladder, value):
    """The symbol of the rung a `with_aside` value's parenthetical is in."""
    size = abs(value)
    for start, symbol in ladder.asides:
        if size >= start:
            return symbol
    return ladder.asides[-1][1]


def format_distance_m(meters):
    """
    Formats a distance (or a non-body radius) in the most meaningful unit on
    the ladder km < AU < mpc < cpc < ly < pc < kpc < Mpc < Gpc: the largest
    unit the value is at least 1 of. Values in a parsec unit add one
    parenthetical: ly when the distance is at least 0.01 ly, else AU when at
    least 0.01 AU, else km (`ladders.LENGTH`).

    Every page and text output passes its distances through this (or its
    `format_distance_km`/`_au`/`_ly`/`_pc` variants) rather than formatting
    them itself; `static/unitladder.js` mirrors it for the maps. Planet,
    moon and star radii are the exception: see `format_body_radius_km`.

    Args:
        meters (float): The distance, in meters.

    Returns:
        str: e.g. "384,400 km", "1.52 AU", "4.2 pc (13.7 ly)",
             "2.4 mpc (495 AU)". `None` gives an en dash.
    """
    return format_quantity(ladders.LENGTH, meters)


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


def format_speed_kms(kms):
    """
    Formats a speed in the most meaningful unit on the ladder km/h < km/s <
    Mm/s < c (`ladders.SPEED`; Boss, 2026-10-01): km/h below 1 km/s, km/s up
    to 1,000 km/s, Mm/s up to a tenth of light speed, then multiples of c.
    One unit, three significant figures, in `format_distance_m`'s number
    style.

    Every page and text output passes its plain speeds through this rather
    than formatting them itself; `static/unitladder.js` mirrors it. Warp and
    fold factors (`navigation.py`) keep their own "x c" columns.

    Args:
        kms (float): The speed, in km/s.

    Returns:
        str: e.g. "36 km/h", "29.8 km/s", "4.5 Mm/s", "0.25 c". `None`
             gives an en dash.
    """
    return format_quantity(ladders.SPEED, None if kms is None else kms * physical_constants.KM_M)


def format_speed_ms(ms):
    """`format_speed_kms` for a value in m/s."""
    return format_quantity(ladders.SPEED, ms)


def format_duration_seconds(seconds):
    """
    Formats a time period or duration in the most meaningful unit on the
    ladder µs < ms < s < minutes < hours < days < years < kyr < Myr < Gyr
    (`ladders.TIME`): the largest unit the value is at least 1 of, as
    `format_distance_m` does. One unit, three significant figures, singular
    when the shown value is exactly 1 ("1 year", "1 day").
    `static/unitladder.js` mirrors it.

    Args:
        seconds (float): The duration, in seconds.

    Returns:
        str: e.g. "12.5 ms", "45 minutes", "27.3 days", "1.88 years",
             "236 Myr". Zero is "0 s"; `None` gives an en dash.
    """
    return format_quantity(ladders.TIME, seconds)


def format_period_years(years):
    """`format_duration_seconds` for a value in (Julian) years -- every
    orbital period the generator stores is in years or Gyr."""
    return format_duration_seconds(None if years is None else years * physical_constants.SECONDS_PER_YEAR)


def _whole_or_tenths(value):
    """`value` to whole units, or to one decimal when under 10 in size
    ("15", "-40", "5.9", "0"), in `format_number`'s style."""
    text = format_number(value, ",.1f" if abs(value) < 10 else ",.0f")
    if text.endswith(".0"):
        text = text[:-2]
    return "0" if text in ("-0", "+0") else text


def format_temperature_k(kelvin, companions=True):
    """
    A surface temperature in kelvin with Celsius and Fahrenheit alongside
    (Boss, 2026-10-01: "temperatures should be reported in K, C, and F for
    surface conditions"). Kelvin to three significant figures in
    `format_distance_m`'s number style (mK, µK and nK under 1 K); °C and °F
    to whole degrees, or one decimal when under 10 in size. A star's
    temperature, or any that is not a surface condition, stays in K alone:
    pass `companions=False`.

    Args:
        kelvin (float): The temperature, in K.
        companions (bool, optional): Show °C and °F too. Defaults to True.

    Returns:
        str: e.g. "288 K (15 °C, 59 °F)", "737 K (464 °C, 867 °F)". `None`
             gives an en dash.
    """
    return format_quantity(ladders.TEMPERATURE, kelvin, companions)


def three_figures_or_tiny(value):
    """`_three_figures`, but scientific below a thousandth, where a
    unit would otherwise run to a row of zeros (UX.81)."""
    if value != 0 and abs(value) < ladders.TINY_NUMBER:
        return scientific_text(value)
    return _three_figures(value)


def format_pressure_pa(pascals, companions=True):
    """
    A pressure with a metric primary on the nPa to TPa ladder and customary
    atm and psi alongside, all to three significant figures (Boss,
    2026-10-01: "Atmospheric pressure and surface conditions should show
    customary units as well as a secondary to help contextualize the metric
    values given"); bar joins them from 10 atm to 10,000 atm and mbar from 1 Pa
    to 0.01 atm.

    Args:
        pascals (float): The pressure, in Pa (the schema's
            `atmospheric_pressure_pa`).
        companions (bool, optional): Show atm, psi and the rest too.
            Defaults to True.

    Returns:
        str: e.g. "101 kPa (1 atm, 14.7 psi)", "9.2 MPa (90.8 atm,
             1,334 psi, 92 bar)". `None` gives an en dash.
    """
    return format_quantity(ladders.PRESSURE, pascals, companions)


def format_pressure_atm(atm):
    """`format_pressure_pa` for a value in standard atmospheres."""
    return format_pressure_pa(None if atm is None else atm * physical_constants.STANDARD_ATMOSPHERE_PA)


def format_gravity_g(gravity_g):
    """
    A surface gravity in m/s² with standard gravity beside it, "3.71 m/s²
    (0.378 g)", for property grids and generator text; a bare "g" would
    collide with the gram, so it never stands alone there. A list row uses
    `format_gravity_g_short`.

    Args:
        gravity_g (float): The gravity, in multiples of standard gravity
            (the schema's `gravity_g`).

    Returns:
        str: `None` gives an en dash.
    """
    return format_quantity(ladders.GRAVITY, None if gravity_g is None else gravity_g * physical_constants.EARTH_GRAVITY)


def format_gravity_g_short(gravity_g):
    """`gravity_g` as "0.378 g" alone, for a list row's chip (UX.22: rows and
    chips show the primary unit; the m/s² is in the detail panel)."""
    if gravity_g is None:
        return "–"
    gravity_g = float(gravity_g)
    if not math.isfinite(gravity_g):
        return "–" if math.isnan(gravity_g) else f"{gravity_g:g} g"
    return f"{three_figures_or_tiny(gravity_g)} g"


def format_mass_kg(kilograms):
    """A mass in kilograms on the mass ladder (UX.22): "5.97 × 10²⁴ kg" is
    shown as "1 Earth mass", "0.299 Jupiter masses", "0.122 solar masses"."""
    return format_quantity(ladders.MASS, kilograms)


def format_mass_solar(solar_masses):
    """`format_mass_kg` for a mass held in solar masses (black holes, neutron
    stars, the quasar's black hole)."""
    return format_quantity(ladders.MASS, None if solar_masses is None else solar_masses * physical_constants.SOLAR_MASS_TO_KG)


def format_power_w(watts):
    """A luminosity or power in watts on the power ladder (UX.22): "47.2 TW",
    "1 L☉", "9.01 × 10⁻³¹ W"."""
    return format_quantity(ladders.POWER, watts)


def format_luminosity_sol(solar_luminosities):
    """`format_power_w` for a luminosity held in solar luminosities."""
    return format_quantity(ladders.POWER, None if solar_luminosities is None
                           else solar_luminosities * physical_constants.SOLAR_LUMINOSITY)


def format_flux_w_m2(flux):
    """A heat flow or flux in W/m² on the flux ladder (UX.22): "87 mW/m²"."""
    return format_quantity(ladders.FLUX, flux)


def format_number_density_cm3(per_cm3):
    """A particle density in particles per cubic centimeter (clouds,
    remnants): "3.5 particles/cm³". `None` gives an en dash."""
    return format_quantity(ladders.NUMBER_DENSITY, None if per_cm3 is None else per_cm3 * 1e6)


def format_magnetic_field_gauss(gauss):
    """A magnetic field held in gauss (the schema's `magnetic_field_gauss`) on
    the field ladder: "50 µT" for a planet, "1.00 × 10¹² G" for a pulsar."""
    return format_quantity(ladders.MAGNETIC_FIELD, None if gauss is None else gauss * 1e-4)


def format_body_radius_km(radius_km, precision=None):
    """
    A planet, moon or star radius: always km in scientific notation (Boss,
    2026-09-30), the one exception to `format_distance_m`'s ladder.

    Args:
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
    return f"{to_scientific_notation(radius_km, precision)} km"


def to_scientific_notation(number, precision=2):
    """
    Converts a number to scientific notation with the specified precision.

    This function is used for formatting large numbers in a compact and
    standardized way, as HTML (`1.23 × 10<sup>4</sup>`), the form the pages embed.

    Args:
        number (float): The number to convert.
        precision (int, optional): The number of decimal places to show. Defaults to 2.

    Returns:
        str: The number in scientific notation as HTML.
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
    return f"{mantissa} × 10<sup>{exponent}</sup>"

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
