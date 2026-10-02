# html/lib/tabledisplay.py

"""
Computes the same "Star Data"/"Planet Data" display strings
`Star`/`Planet`/`BinaryStarProxy.get_table_properties()` used to bake into
the database's now-removed `table_*`/`binary_table_*` columns (see
`schema.sql`'s "v5" header note) -- but on demand, from the raw numeric
columns (`mass_kg`, `radius_km`, `luminosity_w`, `distance_km`, ...) that
were always stored alongside them and still are. Reuses
`stellarObjects.utils`' formatters directly rather than reimplementing the
scientific-notation/Sol-relative math; only the small per-quantity branching
(which unit to show, at what threshold) is mirrored here, operating on plain
numbers instead of a live `Star`/`Planet` object's attributes.

`_HTML_CONFIG.MARKDOWN = True` makes those formatters emit the HTML
"coeff &times; 10<sup>exp</sup>" form (the same form already used for
rendered system Markdown, safely un-escaped by `mdconvert.py`), not the wikitext
"{{Exp|coeff|exp}}" template form -- the wrong one for embedding directly
into an HTML page, which was the whole bug: the interactive HTML viewer used
to read the wikitext form straight out of the database.

That HTML `<sup>` form is only safe to embed where it's actually parsed as
markup (`system.py`'s static table cells, inserted unescaped -- see its own
comment on why). `lib/systemmap.py`'s interactive map instead carries these
same formatted strings through `data-*` attributes that `static/systemmap.js`
reads back with `.textContent` (deliberately never `innerHTML`, to keep every
database-derived value safely un-executable) -- `.textContent` shows tags
literally rather than rendering them, so a raw "10<sup>7</sup>" landed on
screen as the literal text "10<sup>7</sup>" instead of a superscript 7. Any
caller feeding one of these formatters' output into a `data-*` attribute
(never a static HTML page fragment) must run it through `to_plain_text`
first.
"""

import re

from fmt import dash_unless_finite

_SUP_HTML_RE = re.compile(r'<sup>(-?\d+)</sup>')
_SUPERSCRIPT_DIGITS = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def to_plain_text(formatted):
    """
    Converts one of this module's HTML-formatted strings (their one raw-HTML
    pattern, `<sup>exponent</sup>`, from `to_scientific_notation`) into a
    plain-text equivalent using real Unicode superscript digits -- for a
    caller that will hand the result to something that displays it as text
    rather than parsing it as markup (e.g. a `data-*` attribute read back via
    `.textContent`, see this module's own docstring). A no-op for a string
    that never had a `<sup>` in it.

    Args:
        formatted (str): This module's own formatter output.

    Returns:
        str: The same text with every `<sup>N</sup>` replaced by N's
             Unicode-superscript digits.
    """
    return _SUP_HTML_RE.sub(lambda m: m.group(1).translate(_SUPERSCRIPT_DIGITS), formatted)

try:
    from stellarObjects import physical_constants
    from stellarObjects.config import SystemConfig
    from stellarObjects.utils import (
        format_body_radius_km, format_distance_km, format_period_years, format_relative_to_sol,
        to_scientific_notation,
    )

    _HTML_CONFIG = SystemConfig()
    _HTML_CONFIG.MARKDOWN = True
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # fall back to showing the raw number rather than failing outright (see
    # `html/sector.py`'s identical fallback for `milliparsecs_to_ly`).
    _HTML_CONFIG = None
    format_period_years = None


def _dash_unless_positive(formatter):
    """`dash_unless_finite(formatter)`, and a dash for a negative value too:
    a negative mass or luminosity is no real value, and `format_relative_to_sol`
    would spell -1e300 out as a 300-digit percentage (TEST.53)."""
    guarded = dash_unless_finite(formatter)

    def positive(value, *args, **kwargs):
        return "\u2013" if isinstance(value, (int, float)) and value < 0 else guarded(value, *args, **kwargs)
    positive.__name__ = formatter.__name__
    positive.__doc__ = formatter.__doc__
    return positive


def format_star_mass(mass_kg):
    if _HTML_CONFIG is None:
        return f"{mass_kg} kg"
    return format_relative_to_sol(_HTML_CONFIG, mass_kg, physical_constants.SOLAR_MASS_TO_KG, "kg", low_percent_precision=2)


def format_star_luminosity(luminosity_w):
    if _HTML_CONFIG is None:
        return f"{luminosity_w} W"
    return format_relative_to_sol(_HTML_CONFIG, luminosity_w, physical_constants.SOLAR_LUMINOSITY, "W", low_percent_precision=4)


def format_star_radius(radius_km):
    """A planet, moon or star radius: always km in scientific notation
    (Boss, 2026-09-30), unlike every other distance on the site."""
    if _HTML_CONFIG is None:
        return f"{radius_km} km"
    return format_body_radius_km(_HTML_CONFIG, radius_km)


def format_period(period_years):
    """An orbital period, in years, on the shared period ladder
    (`stellarObjects.utils.format_period_years`, UX.14)."""
    if format_period_years is None:
        return f"{period_years} years"
    return format_period_years(period_years)


def format_body_distance(distance_km, is_moon=False):
    """
    A planet's distance from its star or a moon's from its planet, through
    the distance ladder (`utils.format_distance_m`): km under 1 AU, then AU,
    and so on up. `is_moon` is kept for callers; moons and planets now
    follow the same ladder.

    Args:
        distance_km (float): `planets`/`moons`.`distance_km`.
        is_moon (bool): Whether this row came from the `moons` table.

    Returns:
        str: The formatted distance.
    """
    if _HTML_CONFIG is None:
        return f"{distance_km} km"
    return format_distance_km(distance_km)


def _times_reference(ratio):
    """A body's size or mass as a multiple of Earth's or Jupiter's: plain
    digits from 0.001 to a million, scientific notation outside that."""
    if ratio >= 1e6:
        return to_scientific_notation(_HTML_CONFIG, ratio, 2)
    if ratio >= 100:
        return f"{ratio:,.0f}"
    if ratio >= 1:
        return f"{ratio:.2f}"
    if ratio >= 0.001:
        return f"{ratio:.3g}"
    return to_scientific_notation(_HTML_CONFIG, ratio, 2)


def format_body_radius(radius_km):
    """A planet's or moon's radius (MAP.92): km in scientific notation like
    every body radius, then in Earth radii."""
    if _HTML_CONFIG is None:
        return f"{radius_km} km"
    ratio = radius_km / physical_constants.EARTH_RADIUS_KM
    return f"{format_body_radius_km(_HTML_CONFIG, radius_km)} ({_times_reference(ratio)} Earth radii)"


def format_body_mass(mass_kg, gas_giant=False):
    """A planet's or moon's mass (MAP.92): kg in scientific notation, then
    in Earth masses, or Jupiter masses for a gas giant."""
    if _HTML_CONFIG is None:
        return f"{mass_kg} kg"
    reference, label = ((physical_constants.JUPITER_MASS_TO_KG, "Jupiter masses") if gas_giant
                        else (physical_constants.EARTH_MASS_TO_KG, "Earth masses"))
    return f"{to_scientific_notation(_HTML_CONFIG, mass_kg)} kg ({_times_reference(mass_kg / reference)} {label})"


# TEST.53: `None`, NaN, an infinity or a number past float range shows a
# dash, never "nan km" or a traceback.
format_star_mass = _dash_unless_positive(format_star_mass)
format_star_luminosity = _dash_unless_positive(format_star_luminosity)
format_star_radius = dash_unless_finite(format_star_radius)
format_period = dash_unless_finite(format_period)
format_body_distance = dash_unless_finite(format_body_distance)
format_body_radius = _dash_unless_positive(format_body_radius)
format_body_mass = _dash_unless_positive(format_body_mass)
