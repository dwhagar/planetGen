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
`markdown_content`, safely un-escaped by `mdconvert.py`), not the wikitext
"{{Exp|coeff|exp}}" template form -- the wrong one for embedding directly
into an HTML page, which was the whole bug: the interactive HTML viewer used
to read the wikitext form straight out of the database.
"""

try:
    from stellarObjects import physical_constants, program_constants
    from stellarObjects.config import SystemConfig
    from stellarObjects.utils import (
        format_length_km, format_relative_to_sol, to_scientific_notation, years_to_time_string,
    )

    _HTML_CONFIG = SystemConfig()
    _HTML_CONFIG.MARKDOWN = True
except ImportError:
    # The planetGen package isn't on the import path in this deployment --
    # fall back to showing the raw number rather than failing outright (see
    # `html/sector.py`'s identical fallback for `milliparsecs_to_ly`).
    _HTML_CONFIG = None
    years_to_time_string = None


def format_star_mass(mass_kg):
    if _HTML_CONFIG is None:
        return f"{mass_kg} kg"
    return format_relative_to_sol(_HTML_CONFIG, mass_kg, physical_constants.SOLAR_MASS_TO_KG, "kg", low_percent_precision=2)


def format_star_luminosity(luminosity_w):
    if _HTML_CONFIG is None:
        return f"{luminosity_w} W"
    return format_relative_to_sol(_HTML_CONFIG, luminosity_w, physical_constants.SOLAR_LUMINOSITY, "W", low_percent_precision=4)


def format_star_radius(radius_km):
    if _HTML_CONFIG is None:
        return f"{radius_km} km"
    return format_length_km(
        _HTML_CONFIG, radius_km,
        program_constants.RADIUS_KM_SCIENTIFIC_NOTATION_THRESHOLD,
        program_constants.ROUND_RADIUS_KM,
        program_constants.SCIENTIFIC_NOTATION_DECIMAL_PLACES,
    )


def format_period(period_years):
    if years_to_time_string is None:
        return f"{period_years} years"
    return years_to_time_string(period_years)


def format_body_distance(distance_km, is_moon):
    """
    Mirrors the `distance_text` branch in
    `stellarObjects.planetData.Planet.get_table_properties` -- a moon's
    distance (from its parent planet) is always shown in km; a top-level
    planet's distance (from its star) is shown in AU (plus km for context
    under 1 AU), or light-years above `program_constants.LY_THRESHOLD`.

    Args:
        distance_km (float): `planets`/`moons`.`distance_km`.
        is_moon (bool): Whether this row came from the `moons` table.

    Returns:
        str: The formatted distance.
    """
    if _HTML_CONFIG is None:
        return f"{distance_km} km"

    if is_moon:
        return f"{to_scientific_notation(_HTML_CONFIG, distance_km, 4)} km"

    distance_au = distance_km / physical_constants.AU_TO_KM
    distance_ly = distance_au * physical_constants.AU_TO_LY
    if distance_ly < program_constants.LY_THRESHOLD:
        if distance_au < 1:
            return f"{to_scientific_notation(_HTML_CONFIG, distance_km, 1)} km ({distance_au:.3f} AU)"
        return f"{distance_au:.3f} AU"
    return f"{distance_ly:.4f} light-years"
