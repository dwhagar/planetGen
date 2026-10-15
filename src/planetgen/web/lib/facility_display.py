# planetgen/web/lib/facility_display.py

"""
How a facility's type, affiliation, ratings and field values read on a page
(POP.10): the rating chips and the small text cells of the Facilities table.

A rating chip never relies on colour alone (WCAG 1.4.1): each carries a
distinct shape (a circle with a check, a triangle, an octagon with a cross),
the dimension's name and its word ("Elevated"), and an `aria-label` that
reads "Crime: Yellow, elevated". The fills are the Okabe-Ito palette the
design note settles on, with black text (`.rating` in `style.css`).
"""

from planetgen.population import facility_ratings as ratings
from planetgen.web.lib.fmt import esc

DIMENSION_LABELS = {"crime": "Crime", "housing": "Housing", "resources": "Resources",
                    "maintenance": "Maintenance", "health": "Health"}
"""dict: What each dimension is called on a chip."""

COLOR_ICONS = {"green": "check-circle-fill", "yellow": "exclamation-triangle-fill", "red": "x-octagon-fill"}
"""dict: The status shape of each colour: three silhouettes that differ without colour."""

DIMENSION_ICONS = {"crime": "shield-check", "housing": "house-door", "resources": "box-seam",
                   "maintenance": "wrench-adjustable", "health": "heart-pulse"}
"""dict: The dimension glyph shown on a chip."""


def rating_chip(dimension, color, word):
    """One chip: the dimension's glyph, its name, the status shape and the word."""
    label = DIMENSION_LABELS[dimension]
    return (f'<li class="rating rating-{color}" aria-label="{esc(label)}: {color.capitalize()}, {esc(word.lower())}">'
            f'<sl-icon name="{DIMENSION_ICONS[dimension]}" aria-hidden="true"></sl-icon>'
            f'<span class="rating-name">{esc(label)}</span>'
            f'<sl-icon name="{COLOR_ICONS[color]}" aria-hidden="true"></sl-icon>'
            f'<span class="rating-word">{esc(word)}</span></li>')


def rating_chips(facility):
    """
    The five chips of a facility as one list, worst first, or `""` when none is assessed. `facility` is a
    facility as `GET /api/facilities/<id>` returns it.
    """
    chips = []
    for dimension in ratings.RATINGS:
        color = facility["ratings"][dimension]
        if color is not None:
            chips.append((ratings.COLORS.index(color), dimension, color, facility["rating_words"][dimension]))
    if not chips:
        return ""
    chips.sort(key=lambda chip: -chip[0])        # the stable sort keeps the dimensions' order within a colour
    return ('<ul class="rating-strip" aria-label="Ratings">'
            + "".join(rating_chip(dimension, color, word) for _, dimension, color, word in chips) + "</ul>")


def type_html(facility):
    """The facility's type: its icon and name, or `""` for an untyped one."""
    if facility.get("type_id") is None:
        return ""
    icon = f'<sl-icon name="{esc(facility["type_icon"])}" aria-hidden="true"></sl-icon> ' if facility.get("type_icon") else ""
    return f'{icon}{esc(facility["type_name"])}'


def affiliation_text(facility):
    """The affiliation's name (the snapshot kept when the polity is gone), or `None`."""
    return facility.get("affiliation_name")


def fields_text(facility):
    """The field values as `key: value` text for a tooltip or a detail line, in the order the API gave them."""
    parts = []
    for key, value in (facility.get("fields") or {}).items():
        shown = ("yes" if value else "no") if isinstance(value, bool) else value
        parts.append(f"{key}: {shown}")
    return ", ".join(parts)
