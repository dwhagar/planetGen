# planetgen/population/facility_ratings.py

"""
The five ratings of a facility (POP.10): crime, housing, resources,
maintenance and health, each Green, Yellow or Red.

Stored as 1 (Green), 2 (Yellow), 3 (Red) or NULL (not assessed); the API
speaks the colour words. Green is the good end of all five dimensions, so
the worst of the five (the facility's overall rating) is the highest number.
"""

RATINGS = ("crime", "housing", "resources", "maintenance", "health")
"""tuple: The five dimensions, in display order (columns `rating_<name>`)."""

COLORS = ("green", "yellow", "red")
"""tuple: The rating words, best first; a stored 1, 2 or 3."""

WORDS = {
    "crime": ("Low", "Elevated", "High"),
    "housing": ("Ample", "Tight", "Overcrowded"),
    "resources": ("Sufficient", "Strained", "Critical"),
    "maintenance": ("Sound", "Deferred", "Failing"),
    "health": ("Healthy", "At risk", "Crisis"),
}
"""dict: What each colour means for each dimension, as a chip's text."""

SOURCES = ("admin", "generated")
"""tuple: Every `facilities.ratings_source`."""


def color_of(value):
    """The colour word of a stored rating (`None` for not assessed)."""
    return None if value is None else COLORS[value - 1]


def value_of(color):
    """The stored rating of a colour word; `ValueError` for anything else."""
    if color not in COLORS:
        raise ValueError(f"a rating is one of: {', '.join(COLORS)}")
    return COLORS.index(color) + 1


def overall(values):
    """The worst of the ratings given (stored numbers, `None` for not assessed) as a colour word, or `None` when
    none is assessed."""
    assessed = [value for value in values if value is not None]
    return color_of(max(assessed)) if assessed else None


def describe(row):
    """
    The ratings of one `facilities` row as the API returns them.

    Returns:
        dict: `ratings` (`{dimension: colour or None}`), `rating_words` (`{dimension: word or None}`),
            `overall` (the worst colour or `None`), `ratings_source`, `ratings_note` and `ratings_set_at`.
    """
    values = {name: row[f"rating_{name}"] for name in RATINGS}
    return {
        "ratings": {name: color_of(value) for name, value in values.items()},
        "rating_words": {name: None if value is None else WORDS[name][value - 1] for name, value in values.items()},
        "overall": overall(values.values()),
        "ratings_source": row["ratings_source"],
        "ratings_note": row["ratings_note"],
        "ratings_set_at": row["ratings_set_at"],
    }
