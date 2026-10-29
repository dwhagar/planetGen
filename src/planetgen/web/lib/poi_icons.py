# planetgen/web/lib/poi_icons.py

"""
The point-of-interest symbols (MAP.320): the one place that names the seven vector icons of `static/poi/` and builds
their markup, so a page, a legend row and a list row draw the same shape, colour and words for a kind.

Each icon is one SVG file (`static/poi/<id>.svg`, a solid silhouette with its name in its `<title>`); this module puts
one on a page as inline `<svg>` that references the file (`<use href="poi/<id>.svg#poi">`) and carries its words
(`role="img"`, an `aria-label`), so a kind never depends on colour alone. The paint comes from the `poi-*` classes of
`static/style.css`: a world the habitability score colours (the habitable and the inhabited one) takes the class of its
PHI-6 tier, the rest the token of their kind. `static/poiicons.js` lists the same icons for the map, and
`src/tests/test_poi_icons.py` holds the two lists and the files to each other. A page reaches this through the Jinja
global `poi_icon(id, tier)`.

Imports no Flask, so the describe renderer and the list builders can use it.
"""

from collections import namedtuple
from html import escape

from planetgen.physics.habitability_world import TIER_NAMES

PoiIcon = namedtuple("PoiIcon", "id label shape meaning scored")
"""A symbol: its file's name `id`, its accessible name `label`, what its `shape` looks like (the sign that does not
depend on colour), what it `meaning` marks, and whether the habitability score colours it (`scored`)."""

ICONS = (
    PoiIcon("habitable-world", "Habitable world", "a disc with a hill cut out",
            "A planet or moon whose overall habitability (PHI-6) is Blue or Green, coloured by that score.", True),
    PoiIcon("uninhabitable-world", "Uninhabitable world", "a disc pitted with craters",
            "A planet or moon that is not habitable, in brown, never red.", False),
    PoiIcon("inhabited-world", "Inhabited world", "a disc inside a ring",
            "A world with people on it, coloured by its habitability score.", True),
    PoiIcon("gas-giant", "Gas giant", "a disc crossed by three bands",
            "A gas giant.", False),
    PoiIcon("asteroid-field", "Asteroid field", "four scattered rocks",
            "An asteroid belt or field (a belt includes the fields).", False),
    PoiIcon("sector", "Sector", "a square frame with a dot in it",
            "A sector, one cell of the map of space.", False),
    PoiIcon("region", "Region", "a hexagonal frame",
            "A region, a named stretch of charted space.", False),
)
"""tuple: The seven symbols, in the order the legend lists them."""

BY_ID = {icon.id: icon for icon in ICONS}
"""dict: `ICONS` by id."""

SPRITE_ID = "poi"
"""str: The id of the root `<svg>` of every icon file, what a `<use>` points at."""

TIER_CLASSES = tuple(f"poi-tier-{tier}" for tier in range(len(TIER_NAMES)))
"""tuple: The CSS class of each habitability colour by tier number (0 Blue to 3 Red)."""

DEFAULT_TIER = 1
"""int: The tier a score-coloured symbol is drawn in when the caller has none (Green, the habitable colour)."""


def symbol_html(icon_id, file_url, tier=None):
    """
    The `<svg>` of one symbol as markup. `file_url` is a function from a path under `static/` ("poi/region.svg") to its URL
    (the site's `static_url`, which carries the release); `tier` is the PHI-6 tier (0 Blue to 3 Red) that colours a symbol the score colours, and is
    ignored by the others (None draws a scored symbol Green).

    Raises:
        ValueError: `icon_id` is not one of `ICONS`, or `tier` is not 0 to 3.
    """
    if icon_id not in BY_ID:
        raise ValueError(f"no point-of-interest icon named {icon_id!r}")
    icon = BY_ID[icon_id]
    classes = ["poi-icon", f"poi-{icon.id}"]
    if icon.scored:
        chosen = DEFAULT_TIER if tier is None else tier
        if chosen not in range(len(TIER_CLASSES)):
            raise ValueError(f"no habitability tier {tier!r}")
        classes.append(TIER_CLASSES[chosen])
    label = icon.label if not icon.scored or tier is None else f"{icon.label}: {TIER_NAMES[tier]}"
    return (f'<svg class="{" ".join(classes)}" role="img" aria-label="{escape(label)}" data-poi="{icon.id}" '
            f'focusable="false"><use href="{escape(file_url(f"poi/{icon.id}.svg"))}#{SPRITE_ID}"></use></svg>')
