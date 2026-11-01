# planetgen/web/lib/galaxy_outlines.py

"""
The galaxy family outlines (UX.135): the one place that names the six vector drawings of `static/galaxy/` and builds their markup,
so the Generate page's Shape row and the Galaxy Map's key draw the same shape for a family.

Each drawing is one SVG file (`static/galaxy/<id>.svg`, a plain outline with its name in its `<title>`, one stroke weight, no fill,
stroke in `currentColor` so both themes work). They are made-up shapes in the view that tells the family apart, not a copy of any
named galaxy: the arms face-on for the spiral, the arms leaving a bar for the barred spiral, a smooth pointed lens for the
lenticular, an ellipse with an inner contour for the elliptical, a ragged body with clumps for the irregular, a faint dashed round
blob for the dwarf spheroid. The family keys are those of `docs/design/galaxy/galaxy-types.md` section 3; a barred spiral is the
`spiral` family with `bar` set, so it has an outline of its own but no key of its own (`outline_id`).

This module puts one on a page as inline `<svg>` that references the file (`<use href="galaxy/<id>.svg#galaxy">`) and carries its
words (`role="img"`, an `aria-label`). It imports no Flask; a page reaches it through the Jinja global `galaxy_outline(id)`.
"""

from collections import namedtuple
from html import escape

GalaxyOutline = namedtuple("GalaxyOutline", "id label family shape")
"""An outline: its file's name `id`, its accessible name `label`, the family key it draws, and what its `shape` looks like."""

OUTLINES = (
    GalaxyOutline("spiral", "Spiral galaxy", "spiral", "a small core and two arms winding out from it"),
    GalaxyOutline("barred_spiral", "Barred spiral galaxy", "spiral", "a short bar with an arm leaving each end"),
    GalaxyOutline("lenticular", "Lenticular galaxy", "lenticular", "a smooth pointed lens with a bulge in the middle"),
    GalaxyOutline("elliptical", "Elliptical galaxy", "elliptical", "a tilted ellipse with an inner contour"),
    GalaxyOutline("irregular", "Irregular galaxy", "irregular", "a ragged body with three clumps"),
    GalaxyOutline("dwarf_spheroid", "Dwarf spheroidal galaxy", "dwarf_spheroid", "a faint dashed round blob"),
)
"""tuple: The six outlines, in the order the Shape row lists the families."""

BY_ID = {outline.id: outline for outline in OUTLINES}
"""dict: `OUTLINES` by id."""

FAMILY_KEYS = ("spiral", "lenticular", "elliptical", "irregular", "dwarf_spheroid")
"""tuple: The five family keys of the design (`galaxy_shape.family`'s check constraint)."""

SPRITE_ID = "galaxy"
"""str: The id of the root `<svg>` of every outline file, what a `<use>` points at."""


def outline_id(family, bar=False):
    """
    The id of the outline that draws `family`; a spiral with `bar` set takes the barred one.

    Raises:
        ValueError: `family` is not one of `FAMILY_KEYS`.
    """
    if family not in FAMILY_KEYS:
        raise ValueError(f"no galaxy family {family!r}")
    return "barred_spiral" if family == "spiral" and bar else family


def outline_html(outline, file_url):
    """
    The `<svg>` of one outline as markup. `file_url` is a function from a path under `static/` ("galaxy/spiral.svg") to its URL
    (the site's `static_url`, which carries the release).

    Raises:
        ValueError: `outline` is not one of `OUTLINES`.
    """
    if outline not in BY_ID:
        raise ValueError(f"no galaxy outline named {outline!r}")
    return (f'<svg class="galaxy-outline" role="img" aria-label="{escape(BY_ID[outline].label)}" data-galaxy-outline="{outline}" '
            f'focusable="false"><use href="{escape(file_url(f"galaxy/{outline}.svg"))}#{SPRITE_ID}"></use></svg>')
