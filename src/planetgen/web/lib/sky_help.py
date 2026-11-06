# planetgen/web/lib/sky_help.py

"""
The help page about the sky behind a system (`/help/sky`, VIEW.38): what the System View's background is, what its brightness, color and
motion mean, how to read the status line, what hovering a star tells, and what is only a stand-in.

The numbers that live in Python (the limiting magnitude, the range a viewer may ask for, the cell the glow is shared over, the most sprites)
are read from the code that uses them, so the page cannot disagree with the sky it explains. The numbers that live in the page's scripts
(the lower limit on a phone, the star caps) are said in words. Imports no Flask.

`page()` returns the page as `{"title", "description", "intro": [text], "sections": [{"id", "heading", "paragraphs", "table", "facts",
"reference"}]}`, in the shape `templates/sky_help.html` draws.
"""

from planetgen.api import system_sky
from planetgen.galaxy import sky, sky_band
from planetgen.web.lib import help_registry
from planetgen.web.lib.skypack import MAX_LIMIT, MIN_LIMIT

TITLE = "The sky behind a system"
DESCRIPTION = "What the stars, the glow and the clouds behind a system in the System View are, and how to read them."


def _number(value):
    return f"{value:,.0f}" if float(value).is_integer() else f"{value:g}"


def _section(slug, heading, paragraphs, table=None, reference=None):
    return {"id": slug, "heading": heading, "paragraphs": paragraphs, "table": table, "facts": None, "reference": reference}


def page():
    """The help page about the sky of a system, as described at the top of this module."""
    limit = _number(sky.DEFAULT_LIMIT)
    intro = [
        "The dark background of a system in the System View is the sky that system would see. Each star is in the direction it really lies, at its real distance. "
        "It is shown where it was when the light now reaching the system left it.",
        "The sky is for looking at. It does not change what a planet or a ship can do.",
    ]
    sections = [
        _section("what-you-see", "What you see", [
            "Points of light are stars, each with the color of its surface temperature. Brighter stars are bigger and whiter, and a faint star "
            "shows less of its color. A soft band is the light of all the stars too faint to see one by one, and it is brightest toward the "
            "middle of the galaxy. Glowing patches are nebulae and the remnants of supernovae. A dark patch is a cloud that dims what is behind it.",
            "A planet or a moon in front of a star hides it.",
        ]),
        _section("how-faint", "How faint the stars go", [
            f"The sky shows the stars down to a limiting magnitude, {limit} unless you change it. Magnitude counts backward: a star of magnitude 1 "
            f"is brighter than one of magnitude 6, and a difference of 5 is a factor of 100. A phone shows a lower limit and fewer stars so the "
            f"page keeps up. Each step of one magnitude deeper brings in about two and a half to three times as many stars.",
            f"The page can ask for a limit from {_number(MIN_LIMIT)} to {_number(MAX_LIMIT)}. Stars fainter than the limit are not drawn one by one; "
            f"their light is in the soft band.",
        ]),
        _section("where-and-when", "Where a star is, and when", [
            "Light takes time to cross space. A star 100 light years away is seen as it was 100 years ago. The page puts each star where it was "
            "then, so a star that has moved fast shows up a little away from where it is today. For most stars the shift is too small to see; "
            "for a hypervelocity star it can be several minutes of arc.",
            "A star whose light has not had time to reach the system, because the star formed or exploded more recently than the light takes to "
            "travel, is not shown.",
            "Move the orbit clock to its fastest and the nearest stars drift across the sky, while far ones hold still. That is how it would look.",
        ]),
        _section("hover", "Hovering and clicking a star", [
            "Point at a star in the sky to see what it is. The card gives its name if it is a charted star, or its kind. It gives the spectral class, "
            "how far away the star is and how bright it looks. It also says how long ago the light you are seeing left it, as \"seen as it was N years ago\".",
            "Click a star to keep its card open. A charted star links to its page. Stars the survey has not charted yet are described but not linked.",
        ]),
        _section("status-line", "The line under the view", [
            "The line says how many stars the sky holds. If the stars near this system have not all been built yet, it says how far the sky is "
            "built and that it is filling in. The page asks again every half minute, and the new stars fade in. A system with no place in the galaxy "
            "has no sky, and the line says so. The Sky background setting turns the whole sky off, and the browser remembers the choice.",
        ]),
        _section("how-sure", "How sure each part is", [
            "The direction, distance, and brightness of a charted star are the survey's numbers. The stars not yet charted, the soft band, the "
            "dust that dims, and the look of the stars on the screen rest on chosen values where no measurement exists yet. Each chosen value "
            f"is listed with how it was chosen and what would replace it. The soft band is worked out once for each {_number(sky_band.CELL_PC)} "
            f"parsec cell around a system, and up to {_number(system_sky.MAX_SPRITES)} nebulae and remnants are drawn, the widest first.",
            "The galaxy's plane is drawn level with the plane of the orbits, because a system does not yet have a stored tilt to the galaxy.",
        ], reference=("stand-in", None)),
    ]
    return {"title": TITLE, "description": DESCRIPTION, "intro": intro, "sections": sections}


help_registry.register_category("views", "Maps and views", "How to read the maps and the views of the galaxy and its systems.", order=15)
help_registry.register_topic("sky_help_page", "The sky behind a system", "views",
                             "What the stars, the glow and the clouds behind a system are, and how to read them.", order=10)
