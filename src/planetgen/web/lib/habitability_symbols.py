# planetgen/web/lib/habitability_symbols.py

"""
The PHI-6 habitability symbols (UX.101): the one place that names them and builds their markup, so every page draws
the same shape, colour and words for a rating.

The shapes are in the site's icon sprite (`static/icons.svg`, ids `hab-<domain>` and `hab-overall`); this module
puts one on a page as inline `<svg>` that carries its words (`role="img"`, an `aria-label` and a `<title>` tooltip,
"Pressure: Blue"), so the rating never depends on colour alone. The colour comes from the `tier-blue` to `tier-red`
classes of `static/style.css`, which draw in both themes. A page reaches this through the Jinja global
`habitability_symbol(domain, tier)` or, in Python that has the sprite's URL, `symbol_html`; `symbol_row` is the
seven symbols of one body.

Imports no Flask, so the describe renderer and the list builders can use it.
"""

from html import escape

from planetgen.physics import habitability
from planetgen.physics.habitability_world import DOMAINS, TIER_INDEX, TIER_NAMES

OVERALL = "overall"
"""str: The name of the overall symbol, beside the six domains of `DOMAINS`."""

SYMBOL_NAMES = (*DOMAINS, OVERALL)
"""tuple: The seven symbols: pressure, temperature, air, water, radiation, suitability and overall."""

LABELS = {
    "pressure": "Pressure",
    "temperature": "Temperature",
    "air": "Air chemistry",
    "water": "Water chemistry",
    "radiation": "Radiation",
    "suitability": "Suitability",
    OVERALL: "Overall habitability",
}
"""dict: Each symbol's words, as the tooltip and the screen reader say them."""

SPRITE_IDS = {name: f"hab-{name}" for name in SYMBOL_NAMES}
"""dict: The `<symbol>` id of each symbol in `static/icons.svg`."""

SHAPES = {
    "pressure": "a square between two triangles pointing at it",
    "temperature": "a triangle split by a vertical line",
    "air": "a hexagon holding three horizontal lines",
    "water": "a downward triangle holding a circle",
    "radiation": "a dot inside three arcs",
    "suitability": "a circle holding a diamond and a dot",
    OVERALL: "a solid circle in a triangle in a ring",
}
"""dict: What each shape looks like, for the legend on the habitability levels page."""

LIFE_NAMES = ("microbial", "flora", "fauna", "sentient")
"""tuple: The four life scores (Boss, 6:14 PM EDT 2026-10-23), a body's `life_<name>` columns."""

LIFE_MAX = 9
"""int: The top of a life score: 0 is none and 9 abundant."""

LIFE_LABELS = {"microbial": "Microbial life", "flora": "Flora", "fauna": "Fauna", "sentient": "Sentient life"}
"""dict: Each life score's words, as the tooltip and the screen reader say them."""

LIFE_SPRITE_IDS = {name: f"life-{name}" for name in LIFE_NAMES}
"""dict: The `<symbol>` id of each life symbol in `static/icons.svg`."""

LIFE_SHAPES = {
    "microbial": "a cluster of cells",
    "flora": "a stem with two leaves",
    "fauna": "a paw print",
    "sentient": "a head and shoulders",
}
"""dict: What each life symbol looks like, for the legend on the habitability levels page."""

TIER_CLASSES = tuple(f"tier-{name.lower()}" for name in TIER_NAMES)
"""tuple: The CSS class of each rating colour by tier number (0 Blue to 3 Red)."""


def tier_text(tier):
    """The colour name of a stored tier number (0 Blue to 3 Red). `ValueError` for anything else."""
    if isinstance(tier, bool) or not isinstance(tier, int) or not 0 <= tier < len(TIER_NAMES):
        raise ValueError(f"no habitability tier {tier!r}: expected 0 to {len(TIER_NAMES) - 1}")
    return TIER_NAMES[tier]


def symbol_html(name, tier, sprite_url, detail=None):
    """
    One symbol as markup: `name` is one of `SYMBOL_NAMES`, `tier` the stored rating number (0 Blue to 3 Red) and
    `sprite_url` the URL of `static/icons.svg`. `detail` (escaped here) follows the rating in the words, as
    "Pressure: Blue, 101 kPa". `ValueError` for an unknown name or tier.
    """
    if name not in SPRITE_IDS:
        raise ValueError(f"no habitability symbol {name!r}: expected one of {', '.join(SYMBOL_NAMES)}")
    words = f"{LABELS[name]}: {tier_text(tier)}"
    if detail:
        words += f", {detail}"
    words = escape(words)
    return (f'<svg class="icon hab-symbol {TIER_CLASSES[tier]}" role="img" aria-label="{words}" focusable="false">'
            f'<title>{words}</title><use href="{escape(sprite_url)}#{SPRITE_IDS[name]}"></use></svg>')


def overall_tier(body):
    """
    The overall rating of a body's row as a tier number (0 Blue to 3 Red): the colour of its stored PHI-6 score
    (`phi4`, held under the colour gates), not the worst domain. `None` for a body that was not scored (a gas giant,
    or one saved before v76).
    """
    score = body.get("phi4")
    return None if score is None else TIER_INDEX[habitability.tier(score)]


def symbol_row(body, sprite_url, with_overall=True, with_life=True):
    """
    The symbols of one body's row (`phi4` and `tier_<domain>` by name), the overall one first, then those of
    `DOMAINS`, then the cluster of its four life scores (`life_row`) when it has them; an empty string for a body
    that was not scored.
    """
    overall = overall_tier(body)
    if overall is None:
        return ""
    cells = [symbol_html(OVERALL, overall, sprite_url)] if with_overall else []
    cells += [symbol_html(domain, body[f"tier_{domain}"], sprite_url) for domain in DOMAINS]
    if with_life:
        cells.append(life_row(body, sprite_url))
    return f'<span class="hab-symbols">{"".join(cells)}</span>'


def life_text(name, score):
    """The words of a life score: "Flora: 3 of 9", or "Flora: none" for 0. `ValueError` for an unknown name or a
    score that is not a whole number from 0 to `LIFE_MAX`."""
    if name not in LIFE_SPRITE_IDS:
        raise ValueError(f"no life score {name!r}: expected one of {', '.join(LIFE_NAMES)}")
    if isinstance(score, bool) or not isinstance(score, int) or not 0 <= score <= LIFE_MAX:
        raise ValueError(f"no life score {score!r}: expected a whole number from 0 to {LIFE_MAX}")
    return f"{LIFE_LABELS[name]}: none" if score == 0 else f"{LIFE_LABELS[name]}: {score} of {LIFE_MAX}"


def life_html(name, score, sprite_url):
    """
    One life score as markup: its shape and the number 0 to 9 beside it, one `role="img"` whose label is the words
    ("Flora: 3 of 9"); a 0 is drawn muted (`hab-life-none`). `ValueError` for an unknown name or score.
    """
    words = escape(life_text(name, score))
    muted = " hab-life-none" if score == 0 else ""
    return (f'<span class="hab-life{muted}" role="img" aria-label="{words}" title="{words}">'
            f'<svg class="icon hab-life-shape" aria-hidden="true" focusable="false">'
            f'<use href="{escape(sprite_url)}#{LIFE_SPRITE_IDS[name]}"></use></svg>'
            f'<span class="hab-life-score" aria-hidden="true">{score}</span></span>')


def life_row(body, sprite_url):
    """
    The four life scores of a body's row (`life_<name>`, whole numbers 0 to 9) as a cluster, in the order of
    `LIFE_NAMES`; an empty string for a body with no scores.
    """
    scores = [body.get(f"life_{name}") for name in LIFE_NAMES]
    if any(score is None for score in scores):
        return ""
    cells = "".join(life_html(name, score, sprite_url) for name, score in zip(LIFE_NAMES, scores))
    return f'<span class="hab-life-cluster">{cells}</span>'
