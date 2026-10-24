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


# --- Links into the habitability page (UX.181) -------------------------------------------------------------

DOMAIN_ANCHORS = {name: f"domain-{name}" for name in DOMAINS}
"""dict: The id of each domain's table on the habitability page."""

SYMBOL_ANCHORS = {**DOMAIN_ANCHORS, OVERALL: "overall"}
"""dict: Where each of the seven symbols is explained on the habitability page (the overall one is the score rule)."""

LIFE_ANCHORS = {name: f"life-{name}" for name in LIFE_NAMES}
"""dict: The id of each life score's row on the habitability page."""

CHIP_ANCHORS = {
    "habitable": "chip-habitable",
    "habitable moon": "chip-habitable-moon",
    "habitable class": "chip-habitable-class",
    "habitable-class moon": "chip-habitable-class-moon",
    "inhabited": "chip-inhabited",
}
"""dict: The id of each flag chip's entry on the habitability page, by its lower-case words."""

EQUIPMENT_ANCHORS = (
    # (the id of a row of the equipment table, words of an item name that put the item in that row); the first row
    # with a matching word wins, so a "hardened habitat" is the shelter row and a "hardened suit" the radiation row, not the suit it hardens.
    ("equipment-shelter", ("shelter", "habitat")),
    ("equipment-radiation-garment", ("shielded garment", "storm vest", "hardened")),
    ("equipment-powered-frame", ("powered frame",)),
    ("equipment-dust", ("dust",)),
    ("equipment-gas-filter", ("gas filter",)),
    ("equipment-co2-scrubber", ("scrubber",)),
    ("equipment-gas-mix-mask", ("gas-mix",)),
    ("equipment-self-contained-air", ("self-contained",)),
    ("equipment-mixed-gas", ("mixed-gas", "isolated air")),
    ("equipment-pressure-suit", ("pressure suit",)),
    ("equipment-armoured-suit", ("armoured",)),
    ("equipment-vehicle", ("vehicle",)),
    ("equipment-sealed-suit", ("sealed suit", "chemical-barrier", "heat-shielded")),
    ("equipment-thermal-suit", ("insulated", "liquid-cooled", "diving", "full-cover")),
    ("equipment-clothing", ("clothing",)),
)
"""tuple: `(anchor, words)` of the equipment table's rows, in the order `equipment_anchor` tries them."""

EQUIPMENT_FALLBACK = "equipment"
"""str: The equipment section's id, for an item no row of `EQUIPMENT_ANCHORS` names."""

PHI6_CHIP_ANCHOR = "overall"
"""str: Where the "PHI-6 Blue" chip of a body that needs no equipment is explained."""


def equipment_anchor(item):
    """The id on the habitability page of the row that explains an equipment item (`item` is a chip's item name,
    such as "heated sealed suit"), or the equipment section's own id for a name no row claims."""
    name = item.lower()
    for anchor, words in EQUIPMENT_ANCHORS:
        if any(word in name for word in words):
            return anchor
    return EQUIPMENT_FALLBACK


def help_link(inner_html, help_url, anchor):
    """`inner_html` as a link to `anchor` on the habitability page at `help_url`; `inner_html` itself where `help_url`
    is `None` (a page that has no URL for it, or is that page)."""
    if not help_url:
        return inner_html
    return f'<a class="help-link" href="{escape(help_url)}#{anchor}">{inner_html}</a>'


def tier_text(tier):
    """The colour name of a stored tier number (0 Blue to 3 Red). `ValueError` for anything else."""
    if isinstance(tier, bool) or not isinstance(tier, int) or not 0 <= tier < len(TIER_NAMES):
        raise ValueError(f"no habitability tier {tier!r}: expected 0 to {len(TIER_NAMES) - 1}")
    return TIER_NAMES[tier]


def symbol_html(name, tier, sprite_url, detail=None, help_url=None):
    """
    One symbol as markup: `name` is one of `SYMBOL_NAMES`, `tier` the stored rating number (0 Blue to 3 Red) and
    `sprite_url` the URL of `static/icons.svg`. `detail` (escaped here) follows the rating in the words, as
    "Pressure: Blue, 101 kPa". With `help_url` (the habitability page) the symbol links to where it is explained
    there. `ValueError` for an unknown name or tier.
    """
    if name not in SPRITE_IDS:
        raise ValueError(f"no habitability symbol {name!r}: expected one of {', '.join(SYMBOL_NAMES)}")
    words = f"{LABELS[name]}: {tier_text(tier)}"
    if detail:
        words += f", {detail}"
    words = escape(words)
    svg = (f'<svg class="icon hab-symbol {TIER_CLASSES[tier]}" role="img" aria-label="{words}" focusable="false">'
           f'<title>{words}</title><use href="{escape(sprite_url)}#{SPRITE_IDS[name]}"></use></svg>')
    return help_link(svg, help_url, SYMBOL_ANCHORS[name])


def overall_tier(body):
    """
    The overall rating of a body's row as a tier number (0 Blue to 3 Red): the colour of its stored PHI-6 score
    (`phi4`, held under the colour gates), not the worst domain. `None` for a body that was not scored (a gas giant,
    or one saved before v76).
    """
    score = body.get("phi4")
    return None if score is None else TIER_INDEX[habitability.tier(score)]


def symbol_row(body, sprite_url, with_overall=True, with_life=True, help_url=None):
    """
    The symbols of one body's row (`phi4` and `tier_<domain>` by name), the overall one first, then those of
    `DOMAINS`, then the cluster of its four life scores (`life_row`) when it has them; an empty string for a body
    that was not scored. With `help_url` every symbol links into the habitability page.
    """
    overall = overall_tier(body)
    if overall is None:
        return ""
    cells = [symbol_html(OVERALL, overall, sprite_url, help_url=help_url)] if with_overall else []
    cells += [symbol_html(domain, body[f"tier_{domain}"], sprite_url, help_url=help_url) for domain in DOMAINS]
    if with_life:
        cells.append(life_row(body, sprite_url, help_url))
    return f'<span class="hab-symbols">{"".join(cells)}</span>'


def life_text(name, score):
    """The words of a life score: "Flora: 3 of 9", or "Flora: none" for 0. `ValueError` for an unknown name or a
    score that is not a whole number from 0 to `LIFE_MAX`."""
    if name not in LIFE_SPRITE_IDS:
        raise ValueError(f"no life score {name!r}: expected one of {', '.join(LIFE_NAMES)}")
    if isinstance(score, bool) or not isinstance(score, int) or not 0 <= score <= LIFE_MAX:
        raise ValueError(f"no life score {score!r}: expected a whole number from 0 to {LIFE_MAX}")
    return f"{LIFE_LABELS[name]}: none" if score == 0 else f"{LIFE_LABELS[name]}: {score} of {LIFE_MAX}"


def life_html(name, score, sprite_url, help_url=None):
    """
    One life score as markup: its shape and the number 0 to 9 beside it, one `role="img"` whose label is the words
    ("Flora: 3 of 9"); a 0 is drawn muted (`hab-life-none`). With `help_url` it links into the habitability page.
    `ValueError` for an unknown name or score.
    """
    words = escape(life_text(name, score))
    muted = " hab-life-none" if score == 0 else ""
    shape = (f'<span class="hab-life{muted}" role="img" aria-label="{words}" title="{words}">'
             f'<svg class="icon hab-life-shape" aria-hidden="true" focusable="false">'
             f'<use href="{escape(sprite_url)}#{LIFE_SPRITE_IDS[name]}"></use></svg>'
             f'<span class="hab-life-score" aria-hidden="true">{score}</span></span>')
    return help_link(shape, help_url, LIFE_ANCHORS[name])


def life_row(body, sprite_url, help_url=None):
    """
    The four life scores of a body's row (`life_<name>`, whole numbers 0 to 9) as a cluster, in the order of
    `LIFE_NAMES`; an empty string for a body with no scores.
    """
    scores = [body.get(f"life_{name}") for name in LIFE_NAMES]
    if any(score is None for score in scores):
        return ""
    cells = "".join(life_html(name, score, sprite_url, help_url) for name, score in zip(LIFE_NAMES, scores))
    return f'<span class="hab-life-cluster">{cells}</span>'
