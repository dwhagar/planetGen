# planetgen/describe/html.py

"""
A Document as HTML for the web page (docs/design/architecture/system-description.md, section 4, UX.119): the fourth renderer beside
text, markdown and wikitext.

Every number is a `<span class="qty" data-q="length" data-si="...">` (the markup `web/lib/fmt.quantity_html` writes), so
`static/units.js` shows it in the reader's unit system; a quantity with no ladder (a radius, a dose, a count) is plain
text. All text is escaped, whatever the Document holds. The tables carry the classes `doc-facts`, `doc-ratings` and
`doc-composition`, and a rating's tier its colour class (`tier-blue` to `tier-red`), which `static/style.css` draws in
both themes. The six symbols of UX.101 are not drawn yet: a rating names its tier in words.
"""

import html
import math

from planetgen.describe.quantity import LADDER_BY_KIND, _powers, quantity_text
from planetgen.describe.render import share_text
from planetgen.util import chem
from planetgen.util.format import format_quantity

FIRST_HEADING = 2
"""int: The heading level of a top section; a child is one deeper, to `<h6>`."""


def _quantity(quantity):
    kind, value = quantity["kind"], quantity["si"]
    if kind in LADDER_BY_KIND and value is not None and math.isfinite(value):
        companions = quantity.get("companions", True)
        # Plain text inside the span, as `web/lib/fmt.quantity_html` writes it: `static/units.js` rewrites the span's
        # text, with Unicode superscripts, in the reader's units.
        text = html.escape(format_quantity(LADDER_BY_KIND[kind], value, companions=companions))
        return (f'<span class="qty" data-q="{html.escape(kind)}" data-si="{value!r}"'
                f'{"" if companions else " data-nc"}>{text}</span>')
    return _powers(html.escape(quantity_text(quantity)), "markdown")


def inline_html(items, link=None):
    """An inline list as HTML. `link` is a function from a `ref` string (`"planet:31"`) to a URL or `None`; with no
    function, or no URL, a reference is its name."""
    out = []
    for item in items:
        if isinstance(item, str):
            out.append(html.escape(item))
        elif "qty" in item:
            out.append(_quantity(item["qty"]))
        elif "strong" in item:
            out.append(f"<strong>{html.escape(item['strong'])}</strong>")
        elif "swatch" in item:
            colour = item["swatch"]
            word = html.escape(colour["word"])
            out.append(f'<span class="swatch" data-colour="{html.escape(colour["hex"])}"></span> {word}')
        elif "chem" in item:
            out.append(chem.html(item["chem"]))
        else:
            url = link(item["ref"]) if link else None
            name = html.escape(item["name"])
            out.append(f'<a href="{html.escape(url, quote=True)}">{name}</a>' if url else name)
    return "".join(out)


def _table(css, header, rows, caption=None):
    head = "".join(f'<th scope="col">{html.escape(cell)}</th>' for cell in header)
    body = []
    for first, *rest in rows:
        cells = "".join(f"<td>{cell}</td>" for cell in rest)
        body.append(f'<tr><th scope="row">{first}</th>{cells}</tr>')
    title = f"<caption>{html.escape(caption)}</caption>" if caption else ""
    return (f'<div class="table-scroll" tabindex="0"><table class="{css}">{title}'
            f'<thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>')


def _facts(block, link):
    rows = [(html.escape(row["label"]), inline_html(row["value"], link)) for row in block["rows"]]
    return _table("doc-facts", ("Property", "Value"), rows, block.get("caption"))


def _composition(block):
    parts = []
    if block["summary"]:
        parts.append(f'<p class="doc-prose doc-composition-summary"><strong>{html.escape(block["label"])}</strong>: '
                     f'{html.escape(block["summary"])}.</p>')
    if block["parts"]:
        rows = [(html.escape(part["name"].capitalize()), html.escape(share_text(part))) for part in block["parts"]]
        parts.append(_table("doc-composition", ("Component", "Share"), rows, block["label"]))
    return "".join(parts)


def _ratings(block, link):
    overall = block["overall"]
    headline = (f'<p class="doc-prose doc-rating-headline">{html.escape(overall["name"])} habitability is '
                f'{overall["score"]:.2f} ({_tier(overall["tier"])}), and a visitor needs '
                f'{inline_html(block["gear"], link)}.</p>')
    rows = [(html.escape(row["domain"]), _tier(row["tier"]), f'{row["score"]:.2f}', inline_html(row["text"], link))
            for row in block["rows"]]
    return headline + _table("doc-ratings", ("Category", "Rating", "Score", "Notes"), rows)


def _tier(tier):
    return f'<span class="tier tier-{html.escape(tier.lower())}">{html.escape(tier)}</span>'


def _block(block, link):
    if block["type"] == "prose":
        text = " ".join(inline_html(sentence, link) for sentence in block["text"])
        return f'<p class="doc-prose doc-{html.escape(block["role"])}">{text}</p>'
    if block["type"] == "facts":
        return _facts(block, link)
    if block["type"] == "composition":
        return _composition(block)
    return _ratings(block, link)


def blocks_html(blocks, link=None):
    """The blocks of one section as HTML, with no heading and no child sections: what a row of the system page
    shows under the body's own name. `link` is as in `inline_html`."""
    return "".join(_block(block, link) for block in blocks)


def _section(section, depth, link):
    level = min(depth, 6)
    ref = f' data-ref="{html.escape(section["ref"], quote=True)}"' if section["ref"] else ""
    parts = [f'<section class="doc-section doc-{html.escape(section["kind"])}"{ref}>',
             f'<h{level}>{html.escape(section["title"])}</h{level}>']
    parts += [_block(block, link) for block in section["blocks"]]
    parts += [_section(child, depth + 1, link) for child in section["children"]]
    parts.append("</section>")
    return "".join(parts)


def render(document, link=None, first_heading=FIRST_HEADING):
    """The Document as HTML. `link` turns a `ref` into a URL (see `inline_html`); `first_heading` is the level of the
    top sections' headings."""
    return "".join(_section(top, first_heading, link) for top in document["sections"])
