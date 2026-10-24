# planetgen/describe/render.py

"""
The walk the plain-string renderers share (docs/design/architecture/system-description.md, section 4): a
Document's sections in order, each as a heading, its blocks and its child sections. A style
(`text`, `markdown`, `wikitext`) says how a heading, a table, a composition and a rating look.
"""

from planetgen.describe.quantity import inline_text

MAX_DEPTH = 6


def share_text(part):
    """A composition part's share (`"share"`, 0 to 1) or `"level"` word as text."""
    if "level" in part:
        return part["level"]
    share = part["share"]
    return f"{share * 100:.3g}%" if share >= 0.001 else "< 0.1%"


def rating_lines(block, fmt):
    """`(headline, rows)` of a ratings block: the overall line, then `(domain, rating, score, notes)` per row."""
    overall = block["overall"]
    headline = (f"{overall['name']} habitability is {overall['score']:.2f} ({overall['tier']}), "
                f"and a visitor needs {inline_text(block['gear'], fmt)}.")
    rows = [(row["domain"], row["tier"], f"{row['score']:.2f}", inline_text(row["text"], fmt)) for row in block["rows"]]
    return headline, rows


def render_document(document, style, fmt):
    """The Document as one string: `style` is an object with `heading`, `facts`, `composition`, `ratings`
    and `join`; prose is the same everywhere."""
    chunks = []

    def walk(section, depth):
        chunks.append(style.heading(min(depth, MAX_DEPTH), section))
        for block in section["blocks"]:
            if block["type"] == "prose":
                chunks.append(" ".join(inline_text(sentence, fmt) for sentence in block["text"]))
            elif block["type"] == "facts":
                chunks.append(style.facts(section, block, [(row["key"], row["label"], inline_text(row["value"], fmt))
                                                           for row in block["rows"]]))
            elif block["type"] == "composition":
                chunks.append(style.composition(block))
            elif block["type"] == "ratings":
                chunks.append(style.ratings(block))
        for child in section["children"]:
            walk(child, depth + 1)

    for top in document["sections"]:
        walk(top, 1)
    return style.join([chunk for chunk in chunks if chunk])
