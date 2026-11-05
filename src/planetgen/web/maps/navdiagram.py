# planetgen/web/maps/navdiagram.py

"""
The course diagram of the NAV page (NAV.67, NAV.36)
===================================================

The ribbon of `navribbon.layout_ribbon` drawn as SVG, in a panel below the NAV page's map: the route as one line of stops, left to
right, a glyph for each kind of stop, the length of every hop written on it with its course when there is room, and a hop through
unknown space red, dashed and glowing. It is a schematic, not a map: the spacing follows the logarithm of the hop's length
(`navribbon.hop_gap`), so a stop's place on the page says nothing about where it lies in space. The route list below the panel
stays the plain-text form of the same route; each stop and hop here also carries a text alternative (its `<title>`).

The panel holds three drawings of one route, laid out for a narrow, a middle and a wide page (`WIDTHS`); a media query in
`style.css` shows the one that fits, so the text stays near its drawn size (an SVG scaled down to a phone shrinks its text with it)
and a long route wraps on a phone where it runs on in a row on a desktop. The glyphs and colours are classes that `style.css`
colours with the site's theme tokens, in both themes; the glow of an unknown-space hop pulses only where the visitor has not asked
for reduced motion.

The stubs round the stops (NAV.68) and the link to the Galaxy Map (NAV.69) are not drawn here yet; the panel says what a hop is and
that the spacing is not to scale.
"""

import html
import math

from planetgen.galaxy.navigation import format_course
from planetgen.web.maps import navribbon

WIDTHS = (("narrow", 260), ("middle", 520), ("wide", 880))
"""tuple: The drawings of a route, `(size name, layout width in pixels)`: shown below 560 px, from 560 to 959 px and from 960 px
(`.nav-diagram-narrow` and the others in `style.css`)."""

CHAR_WIDTH = 6.4
"""float: What a character of a diagram label takes of its width at `FONT_SIZE`, pixels (to see whether a label fits its gap)."""

FONT_SIZE = 12
"""int: The size of a diagram's text, in the drawing's own pixels."""

NAME_LIMIT = 15
"""int: The longest stop name written under a stop, characters; a longer one is cut with an ellipsis and whole in the title."""

NAME_STAGGER = 15
"""int: How far lower the name of every second stop is written, pixels, so two names beside each other can each take two gaps."""

NODE_RADIUS = 7

KIND_LABELS = {
    "system": "star system", "phenomenon": "phenomenon", "bright-star": "bright star in unknown space",
    "black-hole": "black hole in unknown space", "neutron-star": "neutron star in unknown space",
    "quasar": "quasar in unknown space", "hypervelocity-star": "hypervelocity star in unknown space",
    "stone": "object in unknown space",
}
"""dict: The text for each kind of stop (the legend and the stop's text alternative)."""

_ROLE_TEXT = {"origin": "Start", "destination": "Destination", "stop": "Stop"}


def _esc(text):
    return html.escape(str(text), quote=True)


def _number(value):
    return f"{float(value):.1f}".rstrip("0").rstrip(".")


def glyph(kind, x, y, radius=NODE_RADIUS):
    """
    The SVG shapes of one stop's glyph centred on (`x`, `y`): a disc for a system, a hexagon for a phenomenon, and for each kind of
    unknown-space stop its own outline (a four-point star, a ringed dark disc, a diamond, a double ring, a chevron, a square), so
    the kinds read apart by shape and not by colour alone.
    """
    r = radius
    cx, cy = _number(x), _number(y)
    if kind == "system":
        return f'<circle class="nav-glyph-fill" cx="{cx}" cy="{cy}" r="{_number(r)}"/>'
    if kind == "phenomenon":
        points = " ".join(f"{_number(x + r * 1.15 * math.cos(math.pi / 3 * k))},{_number(y + r * 1.15 * math.sin(math.pi / 3 * k))}"
                          for k in range(6))
        return f'<polygon class="nav-glyph-fill" points="{points}"/>'
    if kind == "bright-star":
        points = []
        for k in range(8):
            reach = r * 1.5 if k % 2 == 0 else r * 0.55
            angle = math.pi / 4 * k - math.pi / 2
            points.append(f"{_number(x + reach * math.cos(angle))},{_number(y + reach * math.sin(angle))}")
        return f'<polygon class="nav-glyph-stone" points="{" ".join(points)}"/>'
    if kind == "black-hole":
        return (f'<circle class="nav-glyph-hole" cx="{cx}" cy="{cy}" r="{_number(r)}"/>'
                f'<circle class="nav-glyph-ring" cx="{cx}" cy="{cy}" r="{_number(r * 1.45)}"/>')
    if kind == "neutron-star":
        points = f"{cx},{_number(y - r * 1.3)} {_number(x + r * 1.3)},{cy} {cx},{_number(y + r * 1.3)} {_number(x - r * 1.3)},{cy}"
        return f'<polygon class="nav-glyph-stone" points="{points}"/>'
    if kind == "quasar":
        return (f'<circle class="nav-glyph-stone" cx="{cx}" cy="{cy}" r="{_number(r * 0.6)}"/>'
                f'<circle class="nav-glyph-ring" cx="{cx}" cy="{cy}" r="{_number(r * 1.35)}"/>')
    if kind == "hypervelocity-star":
        points = f"{_number(x - r)},{_number(y - r)} {_number(x + r * 1.2)},{cy} {_number(x - r)},{_number(y + r)} {_number(x - r * 0.3)},{cy}"
        return f'<polygon class="nav-glyph-stone" points="{points}"/>'
    return f'<rect class="nav-glyph-stone" x="{_number(x - r)}" y="{_number(y - r)}" width="{_number(2 * r)}" height="{_number(2 * r)}"/>'


def _short_name(name):
    return name if len(name) <= NAME_LIMIT else name[:NAME_LIMIT - 1].rstrip() + "…"


def _hop_words(hop, source):
    """The text alternative of a hop: its length, its course when it has one, and that it crosses unknown space."""
    parts = [f"{hop['label']}"]
    if "bearing_deg" in source:
        parts.append(f"course {format_course(source['bearing_deg'], source['mark_deg'])}")
    if hop["unknown_space"]:
        parts.append("through unknown space (uncharted sectors)")
    return ", ".join(parts)


def _hop_gap_pixels(segments):
    return max(abs(s["x2"] - s["x1"]) for s in segments)


def render_drawing(route, names, urls, width, size, ids=None):
    """
    One drawing of the route at layout `width` (`<svg class="nav-diagram-svg nav-diagram-<size>">`).

    Args:
        route (dict): As `navribbon.layout_ribbon` takes it, plus `bearing_deg` and `mark_deg` on a hop for its course label.
        names (dict): `{node: name}` of every stop.
        urls (dict): `{node: url or None}` of every stop (None: no page to link to, such as a stop in unknown space).
        width (int): The layout width, pixels.
        size (str): The size name (a CSS class).
        ids (str, optional): A prefix for the element ids inside the drawing (two drawings of one page must not share ids).
    """
    ribbon = navribbon.layout_ribbon(route, width)
    prefix = ids or f"navdiagram-{size}"
    sources = {index: hop for index, hop in enumerate(route.get("hops") or [])}
    height = ribbon["height"] or navribbon.ROW_HEIGHT
    out = [f'<svg class="nav-diagram-svg nav-diagram-{_esc(size)}" viewBox="0 0 {ribbon["width"]} {_number(height)}" '
           f'width="{ribbon["width"]}" height="{_number(height)}" role="group" aria-labelledby="{prefix}-title" '
           f'font-size="{FONT_SIZE}" focusable="false">',
           f'<title id="{prefix}-title">Course diagram: {len(ribbon["stops"])} stops, from {_esc(names.get(ribbon["stops"][0]["id"], ""))}'
           f' to {_esc(names.get(ribbon["stops"][-1]["id"], ""))}</title>',
           f'<defs><filter id="{prefix}-glow" filterUnits="userSpaceOnUse" x="0" y="0" width="{ribbon["width"]}" height="{_number(height)}">'
           '<feGaussianBlur stdDeviation="2.6" result="blur"/>'
           '<feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge></filter></defs>']

    # The hops under the glyphs.
    for index, hop in enumerate(ribbon["hops"]):
        unknown = hop["unknown_space"]
        classes = "nav-diagram-hop nav-diagram-hop-unknown" if unknown else "nav-diagram-hop"
        filter_attr = f' filter="url(#{prefix}-glow)"' if unknown else ""
        out.append(f'<g class="{classes}" data-from="{_esc(hop["from"])}" data-to="{_esc(hop["to"])}"{filter_attr}>')
        out.append(f'<title>{_esc(_hop_words(hop, sources.get(index, {})))}</title>')
        for segment in hop["segments"]:
            out.append(f'<line class="nav-diagram-line" x1="{_number(segment["x1"])}" y1="{_number(segment["y1"])}" '
                       f'x2="{_number(segment["x2"])}" y2="{_number(segment["y2"])}"/>')
        out.append('</g>')
        gap = _hop_gap_pixels(hop["segments"])
        text = hop["label"] if len(hop["label"]) * CHAR_WIDTH + 6 <= gap else hop["short_label"]
        source = sources.get(index, {})
        label_class = "nav-diagram-label nav-diagram-label-unknown" if unknown else "nav-diagram-label"
        out.append(f'<text class="{label_class}" x="{_number(hop["label_x"])}" y="{_number(hop["label_y"])}" text-anchor="middle" '
                   f'aria-hidden="true">{_esc(text)}</text>')
        if "bearing_deg" in source:
            course = format_course(source["bearing_deg"], source["mark_deg"])
            if len(course) * CHAR_WIDTH + 6 <= gap:
                out.append(f'<text class="nav-diagram-course" x="{_number(hop["label_x"])}" y="{_number(hop["label_y"] - FONT_SIZE - 5)}" '
                           f'text-anchor="middle" aria-hidden="true">{_esc(course)}</text>')

    # The stops over them.
    for stop in ribbon["stops"]:
        name = names.get(stop["id"]) or str(stop["id"])
        words = f"{_ROLE_TEXT[stop['role']]}: {name}, {KIND_LABELS.get(stop['kind'], KIND_LABELS['stone'])}"
        url = urls.get(stop["id"])
        x, y = stop["x"], stop["y"]
        open_tag = (f'<a class="nav-diagram-stop" href="{_esc(url)}" data-stop="{_esc(stop["id"])}">' if url
                    else f'<g class="nav-diagram-stop" data-stop="{_esc(stop["id"])}">')
        out.append(open_tag)
        out.append(f'<title>{_esc(words)}</title>')
        if stop["role"] != "stop":
            out.append(f'<circle class="nav-diagram-end nav-diagram-{stop["role"]}" cx="{_number(x)}" cy="{_number(y)}" r="{NODE_RADIUS + 5}"/>')
            if stop["role"] == "destination":
                out.append(f'<circle class="nav-diagram-end nav-diagram-destination-inner" cx="{_number(x)}" cy="{_number(y)}" r="{NODE_RADIUS + 2}"/>')
        out.append(glyph(stop["kind"], x, y))
        drop = NODE_RADIUS + 15 + (NAME_STAGGER if stop["index"] % 2 else 0)
        shown = _short_name(name)
        # A name that would run past the drawing's edge starts or ends at its stop instead of being centred on it.
        half = len(shown) * CHAR_WIDTH / 2
        anchor, text_x = ("start", x - NODE_RADIUS) if x - half < 2 else ("end", x + NODE_RADIUS) if x + half > ribbon["width"] - 2 else ("middle", x)
        out.append(f'<text class="nav-diagram-name" x="{_number(text_x)}" y="{_number(y + drop)}" text-anchor="{anchor}">{_esc(shown)}</text>')
        out.append('</a>' if url else '</g>')
    out.append('</svg>')
    return "".join(out)


def render_legend(kinds, unknown=False):
    """The legend of the glyphs the route uses: a small drawing and its words for each kind, then the unknown-space hop when `unknown`."""
    items = []
    for kind in KIND_LABELS:
        if kind not in kinds:
            continue
        items.append(f'<li><svg class="nav-legend-glyph" viewBox="0 0 22 22" width="22" height="22" aria-hidden="true" focusable="false">'
                     f'{glyph(kind, 11, 11, 6)}</svg> {_esc(KIND_LABELS[kind].capitalize())}</li>')
    if unknown:
        items.append('<li><svg class="nav-legend-glyph nav-diagram-hop-unknown" viewBox="0 0 34 22" width="34" height="22" aria-hidden="true" '
                     'focusable="false"><line class="nav-diagram-line" x1="1" y1="11" x2="33" y2="11"/></svg> Hop through unknown space</li>')
    return f'<ul class="nav-diagram-legend" role="list">{"".join(items)}</ul>' if items else ""


def render_course_diagram(route, names, urls):
    """
    The NAV page's course diagram panel for `route`, or an empty string when the route has fewer than two stops.

    Args:
        route (dict): The route `GET /api/nav` sends (`path`, `hops`, `stop_places`).
        names (dict): `{node: name}` of every stop.
        urls (dict): `{node: url or None}` of every stop.

    Returns:
        str: A complete `<section class="panel">` block.
    """
    path = (route or {}).get("path") or []
    if len(path) < 2:
        return ""
    places = {place["node"]: place for place in route.get("stop_places") or []}
    kinds = {navribbon._kind(node, places.get(node, {})) for node in path}
    unknown = sum(1 for hop in route.get("hops") or [] if hop.get("unknown_space"))
    drawings = "".join(f'<div class="nav-diagram-{size}-box">{render_drawing(route, names, urls, width, size)}</div>'
                       for size, width in WIDTHS)
    unknown_note = (f' {unknown} hop{"" if unknown == 1 else "s"} cross{"es" if unknown == 1 else ""} unknown space (uncharted '
                    f'sectors); {"it is" if unknown == 1 else "they are"} drawn in red, dashed and glowing, and the route list names '
                    f'{"it" if unknown == 1 else "them"} in words.' if unknown else "")
    return (
        '<section class="panel nav-diagram" aria-labelledby="course-diagram-heading">'
        '<h2 id="course-diagram-heading">Course Diagram</h2>'
        '<p class="hint">The route as a line of stops, each hop with its length'
        ' (and its course where there is room). The spacing follows the logarithm of a hop’s length, so the diagram is a'
        f' schematic, not a map to scale.{_esc(unknown_note)}</p>'
        f'<div class="nav-diagram-scroll">{drawings}</div>'
        f'{render_legend(kinds, bool(unknown))}'
        '</section>'
    )
