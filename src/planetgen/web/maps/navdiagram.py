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

NAV.68: when the route carries its `graph` (NAV.65) the neighbours of each stop are drawn faintly under the route: a short stub from a
stop to each neighbour off the route (as many as the graph lists, which it caps per stop), "and N more" where the cap cut some, and an
arc between two stops of one row that are neighbours but no hop. A link is a neighbour relation (which systems lie near each other),
not a travel lane, and the panel says so; the route's own hops are always drawn, over the neighbours. The link to the Galaxy Map
(NAV.69) is not drawn here yet.
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
    The SVG shapes of one stop's glyph centered on (`x`, `y`): a disk for a system, a hexagon for a phenomenon, and for each kind of
    unknown-space stop its own outline (a four-point star, a ringed dark disk, a diamond, a double ring, a chevron, a square), so
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


def _stub_words(stub, names):
    """The text alternative of a stub: a neighbour of the stop, how far off, and that it is not a lane."""
    return f"A system {stub['label']} from {names.get(stub['stop'], stub['stop'])}: a neighbour, not a travel lane"


def _link_words(link, names):
    return (f"{names.get(link['from'], link['from'])} and {names.get(link['to'], link['to'])} are neighbours {link['label']} apart: "
            "a neighbour relation, not a travel lane")


def _more_words(more, names):
    return f"{more['count']} more neighbours of {names.get(more['stop'], more['stop'])} are not drawn"


def _anchored(x, text, width):
    """`(text-anchor, x)` for `text` at a stop at `x`: centered on it, or starting or ending at it when it would run past the drawing's edge."""
    half = len(text) * CHAR_WIDTH / 2
    if x - half < 2:
        return "start", x - NODE_RADIUS
    if x + half > width - 2:
        return "end", x + NODE_RADIUS
    return "middle", x


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

    # The neighbours, faint, under the route: the arcs between stops, the stubs round them and their "and N more".
    for link in ribbon["links"]:
        control = link["y"] + 2 * (link["peak_y"] - link["y"])
        out.append(f'<g class="nav-diagram-neighbour nav-diagram-link"><title>{_esc(_link_words(link, names))}</title>'
                   f'<path d="M{_number(link["x1"])} {_number(link["y"])} Q{_number((link["x1"] + link["x2"]) / 2)} {_number(control)} '
                   f'{_number(link["x2"])} {_number(link["y"])}"/></g>')
        out.append(f'<text class="nav-diagram-neighbour-label" x="{_number((link["x1"] + link["x2"]) / 2)}" '
                   f'y="{_number(link["peak_y"] + (-3 if link["side"] == "above" else 11))}" text-anchor="middle" aria-hidden="true">'
                   f'{_esc(link["short_label"])}</text>')
    for stub in ribbon["stubs"]:
        out.append(f'<g class="nav-diagram-neighbour nav-diagram-stub"><title>{_esc(_stub_words(stub, names))}</title>'
                   f'<line x1="{_number(stub["x1"])}" y1="{_number(stub["y1"])}" x2="{_number(stub["x2"])}" y2="{_number(stub["y2"])}"/>'
                   f'<circle cx="{_number(stub["x2"])}" cy="{_number(stub["y2"])}" r="2.2"/></g>')
    stop_y = {stop["id"]: stop["y"] for stop in ribbon["stops"]}
    more_spans = []       # where the "and N more" texts lie across the page: a hop's course text keeps clear of them
    for more in ribbon["more"]:
        # Above the stop, over the tips of its stubs (the layout's own place for it is on the side with fewer stubs, where the name is).
        anchor, text_x = _anchored(more["x"], more["text"], ribbon["width"])
        text_width = len(more["text"]) * CHAR_WIDTH
        more_spans.append((text_x - text_width / 2, text_x + text_width / 2) if anchor == "middle"
                          else (text_x, text_x + text_width) if anchor == "start" else (text_x - text_width, text_x))
        out.append(f'<text class="nav-diagram-more" x="{_number(text_x)}" y="{_number(stop_y[more["stop"]] - navribbon.STUB_LENGTH - 4)}" '
                   f'text-anchor="{anchor}"><title>{_esc(_more_words(more, names))}</title>{_esc(more["text"])}</text>')

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
            reach = (hop["label_x"] - len(course) * CHAR_WIDTH / 2, hop["label_x"] + len(course) * CHAR_WIDTH / 2)
            clear = all(reach[1] < left - 4 or reach[0] > right + 4 for left, right in more_spans)
            if len(course) * CHAR_WIDTH + 6 <= gap and clear:
                out.append(f'<text class="nav-diagram-course" x="{_number(hop["label_x"])}" y="{_number(hop["label_y"] - FONT_SIZE - 5)}" '
                           f'text-anchor="middle" aria-hidden="true">{_esc(course)}</text>')

    # The stops over them.
    for stop in ribbon["stops"]:
        name = names.get(stop["id"]) or str(stop["id"])
        words = f"{_ROLE_TEXT[stop['role']]}: {name}, {KIND_LABELS.get(stop['kind'], KIND_LABELS['stone'])}"
        url = urls.get(stop["id"])
        x, y = stop["x"], stop["y"]
        open_tag = (f'<a class="nav-diagram-stop" href="{_esc(url)}" data-stop="{_esc(stop["id"])}" data-index="{stop["index"]}">' if url
                    else f'<g class="nav-diagram-stop" data-stop="{_esc(stop["id"])}" data-index="{stop["index"]}">')
        out.append(open_tag)
        out.append(f'<title>{_esc(words)}</title>')
        # Shown when the stop is the selected one of the course on the Galaxy Map (NAV.69, coursediagram.js).
        out.append(f'<circle class="nav-diagram-select" cx="{_number(x)}" cy="{_number(y)}" r="{NODE_RADIUS + 9}"/>')
        if stop["role"] != "stop":
            out.append(f'<circle class="nav-diagram-end nav-diagram-{stop["role"]}" cx="{_number(x)}" cy="{_number(y)}" r="{NODE_RADIUS + 5}"/>')
            if stop["role"] == "destination":
                out.append(f'<circle class="nav-diagram-end nav-diagram-destination-inner" cx="{_number(x)}" cy="{_number(y)}" r="{NODE_RADIUS + 2}"/>')
        out.append(glyph(stop["kind"], x, y))
        drop = NODE_RADIUS + 15 + (NAME_STAGGER if stop["index"] % 2 else 0)
        shown = _short_name(name)
        anchor, text_x = _anchored(x, shown, ribbon["width"])
        out.append(f'<text class="nav-diagram-name" x="{_number(text_x)}" y="{_number(y + drop)}" text-anchor="{anchor}">{_esc(shown)}</text>')
        out.append('</a>' if url else '</g>')
    out.append('</svg>')
    return "".join(out)


def render_legend(kinds, unknown=False, neighbours=False):
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
    if neighbours:
        items.append('<li><svg class="nav-legend-glyph" viewBox="0 0 34 22" width="34" height="22" aria-hidden="true" focusable="false">'
                     '<g class="nav-diagram-neighbour nav-diagram-stub"><line x1="4" y1="11" x2="29" y2="11"/><circle cx="29" cy="11" r="2.2"/></g></svg>'
                     ' Neighbour (not a travel lane)</li>')
    return f'<ul class="nav-diagram-legend" role="list">{"".join(items)}</ul>' if items else ""


def render_course_diagram(route, names, urls, linked=False):
    """
    The NAV page's course diagram panel for `route`, or an empty string when the route has fewer than two stops.

    Args:
        route (dict): The route `GET /api/nav` sends (`path`, `hops`, `stop_places`).
        names (dict): `{node: name}` of every stop.
        urls (dict): `{node: url or None}` of every stop.
        linked (bool): The diagram sits beside the Galaxy Map's course (NAV.69): says that a click on a stop centers the map on it.

    Returns:
        str: A complete `<section class="panel">` block.
    """
    path = (route or {}).get("path") or []
    if len(path) < 2:
        return ""
    places = {place["node"]: place for place in route.get("stop_places") or []}
    kinds = {navribbon._kind(node, places.get(node, {})) for node in path}
    unknown = sum(1 for hop in route.get("hops") or [] if hop.get("unknown_space"))
    has_graph = bool(route.get("graph"))
    drawings = []
    for size, width in WIDTHS:
        hidden = navribbon.layout_ribbon(route, width)["hidden_links"]
        hidden_note = (f'<p class="hint nav-diagram-hidden">{hidden} neighbour link{"" if hidden == 1 else "s"} between stops on '
                       f'different rows {"is" if hidden == 1 else "are"} not drawn.</p>' if hidden else "")
        drawings.append(f'<div class="nav-diagram-{size}-box">{render_drawing(route, names, urls, width, size)}{hidden_note}</div>')
    drawings = "".join(drawings)
    neighbour_note = (" The faint lines round a stop and the dotted arcs between stops show which systems are its nearest neighbours "
                      "(the few drawn, and \u201cand N more\u201d where there are more); a link is a neighbour relation, not a travel "
                      "lane. The route\u2019s own hops are always drawn." if has_graph else "")
    unknown_note = (f' {unknown} hop{"" if unknown == 1 else "s"} cross{"es" if unknown == 1 else ""} unknown space (uncharted '
                    f'sectors); {"it is" if unknown == 1 else "they are"} drawn in red, dashed and glowing, and the route list names '
                    f'{"it" if unknown == 1 else "them"} in words.' if unknown else "")
    link_note = (" A click on a stop selects it on the map and centers the map on it; a second click, or a click with Ctrl held, opens its page."
                 if linked else "")
    return (
        '<section class="panel nav-diagram" aria-labelledby="course-diagram-heading">'
        '<h2 id="course-diagram-heading">Course Diagram</h2>'
        '<p class="hint">The route as a line of stops, each hop with its length'
        ' (and its course where there is room). The spacing follows the logarithm of a hop’s length, so the diagram is a'
        f' schematic, not a map to scale.{_esc(unknown_note)}{_esc(neighbour_note)}{_esc(link_note)}</p>'
        f'<div class="nav-diagram-scroll">{drawings}</div>'
        f'{render_legend(kinds, bool(unknown), has_graph)}'
        '</section>'
    )
