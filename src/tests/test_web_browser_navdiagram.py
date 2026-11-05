"""
NAV.67, NAV.36: the NAV page's course diagram in a real browser. At a phone and a desktop width one drawing of the route shows and the
page does not scroll sideways, its text is drawn near its set size (about 12 px), no two names or two hop labels overlap, an
unknown-space hop is red and dashed in both themes and glows only where motion is allowed, and a stop is a link to its page.

The panel comes straight from `planetgen.web.maps.navdiagram.render_course_diagram` with the site's own `style.css`, so no database is
needed. Skipped without Playwright or Chromium.
"""

import os
import re

import pytest

sync_api = pytest.importorskip("playwright.sync_api")

from planetgen.web.maps.navdiagram import render_course_diagram  # noqa: E402

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_STYLE = os.path.join(_SRC_DIR, "html", "static", "style.css")


def _route():
    """Eleven stops: short hops, one of thousands of light-years, two through unknown space, a black hole in unknown space among them."""
    path = [1, 2, 3, 4, "stone:black-hole:9", 5, 6, 7, 8, 9, 10]
    lengths = [3.1, 12.0, 0.9, 140.0, 45.0, 6500.0, 8.0, 22.0, 2.4, 61.0]
    hops = [{"from": a, "to": b, "distance_ly": length, "unknown_space": index in (3, 5), "bearing_deg": 10.0 * index, "mark_deg": 2.0}
            for index, (a, b, length) in enumerate(zip(path, path[1:], lengths))]
    places = [{"node": node, "sector_id": None, "sector_position_ly": None,
               **({"kind": "black-hole", "name": "Dark Lantern"} if node == "stone:black-hole:9" else {})} for node in path]
    names = {node: f"Kepler Station {node}" for node in path}
    names["stone:black-hole:9"] = "Dark Lantern"
    urls = {node: (None if isinstance(node, str) else f"/system/{node}") for node in path}
    return {"path": path, "hops": hops, "stop_places": places}, names, urls


@pytest.fixture(scope="module")
def browser():
    try:
        playwright = sync_api.sync_playwright().start()
    except Exception as exc:  # noqa: BLE001 -- any failure here means "no browser"
        pytest.skip(f"Playwright could not start ({exc})")
    try:
        chromium = playwright.chromium.launch()
    except Exception as exc:  # noqa: BLE001
        playwright.stop()
        pytest.skip(f"No Chromium for Playwright ({str(exc).splitlines()[0]})")
    try:
        yield chromium
    finally:
        chromium.close()
        playwright.stop()


def _edge(a, b, length, on_route=False):
    return {"from": a, "to": b, "distance_ly": length, "on_route": on_route}


def _with_neighbours(route, dense):
    """NAV.68: `route` with its graph: a dense case (three neighbours on every stop, 'and N more' on the rest, links between stops) or a
    sparse one (one neighbour on one stop)."""
    path = route["path"]
    edges = [_edge(a, b, hop["distance_ly"], True) for a, b, hop in zip(path, path[1:], route["hops"])]
    more = {}
    if dense:
        for index, stop in enumerate(path):
            edges += [_edge(stop, 9000 + 10 * index + k, 0.8 + 0.4 * k) for k in range(3)]
            if index % 2 == 0:
                more[stop] = 3 + index
        edges += [_edge(path[0], path[2], 15.0), _edge(path[1], path[3], 150.0), _edge(path[6], path[8], 28.0)]
    else:
        edges.append(_edge(path[2], 9000, 1.7))
    return {**route, "graph": {"nodes": [], "edges": edges, "more": more}}


def _open(browser, width, theme="light", motion="no-preference", neighbours=None):
    with open(_STYLE, encoding="utf-8") as handle:
        css = handle.read()
    route, names, urls = _route()
    if neighbours is not None:
        route = _with_neighbours(route, neighbours == "dense")
    panel = render_course_diagram(route, names, urls)
    context = browser.new_context(viewport={"width": width, "height": 900}, device_scale_factor=1, reduced_motion=motion)
    page = context.new_page()
    page.set_content(f"<!doctype html><html data-theme=\"{theme}\"><head><meta name=\"viewport\" content=\"width=device-width\">"
                     f"<style>{css}</style></head><body><main>{panel}</main></body></html>")
    return context, page


SHOWN = """
() => {
  const svg = [...document.querySelectorAll(".nav-diagram-svg")].filter((el) => getComputedStyle(el).display !== "none"
    && getComputedStyle(el.parentElement).display !== "none");
  const only = svg[0];
  const box = only.getBoundingClientRect();
  const scale = box.width / only.viewBox.baseVal.width;
  const rect = (el) => { const r = el.getBoundingClientRect(); return { text: el.textContent, left: r.left, right: r.right, top: r.top, bottom: r.bottom }; };
  return {
    drawings: svg.length, size: only.getAttribute("class"), scale, boxWidth: box.width,
    pageOverflow: document.documentElement.scrollWidth - document.documentElement.clientWidth,
    names: [...only.querySelectorAll(".nav-diagram-name")].map(rect),
    labels: [...only.querySelectorAll(".nav-diagram-label")].map(rect),
    courses: [...only.querySelectorAll(".nav-diagram-course")].map(rect),
    more: [...only.querySelectorAll(".nav-diagram-more")].map(rect),
    neighbourLabels: [...only.querySelectorAll(".nav-diagram-neighbour-label")].map(rect),
    stubs: only.querySelectorAll(".nav-diagram-stub").length,
    arcs: only.querySelectorAll(".nav-diagram-link").length,
    inside: [...only.querySelectorAll(".nav-diagram-stub, .nav-diagram-link, .nav-diagram-more")].every((el) => {
      const r = el.getBoundingClientRect(); return r.left >= box.left - 1 && r.right <= box.right + 1 && r.top >= box.top - 1 && r.bottom <= box.bottom + 1; }),
    fontPx: parseFloat(getComputedStyle(only.querySelector(".nav-diagram-name")).fontSize) * scale,
    links: [...only.querySelectorAll("a.nav-diagram-stop")].map((a) => a.getAttribute("href")),
    stops: only.querySelectorAll(".nav-diagram-stop").length,
  };
}
"""


def _overlaps(a, b):
    return a["left"] < b["right"] - 0.5 and b["left"] < a["right"] - 0.5 and a["top"] < b["bottom"] - 0.5 and b["top"] < a["bottom"] - 0.5


@pytest.mark.parametrize("width,size", [(320, "narrow"), (700, "middle"), (1280, "wide")])
def test_one_drawing_shows_at_each_width_readable_and_without_overlap(browser, width, size):
    context, page = _open(browser, width)
    try:
        result = page.evaluate(SHOWN)
    finally:
        context.close()
    assert result["drawings"] == 1 and f"nav-diagram-{size}" in result["size"], result
    assert result["pageOverflow"] <= 0, f"the page scrolls sideways by {result['pageOverflow']}px at {width}px"
    assert result["stops"] == 11 and len(result["links"]) == 10         # the black hole in unknown space has no page to link to
    assert result["fontPx"] >= 10.0, f"the names are drawn at {result['fontPx']}px at {width}px"
    texts = result["names"] + result["labels"] + result["courses"] + result.get("more", []) + result.get("neighbourLabels", [])
    for index, a in enumerate(texts):
        for b in texts[index + 1:]:
            assert not _overlaps(a, b), f"{a['text']!r} overlaps {b['text']!r} at {width}px"


@pytest.mark.parametrize("theme", ["light", "dark"])
def test_an_unknown_space_hop_is_red_and_dashed_in_both_themes(browser, theme):
    context, page = _open(browser, 1280, theme=theme)
    try:
        line = page.evaluate("""() => {
          const shown = [...document.querySelectorAll(".nav-diagram-wide-box .nav-diagram-hop-unknown .nav-diagram-line")][0];
          const plain = [...document.querySelectorAll(".nav-diagram-wide-box .nav-diagram-hop:not(.nav-diagram-hop-unknown) .nav-diagram-line")][0];
          const error = getComputedStyle(document.documentElement).getPropertyValue("--error").trim();
          const css = (el) => getComputedStyle(el);
          return { stroke: css(shown).stroke, dash: css(shown).strokeDasharray, plain: css(plain).stroke, plainDash: css(plain).strokeDasharray,
                   error, filter: shown.parentElement.getAttribute("filter") };
        }""")
    finally:
        context.close()
    channels = tuple(int(line["error"][i:i + 2], 16) for i in (1, 3, 5))
    assert line["stroke"] == "rgb({}, {}, {})".format(*channels), line
    assert line["dash"] != "none" and line["plainDash"] == "none" and line["stroke"] != line["plain"]
    assert line["filter"].startswith("url(#")


@pytest.mark.parametrize("motion,animated", [("no-preference", True), ("reduce", False)])
def test_the_glow_pulses_only_where_motion_is_allowed(browser, motion, animated):
    context, page = _open(browser, 1280, motion=motion)
    try:
        name = page.evaluate("""() => getComputedStyle(document.querySelector(".nav-diagram-wide-box .nav-diagram-hop-unknown")).animationName""")
        # Red and glowing still, with the pulse off.
        glow = page.evaluate("""() => document.querySelector(".nav-diagram-wide-box .nav-diagram-hop-unknown").getAttribute("filter")""")
    finally:
        context.close()
    assert (name != "none") is animated, name
    assert glow.startswith("url(#")


@pytest.mark.parametrize("width,size", [(320, "narrow"), (700, "middle"), (1280, "wide")])
@pytest.mark.parametrize("neighbours", ["dense", "sparse"])
def test_neighbours_are_drawn_inside_the_panel_without_covering_the_text(browser, width, size, neighbours):
    """NAV.68: faint stubs, 'and N more' and arcs round the stops, in a dense and a sparse case: all inside the drawing, no text over other text."""
    context, page = _open(browser, width, neighbours=neighbours)
    try:
        result = page.evaluate(SHOWN)
    finally:
        context.close()
    assert f"nav-diagram-{size}" in result["size"]
    assert result["stubs"] >= (10 if neighbours == "dense" else 1), result["stubs"]
    if neighbours == "dense":
        assert result["more"], "a stop the cap cut says 'and N more'"
        assert result["arcs"] >= (0 if size == "narrow" else 1)
    assert result["inside"], f"a stub, arc or 'and N more' runs outside the drawing at {width}px"
    texts = result["names"] + result["labels"] + result["courses"] + result["more"] + result["neighbourLabels"]
    for index, a in enumerate(texts):
        for b in texts[index + 1:]:
            assert not _overlaps(a, b), f"{a['text']!r} overlaps {b['text']!r} at {width}px ({neighbours})"
    for sample in result["more"]:
        assert re.search(r"and \d+ more$", sample["text"])


def test_the_page_says_a_link_is_a_neighbour_relation_and_not_a_travel_lane(browser):
    context, page = _open(browser, 700, neighbours="dense")
    try:
        text = page.evaluate("() => document.querySelector('.nav-diagram .hint').textContent")
        legend = page.evaluate("() => document.querySelector('.nav-diagram-legend').textContent")
    finally:
        context.close()
    assert "neighbour relation, not a travel lane" in text and "Neighbour (not a travel lane)" in legend
