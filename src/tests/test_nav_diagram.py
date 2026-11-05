# tests/test_nav_diagram.py

"""
NAV.67: the course diagram of the NAV page (`planetgen.web.maps.navdiagram`) drawn from fixed routes: one drawing per page size, a glyph
and a text alternative for every stop, the length on every hop, an unknown-space hop marked in a way that is not colour alone, and
nothing drawn for a route of fewer than two stops. Pure text, no database.
"""

import re

import pytest

from planetgen.web.maps import navdiagram, navribbon
from planetgen.web.maps.navdiagram import KIND_LABELS, WIDTHS, glyph, render_course_diagram, render_drawing


def _route(path, lengths, unknown=(), kinds=None, courses=True):
    hops = []
    for index, (a, b, length) in enumerate(zip(path, path[1:], lengths)):
        hop = {"from": a, "to": b, "distance_ly": length, "unknown_space": index in unknown}
        if courses:
            hop.update({"bearing_deg": 45.0 + index, "mark_deg": 3.0})
        hops.append(hop)
    return {"path": path, "hops": hops,
            "stop_places": [{"node": node, "sector_id": None, "sector_position_ly": None, **({"kind": kinds[node], "name": f"Stone {node}"}
                                                                                           if kinds and node in kinds else {})}
                            for node in path]}


def _names(path):
    return {node: f"Star {node}" if not isinstance(node, str) else f"Named {node}" for node in path}


def _urls(path):
    return {node: (None if isinstance(node, str) and node.startswith("stone:") else f"/system/{node}") for node in path}


PATH = [11, 12, 13, 14]


def test_a_route_is_drawn_once_for_each_page_size_with_every_stop_and_hop():
    route = _route(PATH, [4.2, 30.0, 650.0])
    html = render_course_diagram(route, _names(PATH), _urls(PATH))
    for size, width in WIDTHS:
        assert f'class="nav-diagram-svg nav-diagram-{size}"' in html
        assert f'viewBox="0 0 {width} ' in html
    assert html.count('class="nav-diagram-stop"') == 4 * len(WIDTHS)
    assert html.count("nav-diagram-label ") + html.count("nav-diagram-label\"") == 3 * len(WIDTHS)
    for name in ("Star 11", "Star 14"):
        assert name in html
    assert "Start: Star 11, star system" in html and "Destination: Star 14, star system" in html and "Stop: Star 12, star system" in html
    assert 'href="/system/12"' in html


def test_a_hop_says_its_length_and_its_course_in_its_title():
    route = _route(PATH, [4.2, 30.0, 650.0])
    html = render_drawing(route, _names(PATH), _urls(PATH), 880, "wide")
    assert "4.2 ly" in html or "pc" in html
    assert re.search(r"<title>[^<]*course 045 mark 003", html)


def test_the_course_is_left_out_of_a_gap_too_narrow_for_it_but_kept_in_the_title():
    route = _route(PATH, [1.0, 1.0, 1.0])
    narrow = render_drawing(route, _names(PATH), _urls(PATH), 260, "narrow")
    assert "nav-diagram-course" not in narrow
    assert "course 045 mark 003" in narrow


def test_an_unknown_space_hop_is_marked_red_dashed_and_in_words_not_by_colour_alone():
    route = _route(PATH, [4.2, 30.0, 650.0], unknown={1})
    html = render_course_diagram(route, _names(PATH), _urls(PATH))
    assert html.count('class="nav-diagram-hop nav-diagram-hop-unknown"') == len(WIDTHS)
    assert html.count('filter="url(#navdiagram-wide-glow)"') == 1
    assert "through unknown space (uncharted sectors)" in html
    assert "1 hop crosses unknown space" in html and "red, dashed and glowing" in html
    assert "Hop through unknown space" in html                           # the legend entry
    assert "nav-diagram-label-unknown" in html
    plain = render_course_diagram(_route(PATH, [4.2, 30.0, 650.0]), _names(PATH), _urls(PATH))
    assert "nav-diagram-hop-unknown" not in plain and "unknown space" not in plain


@pytest.mark.parametrize("kind", ["system", "phenomenon", "bright-star", "black-hole", "neutron-star", "quasar", "hypervelocity-star", "stone"])
def test_every_kind_of_stop_has_its_own_glyph_and_words(kind):
    shape = glyph(kind, 50, 50)
    assert shape and "<" in shape
    assert kind in KIND_LABELS


def test_the_glyphs_differ_by_shape_and_not_only_by_colour():
    # Moved to the origin so that only the outline is compared, not its place.
    outlines = {kind: glyph(kind, 0, 0) for kind in KIND_LABELS}
    assert len(set(outlines.values())) == len(KIND_LABELS)
    assert len({re.sub(r"[-\d.]+", "#", outline) for outline in outlines.values()}) >= 6   # and mostly by kind of shape, not only by size


def test_stops_in_unknown_space_have_no_link_and_a_page_stop_has_one():
    path = [11, "stone:bright-star:3", 14]
    route = _route(path, [4.2, 30.0], kinds={"stone:bright-star:3": "bright-star"})
    html = render_drawing(route, _names(path), _urls(path), 520, "middle")
    assert html.count("<a ") == 2
    assert 'data-stop="stone:bright-star:3"' in html and "bright star in unknown space" in html


def test_a_long_name_is_cut_in_the_drawing_and_whole_in_the_title():
    path = [11, 12]
    html = render_drawing(_route(path, [4.2]), {11: "A Very Long Star Name Indeed", 12: "Short"}, _urls(path), 520, "middle")
    assert "A Very Long St…" in html and "Start: A Very Long Star Name Indeed, star system" in html


def test_names_are_escaped():
    path = [11, 12]
    html = render_drawing(_route(path, [4.2]), {11: '<b>"x"</b>', 12: "ok"}, _urls(path), 520, "middle")
    assert "<b>" not in html and "&lt;b&gt;" in html


def test_a_route_of_fewer_than_two_stops_draws_nothing():
    assert render_course_diagram(None, {}, {}) == ""
    assert render_course_diagram(_route([11], []), {11: "A"}, {11: None}) == ""


def test_a_long_route_wraps_on_the_narrow_drawing_and_runs_on_the_wide_one():
    path = list(range(100, 112))
    route = _route(path, [40.0] * 11)
    narrow = render_drawing(route, _names(path), _urls(path), 260, "narrow")
    wide = render_drawing(route, _names(path), _urls(path), 880, "wide")
    height = lambda svg: float(re.search(r'viewBox="0 0 \d+ ([\d.]+)"', svg).group(1))
    assert height(narrow) > 2 * height(wide) >= navribbon.ROW_HEIGHT
