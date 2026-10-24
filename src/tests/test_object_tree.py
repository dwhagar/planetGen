# tests/test_object_tree.py

"""
UX.172: the summary tree's markup (`web/lib/objecttree.py`): one row per object, linked to its own page, closed rows with
children carry the route that lists them, escaped names, and "Show more" where a level is longer than a page.
"""

import dataclasses
import re

from planetgen.web.lib import objecttree

CONTEXT = objecttree.TreeContext(item_url=lambda ref: "/object/" + ref.replace(":", "-"),
                                 tree_url=lambda ref: "/object/" + ref.replace(":", "-") + "/tree")


def _item(kind, ref, name, child_count=0, counts=None, summary=None):
    return {"ref": ref, "kind": kind, "id": ref.partition(":")[2], "name": name, "child_count": child_count,
            "counts": counts or {}, "summary": summary or {}}


def test_a_row_with_children_has_its_route_a_toggle_and_a_link_for_no_script():
    html = objecttree.nodes_html([_item("planet", "planet:A", "Kepler b", 2, {"moon": 2})], CONTEXT)
    assert 'data-tree-url="/object/planet-A/tree"' in html and "data-open" not in html
    assert '<button type="button" class="tree-toggle" aria-expanded="false" aria-label="Moons of Kepler b" hidden>' in html
    assert 'href="/object/planet-A"' in html
    assert 'class="tree-pending"' in html
    assert "2 moons" in html
    # No link sits inside a button or a summary: a nested control to assistive technology (axe nested-interactive).
    assert "<summary" not in html and "<details" not in html


def test_a_row_without_children_is_a_plain_row_with_no_toggle():
    html = objecttree.nodes_html([_item("moon", "moon:B", "Tethys")], CONTEXT)
    assert "<button" not in html and 'class="tree-node tree-leaf"' in html


def test_a_name_from_the_database_is_escaped():
    html = objecttree.nodes_html([_item("planet", "planet:C", '<script>alert(1)</script>')], CONTEXT)
    assert "<script>" not in html and "&lt;script&gt;" in html


def test_a_loaded_level_is_drawn_open_with_its_children_and_its_show_more_row():
    parent = _item("star", "star:S", "Sol", 1, {"planet": 1})
    child = _item("planet", "planet:P", "Earth")
    html = objecttree.tree_html([parent], CONTEXT, "Bodies", loaded={"star:S": ([child], "TOKEN")})
    assert " data-open>" in html and 'aria-expanded="true"' in html
    assert 'href="/object/planet-P"' in html
    assert 'data-next="TOKEN"' in html


def test_a_tree_longer_than_a_page_ends_in_show_more_with_the_route():
    html = objecttree.tree_html([_item("moon", "moon:B", "Tethys")], CONTEXT, "Moons", token="NEXT", ref="planet:A")
    assert 'data-tree-root' in html and 'data-tree-url="/object/planet-A/tree"' in html
    assert html.count("data-tree-more") == 1


def test_counts_name_what_lies_below_in_order():
    assert objecttree.counts_text({"star": 1, "planet": 5, "moon": 12}) == "1 star, 5 planets, 12 moons"
    assert objecttree.counts_text({}) == ""


AU_KM = 149_597_870.7


def test_a_belt_row_shows_density_and_its_range():
    belt = _item("belt", "belt:5", "Asteroid belt 1", summary={"density": "sparse", "distance_km": 2.7 * AU_KM,
                                                              "lower_limit_km": 2.1 * AU_KM, "upper_limit_km": 3.3 * AU_KM})
    plain = re.sub(r"<[^>]*>", "", objecttree.nodes_html([belt], CONTEXT))
    assert "Sparse" in plain and "2.1 AU to 3.3 AU" in plain and "2.7 AU" not in plain


def test_a_belt_row_without_a_range_falls_back_to_its_distance():
    belt = _item("belt", "belt:5", "Asteroid belt 1", summary={"density": "dense", "distance_km": 2.7 * AU_KM})
    assert "2.7 AU" in objecttree.nodes_html([belt], CONTEXT)


def test_rows_show_no_zone():
    planet = _item("planet", "planet:1", "World", summary={"planet_class": "M", "distance_km": AU_KM, "zone": "Inner Zone"})
    assert "Inner Zone" not in objecttree.nodes_html([planet], CONTEXT)


def _scored(**tiers):
    values = {f"tier_{d}": 0 for d in ("pressure", "temperature", "air", "water", "radiation", "suitability")}
    return {"planet_class": "M", "phi4": 0.8, **values, **{f"tier_{k}": v for k, v in tiers.items()}}


def test_a_scored_body_shows_its_symbols_in_the_row():
    """UX.100: the overall symbol and the six domains, coloured by rating, in the row's summary line."""
    context = dataclasses.replace(CONTEXT, sprite_url="/static/icons.svg")
    body = _item("planet", "planet:1", "Rocky", summary=_scored(pressure=0, temperature=1, air=2, water=0, radiation=3,
                                                                suitability=1))
    html = objecttree.nodes_html([body], context)
    for label in ("Overall habitability: Green", "Pressure: Blue", "Temperature: Green", "Air chemistry: Yellow",
                  "Water chemistry: Blue", "Radiation: Red", "Suitability: Green"):
        assert f'aria-label="{label}"' in html, label
    assert html.count("hab-symbol ") == 7
    assert objecttree.nodes_html([body], CONTEXT).count("hab-symbol") == 0  # no sprite, no symbols
    unscored = _item("planet", "planet:2", "Giant", summary={"planet_class": "J"})
    assert "hab-symbol" not in objecttree.nodes_html([unscored], context)


def test_class_names_are_links_when_the_page_has_class_pages():
    context = dataclasses.replace(CONTEXT, class_url=lambda slug, code: None if code == "Z" else f"/classes/{slug}/{code}")
    planet = _item("planet", "planet:1", "Terra", summary={"planet_class": "M"})
    star = _item("star", "star:1", "Sol", summary={"star_type": "G2V Yellow Main Sequence Star", "role": "single"})
    comet = _item("comet", "comet:1", "Kohoutek", summary={"orbit_type": "elliptical", "period_class": "halley_type"})
    html = objecttree.nodes_html([planet, star, comet], context)
    assert 'Class <a href="/classes/planet/M">M</a>' in html
    assert '<a href="/classes/star-spectral/G">G2</a><a href="/classes/star-luminosity/V">V</a> Yellow Main Sequence Star' in html
    assert "/classes/comet/halley_type" in html
    assert "/classes/" not in objecttree.nodes_html([planet, star, comet], CONTEXT)
    other = _item("planet", "planet:3", "Odd", summary={"planet_class": "Z"})
    assert "/classes/planet/Z" not in objecttree.nodes_html([other], context)
