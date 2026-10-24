# tests/test_object_tree.py

"""
UX.172: the summary tree's markup (`web/lib/objecttree.py`): one row per object, linked to its own page, closed rows with
children carry the route that lists them, escaped names, and "Show more" where a level is longer than a page.
"""

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
