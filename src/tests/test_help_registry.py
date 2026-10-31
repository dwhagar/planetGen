# tests/test_help_registry.py

"""
The help pages are built from the program's own tables, never from a stored list: the registry (`web/lib/help_registry.py`) names each
help page once, the Help menu, the help index and the sitemap read it, and the numbers inside the habitability prose are read from
the limits tables. These tests fail when a page and its source drift apart:

* every route under `/help` is in the registry, and every registry page is a route (the coverage rule of DOC.58 for the help pages);
* the Help menu and the help index list exactly what the registry lists;
* the index lists every class type, every life topic and every life class the program has, so a new one appears with no edit;
* the numbers the habitability page states (gear limits, score edges, life class count, temperature bands) are the ones in the code.
"""

import math
import re

from markupsafe import escape

from planetgen.physics import life_classes, regional_temperature
from planetgen.physics.protection_limits import limit
from planetgen.web import helpers
from planetgen.web.lib import classref, habitability_ref, help_registry, life_help
from tests.test_class_pages import _html, app, client  # noqa: F401


def _help_endpoints(app):  # noqa: F811
    return {rule.endpoint.split(".", 1)[1] for rule in app.url_map.iter_rules() if rule.rule == "/help" or rule.rule.startswith("/help/")}


def test_every_help_route_is_in_the_registry_and_the_other_way_round(app):  # noqa: F811
    registered = {help_registry.INDEX_ENDPOINT}
    for topic in help_registry.topics():
        registered.add(topic.endpoint)
        registered |= {child.endpoint for child in help_registry.children(topic)}
    assert _help_endpoints(app) <= registered | {"classes"}, "a help route is missing from help_registry.TOPICS"
    for endpoint in registered:
        assert f"web.{endpoint}" in {rule.endpoint for rule in app.url_map.iter_rules()}, endpoint


def test_every_topic_has_a_known_category_and_a_summary():
    slugs = {slug for slug, _title, _summary in help_registry.CATEGORIES}
    for topic in help_registry.topics():
        assert topic.category in slugs and topic.summary and topic.label
    assert [t for c in help_registry.categories() for t in c["topics"]] == sorted(
        help_registry.topics(), key=lambda t: [c[0] for c in help_registry.CATEGORIES].index(t.category))


def test_the_menu_is_the_index_then_the_registry(client):  # noqa: F811
    assert helpers.HELP_MENU == (("help_index", "All help"),) + help_registry.menu()
    html = _html(client, "/help")
    bar = html.split('<header class="site-header">')[1].split("</header>")[0]
    menu = re.search(r'<sl-dropdown class="site-help-menu".*?</sl-dropdown>', bar, re.S).group(0)
    assert re.findall(r'data-href="([^"]*)"', menu) == ["/help", "/classes", "/help/habitability", "/help/properties", "/help/life"]


def test_the_index_lists_every_class_type_and_life_page_the_program_has(client):  # noqa: F811
    html = _html(client, "/help")
    for slug, entry in classref.catalog().items():
        assert f'<a href="/classes/{slug}">{escape(entry["name"])}</a>' in html, slug
    for slug, _title, _summary in life_help.TOPICS:
        assert (f'href="/help/life/{slug}"' if slug != "kinds" else 'href="/help/life"') in html, slug
    for number in {row.number for row in life_classes.CLASSES}:
        assert f'href="/help/life/classes/{number}"' in html, number
    for category in help_registry.categories():
        assert escape(category["title"]) in html
        for topic in category["topics"]:
            assert escape(topic.summary) in html


def test_a_class_type_added_to_the_catalog_reaches_the_index(client, monkeypatch):  # noqa: F811
    real = classref.catalog()
    extra = dict(real, invented={"slug": "invented", "name": "Invented Things", "label": "Invented", "summary": "A test type.",
                                 "notes": [], "classes": {"x": {"code": "x", "name": "X", "summary": "x", "facts": []}}})
    monkeypatch.setattr(classref, "catalog", lambda: extra)
    labels = [child.label for child in help_registry.children(help_registry.topics()[0])]
    assert "Invented Things" in labels and len(labels) == len(extra)


def test_the_index_is_in_the_sitemap(app):  # noqa: F811
    from planetgen.web import seo
    with app.test_request_context("/"):
        assert "/help" in seo.top_page_paths()


# --- The numbers on the habitability page are the code's ---------------------------------------------------------------

def _equipment(label):
    return dict(habitability_ref.EQUIPMENT_ITEMS)[label]


def test_the_equipment_limits_on_the_page_are_the_limits_the_score_uses():
    vehicle = _equipment("Vehicle only")
    assert f"heat past {limit('vehicle_heat_c'):g} C" in vehicle and f"{limit('vehicle_gravity_g'):g} g" in vehicle
    assert f"{limit('vehicle_pressure_kpa') / 1000:g} MPa" in vehicle
    assert f"heat past {limit('armoured_heat_c'):g} C" in _equipment("Armoured suit")
    frame = _equipment("Powered frame")
    assert f"{limit('frame_gravity_g'):g} g" in frame and f"{limit('frame_gravity_rec_g'):g} g" in frame
    dust = _equipment("Dust filter or seals")
    assert f"{limit('dust_recommended_mg_m3'):g} mg/m3" in dust and f"{limit('dust_required_mg_m3'):g} mg/m3" in dust


def test_the_score_rule_states_the_edge_scores_of_the_code():
    from planetgen.physics import habitability
    rule = habitability_ref.SCORE_RULE
    for score in habitability._EDGE_SCORES[1:3]:
        assert f"{score:g}" in rule


def test_the_map_row_counts_the_life_classes_of_the_rule_file():
    row = dict(habitability_ref.MAP_ROWS)["Marker: Life classes"]
    assert f"each of the {len(life_classes.CLASSES)} life classes" in row
    assert life_classes.CLASSES[0].id in row


def test_the_temperature_text_states_the_bands_of_regional_temperature():
    text = habitability_ref.temperature_text()
    for name, low, high in regional_temperature.BANDS:
        share = (math.sin(math.radians(high)) - math.sin(math.radians(low))) / 2.0 * (2.0 if name == "polar" else 1.0)
        assert f"{share * 100:.0f} percent" in text, name
    assert {4: "four"}[len(regional_temperature.BANDS)] in text
