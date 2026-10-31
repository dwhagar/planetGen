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

import pytest

from markupsafe import escape

from planetgen.physics import life_classes, regional_temperature
from planetgen.physics.protection_limits import limit
from planetgen.web import helpers
from planetgen.web.lib import classref, habitability_ref, help_registry, life_help, property_help
from tests.test_class_pages import _html, app, client  # noqa: F401


def _help_endpoints(app):  # noqa: F811
    return {rule.endpoint.split(".", 1)[1] for rule in app.url_map.iter_rules() if rule.rule == "/help" or rule.rule.startswith("/help/")}


def test_every_help_route_is_in_the_registry_and_the_other_way_round(app):  # noqa: F811
    registered = {help_registry.INDEX_ENDPOINT}  # "classes" is the registered root of the class pages
    for topic in help_registry.topics():
        registered.add(topic.endpoint)
        registered |= {child.endpoint for child in help_registry.children(topic)}
    assert _help_endpoints(app) <= registered | {"classes"}, "a help route is missing from help_registry.TOPICS"
    for endpoint in registered:
        assert f"web.{endpoint}" in {rule.endpoint for rule in app.url_map.iter_rules()}, endpoint


def test_every_topic_has_a_registered_category_and_a_summary():
    slugs = {entry["slug"] for entry in help_registry.categories()}
    for topic in help_registry.topics():
        assert topic.category in slugs and topic.summary and topic.label
    assert [t for c in help_registry.categories() for t in c["topics"]] != [] and len(help_registry.topics()) >= 4


def test_each_part_of_the_program_registers_its_own_page():
    """The pages come from the modules that own them, not from a list kept here."""
    assert {t.endpoint for t in help_registry.topics()} >= {"classes", "habitability_levels", "property_help_page", "life_help_page"}
    owner = {"classes": classref, "habitability_levels": habitability_ref, "property_help_page": property_help, "life_help_page": life_help}
    for endpoint, module in owner.items():
        source = open(module.__file__, encoding="utf-8").read()
        assert f'register_topic("{endpoint}"' in source, endpoint
    assert "TOPICS" not in open(help_registry.__file__, encoding="utf-8").read().split('"""', 2)[2]


def test_a_page_registered_by_a_new_part_reaches_the_menu_the_index_and_the_terms(client, monkeypatch):  # noqa: F811
    monkeypatch.setattr(help_registry, "_topics", dict(help_registry._topics))
    monkeypatch.setattr(help_registry, "_categories", dict(help_registry._categories))
    monkeypatch.setattr(help_registry, "_term_sources", dict(help_registry._term_sources))
    help_registry.register_category("tides", "Tides and spin", "How tides slow a spin.", order=30)
    help_registry.register_topic("property_help_page", "Tide test", "tides", "Tide summary.", order=5,
                                 children=lambda: [help_registry.Child("classes", {}, "Child page", "")])
    help_registry.register_terms("tides", lambda: [help_registry.Term("tidal lock", "classes", {})])
    assert help_registry.menu()[0] == ("property_help_page", "Tide test")
    assert [c["slug"] for c in help_registry.categories()][-1] == "tides"
    assert help_registry.Term("tidal lock", "classes", {}) in help_registry.terms()
    html = _html(client, "/help")
    assert "Tides and spin" in html and "Tide summary." in html and "Child page" in html


def test_a_topic_in_an_unregistered_category_is_refused(monkeypatch):
    monkeypatch.setattr(help_registry, "_topics", dict(help_registry._topics))
    help_registry.register_topic("classes", "Classes", "nowhere", "x")
    with pytest.raises(KeyError):
        help_registry.categories()


def test_the_life_class_ids_are_link_terms_for_their_help_page():
    found = {term.text: term for term in help_registry.terms()}
    for row in life_classes.CLASSES:
        assert found[row.id].endpoint == "life_help_class" and found[row.id].values == {"number": row.number}


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
