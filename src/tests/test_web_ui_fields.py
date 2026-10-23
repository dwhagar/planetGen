# tests/test_web_ui_fields.py

"""
The one form field (UX.124, `templates/partials/ui_fields.html`): a labelled Shoelace `sl-input`, `sl-select` or
`sl-checkbox` with its name, required flag and an error slot, each with a plain twin in a <noscript> so a form still
posts without script. The macros alone, then the pilot pages (What's nearby, the search page) built from them.
"""

import re

import pytest

from tests.test_web_nearby import _first_system, client  # noqa: F401 -- the nearby fixtures

IMPORT = '{% from "partials/ui_fields.html" import ui_input, ui_select, ui_checkbox %}'


def _render(client, source, **context):
    """Render `source` (after the macro import) in the web app's template environment."""
    app = client.application
    with app.test_request_context("/"):
        return app.jinja_env.from_string(IMPORT + source).render(**context)


def _noscript(html):
    return re.findall(r"<noscript>(.*?)</noscript>", html, re.S)


def test_an_input_carries_its_name_label_value_and_required_flag(client):
    html = _render(client, '{{ ui_input("q", "Name", "Sol", type="search", required=True, placeholder="Any") }}')
    element = re.search(r"<sl-input[^>]*>", html).group(0)
    for part in ('name="q"', 'label="Name"', 'value="Sol"', 'type="search"', "required", 'placeholder="Any"'):
        assert part in element, element
    assert "data-invalid" not in element


def test_a_validation_error_sits_in_the_field_it_is_about_and_marks_it(client):
    html = _render(client, '{{ ui_input("distance", "Distance", "500", error="too far") }}'
                           '{{ ui_input("place", "Place", "") }}')
    _, first, second = html.split('<div class="ui-field">', 2)
    assert re.search(r'<sl-input name="distance"[^>]*data-invalid[^>]*>\s*<span slot="help-text" class="field-error" '
                     r'role="alert">too far</span>', first), first
    assert "field-error" not in second and "data-invalid" not in second


def test_each_field_has_a_plain_twin_that_posts_the_same_name_and_value(client):
    options = [{"value": "a", "label": "Alpha"}, {"value": "b", "label": "Beta"}]
    html = _render(client, '{{ ui_input("q", "Name", "Sol") }}{{ ui_select("pick", "Pick one", options, "b") }}'
                           '{{ ui_checkbox("border", "Also the border", "1", True) }}', options=options)
    twins = _noscript(html)
    assert len(twins) == 3
    assert re.search(r'<input type="text" name="q" value="Sol"', twins[0])
    assert re.search(r'<select name="pick">.*<option value="b" selected>Beta</option>', twins[1], re.S)
    assert re.search(r'<input type="checkbox" name="border" value="1" checked>', twins[2])


def test_a_checkbox_label_can_be_markup_from_a_call_block(client):
    html = _render(client, '{% call ui_checkbox("type", value="black_hole") %}<b>Black Hole</b> (3){% endcall %}')
    assert re.search(r'<sl-checkbox name="type" value="black_hole">\s*<b>Black Hole</b> \(3\)', html)
    assert "<b>Black Hole</b> (3)</label>" in html


def test_a_value_with_markup_is_escaped(client):
    html = _render(client, '{{ ui_input("q", "Name", value) }}', value='"><script>x</script>')
    assert "<script>" not in html and "&#34;&gt;&lt;script&gt;" in html


def test_a_bad_distance_is_shown_on_the_distance_field_not_above_the_form(client, mysql_config):
    system = _first_system(mysql_config)
    html = client.get(f"/nearby?place=system:{system['pid']}&distance=500").get_data(as_text=True)
    assert re.search(r'<sl-input name="distance"[^>]*data-invalid[^>]*>\s*<span slot="help-text" class="field-error" '
                     r'role="alert">the largest distance is 50 pc</span>', html), html
    assert '<p class="error" role="alert">' not in html


def test_a_refusal_that_is_not_about_a_field_stays_above_the_form(client):
    html = client.get("/nearby?place=1,2&distance=5").get_data(as_text=True)
    assert '<p class="error" role="alert">a point is three finite numbers' in html
    assert "data-invalid" not in html


@pytest.mark.parametrize("path, names", [
    ("/nearby?place=0,0,0", {"place": "hidden", "distance": "number", "kinds": "checkbox"}),
    ("/search", {"q": "search", "star_min_radius_km": "number", "planet_max_radius_km": "number"}),
])
def test_the_pilot_forms_post_the_same_fields_without_script(client, path, names):
    """With scripts off the element is never defined; the <noscript> twin carries each field by name and type."""
    html = client.get(path).get_data(as_text=True)
    plain = "".join(_noscript(html)) + re.sub(r"<noscript>.*?</noscript>", "", html, flags=re.S)
    for name, kind in names.items():
        assert re.search(rf'<input[^>]*type="{kind}"[^>]*name="{name}"|<input[^>]*name="{name}"[^>]*type="{kind}"',
                         plain), (path, name, kind)


def test_the_pilot_pages_draw_their_fields_with_the_macro(client):
    nearby = client.get("/nearby?place=0,0,0").get_data(as_text=True)
    assert 'sl-input name="distance"' in nearby and 'sl-checkbox name="kinds"' in nearby
    search = client.get("/search").get_data(as_text=True)
    assert 'sl-input name="q"' in search and 'sl-input name="star_min_radius_km"' in search
    assert 'sl-select name="sector"' in client.get("/nearby").get_data(as_text=True)
