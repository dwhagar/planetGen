# tests/test_habitability_chips.py

"""
The habitability chips of a planet's or moon's page (`planetgen/web/lib/habitability_chips.py`, UX.154, UX.181): the flags,
the equipment chips (solid when required, outlined when recommended, in the colour of their PHI-6 domain) and their links
into the habitability page. Plain dicts in the `queryDb.system_detail` shape; no database.
"""

from planetgen.web.lib.habitability_chips import chips_html

URL = "/help/habitability"


def _scored(**fields):
    values = {f"tier_{domain}": 0 for domain in ("pressure", "temperature", "air", "water", "radiation", "suitability")}
    body = {"habitable": False, "habitable_class": False, "inhabited": False, "moons": [], "phi4": 0.8, **values, **fields}
    body.setdefault("equipment", {"required": ["gas-mix mask"], "recommended": [], "line": "", "code": "S1", "items": [
        {"slot": "supply", "item": "gas-mix mask", "severity": "required", "domain": "Air"}]})
    return body


def test_a_habitable_world_and_a_habitable_class_get_their_own_chips():
    """GEN.643: "Habitable" is a Blue or Green overall rating, "Habitable class" a class that can carry life."""
    both = chips_html(_scored(habitable=True, habitable_class=True))
    assert "Habitable</span>" in both and "Habitable class</span>" in both
    only_class = chips_html(_scored(habitable_class=True))
    assert "Habitable class</span>" in only_class and "Habitable</span>" not in only_class
    assert "Habitable" not in chips_html(_scored())
    assert "Inhabited" in chips_html(_scored(inhabited=True))


def test_a_planet_with_a_habitable_moon_says_so():
    html = chips_html(_scored(moons=[{"habitable": True, "habitable_class": True}]))
    assert "Habitable moon" in html and "Habitable-class moon" in html
    only_class = chips_html(_scored(moons=[{"habitable": False, "habitable_class": True}]))
    assert "Habitable-class moon" in only_class and "Habitable moon" not in only_class
    assert "moon" not in chips_html(_scored(moons=[{"habitable": False, "habitable_class": False}]))


def test_equipment_chips_are_solid_when_required_and_outlined_when_recommended():
    body = _scored(tier_air=2, tier_radiation=3)
    body["equipment"]["items"] = [
        {"slot": "supply", "item": "gas-mix mask", "severity": "required", "domain": "Air"},
        {"slot": "radiation", "item": "hardened habitat", "severity": "recommended", "domain": "Radiation"}]
    html = chips_html(body)
    assert 'flag-eq-required flag-tier flag-tier-2">Gas-mix mask</span>' in html
    assert 'flag-eq-recommended flag-tier flag-tier-3">Hardened habitat <span class="sr-only">(recommended)</span>' in html
    body["equipment"]["items"] = []
    assert "PHI-6 Red" in chips_html(body)


def test_a_body_with_no_scores_has_no_chips():
    assert chips_html({"habitable": False, "equipment": None}) == ""
    assert chips_html({}) == ""


def test_the_chips_link_to_their_sections_of_the_habitability_page():
    html = chips_html(_scored(habitable=True, habitable_class=True), URL)
    assert f'<a class="help-link" href="{URL}#chip-habitable"><span class="flag flag-yes">Habitable</span></a>' in html
    assert f'<a class="help-link" href="{URL}#chip-habitable-class"><span class="flag flag-yes">Habitable class</span></a>' in html
    assert f'<a class="help-link" href="{URL}#equipment-gas-mix-mask"><span class="flag flag-eq' in html
    body = _scored()
    body["equipment"]["items"] = []
    assert f'href="{URL}#overall"><span class="flag flag-tier flag-tier-0">PHI-6 Blue</span>' in chips_html(body, URL)
    assert "help-link" not in chips_html(_scored(habitable=True, habitable_class=True))
