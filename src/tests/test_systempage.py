# tests/test_systempage.py

"""
The system page's body list (`planetgen/web/lib/systempage.py`, UX.7 and
UX.12): the habitable chip, the PHI-6 symbols (UX.100), the habitable-moon chip, moons in their own group,
belt distances, and comets ordered in among the planets by semi-major
axis. Plain dicts in the `queryDb.system_detail` shape; no database.
"""


import re

from planetgen.web.lib.systempage import _planet_row_html, comet_orbit_key_km, system_list_html  # noqa: E402

AU_KM = 149_597_870.7
_SECTIONS = {"overview": [], "stars": {}, "planets": {}, "moons": {}, "belts": {}, "comets": {}}


def _body(id_, name, au, body_type="t", planet_class="M", habitable=False, moons=(), index=0):
    return {
        "id": id_, "name": name, "planet_class": planet_class, "body_type": body_type,
        "habitable": habitable, "inhabited": False, "zone": "Inner Zone",
        "distance_km": au * AU_KM, "period_years": 1.0, "gravity_g": 1.0,
        "orbital_index": index, "moons": list(moons), "star_id": 1,
    }


def _comet(id_, name, perihelion_au, eccentricity, orbit_type="elliptical"):
    return {
        "id": id_, "name": name, "orbit_type": orbit_type, "is_active": False,
        "perihelion_distance_km": perihelion_au * AU_KM, "eccentricity": eccentricity,
        "orbital_period_years": 10.0 if orbit_type == "elliptical" else None, "star_id": 1,
    }


def _system(planets=(), belts=(), comets=()):
    star = {"id": 1, "name": "Sun", "role": "single", "star_type": "G2V"}
    return {"stars": [star], "planets": list(planets), "belts": list(belts), "comets": list(comets),
            "binary_configuration": None}


def test_a_row_names_its_class_and_no_type_label():
    """UX.32: a row shows the class and a Habitable class chip, never a Terrestrial or Gas Giant label."""
    html = system_list_html(_system([
        _body(1, "Rocky", 1.0, habitable=True, index=0),
        _body(2, "Dry", 1.5, index=1),
        _body(3, "Giant", 5.2, body_type="g", planet_class="J", index=2),
    ]), _SECTIONS)
    assert "Habitable class" in html[html.index("Rocky"):html.index("Dry")]
    assert html.count("Habitable class") == 1
    assert "Terrestrial" not in html and "Gas Giant" not in html
    assert "Habitable: " not in html


def test_habitable_moon_chip_and_moon_group():
    moon = _body(9, "Pandora", 0.002, habitable=True)
    html = system_list_html(_system([_body(1, "Polyphemus", 5.0, body_type="g", moons=[moon])]), _SECTIONS)
    assert "Habitable-class moon" in html
    assert '<details class="moon-group"><summary>1 moon of Polyphemus</summary>' in html
    assert '<span class="stat">1 moon</span>' in html  # the count on the planet's row
    # The group follows the planet's own row rather than sitting inside it.
    assert html.index("</details>") < html.index('class="moon-group"')
    no_moon = system_list_html(_system([_body(1, "Plain", 1.0, moons=[_body(9, "Rock", 0.002)])]), _SECTIONS)
    assert "Habitable-class moon" not in no_moon


def test_belt_row_shows_density_range_and_top_minerals():
    composition = [{"component": c, "concentration": k} for c, k in
                   (("iron", "high"), ("nickel", "moderate"), ("olivine", "small"), ("gold", "trace"))]
    belt = {"id": 5, "density": "sparse", "distance_km": 2.7 * AU_KM, "lower_limit_km": 2.1 * AU_KM,
            "upper_limit_km": 3.3 * AU_KM, "orbital_index": 0, "star_id": 1, "composition": composition}
    html = system_list_html(_system(belts=[belt]), _SECTIONS)
    plain = re.sub(r"<[^>]*>", "", html)
    assert "Sparse" in plain and "2.1 AU to 3.3 AU" in plain
    assert "2.7 AU" not in plain  # the nominal distance isn't repeated next to the range
    assert "Iron, nickel, olivine" in html and "gold" not in html


def test_belt_row_without_range_falls_back_to_distance():
    belt = {"id": 5, "density": "dense", "distance_km": 2.7 * AU_KM, "lower_limit_km": None,
            "upper_limit_km": None, "orbital_index": 0, "star_id": 1}
    html = system_list_html(_system(belts=[belt]), _SECTIONS)
    assert "2.7 AU" in html


def test_rows_show_no_zone():
    moon = _body(9, "Rock", 0.002)
    html = system_list_html(_system([_body(1, "World", 1.0, moons=[moon])]), _SECTIONS)
    assert "Inner Zone" not in html


def test_comets_sort_in_by_semi_major_axis():
    planets = [_body(1, "Inner", 1.0, index=0), _body(2, "Outer", 10.0, index=1)]
    comets = [
        _comet(7, "Hyperbolic", 0.5, 1.0, orbit_type="parabolic"),
        # perihelion 1 AU, e = 0.8 -> a = 5 AU: between Inner and Outer.
        _comet(8, "Middling", 1.0, 0.8),
    ]
    html = system_list_html(_system(planets, comets=comets), _SECTIONS)
    order = [html.index(name) for name in ("Inner", "Middling", "Outer", "Hyperbolic")]
    assert order == sorted(order)


def test_comet_orbit_key():
    assert comet_orbit_key_km(_comet(1, "a", 1.0, 0.5)) == (0, 2.0 * AU_KM)
    assert comet_orbit_key_km(_comet(1, "b", 1.0, 1.0, orbit_type="parabolic")) == (1, AU_KM)


_SPRITE = "/static/icons.svg"


def _scored(body, phi4=0.8, **tiers):
    values = {f"tier_{domain}": 0 for domain in ("pressure", "temperature", "air", "water", "radiation", "suitability")}
    body.update(equipment_tier=1, phi4=phi4, **{**values, **{f"tier_{k}": v for k, v in tiers.items()}})
    return body


def test_a_scored_body_shows_its_symbols_and_no_notes():
    """UX.100: the row carries the overall symbol and the six domains, coloured by rating, and none of the old notes."""
    body = _scored(_body(1, "Rocky", 1.0), pressure=0, temperature=1, air=2, water=0, radiation=3, suitability=1)
    html = _planet_row_html(body, {"planets": {}, "moons": {}}, sprite_url=_SPRITE)
    for label in ("Overall habitability: Green", "Pressure: Blue", "Temperature: Green", "Air chemistry: Yellow",
                  "Water chemistry: Blue", "Radiation: Red", "Suitability: Green"):
        assert f'aria-label="{label}"' in html, label
    assert html.count("hab-symbol ") == 7
    assert html.index("Overall habitability") < html.index("Pressure: Blue")
    # The notes the symbols replace: the PHI-6 tooltip, the colour chips with words and the three scores.
    for gone in ("PHI-6", "habitability-detail", "Microbial", "Complex life", "Human operability", "Pressure: Blue</span>"):
        assert gone not in html, gone
    # The symbols sit in the row's summary line with the other stats.
    assert html.index("hab-symbols") < html.index("</summary>")


def test_unscored_bodies_and_pages_without_a_sprite_draw_no_symbols():
    unscored = _body(1, "Giant", 5.0, body_type="g")
    assert "hab-symbol" not in _planet_row_html(unscored, {"planets": {}, "moons": {}}, sprite_url=_SPRITE)
    body = _scored(_body(2, "Rocky", 1.0))
    assert "hab-symbol" not in _planet_row_html(body, {"planets": {}, "moons": {}})


def test_a_moon_row_carries_its_symbols_too():
    moon = _scored(_body(9, "Pandora", 0.002), radiation=2)
    html = system_list_html(_system([_body(1, "Plain", 1.0, moons=[moon])]), _SECTIONS, sprite_url=_SPRITE)
    assert html.count('aria-label="Radiation: Yellow"') == 1


def test_every_body_row_has_an_anchor():
    """NAV.8: a planet, moon, belt and comet each have a row `id` that /system/<id>#<kind>-<id> opens."""
    belt = {"id": 5, "density": "sparse", "distance_km": AU_KM, "lower_limit_km": None, "upper_limit_km": None,
            "orbital_index": 1, "star_id": 1, "composition": []}
    html = system_list_html(_system([_body(7, "Polyphemus", 5.0, moons=[_body(9, "Pandora", 0.002)])], [belt],
                                    [_comet(3, "Halley", 0.6, 0.9)]), _SECTIONS)
    for anchor in ("planet-7", "moon-9", "belt-5", "comet-3"):
        assert html.count(f'<li id="{anchor}"') == 1, anchor
    # A single star has no row of its own: its link opens the system page.
    assert 'id="star-' not in html
