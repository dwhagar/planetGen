# tests/test_systempage.py

"""
The system page's body list (`html/lib/systempage.py`, TODO items 2 and
48): one type chip, the habitable-moon chip, moons in their own group,
belt distances, and comets ordered in among the planets by semi-major
axis. Plain dicts in the `queryDb.system_detail` shape; no database.
"""

import os
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

from systempage import comet_orbit_key_km, system_list_html  # noqa: E402

AU_KM = 149_597_870.7
_SECTIONS = {"overview": "", "stars": {}, "planets": {}, "moons": {}, "belts": {}, "comets": {}}


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


def test_one_type_chip_per_planet():
    html = system_list_html(_system([
        _body(1, "Rocky", 1.0, habitable=True, index=0),
        _body(2, "Dry", 1.5, index=1),
        _body(3, "Giant", 5.2, body_type="g", planet_class="J", index=2),
    ]), _SECTIONS)
    rocky = html[html.index("Rocky"):html.index("Dry")]
    assert "Habitable" in rocky and "Terrestrial" not in rocky
    dry = html[html.index("Dry"):html.index("Giant")]
    assert "Terrestrial" in dry and "Habitable" not in dry
    assert "Gas Giant" in html[html.index("Giant"):]
    assert "Habitable: " not in html


def test_habitable_moon_chip_and_moon_group():
    moon = _body(9, "Pandora", 0.002, habitable=True)
    html = system_list_html(_system([_body(1, "Polyphemus", 5.0, body_type="g", moons=[moon])]), _SECTIONS)
    assert "Habitable moon" in html
    assert '<details class="moon-group"><summary>1 moon of Polyphemus</summary>' in html
    # The group follows the planet's own row rather than sitting inside it.
    assert html.index("</details>") < html.index('class="moon-group"')
    no_moon = system_list_html(_system([_body(1, "Plain", 1.0, moons=[_body(9, "Rock", 0.002)])]), _SECTIONS)
    assert "Habitable moon" not in no_moon


def test_belt_row_shows_its_distance():
    belt = {"id": 5, "density": "sparse", "distance_km": 2.7 * AU_KM, "lower_limit_km": 2.1 * AU_KM,
            "upper_limit_km": 3.3 * AU_KM, "orbital_index": 0, "star_id": 1}
    html = system_list_html(_system(belts=[belt]), _SECTIONS)
    assert "2.7 AU" in html and "2.1 AU to 3.3 AU" in html


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
