# tests/test_class_pages.py

"""
The class reference pages (`web/class_pages.py`, `planetgen/web/lib/classref.py`):
every type and class page renders, the catalog's values come from the
generator's own tables, an unknown type or code is a 404, star types
parse into their spectral and luminosity classes, and the system and
phenomenon pages link their class labels here. No database: the system
and phenomenon pages' data layer is faked like `test_web_system_phen.py`.
"""

import pytest

from planetgen.api.app import create_app
from planetgen.api.config import Config

from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.web.lib import classref  # noqa: E402
from markupsafe import escape  # noqa: E402
from planetgen.physics import constants  # noqa: E402
from planetgen import tuning
from planetgen.generation.comet import PERIOD_CLASS_LABELS  # noqa: E402
from planetgen.generation.phenomena.nebula import NEBULA_CLASS_LETTERS, REMNANT_CLASS_LETTERS  # noqa: E402
from planetgen.physics.stellar_evolution import YERKES_CLASS_NAMES  # noqa: E402
from planetgen.web.lib.systempage import stars_html, system_list_html  # noqa: E402

DB = "planetgen_web_test"


class _Config(Config):
    WEB_DATABASE = DB
    SESSION_COOKIE_SECURE = False
    SECRET_KEY = "test-secret"


@pytest.fixture(scope="module")
def app():
    application = create_app(_Config)
    application.testing = True
    return application


@pytest.fixture
def client(app):
    return app.test_client()


def _html(client, path):
    resp = client.get(path)
    assert resp.status_code == 200, path
    return resp.get_data(as_text=True)


# --- Every page renders ---------------------------------------------------------------

def test_index_lists_every_type_with_its_count(client):
    html = _html(client, "/classes")
    for slug, entry in classref.catalog().items():
        assert f'<a href="/classes/{slug}">{escape(entry["name"])}</a>' in html
        assert f"<td>{len(entry['classes'])}</td>" in html
    assert '<a href="/classes" aria-current="page">Classes</a>' in html


def test_every_type_and_class_page_renders(client):
    for slug, entry in classref.catalog().items():
        html = _html(client, f"/classes/{slug}")
        for code, item in entry["classes"].items():
            assert f'href="/classes/{slug}/{code}"' in html
            page = _html(client, f"/classes/{slug}/{code}")
            # Breadcrumb back to the type, and every fact shown.
            assert f'<li><a href="/classes/{slug}">{escape(entry["name"])}</a></li>' in page
            for label, text in item["facts"]:
                assert f'<th scope="row">{escape(label)}</th><td>{escape(text)}</td>' in page
            assert "<script>" not in page and 'style="' not in page


def test_unknown_type_or_code_is_404(client):
    assert client.get("/classes/galaxy").status_code == 404
    assert client.get("/classes/nebula/Z").status_code == 404
    assert client.get("/classes/nebula/R").status_code == 404  # R is a remnant, not a nebula
    assert client.get("/classes/star-luminosity/D").status_code == 404  # an alias, not its own page
    assert client.get("/classes/asteroid-field/C3").status_code == 404  # the page is the letter's
    assert classref.class_url_parts("nebula", "Z") is None
    assert classref.class_url_parts("galaxy", "A") is None
    assert classref.class_url_parts("planet", None) is None


# --- The catalog is the generator's own tables -------------------------------------------

def test_catalog_types_cover_the_tables():
    cat = classref.catalog()
    assert list(cat) == ["star-spectral", "star-luminosity", "planet", "nebula", "supernova-remnant",
                         "asteroid-field", "black-hole", "rogue-planet", "comet"]
    assert list(cat["star-spectral"]["classes"]) == list(constants.SPECTRAL_CLASS_COLORS)
    assert set(cat["star-luminosity"]["classes"]) == set(YERKES_CLASS_NAMES) - {"D"}
    assert list(cat["planet"]["classes"]) == list(tuning.PLANET_CLASSES)
    assert tuple(cat["nebula"]["classes"]) == NEBULA_CLASS_LETTERS
    assert tuple(cat["supernova-remnant"]["classes"]) == REMNANT_CLASS_LETTERS
    letters = {letter for family in tuning.ASTEROID_FIELD_COMPOSITIONS.values()
               for letter in family["letters"].values()}
    assert set(cat["asteroid-field"]["classes"]) == letters
    assert tuple(cat["black-hole"]["classes"]) == tuning.BLACK_HOLE_MASS_CLASSES
    assert tuple(cat["rogue-planet"]["classes"]) == tuning.ROGUE_PLANET_MASS_BIN_CHOICES
    assert list(cat["comet"]["classes"]) == list(tuning.COMET_PERIOD_CLASSES)
    # Built once and cached.
    assert classref.catalog() is cat


def test_class_pages_show_the_constants(client):
    nebula = tuning.NEBULA_CLASSES["D"]
    html = _html(client, "/classes/nebula/D")
    assert str(escape(nebula["name"])) in html
    assert str(escape(nebula["species"])) in html

    remnant = tuning.NEBULA_CLASSES["T"]
    assert str(escape(remnant["name"])) in _html(client, "/classes/supernova-remnant/T")

    planet = tuning.PLANET_CLASSES["M"]
    html = _html(client, "/classes/planet/M")
    low, high = planet["radius_range"]
    assert f"{classref.number(low)} to {classref.number(high)} km" in html  # UX.20: 1 × 10⁴ km
    assert "<th scope=\"row\">Habitable</th><td>Yes</td>" in html
    assert planet["composition"][1:] in html

    low_k, high_k = constants.TEMP_RANGES["G"]
    html = _html(client, "/classes/star-spectral/G")
    assert f"{classref.number(low_k)} to {classref.number(high_k)} K" in html
    assert constants.SPECTRAL_CLASS_COLORS["G"] in html

    html = _html(client, "/classes/star-luminosity/IA+")
    assert YERKES_CLASS_NAMES["IA+"] in html

    html = _html(client, "/classes/comet/halley_type")
    assert PERIOD_CLASS_LABELS["halley_type"] in html

    family = tuning.ASTEROID_FIELD_COMPOSITIONS["metallic"]
    html = _html(client, f"/classes/asteroid-field/{family['letters']['dense']}")
    assert family["description"][1:] in html
    # One letter can stand for several densities.
    basaltic = tuning.ASTEROID_FIELD_COMPOSITIONS["basaltic"]["letters"]
    assert basaltic["sparse"] == basaltic["typical"]
    facts = dict(classref.class_entry("asteroid-field", basaltic["sparse"])["facts"])
    assert facts["Density"] == "Sparse, typical"
    html = _html(client, "/classes/asteroid-field")
    assert "10^d to 10^(d+1) AU" in html

    low, high = tuning.BLACK_HOLE_INTERMEDIATE_MASS_RANGE_SOLAR
    assert f"{classref.number(low)} to {classref.number(high)} solar masses" in _html(client, "/classes/black-hole/intermediate")

    low, high = tuning.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER
    assert f"{low:.0f} to {high:.0f} Jupiter masses" in _html(client, "/classes/rogue-planet/brown-dwarf")


def test_number_formatting():
    assert classref.number(150) == "150"
    assert classref.number(4131.4) == "4,131.4"
    assert classref.number(0.08) == "0.08"
    assert classref.number(2_000_000) == "2 × 10⁶"
    assert classref.number(1e-5) == "1 × 10⁻⁵"
    assert classref.number(0) == "0"


# --- Star types --------------------------------------------------------------------------

@pytest.mark.parametrize("star_type, expected", [
    ("G2V", ("G", "V")),
    ("G2V Yellow Main Sequence Star", ("G", "V")),
    ("B0IA", ("B", "IA")),
    ("B0IA+ Blue-White Luminous Supergiant Star", ("B", "IA+")),
    ("K5IAB", ("K", "IAB")),
    ("O5VII", ("O", "VII")),
    ("M3D", ("M", "VII")),  # D is a white dwarf, like VII
    ("A1IV", ("A", "IV")),
    ("M9VI", ("M", "VI")),
    ("F0 0", None),
    ("B00", ("B", "0")),
    ("g2v", ("G", "V")),
    ("G10V", None),
    ("X2V", None),
    ("", None),
    (None, None),
])
def test_star_type_classes(star_type, expected):
    assert classref.star_type_classes(star_type) == expected


def test_class_url_parts_resolve_aliases():
    assert classref.class_url_parts("star-luminosity", "D") == ("star-luminosity", "VII")
    assert classref.class_url_parts("asteroid-field", "C3") == ("asteroid-field", "C")
    assert classref.class_url_parts("nebula", "D") == ("nebula", "D")


# --- Links from the system and phenomenon pages ----------------------------------------------

_SECTIONS = {"overview": "", "stars": {}, "planets": {}, "moons": {}, "belts": {}, "comets": {}}


def _planet(planet_class="M"):
    return {"id": 2, "name": "Terra", "planet_class": planet_class, "body_type": "t", "habitable": True,
            "inhabited": False, "zone": "Habitable Zone", "distance_km": 1.5e8, "period_years": 1.0,
            "gravity_g": 1.0, "orbital_index": 0, "moons": [], "star_id": 1, "radius_km": 6371.0}


def _system():
    star = {"id": 1, "name": "Sol", "role": "single", "star_type": "G2V Yellow Main Sequence Star",
            "mass_kg": 1.989e30, "radius_km": 696_000.0, "temperature_k": 5778.0, "luminosity_w": 3.828e26}
    comet = {"id": 3, "name": "Kohoutek", "orbit_type": "elliptical", "period_class": "halley_type",
             "is_active": True, "perihelion_distance_km": 1e8, "eccentricity": 0.9, "orbital_period_years": 70.0,
             "star_id": 1}
    return {"stars": [star], "planets": [_planet()], "belts": [], "comets": [comet], "binary_configuration": None}


def test_system_list_links_classes(app):
    from web.class_pages import class_url
    with app.test_request_context("/system/1"):
        html = system_list_html(_system(), _SECTIONS, class_url)
        table = stars_html(_system()["stars"], class_url)
    assert '<a href="/classes/star-spectral/G">Spectral class G</a>' in html
    assert '<a href="/classes/star-luminosity/V">Luminosity class V</a>' in html
    assert '<a href="/classes/planet/M">Planet class M</a>' in html
    assert '<a href="/classes/comet/halley_type">Halley-type comet</a>' in html
    assert '<a href="/classes/star-spectral/G">G2V Yellow Main Sequence Star</a>' in table
    # Never inside a row's <summary> (a link in the disclosure button).
    for summary in html.split("<summary>")[1:]:
        assert "<a " not in summary.split("</summary>")[0]


def test_system_list_without_hook_or_known_class_stays_text(app):
    html = system_list_html(_system(), _SECTIONS)
    assert "/classes/" not in html and "Class M" in html
    from web.class_pages import class_url
    system = _system()
    system["planets"] = [_planet("Z")]
    system["stars"][0]["star_type"] = "Unknown"
    with app.test_request_context("/system/1"):
        html = system_list_html(system, _SECTIONS, class_url)
        table = stars_html(system["stars"], class_url)
    assert "/classes/planet" not in html and "/classes/star-" not in html + table


class _Fake:
    def __init__(self):
        self.phenomenon = None
        self.system = {
            "id": 5, "name": "Sol", "sector_id": None, "quadrant": None, "location": None, "is_binary": 0,
            "binary_type": None, "binary_configuration": None, "wikijs_url": None, "mediawiki_url": None,
            **{key: value for key, value in _system().items() if key != "binary_configuration"},
        }

    def get_system(self, db, system_id):
        return self.system

    def get_system_sections(self, db, system_id):
        return _SECTIONS

    def get_system_facilities(self, db, system_id):
        return []

    def get_phenomenon(self, db, phenomenon_type, phenomenon_id):
        return self.phenomenon

    def auth_me(self, cookie_header):
        return None


@pytest.fixture
def fake(monkeypatch):
    data = _Fake()
    for name in ("get_system", "get_system_sections", "get_system_facilities", "get_phenomenon", "auth_me"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    return data


def test_system_page_links_classes(client, fake):
    html = _html(client, "/system/5")
    assert '<a href="/classes/star-spectral/G">G2V Yellow Main Sequence Star</a>' in html
    assert '<a href="/classes/planet/M">Planet class M</a>' in html


@pytest.mark.parametrize("kind, detail, label, url", [
    ("nebula", {"nebula_class": "D"}, "D: " + tuning.NEBULA_CLASSES["D"]["name"], "/classes/nebula/D"),
    ("supernova_remnant", {"remnant_class": "T"}, "T: " + tuning.NEBULA_CLASSES["T"]["name"],
     "/classes/supernova-remnant/T"),
    ("asteroid_field", {"field_class": "C3"}, "C3", "/classes/asteroid-field/C"),
    ("black_hole", {"mass_class": "supermassive"}, "Supermassive", "/classes/black-hole/supermassive"),
    ("rogue_planet", {"mass_bin": "brown-dwarf"}, "Brown dwarf", "/classes/rogue-planet/brown-dwarf"),
])
def test_phenomenon_page_links_its_class(client, fake, kind, detail, label, url):
    fake.phenomenon = {"id": 4, "name": "Thing", "sector_id": None, **detail}
    html = _html(client, f"/phenomenon/{kind}/4")
    assert f'<td><a href="{url}">{escape(label)}</a></td>' in html


def test_phenomenon_unknown_class_stays_text(client, fake):
    fake.phenomenon = {"id": 4, "name": "Thing", "sector_id": None, "nebula_class": "Z"}
    html = _html(client, "/phenomenon/nebula/4")
    assert "<td>Z</td>" in html and "/classes/" not in html.split("</header>", 1)[1]
