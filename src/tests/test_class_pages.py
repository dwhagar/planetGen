# tests/test_class_pages.py

"""
The class reference pages (`web/class_pages.py`, `planetgen/web/lib/classref.py`):
every type and class page renders, the catalog's values come from the
generator's own tables, an unknown type or code is a 404, star types
parse into their spectral and luminosity classes, and the system and
phenomenon pages link their class labels here. No database: the system
and phenomenon pages' data layer is faked like `test_web_system_phen.py`.
"""

import re

import pytest

from planetgen.web.app import create_app
from planetgen.api.config import Config

from planetgen.describe import model  # noqa: E402
from planetgen.web.lib import apiclient  # noqa: E402
from planetgen.datasets import mamajek_sequence  # noqa: E402
from planetgen.web.lib import classref  # noqa: E402
from markupsafe import escape  # noqa: E402
from planetgen.physics import constants, planet_classes  # noqa: E402
from planetgen import tuning
from planetgen.util import chem
from planetgen.generation.comet import PERIOD_CLASS_LABELS  # noqa: E402
from planetgen.generation.phenomena.nebula import NEBULA_CLASS_LETTERS, REMNANT_CLASS_LETTERS  # noqa: E402
from planetgen.physics.stellar_evolution import YERKES_CLASS_NAMES  # noqa: E402
from planetgen.web.lib import objecttree  # noqa: E402

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
    # UX.180: Classes is under the Help menu, which is current on it.
    assert '<a href="/classes" aria-current="page">Classes</a>' not in html
    assert 'class="site-help-trigger" aria-haspopup="menu" aria-current="page"><span class="site-help-text">Help</span>' in html
    assert 'data-href="/classes" aria-current="page">Classes</sl-menu-item>' in html


def test_index_lists_the_types_alphabetically(client):
    """UX.77: the class categories are in alphabetical order."""
    html = _html(client, "/classes")
    names = [escape(entry["name"]) for entry in classref.catalog().values()]
    positions = [html.index(f'">{name}</a>') for name in sorted(names, key=str.casefold)]
    assert positions == sorted(positions), "the index is not in alphabetical order"


def test_alphabetical_ignores_case():
    rows = [{"name": "star"}, {"name": "Nebula"}, {"name": "asteroid"}]
    assert [r["name"] for r in classref.alphabetical(rows)] == ["asteroid", "Nebula", "star"]


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
    assert list(cat["planet"]["classes"]) == sorted(planet_classes.LETTERS)
    assert tuple(cat["nebula"]["classes"]) == NEBULA_CLASS_LETTERS
    assert tuple(cat["supernova-remnant"]["classes"]) == REMNANT_CLASS_LETTERS
    letters = {letter for family in tuning.ASTEROID_FIELD_COMPOSITIONS.values()
               for letter in family["letters"].values()}
    assert set(cat["asteroid-field"]["classes"]) == letters
    assert tuple(cat["black-hole"]["classes"]) == tuning.BLACK_HOLE_MASS_CLASSES
    assert tuple(cat["rogue-planet"]["classes"]) == tuning.ROGUE_PLANET_MASS_BIN_CHOICES
    assert list(cat["comet"]["classes"]) == [*tuning.COMET_PERIOD_CLASSES, classref.PARABOLIC_COMET_CLASS]
    # Built once and cached.
    assert classref.catalog() is cat


def test_class_pages_show_the_constants(client):
    nebula = tuning.NEBULA_CLASSES["D"]
    html = _html(client, "/classes/nebula/D")
    assert str(escape(nebula["name"])) in html
    assert chem.mark_html(nebula["species"]) in html     # UX.149: its formulas written out

    remnant = tuning.NEBULA_CLASSES["T"]
    assert str(escape(remnant["name"])) in _html(client, "/classes/supernova-remnant/T")

    html = _html(client, "/classes/planet/M")  # GEN.645: a class page gives the rule
    assert planet_classes.describe("M").sentence in html
    assert "<th scope=\"row\">Core class</th><td>Yes</td>" in html
    html = _html(client, "/classes/planet/S")
    assert "surface pressure below 10 Pa; mass at least 1.5 Earth masses" in html

    low_k, high_k = mamajek_sequence.letter_range_k("G")
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
    assert classref.number_range(low, high, "solar masses") in _html(client, "/classes/black-hole/intermediate")

    low, high = tuning.ROGUE_BROWN_DWARF_MASS_RANGE_JUPITER
    # The upper limit is the mass boundaries table's 83.8 Jupiter masses (`brown_dwarf_upper_limit`), which the page's number rule
    # shows as "13 to 83.8" (rounding it to a whole number, as this line did, read 84 and missed). The page reads the generator's
    # own range, so the test asks `number_range` for the same text as the lines above.
    assert classref.number_range(low, high, "Jupiter masses") in _html(client, "/classes/rogue-planet/brown-dwarf")


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

_SECTIONS = {"overview": [], "stars": {}, "planets": {}, "moons": {}, "belts": {}, "comets": {}}


def _planet(planet_class="M"):
    return {"id": 2, "name": "Terra", "planet_class": planet_class, "body_type": "t", "habitable": True,
            "habitable_class": True, "inhabited": False, "zone": "Habitable Zone", "distance_km": 1.5e8, "period_years": 1.0,
            "gravity_g": 1.0, "orbital_index": 0, "moons": [], "star_id": 1, "radius_km": 6371.0}


def _system():
    star = {"id": 1, "name": "Sol", "role": "single", "star_type": "G2V Yellow Main Sequence Star",
            "mass_kg": 1.989e30, "radius_km": 696_000.0, "temperature_k": 5778.0, "luminosity_w": 3.828e26}
    comet = {"id": 3, "name": "Kohoutek", "orbit_type": "elliptical", "period_class": "halley_type",
             "is_active": True, "perihelion_distance_km": 1e8, "eccentricity": 0.9, "orbital_period_years": 70.0,
             "star_id": 1}
    return {"stars": [star], "planets": [_planet()], "belts": [], "comets": [comet], "binary_configuration": None}


def _tree_context(class_url=None, **more):
    return objecttree.TreeContext(item_url=lambda ref: "/object/" + ref, tree_url=lambda ref: "/object/" + ref + "/tree",
                                  class_url=class_url, **more)


def _items(system):
    """The tree rows (`GET /api/objects/<ref>/children`'s items) of a `_system()`'s star, planets and comets."""
    star = system["stars"][0]
    items = [{"ref": f"star:{star['id']}", "kind": "star", "name": star["name"], "child_count": 0, "counts": {},
              "summary": {"star_type": star["star_type"], "role": star["role"]}}]
    items += [{"ref": f"planet:{p['id']}", "kind": "planet", "name": p["name"], "child_count": 0, "counts": {},
               "summary": {"planet_class": p["planet_class"]}} for p in system["planets"]]
    items += [{"ref": f"comet:{c['id']}", "kind": "comet", "name": c["name"], "child_count": 0, "counts": {},
               "summary": {"orbit_type": c["orbit_type"], "period_class": c["period_class"]}} for c in system["comets"]]
    return items


def test_the_tree_links_classes(app):
    from planetgen.web.class_pages import class_url
    with app.test_request_context("/system/1"):
        html = objecttree.nodes_html(_items(_system()), _tree_context(class_url))
    # UX.182: the class letter or name is the link, with no separate line of class links.
    assert "Planet class M" not in html and "Spectral class" not in html and 'class="class-links"' not in html
    assert 'Class <a href="/classes/planet/M">M</a>' in html
    assert '<a href="/classes/comet/halley_type">Halley-type comet</a>' in html
    assert ('<a href="/classes/star-spectral/G">G2</a><a href="/classes/star-luminosity/V">V</a>'
            " Yellow Main Sequence Star") in html


@pytest.mark.parametrize("orbit_type, period_class", [
    *[("elliptical", code) for code in tuning.COMET_PERIOD_CLASSES], ("parabolic", None),
])
def test_every_comet_links_its_class(app, client, orbit_type, period_class):
    """UX.29: a parabolic comet (no period class) links one too, and
    every link opens a class page."""
    from planetgen.web.class_pages import class_url
    system = _system()
    system["comets"][0].update(orbit_type=orbit_type, period_class=period_class)
    with app.test_request_context("/system/1"):
        html = objecttree.nodes_html(_items(system), _tree_context(class_url))
    code = period_class or classref.PARABOLIC_COMET_CLASS
    name = classref.class_entry("comet", code)["name"]
    assert f'<a href="/classes/comet/{code}">{escape(name)}</a>' in html
    assert name in _html(client, f"/classes/comet/{code}")


def test_the_tree_without_hook_or_known_class_stays_text(app):
    html = objecttree.nodes_html(_items(_system()), _tree_context())
    assert "/classes/" not in html and "Class M" in html
    from planetgen.web.class_pages import class_url
    system = _system()
    system["planets"] = [_planet("Z")]
    system["stars"][0]["star_type"] = "Unknown"
    with app.test_request_context("/system/1"):
        html = objecttree.nodes_html(_items(system), _tree_context(class_url))
    assert "/classes/planet" not in html and "/classes/star-" not in html


class _Fake:
    def __init__(self):
        self.phenomenon = None
        self.system = {
            "id": 5, "name": "Sol", "sector_id": None, "quadrant": None, "location": None, "is_binary": 0,
            "binary_type": None, "binary_configuration": None,
            **{key: value for key, value in _system().items() if key != "binary_configuration"},
        }

    def get_system(self, db, system_id):
        return self.system

    def get_system_description(self, db, system_id):
        return model.document("system", "system:1", "Sol", [model.section("system", "system:1", "Sol")])

    def get_system_facilities(self, db, system_id):
        return []

    def get_object_children(self, db, ref, limit=None, after=None):
        star = "star:FE81000A2B-0000005-000"
        rows = {
            "system:5": [{"ref": star, "kind": "star", "name": "Sol", "child_count": 2, "counts": {"planet": 1, "comet": 1},
                          "summary": {"star_type": "G2V Yellow Main Sequence Star", "role": "single"}}],
            star: [{"ref": "planet:FE81000A2B-0000005-001", "kind": "planet", "name": "Terra", "child_count": 0, "counts": {},
                    "summary": {"planet_class": "M"}},
                   {"ref": "comet:FE81000A2B-0000005-002", "kind": "comet", "name": "Kohoutek", "child_count": 0,
                    "counts": {}, "summary": {"orbit_type": "elliptical", "period_class": "halley_type"}}],
        }.get(ref, [])
        return {"ref": ref, "kind": ref.partition(":")[0], "total": len(rows), "items": rows, "next": None}

    def get_phenomenon(self, db, phenomenon_type, phenomenon_id):
        return self.phenomenon

    def get_phenomenon_contents(self, db, phenomenon_type, phenomenon_id, limit=None):
        return {"total": 0, "counts": {}, "items": []}

    def get_object_sheet(self, db, ref):
        return {"kind": ref.split(":", 1)[0], "groups": []}

    def auth_me(self, cookie_header):
        return None


@pytest.fixture
def fake(monkeypatch):
    data = _Fake()
    for name in ("get_system", "get_system_description", "get_system_facilities", "get_object_children", "get_object_sheet", "get_phenomenon", "get_phenomenon_contents", "auth_me"):
        monkeypatch.setattr(apiclient, name, getattr(data, name))
    return data


def test_system_page_links_classes(client, fake):
    html = _html(client, "/system/5")
    assert '<a href="/classes/star-spectral/G">G2</a><a href="/classes/star-luminosity/V">V</a>' in html
    assert 'Class <a href="/classes/planet/M">M</a>' in html
    assert '<a href="/classes/comet/halley_type">Halley-type comet</a>' in html


@pytest.mark.parametrize("kind, detail, label, url", [
    ("nebula", {"nebula_class": "D", "nebula_type": tuning.NEBULA_CLASSES["D"]["family"]}, "D: " + tuning.NEBULA_CLASSES["D"]["name"], "/classes/nebula/D"),
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


# --- UX.90: the habitability levels page --------------------------------------------

def test_habitability_levels_page_explains_the_five_levels_and_the_six_domains(client):
    from planetgen.web.lib import habitability_ref

    html = _html(client, "/help/habitability")
    chem = client.application.jinja_env.filters["chem"]
    for item, _when in habitability_ref.EQUIPMENT_ITEMS:
        assert str(chem(item)) in html
    assert "Shirtsleeve" not in html
    for domain in ("Pressure", "Temperature", "Air", "Water", "Radiation", "Suitability"):
        assert f'<h3 id="domain-{domain.lower()}">{domain}</h3>' in html
    for colour in ("Blue", "Green", "Yellow", "Red"):
        assert f"<dt>{colour}</dt>" in html
    for chip in ("Habitable", "Habitable moon", "Habitable class", "Habitable-class moon", "Inhabited"):
        assert f'<dt id="chip-{chip.lower().replace(" ", "-")}">{chip}</dt>' in html
    assert "pulsar" in html and "1,000 Sv a year" in html


def test_habitability_levels_page_draws_each_symbol_in_every_rating_colour(client):
    from planetgen.web.lib import habitability_symbols

    html = _html(client, "/help/habitability")
    for name, label in habitability_symbols.LABELS.items():
        for tier, colour in enumerate(("Blue", "Green", "Yellow", "Red")):
            assert f'aria-label="{label}: {colour}"' in html, (name, tier)
            assert f'#hab-{name}"' in html


def test_habitability_levels_page_draws_the_four_life_scores(client):
    html = _html(client, "/help/habitability")
    for label in ("Microbial life", "Flora", "Fauna", "Sentient life"):
        assert f'aria-label="{label}: none"' in html and f'aria-label="{label}: 9 of 9"' in html
    for name in ("microbial", "flora", "fauna", "sentient"):
        assert f'#life-{name}"' in html


def test_habitability_levels_page_says_what_sets_each_life_score_with_the_numbers_from_the_tuning(client):
    from html import unescape

    from planetgen import tuning
    from planetgen.web.lib import habitability_ref, habitability_symbols

    html = _html(client, "/help/habitability")
    assert "What sets it" in html
    for _name, needs in habitability_ref.LIFE_NEEDS:
        assert needs in unescape(html)
    assert f"stays at {tuning.LIFE_SCORE_SIMPLE_FLORA_MAX} or less" in html
    assert f"in full after {tuning.LIFE_SCORE_SENTIENT_FULL_GYR * 1000:g} million years" in html
    assert "Flora together with Fauna is better" in html
    assert [name for name, _ in habitability_ref.LIFE_NEEDS] == list(habitability_symbols.LIFE_NAMES)


def test_habitability_levels_page_says_where_the_symbols_appear(client):
    from planetgen.web.lib import habitability_ref

    html = _html(client, "/help/habitability")
    assert 'id="where-heading"' in html
    for place, _what in habitability_ref.WHERE_ROWS:
        assert place in html
    assert "/lists/planet" in html and "/lists/moon" in html


def test_habitability_levels_page_states_the_habitable_zone_from_the_generators_edges(client):
    from planetgen.physics import orbits
    from planetgen.web.lib import habitability_ref

    inner, outer = orbits.calculate_habitable_zone(constants.SOLAR_LUMINOSITY, orbits.KOPPARAPU_REFERENCE_TEFF_K)
    text = habitability_ref.zone_text()
    assert f"{inner:.2f} and {outer:.2f} AU" in text
    assert (inner, outer) == pytest.approx((0.95, 1.68), abs=0.01)
    low, high = orbits.KOPPARAPU_TEFF_RANGE_K
    assert f"{low:,.0f} to {high:,.0f} K" in text
    light, heavy = orbits.KOPPARAPU_MASS_RANGE_EARTH
    light_inner = orbits.hz_limits(constants.SOLAR_LUMINOSITY, orbits.KOPPARAPU_REFERENCE_TEFF_K, light)[0]
    heavy_inner = orbits.hz_limits(constants.SOLAR_LUMINOSITY, orbits.KOPPARAPU_REFERENCE_TEFF_K, heavy)[0]
    assert f"{light_inner:.2f} AU at {light:g} Earth masses and {heavy_inner:.2f} AU at {heavy:g}" in text
    assert heavy_inner < inner < light_inner
    assert "lost water" in text
    html = _html(client, "/help/habitability")
    assert 'id="zone-heading"' in html and "runaway greenhouse" in html


def test_habitability_levels_page_says_what_sets_each_radiation_term_with_the_scores_limits(client):
    from planetgen.physics.protection_limits import limit
    from planetgen.web.lib import habitability_ref

    text = habitability_ref.radiation_text()
    for word in ("cosmic rays", "Stellar particles", "potassium, thorium and uranium", "ozone"):
        assert word.lower() in text.lower()
    assert f"Blue under {limit('dose_blue_sv_yr') * 1000:g} mSv" in text
    assert f"Yellow to {limit('uv_lethal_index'):g} times" in text
    html = _html(client, "/help/habitability")
    assert 'id="radiation-heading"' in html and "Where the radiation comes from" in html


def test_the_habitability_levels_numbers_are_the_constants_the_score_uses(client):
    """GEN.219: no number on the page is typed in; it is the limit the score reads."""
    from planetgen.physics import habitability
    from planetgen.physics.protection_limits import limit
    from planetgen.web.lib import habitability_ref

    rows = {variable: cells for _domain, table in habitability_ref.domains() for variable, *cells in table}
    assert rows["Air pressure"] == ["50 to 250 kPa", "14.3 to 400 kPa", "6.3 to 10,000 kPa"]
    assert rows["Air temperature"] == ["0 to 31 C", "-20 to 45 C", "-50 to 122 C"]
    assert rows["Carbon dioxide"][0] == f"under {habitability.GAS_LIMITS_KPA['CO2'][0]:g} kPa"
    assert rows["Chaotropicity (how strongly dissolved salts pull the water's structure apart)"][0] == (
        f"under {habitability.CHAOTROPICITY_CRIT_KJ_KG:g} kJ/kg")
    assert rows["Water activity (how freely water is available)"][2] == f"{habitability.WATER_ACTIVITY_CRIT:g} or more"
    assert f"{limit('heat_shielded_c'):g} C" in habitability_ref.EQUIPMENT_RULE
    # A changed limit moves the page.
    wide = habitability.bands({"pressure_no_gear_kpa": 60.0, "dose_blue_sv_yr": 0.03})
    moved = {variable: cells for _domain, table in habitability_ref.domains(wide) for variable, *cells in table}
    assert moved["Air pressure"][0] == "60 to 250 kPa" and moved["Dose per year"][0] == "under 0.03 Sv"
    html = _html(client, "/help/habitability")
    for text in ("Chaotropicity", "too much", "too little", "Hydrogen sulfide"):
        assert text in html


def test_the_help_menu_holds_classes_and_habitability_and_the_old_links_are_gone(client):
    """UX.180: Classes leaves the section row for the Help menu beside the search button, with Habitability."""
    html = _html(client, "/classes")
    bar = html.split('<header class="site-header">')[1].split("</header>")[0]
    assert 'href="/classes"' not in bar.split("site-help-menu")[0]
    help_menu = re.search(r'<sl-dropdown class="site-help-menu".*?</sl-dropdown>', bar, re.S).group(0)
    assert re.findall(r'data-href="([^"]*)"[^>]*>([^<]*)<', help_menu) == [
        ("/help", "All help"), ("/classes", "Classes"), ("/help/habitability", "Habitability"), ("/help/properties", "Properties"), ("/help/life", "Life"),
        ("/references", "References")]
    assert bar.index("site-help-menu") < bar.index("site-search-link") < bar.index("site-gear")
    assert "Habitability levels" not in html
    page = _html(client, "/help/habitability")
    assert 'class="site-help-trigger" aria-haspopup="menu" aria-current="page"><span class="site-help-text">Help</span>' in page
    assert 'data-href="/help/habitability" aria-current="page">Habitability</sl-menu-item>' in page
    assert "<title>Habitability" in page


def _scored_system():
    system = _system()
    system["planets"][0].update({
        "habitable": True, "habitable_class": True, "equipment": {"required": [], "recommended": [], "line": "", "code": "", "items": [
            {"slot": "supply", "item": "gas-mix mask", "severity": "required", "domain": "Air"}]},
        "phi4": 1.0, "tier_pressure": 0, "tier_temperature": 0, "tier_air": 0, "tier_water": 0,
        "tier_radiation": 0, "tier_suitability": 0, "life_microbial": 3, "life_flora": 0, "life_fauna": 0,
        "life_sentient": 0})
    return system


def test_the_tree_symbols_and_the_page_chips_link_to_their_sections_of_the_habitability_page(app, client):
    """UX.181: no separate "Habitability levels" line; each chip is the link (UX.197: the symbols go to the body's page)."""
    from planetgen.web.lib.habitability_chips import chips_html

    url = "/help/habitability"
    planet = _scored_system()["planets"][0]
    item = {"ref": "planet:FE81000A2B-0000005-001", "kind": "planet", "name": "Terra", "child_count": 0, "counts": {},
            "summary": planet}
    html = objecttree.nodes_html([item], _tree_context(sprite_url="/static/icons.svg", help_url=url)) + chips_html(planet, url)
    assert "Habitability levels" not in html
    for anchor in ("chip-habitable", "equipment-gas-mix-mask"):   # UX.224: a habitable world of a habitable class shows the one chip, Habitable
        assert f'<a class="help-link" href="{url}#{anchor}">' in html, anchor
    # UX.197: the tree's symbols and life scores go to the Habitability and Life groups of the body's own page.
    assert html.count('href="/object/planet:FE81000A2B-0000005-001#habitability"') == 7
    assert html.count('href="/object/planet:FE81000A2B-0000005-001#life"') == 4
    page = _html(client, url)
    for anchor in re.findall(r'href="/help/habitability#([^"]+)"', html):
        assert f'id="{anchor}"' in page, anchor
    # A tree built with no help page still links its symbols to the body's page; the page's own legend shows shapes alone.
    assert 'href="/help/habitability' not in objecttree.nodes_html([item], _tree_context(sprite_url="/static/icons.svg"))
    assert "help-link" not in page


def test_every_equipment_item_name_has_a_row_on_the_habitability_page(client):
    """UX.181: a chip's item name reaches a row of the equipment table, never the section's fallback."""
    from planetgen.physics import equipment
    from planetgen.web.lib import habitability_ref, habitability_symbols

    names = {*equipment.BODY_NAMES.values(), *equipment.SUPPLY_NAMES.values(), *equipment.GAS_FILTER_NAMES.values(),
             *equipment.SHELTER_NAMES.values(), *equipment.RAD_SUIT_NAMES.values(), "powered frame", "hardened suit",
             "dust filter", "dust-tight seals", "CO2 scrubber", "hardened sealed suit", "hardened heat-shielded suit"}
    page = _html(client, "/help/habitability")
    for name in names:
        anchor = habitability_symbols.equipment_anchor(name)
        assert anchor != habitability_symbols.EQUIPMENT_FALLBACK, name
        assert f'<tr id="{anchor}">' in page, name
    assert sorted(habitability_ref.EQUIPMENT_ROW_IDS.values()) == sorted(anchor for anchor, _ in habitability_symbols.EQUIPMENT_ANCHORS)
    assert set(habitability_ref.EQUIPMENT_ROW_IDS) == {label for label, _when in habitability_ref.EQUIPMENT_ITEMS}


def test_the_equipment_search_tags_no_longer_carry_a_levels_link(app):
    from planetgen.web import searchpage
    from planetgen.web.searchpage import SearchState

    options = [{"value": "0", "label": "Ideal", "count": 2}]
    with app.test_request_context("/search"):
        groups = searchpage.tag_groups(SearchState(), {"equipment": options})
    assert all("help" not in group for group in groups)


def test_a_range_gives_both_ends_the_one_style():
    """UX.107: a class page range never mixes plain and scientific ("5,000 to 1 × 10⁴ km")."""
    assert classref.number_range(5000, 1e4, "km") == "5 × 10³ to 1 × 10⁴ km"
    assert classref.number_range(100, 900) == "100 to 900"
    assert classref.number_range(0.5, 20000) == "0.5 to 2 × 10⁴"


def test_the_habitability_page_explains_the_map_views_the_temperatures_and_the_chips(client):
    """DOC.75: the Galaxy Map's habitable-world colouring and marker, the regional temperatures and the chips."""
    html = _html(client, "/help/habitability")
    for text in ("On the Galaxy Map", "Color by: Habitable worlds", "show=marker-habitable", "Temperatures across a world",
                 "terminator ring", "solid for a required item, outlined for a recommended one", "0.1 mg/m3"):
        assert text in html, text


def test_a_missing_dataset_does_not_stop_the_site_booting(monkeypatch):
    """update.sh pulls the new code before its downloads finish; a worker that started then died at boot on
    DataMissing (mamajek_dwarf_sequence) and every page was a 500."""
    from planetgen.datasets.fetch import DataMissing

    def missing():
        raise DataMissing("mamajek_dwarf_sequence: not downloaded")

    monkeypatch.setattr(classref, "catalog", missing)
    assert classref.preload() is False
    monkeypatch.setattr(classref, "catalog", lambda: {})
    assert classref.preload() is True
